"""
Prepare the E-HYPE data for SCHISM setup as flux.th and source inputs.

Workflow:
---------
1. Identify SCHISM river boundary segments from `bctides.in`.
2. Match SCHISM rivers with nearest E-HYPE outlets (interactive option available).
3. Export discharge forcing for lateral boundaries (`flux.th`).
4. Select additional E-HYPE outlets within a distance threshold for point sources.
5. Export source/sink files (`source_sink.in`, `vsource.th`, `msource.th`).

Notes:
------
- Lateral boundaries are treated as rivers if `bctides.in` defines them with flag = 1.
- Temperatures and salinity for sources are placeholders (to be replaced with data).
"""

# === Environment setup ===
# conda activate geo_env

# test
import geopandas as gpd

import os
import sys
import datetime as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import num2date, Dataset
from matplotlib import gridspec
from scipy.spatial import cKDTree

# Custom SCHISM tools (make sure to have files)
sys.path.insert(0, '/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import schism_setup
from source_sink_river_functions import create_source_sink, create_vsource, create_msource

# === User inputs ===
outlets = '/work/gg0028/g260114/PROJECTS/FOCCUS/Data/River/EH4_coastal_outlets.txt'
discharge_file = '/work/gg0028/g260114/PROJECTS/FOCCUS/Data/River/netCDF/cout.nc'
rundir = '/work/gg0028/g260114/RUNS/DT_DEMO_RUNS/GermanBight/reference_config/'
catchment_file="/work/gg0028/g260114/PROJECTS/FOCCUS/Data/River/2_EH4_basins_northsea_poly.gpkg" # gpkg format "" if empty
dry_element_file='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Analyse/Veg_HE_dryFlagElement_timestats_second_noq95_3600_31449600.nc' # avoid seeding in dry cells

add_labels = False
start_date = dt.datetime(2017, 1, 2)
end_date = dt.datetime(2018, 1, 1)

method = 'interactive'  # options: 'nn' (nearest neighbour) | 'interactive'

maximum_distance = 0.025  # max distance (deg) for point source inclusion
minimum_depth = 5         # min depth for source insertion (not yet applied)

# === Load E-HYPE data ===

# if catchment file avaialble
catchments = gpd.read_file("2_EH4_basins_northsea_poly.gpkg")

table = pd.read_table(outlets)
dsQ = Dataset(discharge_file)

# Outlet coordinates
lons, lats = table.POURX.values, table.POURY.values

# Time handling
time_var = dsQ.variables['time']
cftime_dates = num2date(time_var[:], units=time_var.units, calendar=time_var.calendar)
time_index = pd.to_datetime([str(d) for d in cftime_dates])
selection = (time_index >= start_date) & (time_index <= end_date)
time_in_seconds = np.asarray(
    (time_index[selection] - time_index[selection][0]) / np.timedelta64(1, 's'),
    int
)

# === Load SCHISM setup ===
cwd = os.getcwd()
os.chdir(rundir)
s = schism_setup()
s.lon, s.lat = np.asarray(s.lon), np.asarray(s.lat)
s.init_node_tree(latlon=True)   # build KD-tree for NN queries

# Identify river boundary indices from bctides.in
with open('bctides.in') as f:
    while True:
        line = f.readline()
        if 'nope' in line:
            nop = int(line.split()[0])
            break

    bds = {}
    bd_nr = 0
    while line:
        parts = line.split()
        if len(parts) > 4 and parts[2].isdigit():
            if int(parts[2]) == 1:
                bds[bd_nr] = 'river'
            else:
                bds[bd_nr] = 'ocean'
            bd_nr += 1
        line = f.readline()

river_indices = [k for k, v in bds.items() if v == 'river']

# Center coords of SCHISM river boundaries
river_segments = [np.asarray(s.bdy_segments[ind]) - 1 for ind in river_indices]
river_coords = [(s.lon[seg].mean(), s.lat[seg].mean()) for seg in river_segments]
riverX, riverY = np.asarray(river_coords)[:, 0], np.asarray(river_coords)[:, 1]

os.chdir(cwd)

# === Match E-HYPE outlets to SCHISM rivers ===
if method == 'interactive':
    # Plot interactive map
    fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines(draw_labels=True)

    # Plot SCHISM rivers and E-HYPE outlets
    ax.plot(riverX, riverY, 'ko', markersize=5, transform=ccrs.PlateCarree(), label='SCHISM rivers')
    ax.plot(lons, lats, 'ro', markersize=5, transform=ccrs.PlateCarree(), label='E-HYPE outlets')

    if add_labels:
        for lon, lat in zip(lons, lats):
            ax.text(lon, lat, f'({lon:.2f}, {lat:.2f})', transform=ccrs.PlateCarree())
    
    # User clicks nearest E-HYPE points
    ax.set_title('RIGHT MB Click nearest red circle for each SCHISM river (in order) to select, click middle mouse to preview values for period')
    #plt.show()
    
    if os.path.exists(catchment_file):
        catchments.plot(ax=ax,linewidth=2,label='catchments')
    
    coords = []
    def onclick(event):
        if not event.inaxes:
            return  # ignore clicks outside axes

        # === Right click: confirm selection ===
        if event.button == 3:
            if len(coords) < len(riverX):
                coords.append((event.xdata, event.ydata))
                ax.plot(event.xdata, event.ydata, 'md')  # mark click
                fig.canvas.draw()
            if len(coords) == len(riverX):
                print("All rivers selected. Close the window.")

        # === Middle click: preview discharge ===
        elif event.button == 2:
            # find nearest E-HYPE outlet
            lon_click, lat_click = event.xdata, event.ydata
            idx = np.argmin((lons - lon_click)**2 + (lats - lat_click)**2)

            # extract discharge series for this outlet
            Q_preview = np.asarray(dsQ['cout'][selection, idx])

            # plot in a new figure
            fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True)
            ax1.plot(time_index[selection], Q_preview)
            ax1.set_title(f"Preview discharge for outlet #{idx}")
            ax1.set_xlabel("Time")
            ax1.set_ylabel("Q [m³/s]")

            ax2.bar([0], [Q_preview.mean()])
            ax2.set_xticks([0])
            ax2.set_xticklabels([f"Outlet {idx}"])
            ax2.set_ylabel("<Q> [m³/s]")

            plt.show(block=False)  # non-blocking, so main figure stays open
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()   # stays open until user closes the window
    plt.close(fig)
else: 
    coords=river_coords    
print("Picked coords:", coords)


    
# Nearest E-HYPE outlets to river boundaries respectivley manual selections
river_nns = np.asarray([
    np.argmin((lons - xq) ** 2 + (lats - yq) ** 2) for xq, yq in coords
])


# === Export flux.th for lateral boundaries ===
Q = np.asarray(dsQ['cout'][selection, river_nns])
M = np.column_stack((time_in_seconds, Q))
np.savetxt(f'flux.th_ehype_{start_date.date()}', M, fmt='%d' + ' %f' * len(riverX))


# plot it
fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(2, 2)
ax0 = plt.subplot(gs[0:1,0:1])
s.plotAtnodes(s.depths,ax=ax0,cmap=plt.cm.gray)
ax0.plot(riverX, riverY,'bo',label='RiverBD')
ax0.plot(lons[river_nns], lats[river_nns],'ro',label='E-hype Ref')
plt.legend()
for nr,coord in enumerate(list(zip(riverX, riverY))):
    plt.text(coord[0],coord[1],nr)
ax1 = plt.subplot(gs[1,:])
ax1.plot(time_index[selection] ,Q)
ax1.legend(['River 1','River 2','River 3','River 4'])
ax1.set_xlabel('T / days')
ax1.set_ylabel('Q / m³/s')
ax2 = plt.subplot(gs[0,1])
ax2.bar( np.arange(Q.shape[1]),Q.mean(axis=0))
ax2.set_xlabel('River Nr')
ax2.set_ylabel('<Q> / m³/s')
plt.tight_layout()
plt.savefig('Lat_BD_rivers_NN')
plt.show()



# === Select point sources (non-boundary outlets) ===
outlet_points = list(zip(lons, lats))
dist, pids = s.node_tree_latlon.query(outlet_points)

positions_in_distance = np.where(dist <= maximum_distance)[0]
positions_in_distance = [i for i in positions_in_distance if i not in river_nns]

table_subset = table.iloc[positions_in_distance]
source_coords = list(zip(table_subset.POURX.values, table_subset.POURY.values))

# Plot selected sources
fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.LAND, edgecolor='black')
ax.plot(lons, lats, 'ro', markersize=5, transform=ccrs.PlateCarree(), label='E-HYPE outlets')
ax.plot(riverX, riverY, 'ko', markersize=5, transform=ccrs.PlateCarree(), label='SCHISM rivers')
ax.plot(table_subset.POURX.values, table_subset.POURY.values, 'md', markersize=8,
        transform=ccrs.PlateCarree(), label='Selected point sources')
ax.set_extent([np.min(s.lon),np.max(s.lon),np.min(s.lat),np.max(s.lat)])        
ax.legend()
plt.show()
plt.savefig('river_selections.png', dpi=300)

# === Write SCHISM source/sink files ===
# Discharge
Q_sources = np.asarray(dsQ['cout'][selection, positions_in_distance])

# Temperature and salinity placeholders
DIM = (len(time_in_seconds), len(source_coords))
T = -9999 * np.ones(DIM)
S = np.zeros(DIM)



# avoid dry cells need precomputed map of dry elements
import xarray as xr
ds=xr.open_dataset(dry_element_file)
dry_elems=np.where(ds.dryFlagElement_max>0)[0]
dry_elems+=1 # make one based counting for depth dict overwrite

for key in  dry_elems:
    print(s.element_depth[key])
    s.element_depth[key]=-9999
    print(s.element_depth[key])
#np.intersect1d(elements,dry_elems)


ds.dryFlagElement_max.values[elements-1]


elements=create_source_sink(s, source_coords, name='source_sink.in',mindepth=minimum_depth,elements=None,numcells=30)
create_vsource(time_in_seconds, Q_sources, name='vsource.th')
create_msource(time_in_seconds, T, S, name='msource.th')

elements=np.asarray(elements)-1


# Plot selected sources
fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.LAND, edgecolor='black')
ax.plot(lons, lats, 'ro', markersize=5, transform=ccrs.PlateCarree(), label='E-HYPE outlets')
ax.plot(np.asarray(list(s.element_lon.values()))[elements],
np.asarray(list(s.element_lat.values()))[elements],
'r+', markersize=5, transform=ccrs.PlateCarree(), label='Applied cells')
ax.plot(riverX, riverY, 'ko', markersize=5, transform=ccrs.PlateCarree(), label='SCHISM rivers')
ax.plot(table_subset.POURX.values, table_subset.POURY.values, 'md', markersize=8,
        transform=ccrs.PlateCarree(), label='Selected point sources')
ax.set_extent([np.min(s.lon),np.max(s.lon),np.min(s.lat),np.max(s.lat)])        
ax.legend()
plt.show()
plt.savefig('river_selections.png', dpi=300)


s.plot_domain_boundaries()
s.plot_mesh(latlon=True)
plt.plot(np.asarray(list(s.element_lon.values()))[elements-1],
np.asarray(list(s.element_lat.values()))[elements-1],
'r+', markersize=5, label='Applied cells')

s.plotAtnodes(ds.dryFlagElement_max.values)
plt.plot(np.asarray(list(s.element_lon.values()))[np.asarray(elements)-1],
np.asarray(list(s.element_lat.values()))[np.asarray(elements)-1],
'r+', markersize=5, label='Applied cells')