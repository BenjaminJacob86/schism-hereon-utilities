"""
Prepare E-HYPE river discharge and sediment data for SCHISM forcing.

Outputs:
--------
1. flux.th              : Lateral river boundary discharge
2. source_sink.in       : Source/sink definition
3. vsource.th           : Volume source discharge
4. msource.th           : Temperature/salinity sources
5. SED_*.th             : Sediment class forcing
6. msource_sed.th       : Combined T/S + sediment source file

Workflow:
---------
1. Read SCHISM setup and identify river boundaries from bctides.in
2. Match SCHISM river boundaries with E-HYPE coastal outlets
   - interactive (manual click)
   - or nearest-neighbour
3. Export lateral boundary discharge (flux.th)
4. Select remaining E-HYPE outlets as point sources
5. Write SCHISM source/sink and forcing files
6. Distribute E-HYPE sediment load into SCHISM sediment classes
"""

# =============================================================================
# Environment & imports
# =============================================================================
# conda activate geo_env

import os
import sys
import datetime as dt

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset, num2date

# SCHISM utilities
sys.path.insert(0, '/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import schism_setup
from source_sink_river_functions import (
    create_source_sink,
    create_vsource,
    create_msource
)

# =============================================================================
# User configuration
# =============================================================================
# --- Input data
outlets_file = '/work/gg0028/g260114/PROJECTS/FOCCUS/Data/E-hype/EH432_coast_outlets.csv'
discharge_file = (
    '/work/gg0028/g260114/PROJECTS/FOCCUS/Data/E-hype/'
    'COUT_eh432_2000-2024_coast.txt'
)
sediment_file = (
    '/work/gg0028/g260114/PROJECTS/FOCCUS/Data/E-hype/'
    'CCTS_eh432_2000-2024_coast.txt'
)

catchment_file = (
    '/work/gg0028/g260114/PROJECTS/FOCCUS/Data/River/'
    '2_EH4_basins_northsea_poly.gpkg'
)

rundir = '/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/NO_Ehype000/'

# --- Time window
start_date = dt.datetime(2017, 1, 1)
end_date   = dt.datetime(2018, 1, 1)

# --- Matching & selection
method = 'interactive'       # 'interactive' | 'nn'
maximum_distance = 0.025     # deg (for point source inclusion)
minimum_depth = 3            # m (NOT YET USED)

# --- Sediment setup
nsed = 8                     # number of sediment classes
add_sed = True

# --- Plotting
add_labels = False

# =============================================================================
# Load E-HYPE outlet locations
# =============================================================================
# Here tweaking might be needed depending on the input file which was not consistent among iterations of ehype on zenodo
try:
    table = pd.read_table(outlets_file, sep=';')
except Exception:
    table = pd.read_table(outlets_file)
    

lons = table.POURX.values
lats = table.POURY.values

# =============================================================================
# Load catchments (optional)
# =============================================================================
add_catchment = False
if os.path.exists(catchment_file):
    try:
        catchments = gpd.read_file(catchment_file)
        add_catchment = True
    except Exception:
        pass

# =============================================================================
# Load E-HYPE discharge time series
# =============================================================================
file_type = discharge_file.split('.')[-1].lower()

if file_type == 'txt':
    # Text-based discharge file
    dfQ = pd.read_csv(discharge_file, sep='\t')
    time_index = pd.to_datetime(dfQ.DATE).dt.tz_localize(None)

    selection = (time_index >= start_date) & (time_index <= end_date)
    time_in_seconds = (
        (time_index[selection] - time_index[selection].iloc[0])
        / np.timedelta64(1, 's')
    ).astype(int).values

    Q_all = dfQ.values[:, 1:]  # shape: time x outlet
    is_nc = False

elif file_type == 'nc':
    dsQ = Dataset(discharge_file)
    time_var = dsQ.variables['time']

    cftime_dates = num2date(
        time_var[:],
        units=time_var.units,
        calendar=time_var.calendar
    )

    time_index = pd.to_datetime([str(d) for d in cftime_dates])
    selection = (time_index >= start_date) & (time_index <= end_date)

    time_in_seconds = (
        (time_index[selection] - time_index[selection][0])
        / np.timedelta64(1, 's')
    ).astype(int).values

    is_nc = True

else:
    raise ValueError("Unsupported discharge file format")

# =============================================================================
# Load SCHISM setup
# =============================================================================
cwd = os.getcwd()
os.chdir(rundir)

s = schism_setup()
s.lon = np.asarray(s.lon)
s.lat = np.asarray(s.lat)
s.init_node_tree(latlon=True)

# =============================================================================
# Identify river boundaries from bctides.in
# =============================================================================
river_indices = []

with open('bctides.in') as f:
    # number of open boundaries
    for line in f:
        if 'nope' in line:
            nop = int(line.split()[0])
            break

    bd_nr = 0
    for line in f:
        parts = line.split()
        if len(parts) > 4 and parts[2].isdigit():
            if int(parts[2]) == 1:   # river boundary flag
                river_indices.append(bd_nr)
            bd_nr += 1

# Compute boundary center coordinates
river_segments = [np.asarray(s.bdy_segments[i]) - 1 for i in river_indices]
river_coords = [
    (s.lon[seg].mean(), s.lat[seg].mean())
    for seg in river_segments
]

riverX = np.array([c[0] for c in river_coords])
riverY = np.array([c[1] for c in river_coords])

os.chdir(cwd)

# =============================================================================
# Match SCHISM rivers to E-HYPE outlets
# =============================================================================
if method == 'interactive':

    fig, ax = plt.subplots(
        figsize=(10, 5),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )

    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.add_feature(cfeature.OCEAN)
    ax.gridlines(draw_labels=True)

    ax.plot(riverX, riverY, 'ko', label='SCHISM rivers')
    ax.plot(lons, lats, 'ro', label='E-HYPE outlets')

    if add_catchment:
        catchments.plot(ax=ax, linewidth=2)

    coords = []

    def onclick(event):
        if not event.inaxes:
            return

        # Right click → select outlet
        if event.button == 3:
            coords.append((event.xdata, event.ydata))
            ax.plot(event.xdata, event.ydata, 'md')
            fig.canvas.draw()

        # Middle click → preview discharge
        elif event.button == 2:
            idx = np.argmin((lons - event.xdata)**2 + (lats - event.ydata)**2)

            if is_nc:
                Q_preview = dsQ['cout'][selection, idx]
            else:
                Q_preview = Q_all[selection, idx]

            plt.figure()
            plt.plot(time_index[selection], Q_preview)
            plt.title(f'Outlet {idx}')
            plt.ylabel('Q [m³/s]')
            plt.show(block=False)

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.title(
        'Right-click: select outlet | Middle-click: preview discharge'
    )
    plt.show()

    river_nns = np.array([
        np.argmin((lons - x)**2 + (lats - y)**2)
        for x, y in coords
    ])

else:
    # Nearest-neighbour matching
    river_nns = np.array([
        np.argmin((lons - x)**2 + (lats - y)**2)
        for x, y in river_coords
    ])

# =============================================================================
# Export flux.th (lateral boundaries)
# =============================================================================
if is_nc:
    Q_river = -dsQ['cout'][selection, river_nns]
else:
    Q_river = -Q_all[selection][:, river_nns]

M = np.column_stack((time_in_seconds, Q_river))
np.savetxt(
    f'flux.th_ehype_{start_date.date()}',
    M,
    fmt='%d' + ' %f' * Q_river.shape[1]
)

# =============================================================================
# Select point sources (non-boundary outlets)
# =============================================================================
outlet_points = list(zip(lons, lats))
dist, node_ids = s.node_tree_latlon.query(outlet_points)

positions_in_distance = np.where(dist <= maximum_distance)[0]
positions_in_distance = [
    i for i in positions_in_distance if i not in river_nns
]

source_coords = list(
    zip(
        table.POURX.values[positions_in_distance],
        table.POURY.values[positions_in_distance]
    )
)

# =============================================================================
# Write source/sink files
# =============================================================================
if is_nc:
    Q_sources = dsQ['cout'][selection, positions_in_distance]
else:
    Q_sources = -Q_all[selection][:, positions_in_distance]

DIM = (len(time_in_seconds), len(source_coords))
T = -9999 * np.ones(DIM)   # placeholder
S = np.zeros(DIM)

create_source_sink(s, source_coords, name='source_sink.in',mindepth=minimum_depth)
create_vsource(time_in_seconds, Q_sources, name='vsource.th')
create_msource(time_in_seconds, T, S, name='msource.th')

# =============================================================================
# Sediment forcing (optional)
# =============================================================================
if add_sed:

    m_sed = pd.read_csv(sediment_file, sep='\t')
    m_sed_in_time = m_sed.values[selection, 1:]

    mgl_to_gl = 0.001  # E-HYPE mg/L → SCHISM g/L

    # Nearest SCHISM nodes
    bd_nodes = s.node_tree_latlon.query(coords)[1]
    src_nodes = s.node_tree_latlon.query(source_coords)[1]

    frac_bd = np.zeros((len(bd_nodes), nsed))
    frac_src = np.zeros((len(src_nodes), nsed))

    os.chdir(rundir)
    for i in range(nsed):
        s.read_gr3(f'bed_frac_{i+1}.ic')
        frac_bd[:, i]  = s.gr3[f'bed_frac_{i+1}'][bd_nodes]
        frac_src[:, i] = s.gr3[f'bed_frac_{i+1}'][src_nodes]
    os.chdir(cwd)

    SEDs = []
    for iclass in range(nsed):
        sed_src = (
            frac_src[:, iclass]
            * m_sed_in_time[:, positions_in_distance]
            * mgl_to_gl
        )
        SEDs.append(sed_src)


        # Write BD sediements
        BD_SED = (
            frac_bd[:, iclass]
            * m_sed_in_time[:, river_nns]
            * mgl_to_gl
        )
        M = np.column_stack((time_in_seconds, BD_SED))
        np.savetxt(
            f'SED_{iclass+1}.th_ehype_{start_date.date()}',
            M,
            fmt='%d' + ' %f' * BD_SED.shape[1]
        )
        
        
        

    create_msource(
        time_in_seconds,
        T,
        S,
        SEDs,
        name='msource_sed.th'
    )
