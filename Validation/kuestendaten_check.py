"""
Validation of Tide Portal files for Gertman Bight.

Assumes the model grid is in UTM coordinates.
If using lon/lat, modify coordinate handling accordingly.
"""
# conda activate geo_env
import os
import sys
#import datetime
from glob import glob
from scipy.spatial import cKDTree
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# === Project imports ===
import dask
sys.path.insert(0, '/work/gg0028/SCHISM/schism-hereon-utilities/')
from schism import *

sys.path.insert(0, '/work/gg0028/SCHISM/schism-hereon-utilities/Validation/')
from kuestendaten import *
import datetime

sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/Lib/')
from TaylorDiagram import * 
# === Settings ===
TP_DIR = '/work/gg0028/g260114/PROJECTS/FOCCUS/downloads'          # Tide portal files
#NCDIR = '/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000/outputs01/'
NCDIR = '/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/NO_Ehype000/outputs_all2/'
NCDIRS=['/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/NO_Ehype000/outputs_all2/','/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000/outputs_all2/']
MAX_DIST = 200  # Maximum distance (m) between station and grid node
MIN_DEPTH = 4 #m
MAX_STACK =  170 #365 the minimum of availalbe days for the comparalbe runs
plt.ion()

# Extract configuration names from paths
def get_config_name(ncdir):
    """Extract configuration name from path."""
    # Extract the directory name after 'RUNS' and before 'outputs'
    parts = ncdir.split('/')
    for i, part in enumerate(parts):
        if 'RUNS' in part and i+1 < len(parts):
            # The next part should be the configuration name
            if i+2 < len(parts) and 'outputs' in parts[i+2]:
                return parts[i+1]
    # Fallback: use last meaningful directory before 'outputs'
    for i, part in enumerate(parts):
        if 'outputs' in part and i > 0:
            return parts[i-1]
    return 'config'

CONFIG_NAMES = [get_config_name(ncdir) for ncdir in NCDIRS]
print(f"Configuration names: {CONFIG_NAMES}")

# ----------------------------------------------------------------------
# Load stations
# ----------------------------------------------------------------------
os.chdir(TP_DIR)

station_files = glob.glob('*.txt')
station_names = [fname.split('!')[0] for fname in station_files]

stations = {}
for fname in station_files:
    key = fname.split('!')[0]
    print(f"Loading {key}")
    stations[key] = KuestenDataFile(fname)


# ----------------------------------------------------------------------
# Extract valid UTM locations
# ----------------------------------------------------------------------
def extract_locations(stations_dict):
    """Extract (x, y, label) for stations with valid UTM coordinates."""
    xs, ys, labels = [], [], []
    for key, obj in stations_dict.items():
        try:
            xs.append(float(obj.East))
            ys.append(float(obj.North))
            labels.append(key)
        except Exception:
            pass  # Skip stations missing coordinates
    return xs, ys, labels


# ----------------------------------------------------------------------
# Plot all station locations
# ----------------------------------------------------------------------
proj_utm32 = ccrs.UTM(zone=32)
proj_geo = ccrs.PlateCarree()

fig = plt.figure(figsize=(8, 10))
ax = plt.axes(projection=proj_utm32)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')

xs, ys, labels = extract_locations(stations)
ax.scatter(xs, ys, transform=proj_utm32, s=30)

for x, y, label in zip(xs, ys, labels):
    ax.text(x, y, label, fontsize=8, transform=proj_utm32)

ax.set_title("Küstenmessstationen — ETRS89 / UTM 32N")
plt.show()


# ----------------------------------------------------------------------
# Filter stations by model domain & proximity
# ----------------------------------------------------------------------
os.chdir('/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000')
s = schism_setup()

x,y=np.asarray(s.x),np.asarray(s.y)
d=np.asarray(s.depths)
MIN_DEPTH = 4 #m
x[d<MIN_DEPTH]=np.inf
tree=cKDTree(list(zip(x,y)))




s.plot_domain_boundaries(append=True, plot_legend=False, latlon=False)
plt.scatter(xs, ys, s=30)
for x, y, label in zip(xs, ys, labels):
    plt.text(x, y, label, fontsize=8)

# Build KD-tree for node search
s.init_node_tree(latlon=False)

station_coords = list(zip(xs, ys))


#dists, nns = s.node_tree_xy.query(station_coords)
# remove drypoins
dists, nns = tree.query(station_coords)


# Remove distant stations
valid_mask = dists <= MAX_DIST
#stations = {
    #k: st for k, keep in zip(stations.keys(), valid_mask) if keep
#}

stations_filtered = {
    k: stations[k]
    for k, remove in zip(stations.keys(), ~valid_mask)
    if not remove
}
stations=stations_filtered

# Update coordinates
xs, ys, labels = extract_locations(stations)
station_coords = list(zip(xs, ys))
dists, nns = s.node_tree_xy.query(station_coords)

# Plot filtered stations
s.plot_domain_boundaries(append=True, plot_legend=False, latlon=False)
plt.scatter(xs, ys, s=30)
for x, y, label in zip(xs, ys, labels):
    plt.text(x, y, label, fontsize=8)
plt.show()


# ----------------------------------------------------------------------
# Initialize storage dictionaries for multiple configurations
# ----------------------------------------------------------------------
config_results = {}
config_stats = {}
config_model_times = {}
config_model_arrays = {}

# ----------------------------------------------------------------------
# Define function to compute error statistics
# ----------------------------------------------------------------------
def cal_error_stats(model_array, station_array, station_names):
    """Compute error statistics between model and observations."""
    residuals = model_array - station_array
    bias = np.nanmean(residuals, axis=0)
    abs_errors = np.abs(residuals)
    rmse_per_station = np.sqrt(np.nanmean(residuals ** 2, axis=0))
    mae_per_station = np.nanmean(abs_errors, axis=0)
    std_data = np.nanstd(station_array, axis=0)
    std_model = np.nanstd(model_array, axis=0)

    # correlation
    R = np.zeros(len(station_names))
    for i in range(len(station_names)):
        obs, model = station_array[:, i], model_array[:, i]
        ivalid = ~np.isnan(obs)
        if np.sum(ivalid) > 1:  # Need at least 2 points for correlation
            R[i] = np.corrcoef(obs[ivalid], model[ivalid])[0, 1]
        else:
            R[i] = np.nan
    # Store metrics
    stats = {
        'bias': bias,
        'rmse': rmse_per_station,
        'cor': R,
        'mae': mae_per_station,
        'std1': std_data,
        'std2': std_model,
        'std_rel': std_model / std_data,
        'residuals': residuals
    }
    return stats

# ----------------------------------------------------------------------
# Loop over model configurations
# ----------------------------------------------------------------------
for iconfig, NCDIR in enumerate(NCDIRS):
    config_name = CONFIG_NAMES[iconfig]
    print(f"\n{'='*60}")
    print(f"Processing configuration: {config_name}")
    print(f"Directory: {NCDIR}")
    print(f"{'='*60}\n")
    
    # ----------------------------------------------------------------------
    # Load SCHISM output
    # ----------------------------------------------------------------------
    acc = schism_outputs_by_variable(ncdir=NCDIR, varlist=['out2d', 'salinity'],use_dask=True,max_stack=MAX_STACK)
    ds = acc.ds['salinity']

    # Build model time vector
    inidate = datetime.datetime(2017, 1, 2, 0, 0, 0)
    findate = inidate + datetime.timedelta(days=float(ds.time[-1].values) / 86400)
    model_times = inidate + np.array([datetime.timedelta(seconds=t) for t in ds.time.values])

    # Store model times for this configuration
    config_model_times[config_name] = model_times

    # ----------------------------------------------------------------------
    # Filter and interpolate station data to model times (only first time)
    # ----------------------------------------------------------------------
    if iconfig == 0:
        plt.figure()
        for key, st in stations.items():
            df = st.data

            # Filter to model time range
            mask = (df['timestamp'] >= inidate) & (df['timestamp'] <= findate)
            df = df.loc[mask].copy()
            df.set_index('timestamp', inplace=True)

            # Interpolate to model times
            df_interp = df.reindex(model_times)
            df_interp['value'] = df_interp['value'].interpolate(method='time')
            df_interp['status'] = df_interp['status'].ffill()

            st.data = df_interp

        plt.show()

        for key, st in stations.items():
            df = st.data
            is_invalid=st.data.value.values<-99
            st.data.value.values[is_invalid]=np.nan

            print(key)
            print(    st.data.value.min())

    # ----------------------------------------------------------------------
    # Compare model vs observations
    # ----------------------------------------------------------------------
    #slow
    subset = ds.sel(nSCHISM_hgrid_node=nns, nSCHISM_vgrid_layers=-1)

    # Save subset for this configuration
    subset.to_netcdf(f'modal_ts_at_stations_{config_name}.nc')

    #model_array = subset.salinity.values
    from dask.diagnostics import ProgressBar
    with ProgressBar():
        model_array = subset.salinity.compute()

    # Store model array for this configuration
    config_model_arrays[config_name] = model_array

    # Station array (same for all configurations)
    if iconfig == 0:
        station_array = np.column_stack([
            stations[k].data['value'].values for k in stations.keys()
        ])
        station_names = list(stations.keys())

if False:
    station_names2 = list(stations.keys())
    # 2. Stack data arrays for each station into a single 2D array (time x station)
    station_array = np.column_stack([stations[k].data['value'].values for k in station_names2])
    # 3. Use the index (assuming all stations have the same datetime index)
    time_index = stations[station_names2[0]].data.index
    # 4. Create an xarray Dataset
    ds = xr.Dataset(
        {
            "value": (["time", "station"], station_array)
        },
        coords={
            "time": time_index,
            "station": station_names
        }
    )
    # 5. Export to NetCDF
    ds.to_netcdf("stations_data.nc")


#########################################] | 100% Completed | 268.34 s    on interactive64

    # Compute and store stats for this configuration
    stats = cal_error_stats(model_array, station_array, station_names)
    config_stats[config_name] = stats
    
    print(f"Completed processing {config_name}")
    print(f"Overall RMSE: {np.sqrt(np.nanmean(stats['residuals'] ** 2)):.3f}")
    print(f"Overall MAE: {np.nanmean(np.abs(stats['residuals'])):.3f}\n")

# ----------------------------------------------------------------------
# Create combined Taylor plot for all configurations
# ----------------------------------------------------------------------
print("Creating combined Taylor plot...")

# Prepare samples for Taylor plot - first configuration
first_config = CONFIG_NAMES[0]
first_stats = config_stats[first_config]
samples = [[first_stats['std2'][i]/first_stats['std1'][i], first_stats['cor'][i], 
            f"{station_names[i]} ({first_config})"] 
           for i in range(len(station_names))]

# Create initial Taylor diagram
dia = plotTaylor(samples=samples)

# Add samples for remaining configurations with different colors
colors = plt.cm.tab10(range(len(NCDIRS)))
for iconfig, config_name in enumerate(CONFIG_NAMES[1:], start=1):
    stats = config_stats[config_name]
    for i in range(len(station_names)):
        stddev = stats['std2'][i] / stats['std1'][i]
        corrcoef = stats['cor'][i]
        name = f"{station_names[i]} ({config_name})"
        # Use different marker style and color for each configuration
        dia.add_sample(stddev, corrcoef,
                      marker='$%d$' % (i+1), ms=10, ls='',
                      mfc=colors[iconfig], mec=colors[iconfig],
                      label=name)

# Update legend to show all configurations
plt.legend(dia.samplePoints,
           [p.get_label() for p in dia.samplePoints],
           loc='upper right', bbox_to_anchor=(1.3, 1.0))
plt.title("Taylor Diagram - Multiple Model Configurations")
plt.tight_layout()
plt.savefig('taylor_plot_combined.png', dpi=300, bbox_inches='tight')
plt.show()

#add data
#for isetup,key in enumerate(setup_names[1:]):
    #compare1=compare[setup_names[isetup+1]]
    #samples.append([ [compare1[tag]['T']['std2'][i]/compare1[tag]['T']['std1'][i],compare1[tag]['T']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'] for tag in names if ~np.isnan(compare1[tag]['T']['bias']).max()])	


#colors=plt.cm.tab10(range(len(setup_names)))  
colors=plt.cm.tab10(range(2))  
nr=1
#for isetup,key in enumerate(setup_names[1:]):
#for isetup,key in enumerate(range(1)):
#    #Add models to Taylor diagram
#    for i,(stddev, corrcoef, name) in enumerate(samples):
#        i,stddev, corrcoef, name
#        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)	



# ----------------------------------------------------------------------
# Print overall statistics for all configurations
# ----------------------------------------------------------------------
print("\n" + "="*60)
print("Overall Statistics Summary")
print("="*60)
for config_name in CONFIG_NAMES:
    stats = config_stats[config_name]
    overall_rmse = np.sqrt(np.nanmean(stats['residuals'] ** 2))
    overall_mae = np.nanmean(np.abs(stats['residuals']))
    print(f"\n{config_name}:")
    print(f"  Overall RMSE: {overall_rmse:.3f}")
    print(f"  Overall MAE: {overall_mae:.3f}")

# ----------------------------------------------------------------------
# Plot RMSE spatial distribution for each configuration
# ----------------------------------------------------------------------
add_val=True
xs,ys=np.asarray(xs),np.asarray(ys)
plt.close('all')

for config_name in CONFIG_NAMES:
    stats = config_stats[config_name]
    for key in stats.keys():
        if key == 'residuals':
            continue
        plt.figure()
        vals=stats[key]
        vmin,vmax=np.nanmin(vals),np.nanmax(vals)
        if vmin < 0:
            vmax=np.max(np.abs(vals))
            vmin=-vmax
            cmap='RdBu_r'
        else:
            cmap='viridis'    
        if key=='cor':
            vmin,vmax=-1,1
        s.plot_domain_boundaries(append=True, plot_legend=False, latlon=False)
        ph = plt.scatter(xs, ys, c=vals, s=30, vmin=vmin, vmax=vmax, cmap=cmap)
        cb = plt.colorbar(ph)
        cb.set_label(key +"[psu]")
        plt.title(f"Station {key} vs Model - {config_name}")
        plt.xlabel("UTM Easting (m)")
        plt.ylabel("UTM Northing (m)")
        if add_val:
            for i in range(len(station_names)):
                xi,yi=xs[i],ys[i]
                val=' {:.2f}'.format(vals[i])
                plt.text(xi,yi,val)
        plt.tight_layout()
        plt.axis([xs.min()-1000,xs.max()+1000,ys.min()-1000,ys.max()+1000])
        plt.show()
        plt.savefig(f'error_stats_{key}_{config_name}.png',dpi=300)

plt.close('all')
####################
# --- Example time series ---subset
#time = pd.date_range("2021-01-01", periods=200)
#salinity = np.sin(np.linspace(0, 10, 200)) + 30


if False:
    # --- Example station coordinates (UTM, meters) ---
    station_x, station_y = xs[istat],ys[istat]
    grid_x = np.linspace(340000, 350000, 50)
    grid_y = np.linspace(6020000, 6030000, 50)

    station_x, station_y = 345000, 6023000
    grid_x = np.linspace(xs[istat]-5000,xs[istat]+5000, 50)
    grid_y = np.linspace(ys[istat]-5000,ys[istat]+5000, 50)





if False:
    for istat in range(len(stations)):
        time,salinity=model_times,model_array[:,istat]
        salinity2=station_array[:,istat]
        # --- Main time series plot ---
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(time, salinity, color='tab:blue',label='sim')
        ax.plot(time, salinity2, 'r.',label='obs')
        ax.set_title("Salinity Time Series at Station")
        ax.set_xlabel("Time")
        ax.set_ylabel("Salinity [psu]")
        ax.legend(loc='lower right')
        plt.title(station_names[istat])
        plt.tight_layout()

        inset_ax = inset_axes(
            ax, width="30%", height="30%", loc='upper left',
            axes_class=plt.Axes
        )
        plt.axes(inset_ax)
        s.plot_domain_boundaries(append=True,latlon=False)
        plt.plot(xs[istat],ys[istat],'ro')
        plt.xlim(xs[istat]-5000,xs[istat]+5000)
        plt.ylim(ys[istat]-5000,ys[istat]+5000)
        plt.xlabel('')
        plt.ylabel('')
        plt.xticks([])
        plt.yticks([])

        plt.savefig('timeseries_'+station_names[istat]+'.png',dpi=300)


# ----------------------------------------------------------------------
# Plot time series for all configurations
# ----------------------------------------------------------------------
for istat in range(len(station_names)):
    salinity2 = station_array[:, istat]
    ivalid = ~np.isnan(salinity2)
    
    # --- Create figure with two subplots (wide left, narrow right) ---
    fig, (ax_ts, ax_box) = plt.subplots(
        1, 2,
        figsize=(12, 5),
        gridspec_kw={'width_ratios': [3, 1]}  # left: wide, right: narrow
    )

    # --- Main time series plot with all configurations ---
    colors_ts = plt.cm.tab10(range(len(CONFIG_NAMES)))
    for iconfig, config_name in enumerate(CONFIG_NAMES):
        time = config_model_times[config_name]
        salinity = config_model_arrays[config_name][:, istat]
        ax_ts.plot(time, salinity, color=colors_ts[iconfig], 
                   label=f'sim ({config_name})', linewidth=1.5, alpha=0.8)
    
    # Plot observations (only once)
    time_first = config_model_times[CONFIG_NAMES[0]]
    ax_ts.plot(time_first, salinity2, 'r.', label='obs', alpha=0.7, markersize=2)
    
    ax_ts.set_title(f"Salinity Time Series at {station_names[istat]}")
    ax_ts.set_xlabel("Time")
    ax_ts.set_ylabel("Salinity [psu]")
    ax_ts.legend(loc='lower right', fontsize=8)

    # --- Inset map inside the time series subplot ---
    inset_ax = inset_axes(ax_ts, width="30%", height="30%", loc='upper left')
    s.plot_domain_boundaries(append=True, latlon=False)
    inset_ax.plot(xs[istat], ys[istat], 'ro')
    inset_ax.set_xlim(xs[istat]-5000, xs[istat]+5000)
    inset_ax.set_ylim(ys[istat]-5000, ys[istat]+5000)
    inset_ax.set_xticks([])
    inset_ax.set_yticks([])

    # --- Right subplot: boxplot of all time series ---
    box_data = [salinity2[ivalid]]
    box_labels = ['Obs']
    for iconfig, config_name in enumerate(CONFIG_NAMES):
        salinity = config_model_arrays[config_name][:, istat]
        box_data.append(salinity[ivalid])
        box_labels.append(f'Model\n({config_name})')
    
    bp = ax_box.boxplot(
        box_data,
        labels=box_labels,
        patch_artist=True,
        medianprops=dict(color='black')
    )
    # Color the boxes
    bp['boxes'][0].set_facecolor('lightcoral')  # Obs
    for i in range(1, len(bp['boxes'])):
        bp['boxes'][i].set_facecolor(colors_ts[i-1])
        bp['boxes'][i].set_alpha(0.7)
    
    ax_box.set_title("Distribution")
    ax_box.set_ylabel("Salinity [psu]")
    ax_box.tick_params(axis='x', rotation=45)

    # --- Layout & save ---
    plt.tight_layout()
    plt.savefig(f"timeseries_{station_names[istat]}.png", dpi=300, bbox_inches='tight')
    #plt.close(fig)

####
