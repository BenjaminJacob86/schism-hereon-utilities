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
import datetime as dt

sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/Lib/')
from TaylorDiagram import * 
# === Settings ===
TP_DIR = '/work/gg0028/g260114/PROJECTS/FOCCUS/downloads'          # Tide portal files
#NCDIR = '/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000/outputs01/'
NCDIR = '/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/NO_Ehype000/outputs_all2/'
NCDIRS=['/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/NO_Ehype000/outputs_all2/','/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000/outputs_all2/']
RUNDIRS=['/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/NO_Ehype000/','/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000/']
MAX_DIST = 200  # Maximum distance (m) between station and grid node
MIN_DEPTH = 4 #m
MAX_STACK =  170 #365 the minimum of availalbe days for the comparalbe runs

# Loading method: 'standard' (load all then select) or 'optimized' (select then load)
# 'optimized' is more memory-efficient and often faster for few stations
LOADING_METHOD = 'optimized'  # Options: 'standard' or 'optimized'

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
# Filter out shallow points before building tree (newer Python versions don't work with np.inf)
# Keep only deep nodes for tree construction
deep_mask = d >= MIN_DEPTH
x_deep = x[deep_mask]
y_deep = y[deep_mask]
deep_node_indices = np.where(deep_mask)[0]  # Original node indices for deep nodes
tree = cKDTree(list(zip(x_deep, y_deep)))


# Build KD-tree for node search
s.init_node_tree(latlon=False)

station_coords = list(zip(xs, ys))


#dists, nns = s.node_tree_xy.query(station_coords)
# Query tree built from deep nodes only
# tree returns indices into the filtered array, need to map back to original indices
tree_indices = tree.query(station_coords)[1]  # Indices in filtered (deep) array
nns = deep_node_indices[tree_indices]  # Map back to original node indices
dists = tree.query(station_coords)[0]  # Distances


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
# Use filtered tree to avoid shallow nodes
# tree returns indices into the filtered array, need to map back to original indices
tree_indices = tree.query(station_coords)[1]  # Indices in filtered (deep) array
nns = deep_node_indices[tree_indices]  # Map back to original node indices
dists = tree.query(station_coords)[0]  # Distances

# Verify that selected nodes have sufficient depth
# If any node is too shallow, it means the tree found a shallow point (shouldn't happen with filtered tree)
# but we check anyway for safety
node_depths = d[nns]
shallow_mask = node_depths < MIN_DEPTH
if np.any(shallow_mask):
    print(f"Warning: {np.sum(shallow_mask)} stations matched to nodes with depth < {MIN_DEPTH}m")
    print("This should not happen with filtered tree. Checking depths...")
    # For shallow matches, try to find deeper nearby nodes
    for i in np.where(shallow_mask)[0]:
        # Find nearest deeper node within reasonable distance
        nearby_mask = (d >= MIN_DEPTH) & (np.sqrt((x - xs[i])**2 + (y - ys[i])**2) <= MAX_DIST * 2)
        if np.any(nearby_mask):
            nearby_indices = np.where(nearby_mask)[0]
            nearby_dists = np.sqrt((x[nearby_indices] - xs[i])**2 + (y[nearby_indices] - ys[i])**2)
            nns[i] = nearby_indices[np.argmin(nearby_dists)]
            dists[i] = nearby_dists[np.argmin(nearby_dists)]
        else:
            print(f"  Station {labels[i]} at ({xs[i]:.0f}, {ys[i]:.0f}) has no nearby deep nodes")

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
config_station_array = None  # Will be set during first iteration
config_station_names = None   # Will be set during first iteration

# ----------------------------------------------------------------------
# Define optimized function to load model data at stations
# ----------------------------------------------------------------------
def load_model_at_stations_optimized(ncdir, nns, varname='salinity', 
                                     layer=-1, max_stack=365, use_dask=True):
    """
    Optimized loading: select nodes first, then load files one by one.
    This is more memory-efficient than loading full dataset then selecting.
    
    Parameters:
    -----------
    ncdir : str
        Directory containing NetCDF output files
    nns : array-like
        Node indices to extract (nearest neighbor indices)
    varname : str
        Variable name to extract (default: 'salinity')
    layer : int
        Vertical layer index (default: -1 for surface)
    max_stack : int
        Maximum number of stack files to load
    use_dask : bool
        Whether to use dask for lazy loading
    
    Returns:
    --------
    model_array : numpy array
        Array of shape (ntime, nstations) with model values at stations
    """
    import glob
    from glob import glob
    
    # Find files
    nr_nc0 = np.sort(glob(f'{ncdir}/out2d*.nc'))[0].split('/')[-1].split('_')[-1]
    varfiles = glob(f'{ncdir}*_{nr_nc0}')
    
    # Find salinity files
    salt_files = [f for f in varfiles if 'salinity' in f or 'salt' in f]
    if not salt_files:
        # Try to find files matching varname
        salt_files = [f for f in varfiles if varname in f]
    
    if not salt_files:
        raise ValueError(f"Could not find {varname} files in {ncdir}")
    
    # Get file pattern
    file_pattern = salt_files[0]
    base_name = file_pattern[file_pattern.rindex('/')+1:file_pattern.rindex('_')]
    
    # Find all stack files
    all_files = []
    for iorder in range(1, 6):  # Support up to 5-digit stack numbers
        pattern = f'{ncdir}{base_name}_{"?"*iorder}.nc'
        files = sorted(glob(pattern))
        all_files.extend(files)
    
    all_files = sorted(set(all_files))  # Remove duplicates and sort
    if max_stack > 0:
        all_files = all_files[:max_stack]
    
    print(f"Loading {len(all_files)} files, selecting {len(nns)} nodes...")
    
    # Load and select nodes from each file, then concatenate
    data_chunks = []
    from dask.diagnostics import ProgressBar
    
    if use_dask:
        # Lazy loading: open each file, select nodes, chunk, then concatenate
        datasets = []
        for f in all_files:
            ds = xr.open_dataset(f, engine='netcdf4')
            # Select nodes immediately (lazy operation)
            ds_subset = ds.sel(nSCHISM_hgrid_node=nns, nSCHISM_vgrid_layers=layer)
            # Chunk for dask
            ds_subset = ds_subset.chunk({'time': -1, 'nSCHISM_hgrid_node': len(nns)})
            datasets.append(ds_subset)
        
        # Concatenate along time dimension (still lazy)
        ds_combined = xr.concat(datasets, dim='time')
        
        # Extract variable and compute
        with ProgressBar():
            model_array = ds_combined[varname].compute().values
        
        # Close datasets
        for ds in datasets:
            ds.close()
    else:
        # Eager loading: load each file, select nodes, store in list
        for f in all_files:
            ds = xr.open_dataset(f)
            # Select nodes immediately
            subset = ds.sel(nSCHISM_hgrid_node=nns, nSCHISM_vgrid_layers=layer)
            data_chunks.append(subset[varname].values)
            ds.close()
        
        # Concatenate along time axis
        model_array = np.concatenate(data_chunks, axis=0)
    
    return model_array


# ----------------------------------------------------------------------
# Define optimized function to load model data at stations
# ----------------------------------------------------------------------
def load_model_at_stations_optimized(ncdir, nns, varname='salinity', 
                                     layer=-1, max_stack=365, use_dask=True):
    """
    Optimized loading: select nodes first, then load files one by one.
    This is more memory-efficient than loading full dataset then selecting.
    
    Parameters:
    -----------
    ncdir : str
        Directory containing NetCDF output files
    nns : array-like
        Node indices to extract (nearest neighbor indices)
    varname : str
        Variable name to extract (default: 'salinity')
    layer : int
        Vertical layer index (default: -1 for surface)
    max_stack : int
        Maximum number of stack files to load
    use_dask : bool
        Whether to use dask for lazy loading
    
    Returns:
    --------
    model_array : numpy array
        Array of shape (ntime, nstations) with model values at stations
    """
    import glob
    from glob import glob
    # xarray should already be imported via schism module, but import here for safety
    try:
        import xarray as xr
    except ImportError:
        # If not available, try to get from schism module context
        import sys
        if 'xarray' in sys.modules:
            xr = sys.modules['xarray']
        else:
            raise ImportError("xarray is required for optimized loading")
    
    # Find files
    nr_nc0 = np.sort(glob(f'{ncdir}/out2d*.nc'))[0].split('/')[-1].split('_')[-1]
    varfiles = glob(f'{ncdir}*_{nr_nc0}')
    
    # Find salinity files
    salt_files = [f for f in varfiles if 'salinity' in f or 'salt' in f]
    if not salt_files:
        # Try to find files matching varname
        salt_files = [f for f in varfiles if varname in f]
    
    if not salt_files:
        raise ValueError(f"Could not find {varname} files in {ncdir}")
    
    # Get file pattern
    file_pattern = salt_files[0]
    base_name = file_pattern[file_pattern.rindex('/')+1:file_pattern.rindex('_')]
    
    # Find all stack files
    all_files = []
    for iorder in range(1, 6):  # Support up to 5-digit stack numbers
        pattern = f'{ncdir}{base_name}_{"?"*iorder}.nc'
        files = sorted(glob(pattern))
        all_files.extend(files)
    
    all_files = sorted(set(all_files))  # Remove duplicates and sort
    if max_stack > 0:
        all_files = all_files[:max_stack]
    
    print(f"Loading {len(all_files)} files, selecting {len(nns)} nodes...")
    
    # Load and select nodes from each file, then concatenate
    from dask.diagnostics import ProgressBar
    
    if use_dask:
        # Lazy loading: open each file, select nodes, chunk, then concatenate
        datasets = []
        for f in all_files:
            ds = xr.open_dataset(f, engine='netcdf4')
            # Select nodes immediately (lazy operation)
            ds_subset = ds.sel(nSCHISM_hgrid_node=nns, nSCHISM_vgrid_layers=layer)
            # Chunk for dask
            ds_subset = ds_subset.chunk({'time': -1, 'nSCHISM_hgrid_node': len(nns)})
            datasets.append(ds_subset)
        
        # Concatenate along time dimension (still lazy)
        ds_combined = xr.concat(datasets, dim='time')
        
        # Extract variable and compute
        with ProgressBar():
            model_array = ds_combined[varname].compute().values
        
        # Close datasets
        for ds in datasets:
            ds.close()
    else:
        # Eager loading: load each file, select nodes, store in list
        data_chunks = []
        for f in all_files:
            ds = xr.open_dataset(f)
            # Select nodes immediately
            subset = ds.sel(nSCHISM_hgrid_node=nns, nSCHISM_vgrid_layers=layer)
            data_chunks.append(subset[varname].values)
            ds.close()
        
        # Concatenate along time axis
        model_array = np.concatenate(data_chunks, axis=0)
    
    return model_array


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

    # correlation - only compute for valid observation time steps
    R = np.zeros(len(station_names))
    for i in range(len(station_names)):
        obs, model = station_array[:, i], model_array[:, i]
        # Find valid time steps where both obs and model have valid data
        # Primary filter: observations must be valid
        ivalid_obs = ~np.isnan(obs)
        # Secondary: also exclude model NaN at those valid obs times
        ivalid = ivalid_obs & ~np.isnan(model)
        
        if np.sum(ivalid) > 1:  # Need at least 2 points for correlation
            R[i] = np.corrcoef(obs[ivalid], model[ivalid])[0, 1]
        else:
            R[i] = np.nan
    # Compute overall metrics (across all stations and times)
    overall_rmse = np.sqrt(np.nanmean(residuals ** 2))
    overall_mae = np.nanmean(np.abs(residuals))
    overall_bias = np.nanmean(bias)
    overall_cor = np.nanmean(R)
    
    # Store metrics
    stats = {
        'bias': bias,
        'rmse': rmse_per_station,
        'cor': R,
        'mae': mae_per_station,
        'std1': std_data,
        'std2': std_model,
        'std_rel': std_model / std_data,
        'residuals': residuals,
        # Overall metrics (precomputed for efficiency)
        'overall_rmse': overall_rmse,
        'overall_mae': overall_mae,
        'overall_bias': overall_bias,
        'overall_cor': overall_cor,
        'mean_rmse': np.nanmean(rmse_per_station),  # Mean of per-station RMSEs
        'mean_mae': np.nanmean(mae_per_station),    # Mean of per-station MAEs
        'mean_cor': overall_cor,  # Same as overall_cor
        'mean_bias': np.nanmean(np.abs(bias))  # Mean absolute bias
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

    # Build model time vector # depends on SCHISM version
    RUNDIR=RUNDIRS[iconfig]          
    p=param(RUNDIR+'/param.nml')
    if type(ds.time[0].values)==np.ndarray:
        inidate=dt.datetime(int(p.get_parameter('start_year')),
                    int(p.get_parameter('start_month')),
                    int(p.get_parameter('start_day')),
                    int(p.get_parameter('start_hour')),0,0)
        #inidate = datetime.datetime(2017, 1, 2, 0, 0, 0)
        findate = inidate + dt.timedelta(days=float(ds.time[-1].values) / 86400)
        model_times = inidate + np.array([dt.timedelta(seconds=t) for t in ds.time.values])
    else:
        #inidate = datetime.datetime(2017, 1, 2, 0, 0, 0)
        #findate = inidate + datetime.timedelta(days=float(ds.time[-1].values) / 86400)
        model_times = ds.time.values
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
    if LOADING_METHOD == 'optimized':
        # Optimized: load files one by one, selecting nodes immediately
        # This is more memory-efficient and often faster for few stations
        # Note: We still need ds for time vector, but don't load full data
        print("Using optimized loading method (file-by-file with node selection)...")
        model_array = load_model_at_stations_optimized(
            NCDIR, nns, varname='salinity', layer=-1, 
            max_stack=MAX_STACK, use_dask=True
        )
    else:
        # Standard: load full dataset, then select nodes
        print("Using standard loading method (load all, then select)...")
        subset = ds.sel(nSCHISM_hgrid_node=nns, nSCHISM_vgrid_layers=-1)

        # Optionally save subset for later use (can be disabled for speed)
        SAVE_SUBSET_NC = False  # Set to True if you need the NetCDF file later
        if SAVE_SUBSET_NC:
            subset.to_netcdf(f'modal_ts_at_stations_{config_name}.nc')

        # Compute only the selected subset (much faster than loading everything)
        from dask.diagnostics import ProgressBar
        with ProgressBar():
            model_array = subset.salinity.compute()

    # Store model array for this configuration
    config_model_arrays[config_name] = model_array
    
    print(f"Completed loading data for {config_name}")

    # Station array (same for all configurations) - create once
    if iconfig == 0:
        station_array = np.column_stack([
            stations[k].data['value'].values for k in stations.keys()
        ])
        station_names = list(stations.keys())
        # Store for use in stats calculation loop
        config_station_array = station_array
        config_station_names = station_names

# ----------------------------------------------------------------------
# Separate loop for statistics calculation
# This allows recalculating stats without reloading expensive data
# ----------------------------------------------------------------------
print("\n" + "="*60)
print("Calculating Statistics")
print("="*60)

for iconfig, config_name in enumerate(CONFIG_NAMES):
    print(f"\nProcessing statistics for: {config_name}")
    
    # Get stored data
    model_array = config_model_arrays[config_name]
    
    # Compute and store stats for this configuration
    stats = cal_error_stats(model_array, config_station_array, config_station_names)
    config_stats[config_name] = stats
    
    print(f"Completed processing {config_name}")
    print(f"Overall RMSE: {stats['overall_rmse']:.3f}")
    print(f"Overall MAE: {stats['overall_mae']:.3f}")
    print(f"Mean Correlation: {stats['overall_cor']:.3f}\n")

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
# Stats are now computed inside the main loop above, so this duplicate block is removed

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
# Define function to compute skill improvement metrics
# ----------------------------------------------------------------------
def compute_skill_improvement(stats_ref, stats_new, config_ref_name, config_new_name):
    """
    Compute skill improvement/deterioration metrics comparing two configurations.
    
    Parameters:
    -----------
    stats_ref : dict
        Statistics dictionary for reference configuration
    stats_new : dict
        Statistics dictionary for new configuration to compare
    config_ref_name : str
        Name of reference configuration
    config_new_name : str
        Name of new configuration
    
    Returns:
    --------
    skill_metrics : dict
        Dictionary containing various skill improvement metrics
    """
    # Use precomputed overall metrics from config_stats (no recalculation needed!)
    rmse_ref = stats_ref['overall_rmse']
    rmse_new = stats_new['overall_rmse']
    mae_ref = stats_ref['overall_mae']
    mae_new = stats_new['overall_mae']
    cor_ref_mean = stats_ref['overall_cor']
    cor_new_mean = stats_new['overall_cor']
    bias_ref_mean = stats_ref['mean_bias']  # Mean absolute bias
    bias_new_mean = stats_new['mean_bias']
    
    # 1. Relative RMSE improvement (%)
    rmse_improvement_pct = (rmse_ref - rmse_new) / rmse_ref * 100
    
    # 2. Skill Score based on RMSE (1 - RMSE_new/RMSE_ref)
    # Range: -∞ to 1, where 1 = perfect, 0 = same skill, negative = worse
    skill_score_rmse = 1 - (rmse_new / rmse_ref)
    
    # 3. Relative MAE improvement (%)
    mae_improvement_pct = (mae_ref - mae_new) / mae_ref * 100
    
    # 4. Correlation improvement (absolute change)
    cor_improvement = cor_new_mean - cor_ref_mean
    
    # 5. Correlation improvement (%)
    cor_improvement_pct = (cor_new_mean - cor_ref_mean) / (1 - cor_ref_mean) * 100 if cor_ref_mean < 1 else 0
    
    # 6. Bias improvement (%)
    bias_improvement_pct = (bias_ref_mean - bias_new_mean) / bias_ref_mean * 100 if bias_ref_mean > 0 else 0
    
    # 7. Combined Skill Index (CSI)
    # Combines normalized improvements in RMSE, correlation, and bias
    # Normalize each component to [0, 1] scale
    rmse_norm = max(0, min(1, skill_score_rmse + 0.5))  # Shift and clip to [0,1]
    cor_norm = (cor_new_mean + 1) / 2  # Map [-1, 1] to [0, 1]
    cor_ref_norm = (cor_ref_mean + 1) / 2
    cor_contrib = (cor_norm - cor_ref_norm) * 0.4  # Weight correlation improvement
    rmse_contrib = rmse_norm * 0.4  # Weight RMSE improvement
    bias_contrib = max(0, (bias_ref_mean - bias_new_mean) / (bias_ref_mean + 0.01)) * 0.2  # Weight bias improvement
    csi = rmse_contrib + cor_contrib + bias_contrib
    
    # 8. Taylor Skill Distance
    # Distance from perfect point (std_ratio=1, cor=1) on Taylor diagram
    std_ratio_ref = np.nanmean(stats_ref['std2'] / stats_ref['std1'])
    std_ratio_new = np.nanmean(stats_new['std2'] / stats_new['std1'])
    # Distance in Taylor space: sqrt((std_ratio-1)^2 + (1-corr)^2)
    taylor_dist_ref = np.sqrt((std_ratio_ref - 1)**2 + (1 - cor_ref_mean)**2)
    taylor_dist_new = np.sqrt((std_ratio_new - 1)**2 + (1 - cor_new_mean)**2)
    taylor_improvement = taylor_dist_ref - taylor_dist_new  # Positive = improvement
    
    skill_metrics = {
        'rmse_improvement_pct': rmse_improvement_pct,
        'skill_score_rmse': skill_score_rmse,
        'mae_improvement_pct': mae_improvement_pct,
        'cor_improvement': cor_improvement,
        'cor_improvement_pct': cor_improvement_pct,
        'bias_improvement_pct': bias_improvement_pct,
        'combined_skill_index': csi,
        'taylor_skill_improvement': taylor_improvement,
        'rmse_ref': rmse_ref,
        'rmse_new': rmse_new,
        'mae_ref': mae_ref,
        'mae_new': mae_new,
        'cor_ref': cor_ref_mean,
        'cor_new': cor_new_mean,
        'bias_ref': bias_ref_mean,
        'bias_new': bias_new_mean,
    }
    
    return skill_metrics

# ----------------------------------------------------------------------
# Print overall statistics for all configurations
# ----------------------------------------------------------------------
print("\n" + "="*60)
print("Overall Statistics Summary")
print("="*60)
for config_name in CONFIG_NAMES:
    stats = config_stats[config_name]
    # Use precomputed overall metrics
    print(f"\n{config_name}:")
    print(f"  Overall RMSE: {stats['overall_rmse']:.3f}")
    print(f"  Overall MAE: {stats['overall_mae']:.3f}")
    print(f"  Mean Correlation: {stats['overall_cor']:.3f}")

# ----------------------------------------------------------------------
# Compute skill improvement between configurations
# ----------------------------------------------------------------------
if len(CONFIG_NAMES) >= 2:
    print("\n" + "="*60)
    print("Skill Improvement Analysis")
    print("="*60)
    print(f"Comparing: {CONFIG_NAMES[0]} (reference) vs {CONFIG_NAMES[1]} (new)")
    
    skill_metrics = compute_skill_improvement(
        config_stats[CONFIG_NAMES[0]], 
        config_stats[CONFIG_NAMES[1]],
        CONFIG_NAMES[0],
        CONFIG_NAMES[1]
    )
    
    print(f"\n{'Metric':<35} {'Value':>15} {'Interpretation':<30}")
    print("-" * 80)
    
    # RMSE metrics
    print(f"{'RMSE Improvement (%)':<35} {skill_metrics['rmse_improvement_pct']:>15.2f} ", end="")
    if skill_metrics['rmse_improvement_pct'] > 0:
        print(f"{'✓ Improvement':<30}")
    else:
        print(f"{'✗ Deterioration':<30}")
    
    print(f"{'Skill Score (RMSE)':<35} {skill_metrics['skill_score_rmse']:>15.3f} ", end="")
    if skill_metrics['skill_score_rmse'] > 0:
        print(f"{'✓ Better than reference':<30}")
    elif skill_metrics['skill_score_rmse'] == 0:
        print(f"{'= Same skill':<30}")
    else:
        print(f"{'✗ Worse than reference':<30}")
    
    # MAE metrics
    print(f"{'MAE Improvement (%)':<35} {skill_metrics['mae_improvement_pct']:>15.2f} ", end="")
    if skill_metrics['mae_improvement_pct'] > 0:
        print(f"{'✓ Improvement':<30}")
    else:
        print(f"{'✗ Deterioration':<30}")
    
    # Correlation metrics
    print(f"{'Correlation Change':<35} {skill_metrics['cor_improvement']:>15.3f} ", end="")
    if skill_metrics['cor_improvement'] > 0:
        print(f"{'✓ Improvement':<30}")
    else:
        print(f"{'✗ Deterioration':<30}")
    
    print(f"{'Correlation Improvement (%)':<35} {skill_metrics['cor_improvement_pct']:>15.2f} ", end="")
    if skill_metrics['cor_improvement_pct'] > 0:
        print(f"{'✓ Improvement':<30}")
    else:
        print(f"{'✗ Deterioration':<30}")
    
    # Bias metrics
    print(f"{'Bias Improvement (%)':<35} {skill_metrics['bias_improvement_pct']:>15.2f} ", end="")
    if skill_metrics['bias_improvement_pct'] > 0:
        print(f"{'✓ Improvement':<30}")
    else:
        print(f"{'✗ Deterioration':<30}")
    
    # Combined metrics
    print(f"{'Combined Skill Index':<35} {skill_metrics['combined_skill_index']:>15.3f} ", end="")
    if skill_metrics['combined_skill_index'] > 0.1:
        print(f"{'✓ Significant improvement':<30}")
    elif skill_metrics['combined_skill_index'] > 0:
        print(f"{'○ Slight improvement':<30}")
    elif skill_metrics['combined_skill_index'] == 0:
        print(f"{'= No change':<30}")
    else:
        print(f"{'✗ Deterioration':<30}")
    
    print(f"{'Taylor Skill Improvement':<35} {skill_metrics['taylor_skill_improvement']:>15.3f} ", end="")
    if skill_metrics['taylor_skill_improvement'] > 0:
        print(f"{'✓ Closer to perfect':<30}")
    else:
        print(f"{'✗ Further from perfect':<30}")
    
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    improvements = []
    deteriorations = []
    
    if skill_metrics['rmse_improvement_pct'] > 0:
        improvements.append(f"RMSE improved by {skill_metrics['rmse_improvement_pct']:.1f}%")
    else:
        deteriorations.append(f"RMSE worsened by {abs(skill_metrics['rmse_improvement_pct']):.1f}%")
    
    if skill_metrics['cor_improvement'] > 0:
        improvements.append(f"Correlation improved by {skill_metrics['cor_improvement']:.3f}")
    else:
        deteriorations.append(f"Correlation decreased by {abs(skill_metrics['cor_improvement']):.3f}")
    
    if skill_metrics['bias_improvement_pct'] > 0:
        improvements.append(f"Bias reduced by {skill_metrics['bias_improvement_pct']:.1f}%")
    else:
        deteriorations.append(f"Bias increased by {abs(skill_metrics['bias_improvement_pct']):.1f}%")
    
    if improvements:
        print("Improvements:")
        for imp in improvements:
            print(f"  ✓ {imp}")
    
    if deteriorations:
        print("Deteriorations:")
        for det in deteriorations:
            print(f"  ✗ {det}")
    
    if skill_metrics['combined_skill_index'] > 0:
        print(f"\nOverall: {CONFIG_NAMES[1]} shows improvement over {CONFIG_NAMES[0]}")
        print(f"Combined Skill Index: {skill_metrics['combined_skill_index']:.3f}")
    elif skill_metrics['combined_skill_index'] < 0:
        print(f"\nOverall: {CONFIG_NAMES[1]} shows deterioration compared to {CONFIG_NAMES[0]}")
        print(f"Combined Skill Index: {skill_metrics['combined_skill_index']:.3f}")
    else:
        print(f"\nOverall: {CONFIG_NAMES[1]} shows similar skill to {CONFIG_NAMES[0]}")
    
    print("="*60)
    
    # ----------------------------------------------------------------------
    # Create skill comparison visualization
    # ----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Skill Comparison: {CONFIG_NAMES[0]} vs {CONFIG_NAMES[1]}', fontsize=14, fontweight='bold')
    
    # 1. RMSE comparison
    ax = axes[0, 0]
    metrics = ['RMSE', 'MAE', '|Bias|']
    ref_vals = [skill_metrics['rmse_ref'], skill_metrics['mae_ref'], skill_metrics['bias_ref']]
    new_vals = [skill_metrics['rmse_new'], skill_metrics['mae_new'], skill_metrics['bias_new']]
    x = np.arange(len(metrics))
    width = 0.35
    ax.bar(x - width/2, ref_vals, width, label=CONFIG_NAMES[0], alpha=0.8, color='tab:blue')
    ax.bar(x + width/2, new_vals, width, label=CONFIG_NAMES[1], alpha=0.8, color='tab:orange')
    ax.set_ylabel('Value [psu]')
    ax.set_title('Error Metrics Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Correlation comparison
    ax = axes[0, 1]
    ax.bar([0, 1], [skill_metrics['cor_ref'], skill_metrics['cor_new']], 
           color=['tab:blue', 'tab:orange'], alpha=0.8, width=0.6)
    ax.set_ylabel('Correlation Coefficient')
    ax.set_title('Correlation Comparison')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(CONFIG_NAMES)
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)
    # Add improvement arrow
    if skill_metrics['cor_improvement'] > 0:
        ax.annotate('', xy=(1, skill_metrics['cor_new']), xytext=(0, skill_metrics['cor_ref']),
                   arrowprops=dict(arrowstyle='->', color='green', lw=2))
        ax.text(0.5, (skill_metrics['cor_ref'] + skill_metrics['cor_new'])/2 + 0.05,
               f"+{skill_metrics['cor_improvement']:.3f}", ha='center', color='green', fontweight='bold')
    
    # 3. Improvement percentages
    ax = axes[1, 0]
    imp_metrics = ['RMSE', 'MAE', 'Bias', 'Correlation']
    imp_vals = [skill_metrics['rmse_improvement_pct'], 
                skill_metrics['mae_improvement_pct'],
                skill_metrics['bias_improvement_pct'],
                skill_metrics['cor_improvement_pct']]
    colors_bar = ['green' if v > 0 else 'red' for v in imp_vals]
    bars = ax.barh(imp_metrics, imp_vals, color=colors_bar, alpha=0.7)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
    ax.set_xlabel('Improvement (%)')
    ax.set_title('Percentage Improvement/Deterioration')
    ax.grid(True, alpha=0.3, axis='x')
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, imp_vals)):
        ax.text(val + (1 if val > 0 else -1), i, f'{val:.1f}%', 
               va='center', ha='left' if val > 0 else 'right', fontweight='bold')
    
    # 4. Combined Skill Index and Taylor Skill
    ax = axes[1, 1]
    skill_metrics_plot = ['Combined\nSkill Index', 'Taylor Skill\nImprovement']
    skill_vals = [skill_metrics['combined_skill_index'], 
                  skill_metrics['taylor_skill_improvement']]
    colors_skill = ['green' if v > 0 else 'red' for v in skill_vals]
    bars = ax.bar(skill_metrics_plot, skill_vals, color=colors_skill, alpha=0.7)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.set_ylabel('Skill Metric Value')
    ax.set_title('Combined Skill Metrics')
    ax.grid(True, alpha=0.3, axis='y')
    # Add value labels
    for bar, val in zip(bars, skill_vals):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + (0.01 if height > 0 else -0.01),
               f'{val:.3f}', ha='center', va='bottom' if height > 0 else 'top', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('skill_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

# ----------------------------------------------------------------------
# Plot RMSE spatial distribution for each configuration
# ----------------------------------------------------------------------
add_val=True
xs,ys=np.asarray(xs),np.asarray(ys)
plt.close('all')

for config_name in CONFIG_NAMES:
    stats = config_stats[config_name]
    # Only plot per-station metrics (exclude overall/mean metrics which are single values)
    for key in stats.keys():
        # Skip overall/mean metrics (single values, not per-station arrays)
        if 'overall' in key or 'mean' in key:
            continue
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







# ----------------------------------------------------------------------
# Plot time series for all configurations (individual station plots)
# ----------------------------------------------------------------------
for istat in range(len(station_names)):
    salinity2 = config_station_array[:, istat]
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

# ----------------------------------------------------------------------
# New: Multi-configuration comparison plots
# These show all simulations together for easier comparison
# ----------------------------------------------------------------------

# 1. Grid plot: All stations in subplots, all configs overlaid
print("\nCreating multi-station comparison plots...")
n_stations = len(station_names)
n_cols = min(3, n_stations)  # 3 columns max
n_rows = int(np.ceil(n_stations / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
if n_stations == 1:
    axes = [axes]
else:
    axes = axes.flatten()

colors_ts = plt.cm.tab10(range(len(CONFIG_NAMES)))

for istat in range(n_stations):
    ax = axes[istat]
    salinity2 = config_station_array[:, istat]
    ivalid = ~np.isnan(salinity2)
    
    # Plot all configurations
    for iconfig, config_name in enumerate(CONFIG_NAMES):
        time = config_model_times[config_name]
        salinity = config_model_arrays[config_name][:, istat]
        ax.plot(time, salinity, color=colors_ts[iconfig], 
               label=f'{config_name}', linewidth=1.5, alpha=0.8)
    
    # Plot observations
    time_first = config_model_times[CONFIG_NAMES[0]]
    ax.plot(time_first, salinity2, 'r.', label='Obs', alpha=0.6, markersize=1.5)
    
    ax.set_title(f"{station_names[istat]}", fontsize=10, fontweight='bold')
    ax.set_xlabel("Time", fontsize=9)
    ax.set_ylabel("Salinity [psu]", fontsize=9)
    ax.grid(True, alpha=0.3)
    if istat == 0:  # Legend only on first subplot
        ax.legend(loc='best', fontsize=7, ncol=1)

# Hide unused subplots
for i in range(n_stations, len(axes)):
    axes[i].axis('off')

plt.suptitle('Multi-Configuration Comparison: All Stations', fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('timeseries_all_stations_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# 2. Single figure with all stations and all configs (zoomed view option)
# Create a figure showing all stations in separate panels
fig, axes = plt.subplots(n_stations, 1, figsize=(14, 3*n_stations), sharex=True)
if n_stations == 1:
    axes = [axes]

for istat in range(n_stations):
    ax = axes[istat]
    salinity2 = config_station_array[:, istat]
    ivalid = ~np.isnan(salinity2)
    
    # Plot all configurations
    for iconfig, config_name in enumerate(CONFIG_NAMES):
        time = config_model_times[config_name]
        salinity = config_model_arrays[config_name][:, istat]
        ax.plot(time, salinity, color=colors_ts[iconfig], 
               label=f'{config_name}', linewidth=1.2, alpha=0.8)
    
    # Plot observations
    time_first = config_model_times[CONFIG_NAMES[0]]
    ax.plot(time_first, salinity2, 'r.', label='Obs', alpha=0.5, markersize=1)
    
    ax.set_ylabel(f"{station_names[istat]}\n[psu]", fontsize=9, rotation=0, ha='right', va='center')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=7, ncol=len(CONFIG_NAMES)+1)

axes[-1].set_xlabel("Time", fontsize=10)
plt.suptitle('Time Series Comparison: All Stations and Configurations', 
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('timeseries_all_stations_vertical.png', dpi=300, bbox_inches='tight')
plt.show()

# 3. Summary statistics plot: Mean and std across all stations for each config
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Calculate mean and std across stations for each time step
for iconfig, config_name in enumerate(CONFIG_NAMES):
    time = config_model_times[config_name]
    model_data = config_model_arrays[config_name]  # shape: (ntime, nstations)
    
    # Mean across stations (spatial mean)
    mean_ts = np.nanmean(model_data, axis=1)
    std_ts = np.nanstd(model_data, axis=1)
    
    ax1.plot(time, mean_ts, color=colors_ts[iconfig], 
            label=f'{config_name} (mean)', linewidth=2, alpha=0.8)
    ax1.fill_between(time, mean_ts - std_ts, mean_ts + std_ts, 
                     color=colors_ts[iconfig], alpha=0.2, label=f'{config_name} (±1 std)')
    
    # Also plot observations mean
    if iconfig == 0:
        obs_mean = np.nanmean(config_station_array, axis=1)
        obs_std = np.nanstd(config_station_array, axis=1)
        ax1.plot(time, obs_mean, 'r-', label='Obs (mean)', linewidth=2, alpha=0.8)
        ax1.fill_between(time, obs_mean - obs_std, obs_mean + obs_std, 
                        color='red', alpha=0.2, label='Obs (±1 std)')

ax1.set_ylabel('Mean Salinity [psu]', fontsize=11)
ax1.set_title('Spatial Mean and Standard Deviation Across All Stations', fontsize=12, fontweight='bold')
ax1.legend(loc='best', fontsize=9, ncol=2)
ax1.grid(True, alpha=0.3)

# Plot difference from observations
for iconfig, config_name in enumerate(CONFIG_NAMES):
    time = config_model_times[config_name]
    model_data = config_model_arrays[config_name]
    obs_data = config_station_array
    
    # Mean difference across stations
    diff = model_data - obs_data
    mean_diff = np.nanmean(diff, axis=1)
    std_diff = np.nanstd(diff, axis=1)
    
    ax2.plot(time, mean_diff, color=colors_ts[iconfig], 
            label=f'{config_name} - Obs', linewidth=2, alpha=0.8)
    ax2.fill_between(time, mean_diff - std_diff, mean_diff + std_diff, 
                     color=colors_ts[iconfig], alpha=0.2)

ax2.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax2.set_ylabel('Model - Obs [psu]', fontsize=11)
ax2.set_xlabel('Time', fontsize=11)
ax2.set_title('Bias: Spatial Mean Difference from Observations', fontsize=12, fontweight='bold')
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('timeseries_summary_statistics.png', dpi=300, bbox_inches='tight')
plt.show()

print("Multi-configuration comparison plots saved.")

####
