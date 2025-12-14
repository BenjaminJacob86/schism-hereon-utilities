"""
Structured grid interpolation script for SCHISM unstructured grid output.

This script interpolates SCHISM model output from unstructured triangular grids to regular
structured grids with variable renaming and CF-compliant metadata.

The script supports two interpolation modes controlled by the interp_surface_only setting:

For 3D variables:
    - If interp_surface_only=True: Only surface layer is extracted (nSCHISM_vgrid_layers=-1)
    - If interp_surface_only=False: Full 3D vertical interpolation is performed using linear 
      interpolation based on SCHISM's hybrid sigma/z-layer coordinates, followed by horizontal 
      interpolation using inverse distance weighting.
      
For 2D variables (e.g., bottom stress):
    - Only horizontal interpolation is applied (always 2D).

Time Format Handling:
    SCHISM outputs can have time in two formats depending on the SCHISM version:
    1. Properly formatted datetime (can be decoded normally by xarray)
    2. Time since model start (e.g., "seconds since 10890000.0" - invalid CF format)
    
    The script automatically detects the time format and handles both cases:
    - For properly formatted time: Uses datetime values directly
    - For time since model start: Reads reference time from param.nml (model start date/time)
      and converts to CF-compliant "seconds since YYYY-MM-DD HH:MM:SS" format
    
    The reference time is read from setup_dir/param.nml using the schism.param module.
    If param.nml is not available, the script will attempt to extract reference time from
    the time units attribute, but this may fail for numeric references.

Variable definitions (names, units, standard names, long names, valid ranges) and global
attributes are loaded from metadata.yaml file.

"""

import sys
import os
import yaml

#sys.path.insert(0, '/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine/')
#sys.path.insert(0, '//work/gg0028/SCHISM/schism-hereon-utilities/')
from schism import *
from matplotlib import pyplot as plt
from glob import glob
import datetime as dt


def get_layer_weights(s, dep, ti):
    """Calculate weights for vertical interpolation."""
    print('calculating weights for vertical interpolation')
    ibelow = np.zeros(s.nnodes, int)
    iabove = np.zeros(s.nnodes, int)
    weights = np.zeros((2, s.nnodes))

    # garbage values below ibtm different type , nan or strange values wrong values
    zcor = s.ds['zCoordinates'][ti, :, :].values
    zcor = np.ma.masked_array(zcor, mask=s.mask3d)

    a = np.sum(zcor <= dep, 1) - 1
    ibelow = a + s.ibbtm - 1
    iabove = np.minimum(ibelow + 1, s.nz - 1)
    inodes = np.where(a > 0)[0]
    ibelow2 = ibelow[inodes]
    iabove2 = iabove[inodes]

    d2 = zcor[inodes, iabove2] - dep
    d1 = dep - zcor[inodes, ibelow2]
    ivalid = d1 > 0.0
    iset = d1 == 0.0
    d1 = d1[ivalid]
    d2 = d2[ivalid]
    weights[0, inodes[ivalid]] = 1 - d1 / (d1 + d2)
    weights[1, inodes[ivalid]] = 1 - weights[0, inodes[ivalid]]
    weights[0, inodes[iset]] = 1
    weights[:, np.sum(weights, 0) == 0.0] = np.nan

    return ibelow, iabove, weights


########## User Settings ############################################

# Set to True to extract only surface layer of 3D variables, False for full 3D interpolation
interp_surface_only = True

########## Settings ############################################

# Load metadata (variable definitions and global attributes) from YAML file
yaml_file = 'metadata.yaml'

# Validate YAML file exists
if not os.path.exists(yaml_file):
    raise FileNotFoundError(f"YAML configuration file not found: {yaml_file}")

with open(yaml_file, 'r') as f:
    var_defs = yaml.safe_load(f)

# Validate required YAML keys exist
required_keys = ['names', 'units', 'std_names_2d', 'std_names_3d', 
                 'long_names_2d', 'long_names_3d', 'valid_ranges', 'global_attributes']
missing_keys = [key for key in required_keys if key not in var_defs]
if missing_keys:
    raise KeyError(f"Missing required keys in YAML file: {missing_keys}")

names = var_defs['names']
units = var_defs['units']
std_names_2d = var_defs['std_names_2d']
std_names_3d = var_defs['std_names_3d']
long_names_2d = var_defs['long_names_2d']
long_names_3d = var_defs['long_names_3d']
valid_ranges = var_defs['valid_ranges']

# Convert None strings to actual None for std_names
for key, value in std_names_2d.items():
    if value == 'null' or value is None:
        std_names_2d[key] = None
for key, value in std_names_3d.items():
    if value == 'null' or value is None:
        std_names_3d[key] = None

# Directory settings - adjust these for your setup
indir = '/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all2/'
outdir = '/work/gg0028/g260114/EDITO/REF_sim_2017/'

# Load setup
setup_dir = '/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
os.chdir(setup_dir)
s = schism_setup()

# Read reference time from param.nml for time format handling
# SCHISM outputs can have time in two formats:
# 1. Properly formatted datetime (can be decoded normally)
# 2. Time since model start (needs reference time from param.nml)
try:
    from schism import param
    p = param(setup_dir + '/param.nml')
    # Create reference time as a scalar datetime64
    reftime_dt = dt.datetime(int(p.get_parameter('start_year')),
                             int(p.get_parameter('start_month')),
                             int(p.get_parameter('start_day')),
                             int(p.get_parameter('start_hour')), 0, 0)
    reftime = np.datetime64(reftime_dt)
    print(f"Reference time from param.nml: {reftime} ({reftime_dt})")
except Exception as e:
    print(f"Warning: Could not read reference time from param.nml: {e}")
    print("  Will attempt to auto-detect time format from dataset")
    reftime = None

os.chdir(outdir)

# Grid and interpolation info
overwrite = False  # overwrite interpolation files if existent
lon, lat = np.asarray(s.lon), np.asarray(s.lat)
dx = 0.005  # interpolation grid spacing in degree
lonmin = 5.11608851
lonmax = 10.40483834
latmin = 53.03857668
latmax = 55.62872055
x = np.arange(lonmin, lonmax, dx)
y = np.arange(latmin, latmax, dx)
reload_weights = True  # reload weights
FillValue = -9999
dThresh = 0.5  # allowed distance for point pre selection

# Time steps per file (will be determined from dataset)
nt = None  # Will be set from dataset

# Vertical levels (only used if interp_surface_only=False)
nz = 21  # number interpolated vertical layers
z_levels = np.linspace(0, 20, nz)  # Example: 20 layers from surface to -30 m

# Global attributes - loaded from YAML, with dynamic fields added
global_atts = var_defs.get('global_attributes', {}).copy()
# Add dynamic fields that change with each run
global_atts['creation_date'] = f'Created on {dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
global_atts['history'] = f"; SCHISM output processed via interpolation_hybrid.py on {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"

# Variable files - make sure out2d is always first in this list
if interp_surface_only:
    varfiles = dict.fromkeys(['out2d', 'temperature', 'salinity', 'horizontalVelX', 'horizontalVelY'])
else:
    varfiles = dict.fromkeys(['out2d', 'temperature', 'salinity', 'horizontalVelX', 'horizontalVelY', 'zCoordinates'])

# Find variable files
for key in varfiles:
    varfiles[key] = list(np.sort(glob(indir + key + '_*.nc')))

# Validate that out2d files exist (required)
if len(varfiles['out2d']) == 0:
    raise FileNotFoundError(f"No out2d files found in {indir}. "
                          f"out2d files are required for interpolation.")

# Check for missing variable files and warn (but continue)
empty_file_lists = [key for key, files in varfiles.items() if len(files) == 0 and key != 'out2d']
if empty_file_lists:
    print(f"Warning: No files found for optional variables: {empty_file_lists}")
    print(f"  These variables will be skipped. Continuing with available variables...")
    # Remove empty variable types from varfiles
    for key in empty_file_lists:
        del varfiles[key]

# Check for inconsistent file counts and warn
file_counts = {key: len(files) for key, files in varfiles.items()}
if len(set(file_counts.values())) > 1:
    print(f"Warning: Inconsistent number of files across variables: {file_counts}")
    print(f"  Using minimum count ({min(file_counts.values())}) to avoid index errors.")
    # Use minimum count to avoid index errors
    min_count = min(file_counts.values())
    for key in varfiles:
        if len(varfiles[key]) > min_count:
            print(f"  Truncating {key} from {len(varfiles[key])} to {min_count} files")
            varfiles[key] = varfiles[key][:min_count]

# Variable names to process
varnames = ['elevation', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY',
            'sigWaveHeight', 'bottomStressX', 'bottomStressY']

# Validate all variables in varnames exist in YAML definitions
missing_vars = [varname for varname in varnames 
                if varname not in names or varname not in units or varname not in valid_ranges]
if missing_vars:
    raise KeyError(f"Variables in varnames list missing from YAML definitions: {missing_vars}")

###################################################################

# Grid and interpolation info
longitude = x
latitude = y
X, Y = np.meshgrid(x, y)
ny, nx = X.shape

# Load or calculate interpolation weights
weight_file_suffix = str(dx).replace('.', '')
parents_file = f'parents_structured_interpolation_{weight_file_suffix}.txt'
nodeweights_file = f'nodeweights_structured_interpolation_{weight_file_suffix}.txt'

if reload_weights:
    try:
        parents = np.loadtxt(parents_file).astype(int)
        ndeweights = np.loadtxt(nodeweights_file)
    except FileNotFoundError:
        print("Weight files not found, calculating weights...")
        reload_weights = False

if not reload_weights:
    # Hybrid finding approach for speed
    xq = X.flatten()
    yq = Y.flatten()
    s.init_element_tree(latlon=True)
    parents, ndeweights = s.find_parent_tri_fast(xq, yq, dThresh=dThresh, latlon=True, k=10)
    
    # Find missing points
    dist, nn = s.element_tree_latlon.query(list(zip(xq, yq)), k=1)
    missing = (parents == -1) & (dist < dThresh)
    
    if np.any(missing):
        print(f"Refining {np.sum(missing)} points with slow fallback...")
        parents_slow, weights_slow = s.find_parent_tri(
            xq[missing], yq[missing], dThresh=dThresh, latlon=True
        )
        parents[missing] = parents_slow
        ndeweights[missing] = weights_slow
    
    np.savetxt(parents_file, parents)
    np.savetxt(nodeweights_file, ndeweights)

############# Begin actual interpolation #########################

print(f"\n{'='*70}")
print(f"Starting interpolation processing")
print(f"Mode: {'Surface-only' if interp_surface_only else 'Full 3D'}")
print(f"Total files to process: {len(varfiles['out2d'])}")
print(f"{'='*70}\n")

outfiles = []
total_files = len(varfiles['out2d'])
for idate in range(total_files):
    file_progress = f"[{idate + 1}/{total_files}]"
    print(f"\n{file_progress} Processing file set {idate + 1} of {total_files}")
    print(f"{'─'*70}")
    
    # Build infiles list, checking for missing files
    infiles = []
    available_keys = []
    for key in varfiles.keys():
        if idate < len(varfiles[key]):
            file_path = varfiles[key][idate]
            if os.path.exists(file_path):
                infiles.append(file_path)
                available_keys.append(key)
            else:
                print(f"{file_progress} Warning: File not found for {key}: {file_path}, skipping this variable type")
        else:
            print(f"{file_progress} Warning: No file available for {key} at index {idate}, skipping this variable type")
    
    # Validate that out2d file exists (required)
    if len(infiles) == 0 or available_keys[0] != 'out2d':
        print(f"{file_progress} Error: out2d file is missing for date {idate + 1}. Skipping this date.")
        continue

    print(f"{file_progress} Opening and merging datasets...")
    print(f"{file_progress} Available variable types: {', '.join(available_keys)}")
    
    # Open datasets and detect time format - SCHISM outputs can have time in two formats:
    # 1. Properly formatted datetime (can be decoded normally)
    # 2. Time since model start (numeric values, needs reference time from param.nml)
    # Note: xarray may open files with decode_times=True even with malformed units,
    # leaving time as numeric values, so we need to check the actual time values
    try:
        dsout2d = xr.open_dataset(infiles[0], decode_times=True)
        # Check if time was actually decoded to datetime or is still numeric
        time_values = dsout2d.time.values
        time_units = dsout2d.time.attrs.get('units', '')
        
        # Check if time values are datetime-like or numeric
        is_datetime = False
        if len(time_values) > 0:
            test_val = time_values[0]
            # Check if it's a datetime type
            if isinstance(test_val, (np.datetime64, dt.datetime)):
                is_datetime = True
            elif hasattr(test_val, 'dtype'):
                # Check dtype - datetime64 or numeric
                if np.issubdtype(test_val.dtype, np.datetime64):
                    is_datetime = True
                elif np.issubdtype(test_val.dtype, np.number):
                    # Values are numeric - check if units suggest malformed time
                    # Malformed units like "seconds since 10890000.0" indicate time since model start
                    if 'since' in time_units:
                        ref_part = time_units.split('since')[1].strip()
                        # If reference part is just a number (not a date string), it's malformed
                        # Check if it's a pure number (possibly with decimal point)
                        ref_clean = ref_part.replace('.', '').replace('-', '').replace(':', '').replace(' ', '').replace('T', '')
                        if ref_clean.isdigit() or ('.' in ref_part and ref_clean.replace('e', '').replace('E', '').replace('+', '').isdigit()):
                            # Numeric reference (e.g., "10890000.0") - malformed, needs reference time
                            is_datetime = False
                        else:
                            # Date string in units, but values are numeric - not decoded properly
                            is_datetime = False
                    else:
                        # No 'since' in units and values are numeric - not decoded
                        is_datetime = False
        
        if not is_datetime:
            # Time is numeric, need to close and reopen with decode_times=False
            dsout2d.close()
            print(f"{file_progress} Time format detected: numeric values (seconds since model start)")
            print(f"{file_progress}   Time units: {time_units}")
            print(f"{file_progress}   Opening with decode_times=False and will use reference time from param.nml")
            dsout2d = xr.open_dataset(infiles[0], decode_times=False)
            time_format = 'seconds_since_start'  # Needs reference time
        else:
            time_format = 'datetime'  # Properly formatted
            print(f"{file_progress} Time format detected: properly formatted datetime")
    except (ValueError, TypeError) as e:
        if 'unable to decode time units' in str(e) or 'decode_times' in str(e) or 'time' in str(e).lower():
            print(f"{file_progress} Time format issue detected (exception during opening)")
            print(f"{file_progress}   Error: {str(e)}")
            print(f"{file_progress}   Opening with decode_times=False")
            dsout2d = xr.open_dataset(infiles[0], decode_times=False)
            time_format = 'seconds_since_start'  # Needs reference time
        else:
            raise
    
    # Validate dataset structure - check for required variables/dimensions
    required_dims = ['time']
    missing_dims = [dim for dim in required_dims if dim not in dsout2d.dims]
    if missing_dims:
        print(f"{file_progress} Error: Dataset {infiles[0]} missing required dimensions: {missing_dims}. Skipping this date.")
        dsout2d.close()
        continue
    
    if 'dryFlagElement' not in dsout2d.variables:
        print(f"{file_progress} Error: Dataset {infiles[0]} missing required variable: 'dryFlagElement'. Skipping this date.")
        dsout2d.close()
        continue
    
    # Open additional datasets (skip out2d which is already open)
    # Use same time decoding approach as out2d
    ds_list = []
    for i in range(1, len(infiles)):
        try:
            if time_format == 'seconds_since_start':
                # Open with decode_times=False to match out2d
                ds = xr.open_dataset(infiles[i], decode_times=False)
            else:
                ds = xr.open_dataset(infiles[i], decode_times=True)
            
            if interp_surface_only:
                # Extract surface layer only for 3D variables (still lazy)
                if 'nSCHISM_vgrid_layers' in ds.dims:
                    ds = ds.sel(nSCHISM_vgrid_layers=-1)
            ds_list.append(ds)
        except Exception as e:
            print(f"{file_progress} Warning: Could not open {infiles[i]}: {e}. Skipping this variable type.")
            continue

    # Align all datasets to avoid NaNs (works with lazy arrays, files still open)
    dsout2d, *ds_list = xr.align(dsout2d, *ds_list, join="outer")

    # Merge datasets (works with lazy arrays, files still open)
    dsin = xr.merge([dsout2d] + ds_list)
    
    print(f"{file_progress} Loading data into memory...")
    # Now load into memory and close files - we need data for numpy operations
    # This is done after lazy operations to benefit from xarray's optimizations
    dsin = dsin.load()
    dsout2d.close()
    for ds in ds_list:
        ds.close()

    # Determine number of time steps from dataset
    nt = len(dsin.time)
    print(f"{file_progress} Dataset loaded: {nt} time steps, {len(dsin.data_vars)} variables")

    # Output file name
    if interp_surface_only:
        fname = infiles[0].split('/')[-1].replace('out2d', 'schism-wwm_interp')
    else:
        fname = infiles[0].split('/')[-1].replace('out2d', 'schism-wwm_3Dinterp')
    
    outfile = outdir + fname
    outfiles.append(outfile)
    print(f"{file_progress} Output file: {outfile}")
    
    if (not os.path.exists(outfile)) or overwrite:
        if os.path.exists(outfile) and overwrite:
            print(f"{file_progress} Overwriting existing file...")
        print(f"{file_progress} Starting interpolation...")
        # Calculate velocity magnitude if velocity components exist
        if "horizontalVelX" in dsin.keys():
            dsin['velocity_magnitude'] = np.sqrt(dsin['horizontalVelX'] ** 2 + dsin['horizontalVelY'] ** 2)

        # Setup for 3D interpolation if needed
        if not interp_surface_only:
            s.ds = dsin
            s.nz = dsin['zCoordinates'].shape[-1]
            s.ibbtm = s.ds.bottom_index_node[:].values - 1
            s.mask3d = np.zeros((s.nnodes, s.nz), bool)
            for inode in range(s.nnodes):
                s.mask3d[inode, :s.ibbtm[inode]] = True

        # Reformat time - handle both SCHISM time formats
        # time_format was already determined when opening the dataset
        time = dsin.time.values
        attrs = dsin.time.attrs.copy()
        
        if time_format == 'seconds_since_start':
            # Time is in seconds since model start - need reference time from param.nml
            if reftime is None:
                print(f"{file_progress} Error: Time format is 'seconds_since_start' but reference time (reftime) is None.")
                print(f"{file_progress}   This should not happen if param.nml was read correctly.")
                raise ValueError("Reference time (reftime) is required for time conversion but is None")
            
            # Use the reference time from param.nml
            t0 = reftime
            
            # Time values are already in seconds since model start
            # They represent elapsed time, so we keep them as-is but set the reference
            seconds = time.astype(float)  # Ensure numeric
            
            print(f"{file_progress} Converting time: {len(seconds)} time steps")
            print(f"{file_progress}   Time range: {seconds[0]:.0f} to {seconds[-1]:.0f} seconds since model start")
            print(f"{file_progress}   Reference time: {t0}")
            
        else:
            # Time is already in datetime format
            t0 = time[0]
            try:
                seconds = (time - t0) / np.timedelta64(1, 's')
            except:
                # Fallback if timedelta conversion fails
                if hasattr(time, 'astype'):
                    # Try converting nanoseconds to seconds
                    time_diff = time - t0
                    if hasattr(time_diff, 'astype'):
                        seconds = time_diff.astype('timedelta64[s]').astype(float)
                    else:
                        seconds = np.array([(t - t0).total_seconds() for t in time])
                else:
                    seconds = np.array([(t - t0).total_seconds() for t in time])
            
            print(f"{file_progress} Converting time: {len(seconds)} time steps from datetime format")
        
        # Format reference time for CF-compliant units string
        # Ensure we have a scalar datetime64 value
        if hasattr(t0, 'item') and not isinstance(t0, (np.datetime64, dt.datetime)):
            # numpy array - extract scalar
            t0 = t0.item()
        elif isinstance(t0, np.ndarray) and t0.size == 1:
            # numpy array with one element
            t0 = t0.item()
        
        # Now format as string
        if isinstance(t0, np.datetime64):
            t0_str = str(t0)[:19].replace('T', ' ')
        elif isinstance(t0, dt.datetime):
            t0_str = t0.strftime('%Y-%m-%d %H:%M:%S')
        else:
            # Try to convert to datetime64 and format
            try:
                t0_dt64 = np.datetime64(t0)
                t0_str = str(t0_dt64)[:19].replace('T', ' ')
            except:
                t0_str = str(t0)[:19].replace('T', ' ')
        
        tunit = "seconds since {:s}".format(t0_str)
        attrs['units'] = tunit
        attrs['_FillValue'] = False
        
        print(f"{file_progress} Time units set to: {tunit}")

        # Process each variable
        total_vars = len(varnames)
        for var_idx, varname in enumerate(varnames):
            var_progress = f"[{var_idx + 1}/{total_vars}]"
            if varname not in dsin.keys():
                print(f"{file_progress} {var_progress} Warning: {varname} not found in dataset, skipping...")
                continue
            
            # Validate variable has required metadata in YAML
            if varname not in names:
                print(f"{file_progress} {var_progress} Warning: {varname} missing from 'names' in YAML, using variable name as-is")
            if varname not in units:
                raise KeyError(f"Variable {varname} missing 'units' definition in YAML")
            if varname not in valid_ranges:
                raise KeyError(f"Variable {varname} missing 'valid_ranges' definition in YAML")
                
            # Determine if variable is 3D
            if 'nSCHISM_vgrid_layers' in dsin[varname].dims:
                dim3D = True
                var_type = "3D"
            else:
                dim3D = False
                var_type = "2D"
            
            print(f"{file_progress} {var_progress} Processing {varname} ({var_type})...")
            
            # Select appropriate standard name and long name based on 2D/3D and interpolation mode
            # For 2D variables or surface-only extraction: use 2D names
            # For full 3D interpolation: use 3D names
            if dim3D and not interp_surface_only:
                std_name = std_names_3d.get(varname)
                long_name = long_names_3d.get(varname, long_names_2d.get(varname, varname))
            else:
                std_name = std_names_2d.get(varname)
                long_name = long_names_2d.get(varname, varname)

            # Initialize data array
            if dim3D and not interp_surface_only:
                # Full 3D interpolation
                data = np.ones((nt, nz, ny, nx)) * FillValue
            else:
                # 2D or surface-only
                data = np.ones((nt,) + X.shape) * FillValue

            iuse = (parents != -1)
            target_inds_all = np.where(iuse)[0]
            ii_all, jj_all = np.unravel_index(target_inds_all, X.shape)

            # Progress tracking for time steps
            if nt > 10:  # Only show progress for files with many time steps
                progress_interval = max(1, nt // 10)  # Update every 10%
            else:
                progress_interval = 1
            
            for t in range(nt):
                if nt > 10 and t % progress_interval == 0:
                    progress_pct = int((t + 1) / nt * 100)
                    print(f"{file_progress} {var_progress}   Time step {t + 1}/{nt} ({progress_pct}%)", end='\r')
                wet = dsin.dryFlagElement[t, :].values == 0
                iuse = (parents != -1)
                valid_parents = parents[iuse]
                ivalid = wet[valid_parents]
                valid_parents = valid_parents[ivalid]
                weights = ndeweights[iuse, :][ivalid, :]
                target_inds = np.where(iuse)[0][ivalid]
                ii, jj = np.unravel_index(target_inds, X.shape)

                if dim3D and not interp_surface_only:
                    # Full 3D interpolation
                    for k in range(nz):
                        ibelow, iabove, weights_vert = get_layer_weights(s, -z_levels[k], t)
                        s.nodeinds = np.arange(s.nnodes)
                        field3D = dsin[varname][t, :].values
                        interp_slice = weights_vert[0, :] * field3D[s.nodeinds, ibelow] + \
                                      weights_vert[1, :] * field3D[:, :][s.nodeinds, iabove]
                        interp_slice[np.isnan(interp_slice)] = FillValue
                        data[t, k, ii, jj] = (interp_slice[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
                else:
                    # 2D or surface-only interpolation
                    interp_slice = dsin[varname][t, :].values
                    data[t, ii, jj] = (interp_slice[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
            
            if nt > 10:
                print()  # New line after progress updates

            valid_data = data.flatten()[data.flatten() != FillValue]
            
            if len(valid_data) == 0:
                print(f"{file_progress} {var_progress} Warning: No valid data for {varname}, skipping...")
                continue
            
            print(f"{file_progress} {var_progress} Writing {varname} to output file...")

            # Create coordinate arrays
            if varname == 'elevation':
                dslat = xr.DataArray(name="latitude", data=latitude, dims=["latitude"], 
                                     coords=dict(latitude=(["latitude"], latitude)),
                                     attrs=dict(description="latitude", units="degrees_north",
                                               standard_name="latitude", _FillValue=False))

                dslon = xr.DataArray(name="longitude", data=longitude, dims=["longitude"],
                                     coords=dict(longitude=(["longitude"], longitude)),
                                     attrs=dict(description="longitude", units="degrees_east",
                                               standard_name="longitude", _FillValue=False))

                dstime = xr.DataArray(name="time", data=seconds, dims=["time"],
                                     coords=dict(time=(["time"], seconds)), attrs=attrs)

                # Create depth coordinate if doing 3D
                if not interp_surface_only:
                    dsdepth = xr.DataArray(name="depth", data=z_levels, dims=["depth"],
                                          coords={"depth": (["depth"], z_levels)},
                                          attrs={"standard_name": "depth",
                                                "long_name": "Depth below sea surface",
                                                "units": "m", "positive": "down", "axis": "Z",
                                                "_FillValue": FillValue})

                # Create data array
                da = xr.DataArray(name=names[varname], data=data,
                                 dims=["time", "latitude", "longitude"],
                                 coords=dict(latitude=dslat, longitude=dslon, time=dstime),
                                 attrs=dict(description="ssh", units=units[varname],
                                           long_name=long_name, _FillValue=FillValue,
                                           actual_range=[np.min(valid_data), np.max(valid_data)],
                                           valid_range=valid_ranges[varname],
                                           **({ 'standard_name': std_name } 
                                              if std_name is not None else {})))

                # Create dataset
                if not interp_surface_only:
                    # Name depth coordinate as "depth" so 3D variables can reference it
                    ds = xr.Dataset({da.name: da, "depth": dsdepth})
                else:
                    ds = xr.Dataset({da.name: da})
                ds.attrs.update(global_atts)
                ds.to_netcdf(outfile, mode="w")

            elif dim3D and not interp_surface_only:
                # 3D variable - need to reference existing depth coordinate
                # Open existing file to get depth coordinate using context manager
                with xr.open_dataset(outfile) as ds_existing:
                    depth_coord = ds_existing.depth.copy()
                
                # Create data array with depth coordinate reference
                da = xr.DataArray(name=names[varname], data=data,
                                 dims=["time", "depth", "latitude", "longitude"],
                                 coords={"depth": depth_coord},
                                 attrs=dict(description=varname, units=units[varname],
                                           long_name=long_name, _FillValue=FillValue,
                                           actual_range=[np.min(valid_data), np.max(valid_data)],
                                           valid_range=valid_ranges[varname],
                                           **({ 'standard_name': std_name } 
                                              if std_name is not None else {})))
                da.to_netcdf(outfile, mode='a')

            else:
                # 2D variable or surface-only
                da = xr.DataArray(name=names[varname], data=data,
                                 dims=["time", "latitude", "longitude"],
                                 attrs=dict(description=varname, units=units[varname],
                                           long_name=long_name, _FillValue=FillValue,
                                           actual_range=[np.min(valid_data), np.max(valid_data)],
                                           valid_range=valid_ranges[varname],
                                           **({ 'standard_name': std_name } 
                                              if std_name is not None else {})))
                da.to_netcdf(outfile, mode='a')

        # Add depth field (bathymetry) if doing surface-only interpolation
        if interp_surface_only:
            Ddata = np.ones(X.shape) * FillValue
            iuse = (parents != -1)
            valid_parents = parents[iuse]
            weights = ndeweights[iuse, :]
            target_inds = np.where(iuse)[0]
            ii, jj = np.unravel_index(target_inds, X.shape)
            Ddata[ii, jj] = (np.asarray(s.depths)[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
            valid_data = Ddata[ii, jj].flatten()
            varname = 'depth'
            
            if 'depth' in names:
                # Depth is always 2D, so use 2D standard name
                depth_std_name = std_names_2d.get('depth')
                da_depth = xr.DataArray(name=names[varname], data=Ddata,
                                       dims=["latitude", "longitude"],
                                       coords=dict(latitude=dslat, longitude=dslon),
                                       attrs=dict(description="bathymetry", units=units[varname],
                                                 long_name=long_names_2d.get('depth', 'depth'), _FillValue=FillValue,
                                                 actual_range=[np.min(valid_data), np.max(valid_data)],
                                                 valid_range=valid_ranges[varname],
                                                 **({ 'standard_name': depth_std_name } 
                                                    if depth_std_name is not None else {})))
                da_depth.to_netcdf(outfile, mode='a')
        
        print(f"{file_progress} ✓ File processing complete: {outfile}")
    else:
        print(f"{file_progress} ⏭ Skipping (file exists and overwrite=False): {outfile}")

print(f"\n{'='*70}")
print(f"Interpolation complete! Processed {len(outfiles)} file(s)")
print(f"Output directory: {outdir}")
print(f"{'='*70}\n")

