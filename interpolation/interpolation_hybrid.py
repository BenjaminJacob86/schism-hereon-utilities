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

# Validate file lists are not empty
empty_file_lists = [key for key, files in varfiles.items() if len(files) == 0]
if empty_file_lists:
    raise FileNotFoundError(f"No files found for variables: {empty_file_lists}. "
                          f"Check input directory: {indir}")

# Validate all variable files have the same number of files
file_counts = {key: len(files) for key, files in varfiles.items()}
if len(set(file_counts.values())) > 1:
    raise ValueError(f"Inconsistent number of files across variables: {file_counts}")

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
    
    infiles = [varfiles[key][idate] for key in varfiles.keys()]

    # Validate input files exist
    missing_files = [f for f in infiles if not os.path.exists(f)]
    if missing_files:
        raise FileNotFoundError(f"Input files not found: {missing_files}")

    print(f"{file_progress} Opening and merging datasets...")
    # Open datasets - keep lazy for efficient memory usage during align/merge
    # xarray can optimize operations on lazy arrays, only loading what's needed
    dsout2d = xr.open_dataset(infiles[0])
    
    # Validate dataset structure - check for required variables/dimensions
    required_dims = ['time']
    missing_dims = [dim for dim in required_dims if dim not in dsout2d.dims]
    if missing_dims:
        raise ValueError(f"Dataset {infiles[0]} missing required dimensions: {missing_dims}")
    
    if 'dryFlagElement' not in dsout2d.variables:
        raise ValueError(f"Dataset {infiles[0]} missing required variable: 'dryFlagElement'")
    
    if interp_surface_only:
        # Extract surface layer only for 3D variables (still lazy)
        ds_list = [xr.open_dataset(infiles[i]).sel(nSCHISM_vgrid_layers=-1) 
                   for i in range(1, len(infiles))]
    else:
        # Load full 3D fields (still lazy)
        ds_list = [xr.open_dataset(infiles[i]) for i in range(1, len(infiles))]

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

        # Reformat time
        time = dsin.time.values
        attrs = dsin.time.attrs
        t0 = time[0]
        try:
            seconds = (time - t0) / np.timedelta64(1, 's')
        except:
            seconds = (time - 0)
        t0 = str(t0)[:19].replace('T', ' ')
        tunit = "seconds since {:s}".format(t0)
        attrs['units'] = tunit
        attrs['_FillValue'] = False

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

