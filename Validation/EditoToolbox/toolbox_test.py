import os
import sys
import logging
import json
import yaml
from pathlib import Path
import xarray as xr
import time

# insert path of config file:
sys.path.append('/work/gg0028/Validation_Toolbox/wp3-validation-toolbox/toolbox/')
#os.chdir('/work/gg0028/Validation_Toolbox/wp3-validation-toolbox/toolbox/')
os.chdir('/work/gg0028/Validation_Toolbox/wp3-validation-toolbox/toolbox/')


# --- Imports from your provided scripts ---
try:
    # Configuration loading
    from toolbox.config import ValidationConfig
    
    # Core validation logic from gridded_validation_EDITO_Hackathon.py
    # (Import necessary functions directly or adapt the main logic)
    from toolbox.insitu_mooring_validation import (
        setup_environment, 
        load_datasets, 
        gather_validation_info, 
        display_validation_info,
        process_data, 
        compute_validation_metrics,
        save_validation_results)
#        run_validation # You can call this directly or replicate its steps
#    )
    
    # Plotting logic (optional, can be run after validation)
    from toolbox.insitu_plots import process_validation_results as run_plotting
    
    # Dataclasses (ensure this is importable)
    from toolbox.core.dataclasses import InSituValidationData, DomainBounds

    
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Please ensure configs.py, insitu_mooring_validation.py, insitu_mooring_plots.py, and dataclasses.py are in the Python path or the current directory.")
    # You might need to install dependencies like xarray, click, pydantic, matplotlib, etc.
    # Example: !pip install xarray pandas matplotlib click pydantic netcdf4 dask

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout # Log to notebook output
)
logger = logging.getLogger("validation_notebook")

print("Setup and imports complete.")

#########################################

# --- Configuration Setup ---
#config_file_path = '../configs/insitu_TEMP_NorthSea-CLIMAT-zarr.yaml'
#config_file_path = '../configs/insitu_SALT_NorthSea-CLIMAT-zarr.yaml'
#config_file_path = '../configs/insitu_SALT_NorthSea-CLIMAT.yaml'
config_file_path = '../configs/insitu_SALT_NorthSea-CLIMAT.yaml'

try:
    #with open(config_file_path, 'r') as f:
    #    config_dict = json.load(f)
    with open(config_file_path, 'r') as f:
        config_dict = yaml.safe_load(f)
    logger.info(f"Loaded configuration from: {config_file_path}")
    
    # Update the global configs using the function from configs.py
    #update_config(config_dict)
    
    # Create a ValidationConfig instance from the 'validation' part of the configs
    # This mimics how gridded_validation.py uses the configs
    if 'validation' in config_dict:
         # Use ** to unpack the dictionary into keyword arguments
        validation_settings = config_dict['validation']
      
        # Add other top-level settings if needed by ValidationConfig Pydantic model
        # Example: Add interpolation settings if they are direct fields in ValidationConfig
        if 'interpolation' in config_dict:
            #print('Yes, interpolation is in dict')
            if 'spatial_tolerance' in config_dict['interpolation']:
                validation_settings['spatial_tolerance'] = config_dict['interpolation']['spatial_tolerance']
            if 'max_temporal_tolerance' in config_dict['interpolation'] and 'days' in config_dict['interpolation']['max_temporal_tolerance']:
                validation_settings['temporal_tolerance_days'] = config_dict['interpolation']['max_temporal_tolerance']['days']
        #print(validation_settings)
        
        # Instantiate the configs object
        #config = ValidationConfig(**validation_settings)
        config = ValidationConfig.load_file(config_file_path)
    else:
        # If no specific 'validation' section, try using the whole dict
        # This might need adjustments based on your ValidationConfig model in configs.py
        logger.warning("No 'validation' section in configs file, attempting to use the full configs.")
        config = ValidationConfig(**config_dict)

    # --- Modify specific configurations here if needed ---
    # Example: Override the output directory
    # configs.output_dir = './validation_output_notebook'
    # set_config('validation.output_dir', './validation_output_notebook') # Update global configs too
    
    # Set skip_confirmation to True for non-interactive execution in the notebook
    config.skip_confirmation = True 
    config.yes = True # Corresponds to -y flag
    # set_config('validation.skip_confirmation', True)
    # set_config('validation.yes', True)


    logger.info(f"Using Model URL: {config.model_url}")
    logger.info(f"Using Source Data URI: {config.observations_uri}")
    logger.info(f"Using Output Directory: {config.output_dir}")
    logger.info(f"Skip confirmation: {config.skip_confirmation}")

except FileNotFoundError:
    logger.error(f"Configuration file not found: {config_file_path}")
    config = None
except Exception as e:
    logger.error(f"Error loading or processing configuration: {e}")
    config = None



################ RUN Validation ################

validation_successful = False
validation_results_dir = None

logger.info(f"Starting {config.dataset_name.upper()} validation workflow")

# Step 1: Setup environment (create output dir)
config.output_dir = setup_environment(config.output_dir)
validation_results_dir = config.output_dir # Store for plotting step
logger.info(f"Output directory set to: {validation_results_dir}")

# Step 2-3: Load datasets
logger.info("Loading datasets...")
# Pass the configs *object* to load_datasets --> ~300 seconds
model_dataset, obs_data, model_domain, variable_units = load_datasets(config)
logger.info("Datasets loaded successfully.")

# (Optional) Step 4-5: Gather and display validation info (skipped in notebook)
# validation_info = gather_validation_info(model_dataset, obs_data, model_domain, obs_domain)
# Requires user interaction
logger.info("Skipping interactive validation info display.")


# Step 6-7: Process and align datasets
logger.info("Processing and aligning data...")
model_insitu, obs_data_arr, obs_data_resample = (model_dataset, obs_data, config)
logger.info("Data processing complete.")



 surface=_model_surface_2d(model_dataset)

def _model_surface_2d(model_data: xr.DataArray) -> xr.DataArray:
    """time=0 and shallowest vertical level so mask/extraction stays 2D horizontal."""
    logger.debug("[spatial_temporal_alignment] _model_surface_2d: start")
    t0 = time.perf_counter()
    da = model_data.isel(time=0)
    for dim in ("depth", "z", "level"):
        if dim in da.dims:
            da = da.isel({dim: 0})
    if len(da.dims) != 2:
        da = da.squeeze(drop=True)
    if len(da.dims) != 2:
        raise ValueError(
            f"Expected 2D (latitude, longitude) after time=0 and surface slice; got dims {da.dims}"
        )
    logger.debug(
        f"[spatial_temporal_alignment] _model_surface_2d: done in {time.perf_counter() - t0:.3f}s, "
        f"out_dims={dict(da.sizes)}"
    )
    return da





def process_data(model_dataset: xr.Dataset, obs_data: xr.Dataset,
                 config: ValidationConfig) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Process and align datasets for validation.

    Args:
        model_dataset: Model dataset
        obs_data: Observation dataset
        config: Validation configuration

    Returns:
        Tuple of aligned model and observation data arrays
    """
    # Get variable-specific preparers
    if not config.variable_type:
        raise ValueError(
            "variable_type must be specified in config. "
            "Add 'variable_type: temperature' (or 'chlorophyll', 'ssh', etc.) to your config file."
        )

    logger.info(f"Using variable type: {config.variable_type}")
    model_preparer, obs_preparer = get_preparers(config.variable_type)
    # Prepare datasets using variable-specific preparers
    logger.info(f"Preparing model data with variable: {config.model_var}")
    model_data = model_preparer(model_dataset, config.model_var)
    #

    logger.info(f"Preparing observation data with variable: {config.obs_var}")
    obs_data_arr = {}
    for site in obs_data.keys():
        obs_data_arr[site] = obs_preparer(obs_data[site], config.obs_var)

    # Align datasets through regridding
    logger.info("Selecting model data to match observation locations")
    model_insitu, excluded_sites = select_insitu_gridpoint(model_data, obs_data_arr)

    # Check if there are excluded sites to then remove them for the insitu dict
    if excluded_sites:  # True if not empty
        for site in excluded_sites.keys():
            obs_data_arr.pop(site, None)
        logger.info(f"Removed {len(excluded_sites)} excluded sites.")
    else:
        logger.info("No excluded sites to remove.")

    #Resample each in-situ DataArray to match the temporal resolution of the model_data
    logger.info("Resampling insitu observations to model temporal resolution")
    obs_data_resample = resample_insitu_obs(model_data, obs_data_arr)

    return model_insitu, obs_data_arr, obs_data_resample












def compute_station_grid_indices(
    model_data: xr.DataArray,
    obs_data: dict[str, xr.DataArray],
    *,
    max_distance_km: float = 10.0,
    max_window_radius: int = 32,
    k_neighbors: int = 10,
    mask_netcdf_path: Optional[str] = None,
    mask_file_engine: str = "netcdf4",
    mask_file_data_var: Optional[str] = None,
) -> tuple[dict[str, tuple[int, int]], dict[str, str]]:
    """
    Compute (iy, ix) indices of nearest valid model grid point for each station.

    If latitude and longitude are uniformly spaced (regular grid), uses a cheap
    index + expanding-window search. Otherwise builds a cKDTree once (requires scipy).

    Intended for the sequential_files strategy: compute once, then reuse indices
    while iterating through a NetCDF file list.

    ``mask_netcdf_path`` / ``mask_file_*``: read the 2D surface mask from one file so
    ``open_mfdataset`` dask arrays are not computed for this step (speed + avoids OOM).
    """
    if not isinstance(obs_data, dict):
        raise TypeError(f"obs_data must be dict[str, DataArray], got {type(obs_data)}")

    logger.info(
        f"[spatial_temporal_alignment] compute_station_grid_indices: START n_sites={len(obs_data)}"
    )
    _log_xarray_step("compute_station_grid_indices model_data", model_data)

    t0 = time.perf_counter()
    lats = np.asarray(model_data["latitude"].values)
    lons = np.asarray(model_data["longitude"].values)
    logger.info(
        f"[spatial_temporal_alignment] compute_station_grid_indices: lat/lon coords loaded in "
        f"{time.perf_counter() - t0:.3f}s"
    )
    regular = _is_uniform_spacing_1d(lats) and _is_uniform_spacing_1d(lons)
    logger.info(
        f"[spatial_temporal_alignment] compute_station_grid_indices: uniform_spacing_regular={regular}"
    )

    indices: dict[str, tuple[int, int]] = {}
    excluded: dict[str, str] = {}

    if regular:
        logger.info(
            "Regular lat/lon spacing: using index + local window nearest valid point"
        )
        t_loop = time.perf_counter()
        # One materialization for all stations: repeated dask .values in the window loop
        # is catastrophically slow on multi-file mfdataset (can look "stuck" for >1h).
        t_mat = time.perf_counter()
        surface_np = _materialize_mask_surface_np(
            model_data,
            mask_netcdf_path=mask_netcdf_path,
            mask_file_engine=mask_file_engine,
            mask_file_data_var=mask_file_data_var,
        )
        logger.info(
            f"[spatial_temporal_alignment] regular path: surface array ready in "
            f"{time.perf_counter() - t_mat:.3f}s"
        )
        _log_array_bytes("regular path surface_np (shared by all stations)", surface_np)
        if surface_np.shape != (lats.size, lons.size):
            logger.warning(
                f"surface_np shape {surface_np.shape} vs expected ({lats.size}, {lons.size}) "
                "(check latitude/longitude dimension order)"
            )

        for site, da in obs_data.items():
            logger.info(f"Selecting nearest valid grid point for {site}")
            iy, ix, invalid_site = nearest_valid_point_indexed(
                model_data,
                da,
                lats=lats,
                lons=lons,
                max_distance_km=max_distance_km,
                max_window_radius=max_window_radius,
                surface_np=surface_np,
            )
            if invalid_site is not None or iy is None or ix is None:
                excluded[invalid_site or site] = "Too far from model grid"
            else:
                indices[site] = (int(iy), int(ix))
        logger.info(
            f"[spatial_temporal_alignment] compute_station_grid_indices: regular-path station loop "
            f"done in {time.perf_counter() - t_loop:.3f}s"
        )
    else:
        logger.info("Irregular lat/lon spacing: using cKDTree (one build, all stations)")
        t_tree = time.perf_counter()
        kdtree_index = _build_kdtree_index(
            model_data,
            mask_netcdf_path=mask_netcdf_path,
            mask_file_engine=mask_file_engine,
            mask_file_data_var=mask_file_data_var,
        )
        logger.info(
            f"[spatial_temporal_alignment] compute_station_grid_indices: KDTree index ready in "
            f"{time.perf_counter() - t_tree:.3f}s"
        )
        t_loop = time.perf_counter()
        for site, da in obs_data.items():
            logger.info(f"Selecting nearest valid grid point for {site}")
            iy, ix, invalid_site = _nearest_valid_point_kdtree(
                kdtree_index,
                da,
                max_distance_km=max_distance_km,
                k_neighbors=k_neighbors,
            )
            if invalid_site is not None or iy is None or ix is None:
                excluded[invalid_site or site] = "Too far from model grid"
            else:
                indices[site] = (int(iy), int(ix))
        logger.info(
            f"[spatial_temporal_alignment] compute_station_grid_indices: kdtree-path station loop "
            f"done in {time.perf_counter() - t_loop:.3f}s"
        )

    logger.info(
        f"[spatial_temporal_alignment] compute_station_grid_indices: DONE n_ok={len(indices)} "
        f"n_excluded={len(excluded)}"
    )
    return indices, excluded


def _nearest_1d_index(values: np.ndarray, target: float) -> int:
    # values are assumed monotonic; argmin works fine either way
    return int(np.abs(values - target).argmin())







#just nn lon lat manual
lons=[]
lats=[]
for key in obs_data.keys():
    lons.append(obs_data[key].longitude.values)
    lats.append(obs_data[key].latitude.values)
model_insitu=model_dataset.interp(longitude=lons,latitude=lats)

    
model_insitu    obs_data

 





# Step 8: Compute validation metrics
logger.info("Computing validation metrics...")
validation_data = compute_validation_metrics(model_insitu, obs_data_arr, obs_data_resample, variable_units)
logger.info("Metrics computation complete.")

# Step 9: Save results
logger.info("Saving validation results...")
metrics_file, insitu_file, summary_file = save_validation_results(validation_data, config, variable_units)
logger.info(f"Results saved:")
logger.info(f"  Metrics: {metrics_file}")
logger.info(f"  Datasets: {insitu_file}")
logger.info(f"  Summary: {summary_file}")

logger.info(f"Validation workflow completed successfully for dataset '{config.dataset_name}'")
print(f"\nValidation complete. Results saved to: {config.output_dir}")
validation_successful = True

