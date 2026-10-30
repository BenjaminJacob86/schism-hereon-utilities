#!/usr/bin/env python3
"""
Create a combined Zarr dataset for the validation toolbox
from SCHISM ugrid NetCDF outputs.
"""

import os
import sys
from pathlib import Path
import xarray as xr
import numpy as np
import shutil

# -----------------------------------------------------------------------------
# USER CONFIG
# -----------------------------------------------------------------------------
SCHISM_UTILS = '/gpfs/work/ksddata/code/schism/scripts/schism-hereon-utilities/'
SCHISM_UTILS2 = '/work/gg0028/SCHISM/schism-hereon-utilities/'

RUN_DIR = Path('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/')

NC_DIR = RUN_DIR / 'outputs01/'

ZARR_OUT = Path('/scratch/g/g260114/GB/Veg_CNTRL/schismGB4valid.zarr')

ZARR_OUT.parent.mkdir(parents=True, exist_ok=True)

VARIABLE_FILES = ['out2d'] #, 'salinity', 'temperature']
STACKS = (1, 20)   # inclusive start, exclusive end

# -----------------------------------------------------------------------------
# IMPORT SCHISM SETUP
# -----------------------------------------------------------------------------
sys.path.insert(0, SCHISM_UTILS)
sys.path.insert(0, SCHISM_UTILS2)
from schism import schism_setup  # noqa

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
def build_file_list(prefix, stacks):
    return [
        f"{prefix}_{i}.nc"
        for i in range(stacks[0], stacks[1])
    ]


def open_dataset(files):
    return xr.open_mfdataset(
        files,
        combine="by_coords",
        parallel=False,
        chunks="auto",
    )


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
def main():

    # --- Load SCHISM grid ---
    os.chdir(RUN_DIR)
    s = schism_setup()

    os.chdir(NC_DIR)

    print("Opening out2d dataset...")
    
    # Load out2d once
    files2d = build_file_list("out2d", STACKS)
    ds2d = open_dataset(files2d)

    ds2d = ds2d.drop_vars('dryFlagSide') #not compatiblehere	

    # Collect SCHISM coordinate variables
    coord_vars = {
        var: ds2d[var]
        for var in ds2d.variables
        if "SCHISM" in var
    }

    datasets = []

    # -------------------------------------------------------------------------
    # LOAD ALL VARIABLE FILES
    # -------------------------------------------------------------------------
    for variable_file in VARIABLE_FILES:

        if  variable_file == 'out2d':
            ds=ds2d.copy()
        
        else:
            print(f"Opening dataset: {variable_file}")

            files = build_file_list(variable_file, STACKS)
            ds = open_dataset(files)

            # Attach SCHISM coordinate variables
            ds = ds.assign_coords(coord_vars)

        datasets.append(ds)

    # -------------------------------------------------------------------------
    # MERGE ALL VARIABLES
    # -------------------------------------------------------------------------
    print("Merging datasets...")

    ds_merged = xr.merge(datasets)
    
    
    # Overwrite time if time is not a proper format like in the older schism versions
    if type(ds_merged.time[0].values) != np.datetime64:
        from schism import param
        
        p=param(RUN_DIR/'param.nml')
        
        

        year = int(p.get_parameter('start_year'))
        month = int(p.get_parameter('start_month'))
        day = int(p.get_parameter('start_day'))
        hour = int(p.get_parameter('start_hour'))
        dt = np.datetime64(f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}")
        new_time= dt + ds_merged.time.values * np.timedelta64(1,'s')

        ds_merged = ds_merged.assign_coords(time=new_time)

    # -------------------------------------------------------------------------
    # ADD LON/LAT
    # -------------------------------------------------------------------------
    print("Assigning lon/lat from SCHISM setup...")

    ds_merged = ds_merged.assign_coords({
        "SCHISM_hgrid_node_x": (
            "nSCHISM_hgrid_node",
            np.asarray(s.lon)
        ),
        "SCHISM_hgrid_node_y": (
            "nSCHISM_hgrid_node",
            np.asarray(s.lat)
        ),
    })

    # -------------------------------------------------------------------------
    # WRITE ZARR
    # -------------------------------------------------------------------------
    if ZARR_OUT.exists():
        print(f"Removing existing Zarr store: {ZARR_OUT}")
        shutil.rmtree(ZARR_OUT)

    print(f"Writing Zarr to: {ZARR_OUT}")

    ds_merged.to_zarr(
        ZARR_OUT,
        mode="w",
        consolidated=True,
    )

    print("Done.")


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
