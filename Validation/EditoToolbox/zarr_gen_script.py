#!/usr/bin/env python3
"""
Create a combined Zarr dataset for the validation toolbox
from SCHISM NetCDF outputs.
"""

import os
import sys
from pathlib import Path
import xarray as xr
import numpy as np
import shutil


#crashess for parallel=True writing

# -----------------------------------------------------------------------------
# USER CONFIG
# -----------------------------------------------------------------------------
SCHISM_UTILS = '/gpfs/work/ksddata/code/schism/scripts/schism-hereon-utilities/'
RUN_DIR = Path('/gpfs/work/bergm/data/Seyed_Setup/Run01_wpile')
NC_DIR = RUN_DIR / 'outputs'

ZARR_OUT = Path('/gpfs/work/work_scratch/jacobb/zarrSNS/schism4valid.zarr')
ZARR_OUT.parent.mkdir(parents=True, exist_ok=True)

VARIABLE_FILE = 'out2d'   # e.g. 'out2d', 'salinity', etc.
STACKS = (1, 7)           # inclusive start, exclusive end

# -----------------------------------------------------------------------------
# IMPORT SCHISM SETUP
# -----------------------------------------------------------------------------
sys.path.insert(0, SCHISM_UTILS)
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
        parallel=True,
        chunks="auto",   # 👈 your requested setting
    )


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
def main():

    # --- Load SCHISM grid ---
    os.chdir(RUN_DIR)
    s = schism_setup()

    os.chdir(NC_DIR)

    print("Opening datasets...")

    if VARIABLE_FILE != "out2d":
        # Load coordinate dataset
        files2d = build_file_list("out2d", STACKS)
        ds2d = open_dataset(files2d)

        # Load variable dataset
        files = build_file_list(VARIABLE_FILE, STACKS)
        ds = open_dataset(files)

        print("Merging coordinate variables from out2d...")

        coord_vars = {
            var: ds2d[var]
            for var in ds2d.variables
            if "SCHISM" in var
        }

        ds = ds.assign_coords(coord_vars)

    else:
        files = build_file_list(VARIABLE_FILE, STACKS)
        ds = open_dataset(files)

        print("Assigning lon/lat from SCHISM setup...")

        ds = ds.assign_coords({
            "SCHISM_hgrid_node_x": ("nSCHISM_hgrid_node", np.asarray(s.lon)),
            "SCHISM_hgrid_node_y": ("nSCHISM_hgrid_node", np.asarray(s.lat)),
        })

    # -----------------------------------------------------------------------------
    # WRITE ZARR
    # -----------------------------------------------------------------------------
    if ZARR_OUT.exists():
        print(f"Removing existing Zarr store: {ZARR_OUT}")
        shutil.rmtree(ZARR_OUT)

    print(f"Writing Zarr to: {ZARR_OUT}")

    ds.to_zarr(
        ZARR_OUT,
        mode="w",
        consolidated=True,   # important for fast reads
    )

    print("Done.")


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()