
# write a combined zarr file understandable  for the validation toolbox
# based on schism netcdf

Write_Combined_Zarr_For_Validation_Toolbox.py

#conda activate validation-toolbox

import sys
sys.path.insert(0, '/gpfs/work/ksddata/code/schism/scripts/schism-hereon-utilities/')
from schism import *


run_dir='/gpfs/work/bergm/data/Seyed_Setup/Run01_wpile'
nc_dir='/gpfs/work/bergm/data/Seyed_Setup/Run01_wpile/outputs/'

zarr_out_dir='/gpfs/work/work_scratch/jacobb/zarr/' #  merged zarr output folder
zarr_name='schism4valid.zarr'

os.makedirs(zarr_out_dir)
outfile=os.path.join(zarr_out_dir,zarr_name)

variable_file='out2d'
stacks=(1,7)


# load schism setup
os.chdir(run_dir)
s=schism_setup()


## we always will need the coordinate information which is only in the out2d files
# so for other variable files we need to append




os.chdir(nc_dir)

if variable_file!= 'out2d': #load also varaible fiel for coordiantes

    files2d=[f'out2d_2025120{i}.nc' for i in range(1,3)]
    files=[f'variable_file_2025120{i}.nc' for i in range(1,3)]


    ds2d=xr.open_mfdataset(files2d)
    ds=xr.open_mfdataset(files)

    #for var in ds2d.variables:
    #    if 'SCHISM' in var:
    #        print(var)
            #ds[var]=ds2d[var]

    # add the variables only in out2d to the variable for valid        
    for var in ds2d.variables:
        if 'SCHISM' in var:
            print(var)
            ds = ds.assign_coords({var: ds2d[var]})        

else:
    files=[f'{variable_file}_{i}.nc' for i in range(stacks[0],stacks[1])]

    ds=xr.open_mfdataset(files)
    # need to enforce lon lat coordinates
    ds.SCHISM_hgrid_node_x[:]=np.asarray(s.lon)
    ds.SCHISM_hgrid_node_y[:]=np.asarray(s.lat)

#open 2ds varaible to append cooridnate varialbes

ds.to_zarr(outfile)

