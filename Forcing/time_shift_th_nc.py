# interp 12 to hourly	
import xarray as xr
from glob import glob
import numpy as np
import sys
import os
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities')
sys.path.insert(0,'/work/gg0028/g260114/RUNS/GermanBight/GB_HR_Ballje')
sys.path.insert(0,'/gpfs/work/jacobb/data/RUNS/GB_template/schism-hzg-utilities/')
from schism import *
s=schism_setup()
nz= len(s.vgrid[1])


files=['elev2D.th.nc','SAL_3D.th.nc','TEM_3D.th.nc','uv3D.th.nc']

nshift=2 # shift forcing by nshift timesteps

for fname in files:
	print(fname)
	fnameout=fname[:-2]+'2hour_shift.nc'
	ds=xr.open_dataset(fname)
	time=ds['time'].values
	out=ds['time_series'].values

	s.write_bdy_netcdf(fnameout,time[:-nshift],out[nshift:])	
	os.rename(fname,fname+'original')
	os.symlink(fnameout,fname)
	
