"""
merge DWD wind files that were
downloaded and preprocessed in ksd routine
for application for WW3 forecast
"""

import os
from netCDF4 import Dataset,MFDataset
from glob import glob
import numpy as np
import shutil
import xarray as xr
import netCDF4

dwd_routine='/h/ksddata06/routine-ksd/dwd1/' # directory of routinely downloaded dwd files maped to icon grid
files=np.sort(glob(dwd_routine+'*.nc'))[1:] # fist is in complete
ncdwd=xr.open_mfdataset(files)

outdir='/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd/'
timerange=files[0][-11:-3]+'_'+files[-1][-11:-3]
ncdwd.to_netcdf(outdir+'ICONglob_'+timerange+'.nc')


## merge individual (u,v) files
#outnames=[]
#for var in 'u','v':
#	files=np.sort(glob(datadir+'*'+var.upper()+'_10M*nc'))
#	mon0=files[0]
#	mon0=mon0[mon0.index('level_')+6:mon0.index('_'+var.upper()+'_')-6]
#	mon1=files[-1]
#	mon1=mon1[mon1.index('level_')+6:mon1.index('_'+var.upper()+'_')-6]
#	year=files[0]
#	year=mon0[:4]
#	mon0=mon0[4:6]
#	mon1=mon1[4:6]
#	fnameout=''.join((var+'DWDwnd_',year,mon0,'to',mon1,'.nc'))
#	outnames.append(fnameout)
#
#	shutil.copy(files[0], './'+fnameout)
#	files=files[1:]
#	nca=Dataset(fnameout, "a")
#	varname='10'+var
#	print('merging ' +var+' files')
#	for file in files:
#			print('appending '+ file)
#			ncb=Dataset(file, "r")
#			lena=len(nca['time'])
#			lenb=len(ncb['time'])
#
#			nca['time'][lena:lena+lenb]=ncb['time'][:] # zero start
#			append_data=ncb[varname][:]
#			growdata=nca[varname]
#			growdata[lena:lena+lenb,:]=append_data
#			ncb.close()
#	nca.sync()
#	nca.close()
#
## merge u and v files	
#ncu=Dataset(outnames[0],"a")
#ncv=Dataset(outnames[1])
#
#
#
#outname='u'+outnames[1]
#
##src=netCDF4.Dataset(outnames[0]) 
##srcv=netCDF4.Dataset(outnames[1]) 
##dst=netCDF4.Dataset(outname, "w")
##ncb=Dataset('/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd/ERA5wnd_2017-01.nc')
#
#with netCDF4.Dataset(outnames[0]) as src, netCDF4.Dataset(outnames[1]) as srcv ,netCDF4.Dataset(outname, "w") as dst:
#	# copy global attributes all at once via dictionary
#	dst.setncatts(src.__dict__)
#	# copy dimensions
#	for name, dimension in src.dimensions.items():
#		dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
#    # copy all file data except for the excluded
#	for name, variable in src.variables.items():
#		x = dst.createVariable(name.replace('10u','u10'), variable.datatype, variable.dimensions,fill_value=9.969209968386869e+36)
#		dst[name.replace('10u','u10')][:] = src[name][:]
#		# copy variable attributes all at once via dictionary
#		dst[name.replace('10u','u10')].setncatts(src[name].__dict__)
#	variable=srcv['10v']	
#	x = dst.createVariable('v10', variable.datatype, variable.dimensions,fill_value=9.969209968386869e+36)
#	dst['v10'][:] = srcv['10v'][:]
#	# copy variable attributes all at once via dictionary
#	dst['v10'].setncatts(srcv['10v'].__dict__)
#		
#
#for infile in *uvDWD*
#do
#        outfile=${infile}[-1]_proj.nc
#done
#