from netCDF4 import Dataset
import netcdftime
import numpy as np
from datetime import datetime as dt
import sys, os
import shutil

####### 
# Script to create atmospheric forcing files for SCHISM from coastDAT3 data.
# cD3 data needs to be remapped from COSMOs rotated grid to regular lonlat grid in base_dir
# For downwelling SW radiation cD3: ASWDIFD_S + ASWDIR_S --> full downwelling SW rad,
# to be used instead of ASOB_S (net downwelling SW rad) as SCHISM considers albedo and SW then
# then has more tuning potential
#
# Author: jan.kossack@hereon.de, based on script by Richard Hofmeister
# Date: 15.04.2021

# TODO:
########  

#### USER INPUT ####
years=range(2011,2017)
months=range(1,13)
base_dir='/work/gg0877/KST/cD3/cD3_011_ERAi'
target_dir='/work/gg0877/KST/cD3/schism-sflux_cD3_011_ERAi'
#####

# Remove complete (!) out dir
#if os.path.exists(target_dir) and os.path.isdir(target_dir):
#  shutil.rmtree(target_dir)
#os.makedirs(target_dir)

def create_grid(nc,inc):


    nc.createDimension('time',None)
    tv = nc.createVariable('time','f8',('time'))
    tv.long_name = 'Time'
    tv.standard_name = 'time'
    tv.units = 'days since 1948-01-01 00:00:00'
    tv.base_date = [1948,1,1,0]
    ut = netcdftime.utime(tv.units)

    incv = inc.variables
    
    # copy some global attributes
    for attr in ['experiment_id','references']:
      nc.setncattr(attr,inc.getncattr(attr))

    hstr = dt.strftime(dt.now(),'%a %b %d %H:%M:%S %Y')+': create_schism_sflux.py\n'
    nc.setncattr('history',unicode(hstr+inc.getncattr('history')))

    # write time
    iut = netcdftime.utime(incv['time'].units)
    tv[0:len(inc.dimensions['time'])] = ut.date2num(iut.num2date(incv['time'][:]))
    # write grid
    nc.createDimension('nx_grid',len(inc.dimensions['lon']))
    nc.createDimension('ny_grid',len(inc.dimensions['lat']))
    lon = incv['lon'][:]
    lat = incv['lat'][:]
    gridlon,gridlat = np.meshgrid(lon,lat)

    lv = nc.createVariable('lon','f4',('ny_grid','nx_grid'))
    lv.long_name = 'Longitude'
    lv.standard_name = 'longitude'
    lv.units = 'degrees_east'
    lv[:] = gridlon

    lv = nc.createVariable('lat','f4',('ny_grid','nx_grid'))
    lv.long_name = 'Latitude'
    lv.standard_name = 'latitude'
    lv.units = 'degrees_north'
    lv[:] = gridlat

    nc.sync()


for year in years:
  for month in months:
    # create output file
    ncfile='%s/cD3.air.%04d_%02d.nc'%(target_dir,year,month)
    nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

    # open wind input file
    inncfile='%s/UVlat_10M/cD3_011_ERAi.UVlat_10M.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc = Dataset(inncfile)
    incv=inc.variables

    # create grid from input file
    create_grid(nc,inc)
    
    # copy wind speeds
    vv = nc.createVariable('uwind','f4',('time','ny_grid','nx_grid'))
    vv.units = 'm/s'
    vv.standard_name = 'eastward_wind'
    vv.coordinates = 'lat lon'
    vv[:] = incv['U_10M'][:].squeeze()

    vv = nc.createVariable('vwind','f4',('time','ny_grid','nx_grid'))
    vv.units = 'm/s'
    vv.standard_name = 'northward_wind'
    vv.coordinates = 'lat lon'
    vv[:] = incv['V_10M'][:].squeeze()
    inc.close()

    # open temp file, copy temp
    inncfile='%s/T_2M/cD3_011_ERAi.T_2M.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc = Dataset(inncfile)
    incv=inc.variables
    
    vv = nc.createVariable('stmp','f4',('time','ny_grid','nx_grid'))
    vv.units = 'K'
    vv.standard_name = 'air_temperature'
    vv.coordinates = 'lat lon'
    vv[:] = incv['T_2M'][:].squeeze()
    inc.close()

    # open pressure file, copy pressure
    inncfile='%s/PMSL/cD3_011_ERAi.PMSL.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc = Dataset(inncfile)
    incv=inc.variables
    
    vv = nc.createVariable('prmsl','f4',('time','ny_grid','nx_grid'))
    vv.units = 'Pa'
    vv.standard_name = 'air_pressure_at_sea_level'
    vv.coordinates = 'lat lon'
    vv[:] = incv['PMSL'][:].squeeze()
    inc.close()

    # open spec. hum file, copy data
    inncfile='%s/QV_2M/cD3_011_ERAi.QV_2M.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc = Dataset(inncfile)
    incv=inc.variables
    
    vv = nc.createVariable('spfh','f4',('time','ny_grid','nx_grid'))
    vv.units = '1'
    vv.standard_name = 'specific_humidity'
    vv.coordinates = 'lat lon'
    vv[:] = incv['QV_2M'][:].squeeze()
    inc.close()

    # close air file
    nc.close()

    # ===============================================
    # write rad file
    ncfile='%s/cD3.rad.%04d_%02d.nc'%(target_dir,year,month)
    nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

    # short wave rad input: sum of direct + diffuse downwelling SW
    inncfile1='%s/ASWDIR_S/cD3_011_ERAi.ASWDIR_S.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inncfile2='%s/ASWDIFD_S/cD3_011_ERAi.ASWDIFD_S.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc1 = Dataset(inncfile1)
    inc2 = Dataset(inncfile2)
    incv1=inc1.variables
    incv2=inc2.variables

    # create grid from input file
    create_grid(nc,inc1)
    
    # copy short wave flux
    vv = nc.createVariable('dswrf','f4',('time','ny_grid','nx_grid'))
    vv.units = 'W/m^2'
    vv.standard_name = 'surface_downwelling_shortwave_flux_in_air'
    vv.coordinates = 'lat lon'
    vv[:] = incv1['ASWDIR_S'][:].squeeze() + incv2['ASWDIFD_S'][:].squeeze()

    inc1.close()
    inc2.close()

    # open long wave rad input
    inncfile='%s/ALWD_S/cD3_011_ERAi.ALWD_S.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc = Dataset(inncfile)
    incv=inc.variables

    vv = nc.createVariable('dlwrf','f4',('time','ny_grid','nx_grid'))
    vv.units = 'W/m^2'
    vv.standard_name = 'surface_downwelling_longwave_flux_in_air'
    vv.coordinates = 'lat lon'
    vv[:] = incv['ALWD_S'][:].squeeze()

    inc.close()
    nc.close()

    # ===============================================
    # write precipitation file
    ncfile='%s/cD3.prec.%04d_%02d.nc'%(target_dir,year,month)
    nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

    # open short wave rad input
    inncfile='%s/TOT_PREC/cD3_011_ERAi.TOT_PREC.lonlat.%04d%02d.nc'%(base_dir,year,month)
    inc = Dataset(inncfile)
    incv=inc.variables

    # create grid from input file
    create_grid(nc,inc)
    
    # copy precipitation flux
    t = nc.variables['time'][:2]
    timestep = (t[1]-t[0])*86400.
    print('TOT_PREC  timestep: %0.2f s'%timestep)

    vv = nc.createVariable('prate','f4',('time','ny_grid','nx_grid'))
    vv.units = 'kg/m^2/s'
    vv.standard_name = 'precipitation_flux'
    vv.coordinates = 'lat lon'
    vv[:] = incv['TOT_PREC'][:].squeeze()/timestep

    inc.close()
    nc.close()

