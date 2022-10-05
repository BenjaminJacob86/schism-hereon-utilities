from netCDF4 import Dataset
import netcdftime
import numpy as np
from datetime import datetime as dt


years=range(2000,2015)
years=range(2013,2015)
months=range(1,13)
base_dir='original'
target_dir='schism-sflux'

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
    ncfile='%s/cDII.air.%04d_%02d.nc'%(target_dir,year,month)
    nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

    # open input file
    inncfile='%s/cDII.00.KS.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
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
    #inc.close()

    # open temp file, copy temp
    #inncfile='%s/cDII.00.T_2M.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
    #inc = Dataset(inncfile)
    #incv=inc.variables
    
    vv = nc.createVariable('stmp','f4',('time','ny_grid','nx_grid'))
    vv.units = 'K'
    vv.standard_name = 'air_temperature'
    vv.coordinates = 'lat lon'
    vv[:] = incv['T_2M'][:].squeeze()
    
    #inc.close()

    # open pressure file, copy pressure
    #inncfile='%s/cDII.00.PMSL.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
    #inc = Dataset(inncfile)
    #incv=inc.variables
    
    vv = nc.createVariable('prmsl','f4',('time','ny_grid','nx_grid'))
    vv.units = 'Pa'
    vv.standard_name = 'air_pressure_at_sea_level'
    vv.coordinates = 'lat lon'
    vv[:] = incv['PMSL'][:].squeeze()
    
    #inc.close()

    # open spec. hum file, copy data
    #inncfile='%s/cDII.00.QV_2M.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
    #inc = Dataset(inncfile)
    #incv=inc.variables
    
    vv = nc.createVariable('spfh','f4',('time','ny_grid','nx_grid'))
    vv.units = '1'
    vv.standard_name = 'specific_humidity'
    vv.coordinates = 'lat lon'
    vv[:] = incv['QV_2M'][:].squeeze()
    
    #inc.close()

    # close air file
    nc.close()

    # ===============================================
    # write rad file
    ncfile='%s/cDII.rad.%04d_%02d.nc'%(target_dir,year,month)
    nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

    # open short wave rad input
    #inncfile='%s/cDII.00.ASOB_S.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
    #inc = Dataset(inncfile)
    #incv=inc.variables

    # create grid from input file
    create_grid(nc,inc)
    
    # copy short wave flux
    vv = nc.createVariable('dswrf','f4',('time','ny_grid','nx_grid'))
    vv.units = 'W/m^2'
    vv.standard_name = 'surface_downwelling_shortwave_flux_in_air'
    vv.coordinates = 'lat lon'
    vv[:] = incv['ASOB_S'][:].squeeze()

    #inc.close()

    # open long wave rad input
    #inncfile='%s/cDII.00.ALWD_S.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
    #inc = Dataset(inncfile)
    #incv=inc.variables

    vv = nc.createVariable('dlwrf','f4',('time','ny_grid','nx_grid'))
    vv.units = 'W/m^2'
    vv.standard_name = 'surface_downwelling_longwave_flux_in_air'
    vv.coordinates = 'lat lon'
    vv[:] = incv['ALWD_S'][:].squeeze()

    #inc.close()
    nc.close()

    # ===============================================
    # write precipitation file
    ncfile='%s/cDII.prec.%04d_%02d.nc'%(target_dir,year,month)
    nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

    # open short wave rad input
    #inncfile='%s/cDII.00.TOT_PREC.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
    #inc = Dataset(inncfile)
    #incv=inc.variables

    # create grid from input file
    create_grid(nc,inc)
    
    # copy precipitation flux
    t = nc.variables['time'][:2]
    timestep = (t[1]-t[0])*86400.
    print('  timestep: %0.2f s'%timestep)

    vv = nc.createVariable('prate','f4',('time','ny_grid','nx_grid'))
    vv.units = 'kg/m^2/s'
    vv.standard_name = 'precipitation_flux'
    vv.coordinates = 'lat lon'
    vv[:] = incv['TOT_PREC'][:].squeeze()/timestep

    inc.close()
    nc.close()

