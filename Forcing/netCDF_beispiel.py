import numpy as np
from netCDF4 import Dataset
import datetime as dt


#m=np.loadtxt # would be the analogue of load in matlab

ncfile='example.nc'
nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

# example Data
lon=np.arange(-6,6,0.1)
lat=np.arange(-4,4,0.1)
time=np.arange(0,30)
lon2d,lat2d=np.meshgrid(lon,lat)
lon3d,lat3d,time3d=np.meshgrid(lon,lat,time)
#lon3d.shape
# 80, 120, 30    := nlat, lon, ntime
 
basedate=(2012,6,1,0,0,0)
variable=np.sin(time3d)+np.cos(lat3d)+np.sin(lon3d)



# meta infow
hstr = dt.datetime.strftime(dt.datetime.now(),'%a %b %d %H:%M:%S %Y')+': created by some script\n'
nc.setncattr('history',(hstr))

# dimensions
nc.createDimension('time',None)  # nlimiteted
tv = nc.createVariable('time','f8',('time'))
tv.long_name = 'Time'
tv.standard_name = 'time'
tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
tv.base_date =  list(basedate)
tv[:]=time

nc.createDimension('lon',len(lon))  # 
lv = nc.createVariable('lon','f4',('lon'))
lv.long_name = 'Longitude'
lv.standard_name = 'longitude'
lv.units = 'degrees_east'
lv[:] = lon

nc.createDimension('lat',len(lat))  # 
lv = nc.createVariable('lat','f4',('lat'))
lv.long_name = 'Latitude'
lv.standard_name = 'latitude'
lv.units = 'degrees_north'
lv[:] = lat

# copy wind speeds
vv = nc.createVariable('some varible','f4',('time','lat','lon'))
vv.units = 'm/s'
vv.standard_name = 'eastward_wind'
vv.coordinates = 'lat lon'
for ti in range(len(time)):
	vv[ti,:]=variable[:,:,ti]	# set variable with timea as first dimension
nc.sync()