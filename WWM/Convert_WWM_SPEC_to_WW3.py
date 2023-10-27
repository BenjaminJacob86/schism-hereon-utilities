# WWM WW3 convert


import xarray as xr


dsww3=xr.open_dataset('ww3_2013_11_12_spec.nc')


dswwm=xr.open_dataset('wwm_station_0001.nc')


# i dont have frequency 


Dimensions:       (time: 1465, station: 957, string40: 40, frequency: 36,
                   direction: 24)
Coordinates:
  * time          (time) datetime64[ns] 2013-11-01 ... 2014-01-01
  * station       (station) float64 1.0 2.0 3.0 4.0 ... 954.0 955.0 956.0 957.0
  * string40      (string40) float64 nan nan nan nan nan ... nan nan nan nan nan
  * frequency     (frequency) float32 0.04 0.044 0.0484 ... 0.929 1.022 1.124
  * direction     (direction) float32 90.0 75.0 60.0 45.0 ... 135.0 120.0 105.0
Data variables:
    station_name  (station, string40) |S1 ...
    longitude     (time, station) float32 ...
    latitude      (time, station) float32 ...
    frequency1    (frequency) float32 ...
    frequency2    (frequency) float32 ...
    efth          (time, station, frequency, direction) float32 ...
    dpt           (time, station) float32 ...
    wnd           (time, station) float32 ...
    wnddir        (time, station) float32 ...
    cur           (time, station) float32 ...
    curdir        (time, station) float32 ...
Attributes: (12/19)
    product_name:              ww3.201311_spec.nc
    area:                      NBSext
    data_type:                 OCO spectra 2D
    format_version:            1.1
    southernmost_latitude:     n/a
    northernmost_latitude:     n/a
    ...                        ...
    start_date:                2013-11-01 00:00:00
    stop_date:                 2013-12-01 00:00:00
    field_type:                hourly
    history:                   Thu Apr 13 10:08:40 2023: ncrcat ww3.201311_sp...
    NCO:                       netCDF Operators version 4.7.5 (Homepage = htt...
    nco_openmp_thread_number:  1

	
time=dswwm.ocean_time.rename({'time'})

dswwm['ACOUT_2D']

ds=xr.open_dataset('../wwm_hist_0001.nc')

# HS: hs  significant wave height
#'DM': dir  NETCDF var. name for the mean wave direction 
# 'DSPR' : spr? #
# 1/TP ''fp
#  1/TM=2      t02



#sea surface wave directional variance spectral density

#wo-dimensional wave variance spectrum 

#2-dimensional spectrum


ds1=xr.open_dataset('../outputs/out2d_1.nc')
ds=xr.open_dataset('../wwm_hist_0001.nc')

ds['ocean_time']=ds['ocean_time'][0].values +np.timedelta64(3600,'s')*np.arange(24)



# 'sigWaveHeight'
timeout=xr.DataArray(ds['ocean_time'][0].values +np.timedelta64(3600,'s')*np.arange(24),dims=('ocean_time'))
hsout=xr.DataArray(ds1['sigWaveHeight'].values,dims=('ocean_time','mnp'))
hsout['mnp']=ds['mnp']
hsout=hsout.rename('hs')
#timeout.to_dataset(name='ocean_time')
timeout=timeout.rename('ocean_time')

# 'peakPeriod', -> frequency
FRHIGH         = 0.6345
vals=1/ds1['peakPeriod'].values
vals[np.isinf(vals)]=FRHIGH
fpout=xr.DataArray(vals,dims=('ocean_time','mnp')) #sec
fpout['mnp']=ds['mnp']#
fpout=fpout.rename('fp')
# can have inf

#'meanDirSpreading'
#'spr' 
sprout=xr.DataArray(ds1['meanDirSpreading'].values,dims=('ocean_time','mnp'))
sprout['mnp']=ds['mnp']
sprout=sprout.rename('spr')

# 'zeroDowncrossPeriod' -> frequebcy
vals=1/ds1['zeroDowncrossPeriod'].values
vals[np.isinf(vals)]=FRHIGH
t02out=xr.DataArray(vals,dims=('ocean_time','mnp')) #sec
t02out['mnp']=ds['mnp']#
t02out=t02out.rename('t02')

dsout=timeout.to_dataset(name='time')
dsout['hs']=hsout
dsout['fp']=fpout
dsout['spr']=sprout
dsout['t02']=t02out

dsout.to_netcdf('wwm_bd_test.nc')




#dsout=ds[['HS','DSPR']]

# output in  history is not hourly check 

# use schout outputs




# rename and convert perods to frequency
hs=ds['HS']
spr=ds['DSPR']
t02=xr.zeros_like(ds['TM02'])
t02.attrs['full-name']='Zero-crossing wave frequency'
t02=t02.rename('t02')


import x
 
 




# modifiy hist output



 NCDF_HS_NAME   = 'hs'               ! NETCDF var. name for the significant wave height (normally it is just 'hs')
 NCDF_DIR_NAME  = 'dir'              ! NETCDF var. name for the mean wave direction (normally it is just 'dir')
 NCDF_SPR_NAME  = 'spr'              ! NETCDF var. name for the mean directional spreading (normally it is just 'spr')
 NCDF_FP_NAME   = 'fp'               ! NETCDF var. name for the peak freq. (normally it is just 'fp')
 NCDF_F02_NAME  = 't02'              ! NETCDF var. name for the zero down crossing freq. (normally it is just 't02')
