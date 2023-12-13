from netCDF4 import Dataset, date2num
#import netcdftime
#from cftime import utime
import cftime
import numpy as np
from datetime import datetime as dt
import os

# use within a year
year0=2023
startyear=year0 #07
endyear=2023
years=range(startyear,endyear+1)

from glob import glob



erafiles=np.sort(glob('ERA5*'))
file0=erafiles[0]
month0_year0=int(file0[file0.index('.')-2:-3])
nmonth=len(erafiles)

base_dir='./'
evap_dir=''
lat_sens_dir=''
# daten jetzt auf mistral: /work/gg0028/g260099/ERA_INTERIM/RAW_ERA5_CDS
target_dir='./era5/'

do_prec=1       # do precipiation files
do_non_prec=1   # do other files


load_spfh=False # True load specific humidity from era forcases. Else compute from 2m dew point 
spfh_dir=''


if not os.path.exists(target_dir):
	os.mkdir(target_dir)

netflux=False	
#therm='str' #'strd'  # str: longwave strd net longave downwelling (= Ldown - Lup)  strd: surface thermal radiation downward  
therm='strd'
#! net_sfc_flux_d = - (sen_flux + lat_flux + (longwave_u - longwave_d))


# set start of simulation as reference time
basedate=(year0,month0_year0,1,0,0,0) # fixed start
basedate=[]   #  if basedate =[] use beginnig of each month this is how Joseph does it

# 30.11.2018 with forecaste and Ana in one file process of averaging of forecast becomse obsolete
# as timesteps are centralised now

# reverse the latitute order and associated matrices to have ascending order for schism forcing


def TEMP2WVP(tkelvin):
	#! --- compute the saturation water vapor pressure from
	#! --- temperature in kelvin  (Lowe,1977; JAM,16,100,1977)
	a=[6984.505294e0,-188.9039310e0,2.133357675e0,-1.288580973e-2,4.393587233e-5,-8.023923082e-8,6.136820929e-11]
	WVP = a[0]+tkelvin*(a[1]+tkelvin*(a[2]+tkelvin*(a[3]+tkelvin*(a[4]+tkelvin*(a[5]+tkelvin*a[6])))))
	return WVP


def create_grid(nc,inc,basedate):
	nc.createDimension('time',None)
	tv = nc.createVariable('time','f8',('time'))
	tv.long_name = 'Time'
	tv.standard_name = 'time'
	tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
	tv.base_date =  list(basedate)
	#ut = netcdftime.utime(tv.units)
	#ut = utime(tv.units)
	#
	incv = inc.variables
	#
	# copy some global attributes
	#for attr in ['experiment_id','references']:
	#  nc.setncattr(attr,inc.getncattr(attr))
	#
	hstr = dt.strftime(dt.now(),'%a %b %d %H:%M:%S %Y')+': create_schism_sflux.py\n'
	#nc.setncattr('history',unicode(hstr+inc.getncattr('history')))
	nc.setncattr('history',(hstr+inc.getncattr('history')))
	#
	# write time
	#iut = netcdftime.utime(incv['time'].units)
	#iut = utime(incv['time'].units)
	#iut = cftime(incv['time'].units)
	tv[0:len(inc.dimensions['time'])]=cftime.date2num(cftime.num2date(incv['time'][:],units=incv['time'].units),units=tv.units)
	#tv[0:len(inc.dimensions['time'])] = cftime.date2num(iut.num2date(incv['time'][:]))
	# write grid
	nc.createDimension('nx_grid',len(inc.dimensions['longitude']))
	nc.createDimension('ny_grid',len(inc.dimensions['latitude']))
	lon = incv['longitude'][:]
	lat = incv['latitude'][::-1]
	gridlon,gridlat = np.meshgrid(lon,lat)
	#
	lv = nc.createVariable('lon','f4',('ny_grid','nx_grid'))
	lv.long_name = 'Longitude'
	lv.standard_name = 'longitude'
	lv.units = 'degrees_east'
	lv[:] = gridlon
	#
	lv = nc.createVariable('lat','f4',('ny_grid','nx_grid'))
	lv.long_name = 'Latitude'
	lv.standard_name = 'latitude'
	lv.units = 'degrees_north'
	lv[:] = gridlat
	#
	nc.sync()


	
rho_water=1000 # kg/m3	
	
for year in years:

	
	if year==year0:
		months=range(month0_year0,12+1)
	else:
		 months=range(1,12+1) 	

	for month in months(:nmonth):
		print('month:' +str(month))
		# open input file
		inncfile='%s/ERA5CDS%04d%02d.nc'%(base_dir,year,month)
		incan = Dataset(inncfile)
		incv=incan.variables
		t=incv['time'][:]
		timestep=np.double(np.diff(t[0:2]))*3600.0
	
		if do_non_prec:

			# create output file
			ncfile='%s/cDII.air.%04d_%02d.nc'%(target_dir,year,month)
			nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')


			# create grid from input file
			if basedate==[]:
				create_grid(nc,incan,(year,month,1,0,0,0))
			else:
				create_grid(nc,incan,basedate)


			# copy wind speeds
			vv = nc.createVariable('uwind','f4',('time','ny_grid','nx_grid'))
			vv.units = 'm/s'
			vv.standard_name = 'eastward_wind'
			vv.coordinates = 'lat lon'

			vv[:] = incv['u10'][:,::-1,:].squeeze()

			vv = nc.createVariable('vwind','f4',('time','ny_grid','nx_grid'))
			vv.units = 'm/s'
			vv.standard_name = 'northward_wind'
			vv.coordinates = 'lat lon'
			vv[:] = incv['v10'][:,::-1,:].squeeze()
			#inc.close()

			# open temp file, copy temp  - same inncfile
			#inncfile='%s/cDII.00.T_2M.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
			#inc = Dataset(inncfile)
			#incv=inc.variables

			vv = nc.createVariable('stmp','f4',('time','ny_grid','nx_grid'))
			vv.units = 'K'
			vv.standard_name = 'air_temperature'
			vv.coordinates = 'lat lon'
			#vv[:] = incv['sst'][:,::-1,:].squeeze()
			#vv[:] = incv['t2m'][:].squeeze()
			sst = incv['sst'][:,::-1,:].squeeze()
			data = incv['t2m'][:,::-1,:].squeeze()
			data.mask=sst.mask
			vv[:]=data


			#inc.close()

			# open pressure file, copy pressure
			#inncfile='%s/cDII.00.PMSL.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
			#inc = Dataset(inncfile)
			#incv=inc.variables

			vv = nc.createVariable('prmsl','f4',('time','ny_grid','nx_grid'))
			vv.units = 'Pa'
			vv.standard_name = 'air_pressure_at_sea_level'
			vv.coordinates = 'lat lon'
			vv[:] = incv['msl'][:,::-1,:].squeeze()

			#inc.close()



			#inncfile='%s/cDII.00.QV_2M.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
			#inc = Dataset(inncfile)
			#incv=inc.variables

			vv = nc.createVariable('spfh','f4',('time','ny_grid','nx_grid'))
			vv.units = '1'
			vv.standard_name = 'specific_humidity'
			vv.coordinates = 'lat lon'

			if load_spfh:
					# open input file
					spfh_file='%s/ERA5CDS%04d%02d.nc'%(spfh_dir,year,month)
					incspfh = Dataset(spfh_file)
					spfhv=incspfh.variables
					tspfh=incv['time'][:] # time is compatible ensure via download
					if (tspfh==t).min()==False:
						print('specific humidity timevectors do not correspondent to analysis one')
					#timestepspfh=np.double(np.diff(tspfh[0:2]))*3600.0
					incv['d2m'][:,::-1,:]
					vv[:]=spfhv['q'][:,0,::-1,:].squeeze()
			
			else: # compute spfh
				#vv[:] = incv['QV_2M'][:].squeeze()
				# derive specific humidity from dew point temperature
				SWV=TEMP2WVP(incv['d2m'][:,::-1,:]).squeeze()
				vv[:]=(0.622*SWV/(1013.-((1.-0.622)*SWV)))
				#vv[:]=(0.622*TEMP2WVP(incv['d2m'][:])/(1013.-((1.-0.622)*TEMP2WVP(incv['d2m'][:]))))

			#inc.close()
	
			# close air file
			nc.close()

			# ===============================================
			# write rad file
			ncfile='%s/cDII.rad.%04d_%02d.nc'%(target_dir,year,month)
			nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

			# open short wave rad input  - forecast files are integrals and time shifted 1/2 dt 
			inncfile='%s/ERA5CDS%04d%02d.nc'%(base_dir,year,month)
			#inncfileP='%s/ERA5CDS_%04d%02d.nc'%(base_dir,year-(month==1),(month-1)+12*(month==1)) #previous file
			#inncfileN='%s/ERA5CDS_%04d%02d.nc'%(base_dir,year+(month==12),(month+1)%13+int(month==12)) #following file
			#create grid from input file
			# use input file from wind file (timestep error in radiation)
			# create grid from input file
			if basedate==[]:
				create_grid(nc,incan,(year,month,1,0,0,0))
			else:
				create_grid(nc,incan,basedate)
			#create_grid(nc,incan,basedate)
			#inc.close()


			# files are diference formed need to devide by time step	
			inc = Dataset(inncfile)
			incv=inc.variables

			#incp = Dataset(inncfileP)
			#incpv=incp.variables
	
			#incn = Dataset(inncfileN)
			#incnv=incn.variables




			t=incv['time'][:]
			timestep=np.double(np.diff(t[0:2]))*3600.0



			# open long wave rad input
			#inncfile='%s/%04d%02dfc.nc'%(base_dir,year,month)
			#    inc = Dataset(inncfile)
			#    incv=inc.variables

			vv = nc.createVariable('dlwrf','f4',('time','ny_grid','nx_grid'))
			vv.units = 'W/m^2'
			vv.standard_name = 'surface_downwelling_longwave_flux_in_air'
			vv.coordinates = 'lat lon'


			#tempP=incpv[therm][-1,::-1,:] # surface net thermal radiation
			#temp=incv[therm][:,::-1,:]
			#tempN=incnv[therm][0,::-1,:]
			# average
			#temp[1:-1,:,:]=temp[0:-2,:,:]+temp[2:,:,:]
			#temp[0,:,:]=tempP+temp[1,:,:] 	
			#temp[-1,:,:]=temp[-1,:,:]+tempN 	
			#vv[:] = temp[:].squeeze()/2/timestep # downward short wave radiation net?
			
			if netflux==True:
	
				# open latent_sensible file
				lat_sens_ncfile='%s/ERA5CDS_lsheat_%04d%02d.nc'%(lat_sens_dir,year,month)
				inc_lat_sens = Dataset(lat_sens_ncfile)
				inc_lat_sens_v=inc_lat_sens.variables
				
				
				# net_sfc_flux_d = - (sen_flux + lat_flux + (longwave_u - longwave_d)) # units all check out to [J/m**2] -> /s =[W/m2]
				# The ECMWF convention for vertical fluxes is positive downwards -> add values 
				vv[:] = ( incv['str'][:,::-1,:].squeeze() +inc_lat_sens_v['slhf'][:,::-1,:].squeeze() + inc_lat_sens_v['sshf'][:,::-1,:].squeeze() ) /timestep # downward short wave radiation net?

			else:	
				vv[:] = incv[therm][:,::-1,:].squeeze()/timestep # downward short wave radiation net?

			#    inc.close()




			# copy short wave flux
			vv = nc.createVariable('dswrf','f4',('time','ny_grid','nx_grid'))
			vv.units = 'W/m^2'
			vv.standard_name = 'surface_downwelling_shortwave_flux_in_air'
			vv.coordinates = 'lat lon'

			#tempP=incpv['ssrd'][-1,::-1,:] # [J/m^2]
			#temp=incv['ssrd'][:,::-1,:]
			#t1empN=incnv['ssrd'][0,::-1,:]
			# average
			#temp[1:-1,:,:]=temp[0:-2,:,:]+temp[2:,:,:]
			#temp[0,:,:]=tempP+temp[1,:,:] 	
			#temp[-1,:,:]=temp[-1,:,:]+tempN 	
			#vv[:] = temp[:].squeeze()/2/timestep # downward short wave radiation net?
			vv[:] = incv['ssrd'][:,::-1,:].squeeze()/timestep # downward short wave radiation net?
			#vv[:] = incv['ASOB_S'][:].squeeze()

			inc.close()
			nc.close()
			
		# ===============================================
		if do_prec:
			# write precipitation file
			ncfile='%s/cDII.prec.%04d_%02d.nc'%(target_dir,year,month)
			nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

			

			# open short wave rad input
			#inncfile='%s/cDII.00.TOT_PREC.%04d_%02d.lonlat_gridCE.nc'%(base_dir,year,month)
			inc = Dataset(inncfile)
			incv=inc.variables

			# create grid from input file
			#create_grid(nc,incan,basedate)
			# create grid from input file
			if basedate==[]:
				create_grid(nc,incan,(year,month,1,0,0,0))
			else:
				create_grid(nc,incan,basedate)
			incan.close()
			# copy precipitation flux
			#t = nc.variables['time'][:2]
			#timestep = (t[1]-t[0])*86400.

			#	tempP=incpv['sf'][-1,::-1,:]+incpv['tp'][-1,::-1,:]  # [m] snow fall + total precipitation as water equivalents
			#	temp=incv['sf'][:,::-1,:]+incv['tp'][:,::-1,:]
			#	tempN=incpv['sf'][0,::-1,:]+incpv['tp'][0,::-1,:]

			# wo snow fall
			#tempP=incpv['tp'][-1,::-1,:]  # [m] snow fall + total precipitation as water equivalents
			#temp=incv['tp'][:,::-1,:]
			#tempN=incnv['tp'][0,::-1,:]

			# average
			#temp[1:-1,:,:]=temp[0:-2,:,:]+temp[2:,:,:] # total 
			#temp[0,:,:]=tempP+temp[1,:,:] 	
			#temp[-1,:,:]=temp[-1,:,:]+tempN  # next step 	

			#temp=temp/2/timestep #convert to rate
			# do i have to convert ?  m to kg/m2 or are values 
			#rho_water=1000 # kg/m3

			#  build difference ?
			#temp=temp*rho_water


			print('  timestep: %0.2f s'%timestep)
			vv = nc.createVariable('prate','f4',('time','ny_grid','nx_grid'))
			vv.units = 'kg/m^2/s'
			vv.standard_name = 'precipitation_flux'
			vv.coordinates = 'lat lon'
			#vv[:] = incv['TOT_PREC'][:].squeeze()/timestep
			#vv[:] = (temp[:]/2/timestep).squeeze() # downward long wave radiation net?
			##vv[:] = (temp[:]).squeeze() 
			if netflux==True:
			
				# open evap input file
				evpncfile='%s/ERA5CDS_evap_%04d%02d.nc'%(evap_dir,year,month)
				incevp = Dataset(evpncfile)
				incevpv=incevp.variables

				# in Era fluxes postive Downwards (negative values are evaporation and positive condensation)
				vv[:] = (incv['tp'][:,::-1,:].squeeze() +  (incevpv['e'][:,::-1,:].squeeze() + incevpv['es'][:,::-1,:].squeeze()  ) ) *rho_water/timestep # downward short wave radiation net?		
				# sign correct ?
			else:
				vv[:] = (incv['tp'][:,::-1,:].squeeze()  ) *rho_water/timestep # downward short wave radiation net?		
			inc.close()
			nc.close()

bd=[]			
for o,l in zip(s.bdy_segments,s.land_segments):
	bd+=l+o
bd=np.asarray(bd)-1
lon,lat=np.asarray(s.lon),np.asarray(s.lat)
plt.plot(lon[bd],lat[bd])




import os

""" link date carrayin atmospheric file name fron dwd_script.m to alphanumerical format expexcted by schism  """

indir=target_dir #'/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/sflux_new/out_dwd_script/'
linkdir=os.getcwd()+'/' #indir #'/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/sflux_new/'


air=list(np.sort(glob(indir+'*air*')))
rad=list(np.sort(glob(indir+'*rad*')))
prc=list(np.sort(glob(indir+'*prec*')))
n=0
#for i in range(len(air)): 
for i,files in enumerate(zip(air,rad,prc)):
	#fair,frad,fprc=files
	fair,frad,fprc=air[i],rad[i],prc[i]
	n=i+1
	#cmd='ln -s {:s} {:s}sflux_air_1.{:04d}.nc'.format(fair,linkdir,n)
	os.system('ln -s {:s} {:s}sflux_air_1.{:04d}.nc'.format(fair,linkdir,n))
	os.system('ln -s {:s} {:s}sflux_rad_1.{:04d}.nc'.format(frad,linkdir,n))
	os.system('ln -s {:s} {:s}sflux_prc_1.{:04d}.nc'.format(fprc,linkdir,n))
print('linked ' +str(n) +' files')




