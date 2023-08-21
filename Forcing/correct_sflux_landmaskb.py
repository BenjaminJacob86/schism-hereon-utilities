from netCDF4 import Dataset
import cftime
import numpy as np

import sys
import os
sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/home/g/g260114/git/schism-hereon-utilities/')
from schism import *
from datetime import datetime as dt
import matplotlib
matplotlib.use('agg')
#matplotlib.use('GtkAgg')
#matplotlib.use('GtkAgg')
from matplotlib import pyplot as plt
from scipy.spatial import KDTree


cwd=os.getcwd()



#########  S  E T  T  I N G S ##########################################
# schism dir
setupdir='/work/gg0028/g260114/SETUPS/Ghana/'
atmodir='/work/gg0028/g260114/SETUPS/Ghana/sflux/era5/'
target_dir='/work/gg0028/g260114/SETUPS/Ghana/sflux/era5_mask_reduce/'#'schism-sflux_noextra_snow'





year0=2022
endyear=2022
years=range(year0,endyear+1)
month0_year0=1

all_land2ocean=True


if not os.path.exists(target_dir):
        os.mkdir(target_dir)


# set start of simulation as reference time
basedate=(year0,month0_year0,1,0,0,0)



#########################################################################

if not  os.path.isdir(target_dir):
	os.mkdir(target_dir)


# load schism setup
os.chdir(setupdir)
europe=schism_setup()
os.chdir(cwd)

# removed adding snowfall to totl precip




def create_reduced_grid(nc,inc,basedate,lonrange,latrange):


	# get sub index range for input inc
	lon=inc['lon']
	lat=inc['lat']

	loninds=np.arange(np.sum(np.asarray(lon)[0,:]<lonrange[0])-1,np.sum(np.asarray(lon)[0,:]<=lonrange[1]) )
	latinds=np.arange(np.sum(np.asarray(lat)[:,0]<latrange[0])-1,np.sum(np.asarray(lat)[:,0]<=latrange[1]) )

	lon=lon[latinds,loninds]
	lat=lat[latinds,loninds]


	nc.createDimension('time',None)
	tv = nc.createVariable('time','f8',('time'))
	tv.long_name = 'Time'
	tv.standard_name = 'time'
	tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
	tv.base_date =  list(basedate)
	#ut = cftime.utime(tv.units)

	incv = inc.variables

	# copy some global attributes
	#for attr in ['experiment_id','references']:
	#  nc.setncattr(attr,inc.getncattr(attr))

	hstr = dt.strftime(dt.now(),'%a %b %d %H:%M:%S %Y')+': create_schism_sflux.py\n'
	nc.setncattr('history',(hstr+inc.getncattr('history')))

	# write time
	#iut = cftime.utime(incv['time'].units)
	#tv[0:len(inc.dimensions['time'])] = ut.date2num(iut.num2date(incv['time'][:]))
	tv[0:len(inc.dimensions['time'])]=cftime.date2num(cftime.num2date(incv['time'][:],units=incv['time'].units),units=tv.units)
	# write grid
	nc.createDimension('nx_grid',len(loninds))
	nc.createDimension('ny_grid',len(latinds))
	#lon = incv['longitude'][:]
	#lat = incv['latitude'][::-1]
	gridlon,gridlat = lon,lat

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





def create_grid(nc,inc,basedate):


	nc.createDimension('time',None)
	tv = nc.createVariable('time','f8',('time'))
	tv.long_name = 'Time'
	tv.standard_name = 'time'
	tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
	tv.base_date =  list(basedate)
	ut = cftime.utime(tv.units)

	incv = inc.variables

	# copy some global attributes
	#for attr in ['experiment_id','references']:
	#  nc.setncattr(attr,inc.getncattr(attr))

	hstr = dt.strftime(dt.now(),'%a %b %d %H:%M:%S %Y')+': create_schism_sflux.py\n'
	nc.setncattr('history',(hstr+inc.getncattr('history')))

	# write time
	iut = cftime.utime(incv['time'].units)
	tv[0:len(inc.dimensions['time'])] = ut.date2num(iut.num2date(incv['time'][:]))
	# write grid
	nc.createDimension('nx_grid',len(inc.dimensions['longitude']))
	nc.createDimension('ny_grid',len(inc.dimensions['latitude']))
	lon = incv['longitude'][:]
	lat = incv['latitude'][::-1]
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




def plotbound(s):
        for oseg,lseg in zip(s.bdy_segments,s.land_segments):
                plt.plot(np.asarray(s.lon)[np.asarray(oseg)-1],np.asarray(s.lat)[np.asarray(oseg)-1],'k')
                plt.plot(np.asarray(s.lon)[np.asarray(lseg)-1],np.asarray(s.lat)[np.asarray(lseg)-1],'r')




# get indices of forcing area just covering model area
# get land mask
# get indices of land_mask values to be replaced with next neighbour ocean values



# open input file
year=year0
month=month0_year0
inncfile='%s/cDII.air.%04d_%02d.nc'%(atmodir,year,month)
incan = Dataset(inncfile)
incv=incan.variables


####### reduce grid to model coverage ##############################	

# get rectangle arround model domain
lonrange=(np.min(europe.lon),np.max(europe.lon))
latrange=(np.min(europe.lat),np.max(europe.lat))	

# get sub index range for input inc
lon=incv['lon']
lat=incv['lat']

loninds=np.arange(np.sum(np.asarray(lon)[0,:]<lonrange[0])-1,np.sum(np.asarray(lon)[0,:]<=lonrange[1]) )
latinds=np.arange(np.sum(np.asarray(lat)[:,0]<latrange[0])-1,np.sum(np.asarray(lat)[:,0]<=latrange[1]) )

# lon lat of atmo model domain part necessary to force
lon=lon[latinds,loninds]
lat=lat[latinds,loninds]
###############################################################

# plot variable making clear land mask
data=incv['stmp']
pdata=data[0,latinds,loninds]
plt.pcolormesh(lon,lat,pdata)
ch=plt.colorbar()
plotbound(europe)
ch.set_label('stmp')
plt.savefig('landmask_check',dpi=600)
plt.close()
# scan along boundarie nodes 
# check within distance if land values
# set those landvalues to next neighbour ocean values



########  determined nodes along boundary #################################
#####  needing to be sorounded by ocean values ############################

# distance is lon lat resolution
dx=1.5*np.diff(lon[0,:2]) # time 1.5 to ensure domainboundary within selection
dy=1.5*np.diff(lat[:2,0])

land_nodes=np.asarray(europe.land_nodes)-1
bdy_nodes=np.asarray(europe.bdy_nodes)-1
model_lon=np.asarray(europe.lon)
model_lat=np.asarray(europe.lat)
#plt.plot(model_lon[bdy_nodes],model_lat[bdy_nodes],'b.')

check_nodes=np.concatenate((land_nodes,bdy_nodes)) # nodes that need ocean forcing
						   # to avoid land mask artefacts
	
# find atmo model (reduced) grid nodes which need to be ocean values (land mask correction)
check_grid=np.zeros(np.shape(lon))
for node in check_nodes:
	check_grid[
	 (model_lon[node]-dx <= lon)  & (lon <= model_lon[node]+dx) & \
	 (model_lat[node]-dy <= lat)  & (lat <= model_lat[node]+dy)]=1 

plt.pcolormesh(lon,lat,check_grid)
plotbound(europe)
#plt.show()
plt.savefig('grid_check_replace_ocean_values',dpi=600)

# atmo grid indices to replace with nn ocean values if land
#replace_dest=(check_grid==1) & (pdata==np.max(pdata))

if all_land2ocean:
	replace_dest= pdata.mask==True
else:
	replace_dest=(check_grid==1) & (pdata.mask==True)

#isocean=pdata<np.max(pdata)
isocean=pdata.mask==False

#plt.pcolormesh(lon,lat,replace_dest)

# finde next ocean point neighbours
# use indices not lon lat
# draw back distance are treated equally
# does not work for lon lat ?


# prepare indices for kde quick nn look up 
#(use array indices insetad long lat , but distances not euclidc)

x,y=np.mgrid[0:pdata.shape[0],0:pdata.shape[1]]

# transpose coordinates to  column vectors ??
xygood=np.array((x[isocean],y[isocean])).T	# pool of ocean points
xyreplace=np.array((x[replace_dest],y[replace_dest])).T # nodes at boundary to be exchahge with nearest ocean nodes
							
 
tree=KDTree(xygood) # lookup tree of good(ocean) points   #distances,indices=tree.query(xyreplace)

# test it
isocean2=isocean.copy()
isocean2[replace_dest]=isocean[isocean][tree.query(xyreplace)[1]] # exchange coomand

nn_ocean_index=tree.query(xyreplace)[1]

plt.subplot(2,2,1)
plt.pcolormesh(isocean*1.0)
plt.colorbar()
plt.title('initial ocean area')
plt.subplot(2,2,2)
plt.pcolormesh(isocean2*1.0)
plt.colorbar()
plt.title('reser ocean area')
plt.subplot(2,2,3)
plt.pcolormesh((isocean2)*1.0-(isocean)*1.0)
plt.colorbar()
plt.title('added ocean area')
plt.savefig('nn_copied_ocean_area',dpi=600)


filetypes=('air','rad','prec')

for year in years:

	if year==year0:
		months=range(month0_year0,12+1)
	else:
		months=range(1,12+1)

	#month=6

	for month in months:


		print('doing month ' + str(month))

		for filetype in filetypes:

			print('doing ' + filetype + ' file')

			# create output file
			ncfile='%s/cDII.%s.%04d_%02d.nc'%(target_dir,filetype,year,month)
			nc = Dataset(ncfile,'w',format='NETCDF3_CLASSIC')

			# open input file
			inncfile='%s/cDII.%s.%04d_%02d.nc'%(atmodir,filetype,year,month)
			incan = Dataset(inncfile)
			incv=incan.variables


			# create grid from input file
			create_reduced_grid(nc,incan,basedate,lonrange,latrange)


			# loop over variables | reduces to minimun necessary area | exchange critical land points with ocean points

			for key in incv.keys():	

				if not key in ('time','lon','lat'):	

					# copy wind speeds
					vv = nc.createVariable(key,'f4',incv[key].dimensions)
					vv.units = incv[key].units
					vv.standard_name = incv[key].standard_name
					vv.coordinates = incv[key].coordinates

					# exchange points with ocean nn
					pdata=incv[key][:,latinds,loninds] # load reduced grid
					pdata[:,replace_dest]=pdata[:,isocean][:,nn_ocean_index] # exchange coomand

					vv[:]=pdata[:]
			incan.close()
			nc.close()


