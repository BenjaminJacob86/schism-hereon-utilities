## interpolate schism output to structured grid output
# for webpage presentation
import numpy as np
import cmocean # import ocean colormaps (before matplotlib)
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable # allow axis adjustment of colorbar
from glob import glob
import xarray as xr  # xarray netcdf access
from netCDF4 import Dataset,MFDataset, date2num
from cftime import utime
import sys
import os
# import hereon schism class schism.py
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/') 
from schism import *
from matplotlib.path import Path

plt.ion()
rundir='/gpfs/work/jacobb/data/RUNS/routine_linked/'
netcdfdir=rundir+'outputs/'
wavedir='/gpfs/work/jacobb/data/RUNS/routine_GB_wave/old/'
##############


os.chdir(rundir)
s=schism_setup()
lon,lat=np.asarray(s.lon),np.asarray(s.lat)


###### load nectf data #####  using netcdf4 instead
files=np.sort(glob(netcdfdir +'schout_????????.nc')) # last files to have most uptodate day all for long timesiers
wavefiles=np.sort(glob(wavedir +'wwm_hist_????.nc'))

#ds=xr.open_mfdataset(files)
#ds=ds.sel(time=ds['time'][11]) # reduce time dimension to 1 to plot one value at noon
#dsw=xr.open_mfdataset(wavefiles)
#dsw=ds.sel(time=ds['time'][11]) # reduce time dimension to 1 to plot one value at noon
#############
def create_grid(nc,vlevels=2):
	nc.createDimension('time',None)
	tv = nc.createVariable('time','f8',('time'))
	tv.long_name = 'Time'
	tv.standard_name = 'time'
	tv.units = ncin['time'].units
	tv.base_date = ncin['time'].base_date

	# level
	nc.createDimension('level',vlevels)
	lv = nc.createVariable('lvl','f4',('level'))
	lv.long_name = 'sigma_level'
	lv.standard_name = 'sigma_level'
	lv.units = '#'
	lv[:] = np.arange(vlevels)


	# lon
	nc.createDimension('lon',len(xq))
	#lv = nc.createVariable('lon','f4',('ny_grid','nx_grid'))
	lv = nc.createVariable('lon','f4',('lon'))
	lv.long_name = 'Longitude'
	lv.standard_name = 'longitude'
	lv.units = 'degrees_east'
	lv[:] = xq

	# lat
	nc.createDimension('lat',len(yq))
	lv = nc.createVariable('lat','f4',('lat'))
	lv.long_name = 'Latitude'
	lv.standard_name = 'latitude'
	lv.units = 'degrees_north'
	lv[:] = yq

def get_bdy_latlon(s):
	bdnodes=[]
	for land,ocean in list(zip(s.land_segments,s.bdy_segments)):		
		bdnodes.append(ocean)
		bdnodes.append(land[1:])
	bdnodes=np.hstack(bdnodes)		
	bdylon=np.asarray(s.lon)[bdnodes-1]
	bdylat=np.asarray(s.lat)[bdnodes-1]
	return (bdylon,bdylat)
bds=get_bdy_latlon(s)
Poly=Path(list(zip(bds[0],bds[1])))

# reduce grid to non nudging zone
cutout_nudging_zone=True
if cutout_nudging_zone:
	nu_nodes=np.loadtxt('cutoff_smoothed2.gr3',skiprows=2,max_rows=s.nnodes)
	nu_elems=np.asarray(np.loadtxt('cutoff_smoothed2.gr3',skiprows=2+s.nnodes),int)
	nu_elems=nu_elems[:,2:]-1
	isnu=np.asarray(nu_nodes[:,-1],bool)
#	keep=isnu[nu_elems].mean(axis=1)<1
#
#	s.nvplt=s.nvplt[keep,:]
	# causes gaps
	#s.lon[isnu]=np.nan
	#s.lat[isnu]=np.nan
	#s.lon[isnu]=-99999
	#s.lat[isnu]=-99999
#
#	A=A[keep]
#######



s.init_node_tree(latlon=True)


# Destination grid
ny=650 # nr of nodes on y driection
nx=np.int(np.ceil(ny*np.cos(54.5*np.pi/180)))
xq=np.linspace(np.nanmin(lon),np.nanmax(lon),nx)
yq=np.linspace(np.nanmin(lat),np.nanmax(lat),ny)
X, Y = np.meshgrid(xq, yq)

#invalid=0
invalid=~(Poly.contains_points(list(zip(X.flatten(),Y.flatten())))) 
for island in s.island_segments:
	islandPoly=Path(list(zip(lon[np.asarray(island)-1],lat[np.asarray(island)-1])))
	#plt.plot(islandPoly.vertices[:,0],islandPoly.vertices[:,1],'r')
	invalid+=islandPoly.contains_points(list(zip(X.flatten(),Y.flatten())))

lon[isnu]=np.nan
lat[isnu]=np.nan

# mask islands
#mask too far points	
dd,nn=s.node_tree_latlon.query(list(zip(X.flatten(),Y.flatten())))
invalid += (dd >0.01)

check_plot=True
if check_plot:
	D=np.asarray(s.depths)
	D[isnu]=np.nan
	vintp=D[nn]
	#D[isnu]=np.nan
	vintp[invalid]=np.nan
	plt.figure()
	plt.pcolormesh(X,Y,vintp.reshape(X.shape))










variables={'sea_surface_height':{'var':'zcor','level':-1,'units':'m','time_step':1},'sea_surface_temperature':{'var':'temp','level':-1,'units':'degC','time_step':1},'sea_surface_salinity':{'var':'salt','level':-1,'units':'psu','time_step':1},'sea_floor_temperature':{'var':'temp','level':0,'units':'degC','time_step':1},'sea_floor_salinity':{'var':'salt','level':0,'units':'psu','time_step':1}}

wave_variables={'significant_wave_height':{'var':'HS','units':'m','time_step':1},'zero_crossing_wave_period':{'var':'TM02','units':'s','time_step':1}}


#ncsin=MFDataset(files)
#ncsin['time'].units

step=1
nodes=np.arange(s.nnodes)
for file,filew in zip(files[:1],wavefiles[:1]):
	
	nc=Dataset('schism_routine_structured.nc','w',data_model='NETCDF3_CLASSIC')
	ncin=MFDataset(file)
	create_grid(nc,vlevels=1)
	

	# unterpolate from schism input	
	
	dryelems2d=np.asarray(ncin['wetdry_elem'][:,:],bool)
	drynodes2d=[np.isin(nodes,np.unique(s.nvplt[dryelems2d[ti,:],:])) for ti in range(0,24,step)]
	#drynodes=nodes in np.unique(s.nvplt[dryelems2d[ti,:],:]) 
	
	#tv = nc.createVariable('time','f8',('time'))
	nc['time'][:]=ncin['time'][:]
	#tv=ncin['time']
	#tv.long_name = variable
	#tv.standard_name = variable
	#tv.units = variables[variable]['units']
	
	
	for variable in variables:
		variable

		tv = nc.createVariable(variable,'f8',('time','lat','lon'))
		tv.long_name = variable
		tv.standard_name = variable
		tv.units = variables[variable]['units']


		
		for ti in range(0,24,step):#len(ncsin['elev'])
		
			datain=ncin[variables[variable]['var']][ti,:,variables[variable]['level']]
			# mask dryelements
			drynodes=drynodes2d[ti]
			datain[drynodes]=np.nan
			
			datain[isnu]=np.nan
			dataout=datain[nn]
			dataout[invalid]=np.nan
			dataout.reshape(X.shape)
			tv[ti,:]=dataout
		nc.sync()	
		
	# unterpolate from wwm input	
	ncin.close()
	ncin=Dataset(filew)

	# unterpolate from schism input	
	for variable in wave_variables:
		
		tv = nc.createVariable(variable,'f8',('time','lat','lon'))
		tv.long_name = variable
		tv.standard_name = variable
		tv.units = wave_variables[variable]['units']
		
		for ti in range(0,24,step):
		
			datain=ncin[wave_variables[variable]['var']][ti,:]
			# mask dryelements
			drynodes=drynodes2d[ti]
			datain[drynodes]=np.nan
			
			datain[isnu]=np.nan
			dataout=datain[nn]
			dataout[invalid]=np.nan
			dataout.reshape(X.shape)
			tv[ti,:]=dataout
		nc.sync()	
	ncin.close()	
	# unterpolate from wwm input	
	nc.close()
