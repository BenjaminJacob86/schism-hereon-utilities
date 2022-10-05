#export OMP_NUM_THREADS=1 # call before python
"""
Script for validation of schism output against cmems tide gauge data and optional intercomparison with amm15
validates SCHISM vs tide gauge stations (for 1 year) from EMOD net 
stored by SG in /work/gg0028/g260099/OSTIA/VALIDATION/VAL_TIDE/EMOD_TIDEG_DATA/
as monthly files.
1 - search data in subfolders by name, which is in specified <year>
2 - second sub select tide gauge stations within <dtol> distance from
closest schism grid nodes. 
3 - plot 
4 - derive skill parameter (bias,rmse,correlation) after interpolation of model data
to ovebservational time steps - the temporal resolution is reduced to the model time step
Comparison is done with zcor/elev data from SCHISMs closest grid point
Uses schism class stored in /pf/g/g260114/Programs/python/scripts/schism.py

call after 
module load python3
if mistralpp busy, performance better when called from run_sript

Output:
 pdf with validation images
'errorstats.csv' - csv file with error statistics
"""

__author__ = "Benjamin Jacob"
__copyright__   = "Copyright 2020 - 03\2021 Helmholtz-Zentrum Geesthacht GmbH"
__copyright__   = "Copyright 04\2021 - Helmholtz-Zentrum Hereon GmbH"

__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"

import os
import netCDF4
import sys
import csv
import matplotlib
from matplotlib import pyplot as plt
background=False
if background:
	matplotlib.use('Agg') # backend
else:
	plt.ion()	
import datetime as dt
import glob
from scipy.spatial import cKDTree
import datetime as dt
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')
from schism import * # import schism functions
from techit2 import * # import latex script
from TaylorDiagram import * # import taylordiagram
from data_and_model_classes import cmems
import pandas as pd
import utide
from matplotlib.dates import date2num
import xarray as xr
from numba import jit
import time
from matplotlib import path
plt.ion()

########## settings #################################
osafdir='/gpfs/work/ksddata/observation/remote/AVHRR_unzip/'
#'/gpfs/work/ksddata/observation/remote/AVHRR/' #tgdir # directories (have to end with '/')
oceandir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/' #amm15 directory

## schism setup(s) for cross validation TG stations are selected based on coverage of first setup
setupdir=['/gpfs/work/jacobb/data/SETUPS/check_wei/DKRZrun/sns_combined/']  # schism run directory, containts hgrid.* etc.
#ncdir=[setupdir[0] + 'outputs/Combined/'] 		  #   directory of schism nc output 
					  # True: use station output False: extract from netcdf (slow)

					  # schism setup(s) for cross validation TG stations are selected based on coverage of first setup
setupdir=['/gpfs/work/chenw1/SCHISM/cfgs/Control_Run/SNS_GB_Combined/']
#setupdir+=['/gpfs/work/xu/data/NorthSea/SNS_Helgoland200m_test/RUN_cosdatIII/']
ncdir=[setupdir[0] + 'outputs/Combined_drag_turned_1/']
					  
					  
#setupdir+=['/gpfs/work/jacobb/data/SETUPS/NorthSea/NorthSea200mHelgoland/RERUN2/']  # schism run directory, containts hgrid.* etc.
#ncdir+=[setupdir[1] + 'combined_no_nudge/'] 		  #   directory of schism nc output 
setup_names=['schism_sns'] #,'schism_era']
######################
 
#outdir='/gpfs/work/jacobb/data/SETUPS/check_wei/DKRZrun/'+'/OSISAF2/'	  # output directory where images will be stored
outdir='/gpfs/work/jacobb/data/ValidOsisafWei/'

year=2018							  # year to be analyesed 	 
dtol=0.05           				  # distance tolerance in degree lon/lat towards tg stations 

use_amm=True						  # compare against amm15
# images ############################
pic_format='.eps'
dpivalue=300
cmap='jet'
Rcircle=150 							# radius of colorcoded cricles used in scatter maps
limit_to_data=True   					# limit map plots to bounding box of data.
cmap=plt.cm.viridis
cmap2=plt.cm.jet
cmap.set_over('m')
cmap.set_under('k')
cmap2.set_over('m')
cmap2.set_under('k')



### Latex output
put_pics_to_texdoc=True    										# images will be put in tex document
latexname='OSISAF_valid.tex'									
latextitle='SNS 2018'
latextext='Tide Gauge validation of SNS ' + str(year) +'n against Emodnet data.'
#############################################################

######### load SCHISM setup   ##################################
cwd=os.getcwd()
setups={}
for i,folder in enumerate(setupdir):
	os.chdir(folder)
	s=schism_setup()
	s.nntree = cKDTree(list(zip(s.lon,s.lat))) 
	setups[i]=s
	
	schismfiles=[] 
	for iorder in range(6): # check for schout_nc files until 99999
		schismfiles+=glob.glob(ncdir[i]+'schout_'+'?'*iorder+'.nc')
	nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
	schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
	nrs=list(np.asarray(nrs)[np.argsort(nrs)])
	s.nc=xr.open_mfdataset(schismfiles)
	
	# bounadary path
	bdnodes=[]
	x,y=np.asarray(s.lon),np.asarray(s.lat)
	for ocean, land in zip(s.bdy_segments,s.land_segments):
		bdnodes+=ocean #+land[1:]
		bdnodes+=land[1:]
	bdnodes=np.asarray(bdnodes)-1
	s.p=path.Path([(x[bdnodes][i],y[bdnodes][i]) for i in range(len(bdnodes))])
	setups[i]=s

# select surface laywrs	
for s in setups.values():
	s.nc=s.nc.sel(nSCHISM_vgrid_layers=-1)

# time range
t0=np.max([(s.nc['time'][0].values) for s in setups.values() ])
t1=np.min([(s.nc['time'][-1].values) for s in setups.values() ])
year0=np.int(str(t0)[:4])
year1=np.int(str(t1)[:4])
years=np.unique(range(year0,year1+1))

# geo range
lonmin=np.min([np.min(s.lon) for s in setups.values() ])
lonmax=np.max([np.max(s.lon) for s in setups.values() ])
latmin=np.min([np.min(s.lat) for s in setups.values() ])
latmax=np.max([np.max(s.lat) for s in setups.values() ])


# corresponding osisaf	
print('loading osisaf data')
import dask
dask.config.set({"array.slicing.split_large_chunks": True})
osisaffiles=np.sort(np.concatenate([glob.glob(osafdir+'*{:d}*.nc'.format(year)) for year in years]))
#osisaf=xr.open_mfdataset(osisaffiles) #memory laoding all
osisaf=xr.open_dataset(osisaffiles[0])

# reduce extend
#osisaf=osisaf.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
osisaf=osisaf.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax))
LON,LAT=np.meshgrid(osisaf['lon'],osisaf['lat'])
SST=osisaf['sea_surface_temperature'][0,:].values
time=osisaf['time'].values
#plt.pcolormesh(LON,LAT,SST)



# mask for Satelite outside domain
restrict_to_model=True
nnTree=cKDTree(list(zip(LON.flatten(),LAT.flatten())))
for s in setups.values():
	s.init_node_tree()
	if restrict_to_model:
		mask=s.p.contains_points(list(zip(LON.flatten(),LAT.flatten()))).reshape(LON.shape)
		s.mask=~mask
	else:
		s.mask=np.zeros(LON.shape,bool)
		
	s.nn=s.node_tree_latlon.query(list(zip(LON.flatten(),LAT.flatten())))[1]
	#s.nn=nnTree.query(list(zip(s.lon,s.lat)))[1]  # nn of schism to osisaf
varnames=['temp']*len(setups)	

# corresponding cmems	
if use_amm:	
	pattern='TEM'   
	cmemsfiles=np.sort(np.concatenate([glob.glob(oceandir+'*'+pattern+'*{:d}*.nc'.format(year)) for year in years]))
	extend=[lonmin,lonmax,latmin,latmax]
	varnames+=['thetao']	
	class cmems():
		def __init__(self,cmemsfiles,extend,t0,t1):
			self.nc = xr.open_mfdataset(cmemsfiles)
			self.nc = self.nc.sel(lon=slice(extend[0],extend[1]),lat=slice(extend[2],extend[3]),time=slice(str(t0)[:19],str(t1)[:19]),depth=0)
			self.LON,self.LAT=np.meshgrid(self.nc['lon'],self.nc['lat'])
			self.tree=cKDTree(list(zip(self.LON.flatten(),self.LAT.flatten())))
	amm15=cmems(cmemsfiles,extend,t0,t1)			
	amm15.nn=amm15.tree.query(list(zip(LON.flatten(),LAT.flatten())))[1]
	amm15.mask=s.mask
	setups[len(setups)]=amm15
	setup_names+=['amm15']
	

vmin,vmax=0,16
dmin,dmax=-2.5,2.5







ti=0
timei=time[ti]
Tmean={}
SST=osisaf['sea_surface_temperature'][ti,:].values-273.15
SST=np.ma.masked_array(SST,mask=s.mask | np.isnan(SST))
ivalid=SST.mask==False
Tmean[0]=np.nanmean(SST)


mean=np.zeros(SST.shape)
var=np.zeros(SST.shape)
means={i:np.zeros(SST.shape) for i in range(len(setups))}
vars={i:np.zeros(SST.shape) for i in range(len(setups))}
count={i:np.zeros(SST.shape) for i in range(len(setups))}
bias={i:np.zeros(SST.shape) for i in range(len(setups))}
rmse={i:np.zeros(SST.shape) for i in range(len(setups))}
diffsqr={i:np.zeros(SST.shape) for i in range(len(setups))}
cor={i:np.zeros(SST.shape) for i in range(len(setups))}
stds={i:np.zeros(SST.shape) for i in range(len(setups))}
cov={i:np.zeros(SST.shape) for i in range(len(setups))}

names=['osisaf']+setup_names	
for i,s in enumerate(setups.values()):
	nci=s.nc.interp(time=str(timei)[:19])
	#SSTi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
	SSTi=np.ma.masked_array(nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape),mask=s.mask)
	Tmean[i+1]=np.nanmean(SSTi[ivalid])
	
day=(timei-t0)/np.timedelta64(1,'D')
t_ts=[day+(-np.timedelta64(1,'s'))/np.timedelta64(1,'D')	]	



# plotting
if not os.path.exists(outdir): os.mkdir(outdir)
n=1+len(setups)
m=3
fig, axes = plt.subplots(nrows=m, ncols=n)
nt=len(time)

nt=len(osisaffiles)
#nt=10
for file in osisaffiles[:nt]:
	osisaf=xr.open_dataset(file)
	osisaf=osisaf.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
	time=osisaf['time'].values

	for ti,timei in enumerate(time):
		#if ti%10 == 0:
		#	print('{:d}/{:d}'.format(ti,nt))
		SST=osisaf['sea_surface_temperature'][ti,:].values-273.15
		SST=np.ma.masked_array(SST,mask=s.mask | np.isnan(SST))
		iuse=SST.mask==False
		mean[iuse]+=SST[iuse]
		ivalid=SST.mask==False
		vmin=np.nanmin(SST[ivalid])
		vmax=np.nanquantile(SST[ivalid],0.99)
		Tmean[0]=np.hstack((Tmean[0],SST.mean()))
		
		plt.clf()
		phs={}
		plt.subplot(m,n,1)
		phs[0]=plt.pcolormesh(LON,LAT,SST,vmin=vmin,vmax=vmax,shading='gouraud')
		plt.xticks([])
		plt.yticks([])

		plt.title('osisaf')
		ch=plt.colorbar()
		SSTis=[]
		plt.suptitle(str(timei)[:19])

		for i,s in enumerate(setups.values()):
			plt.subplot(m,n,2+i)
			nci=s.nc.interp(time=str(timei)[:19])
			#SSTi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
			SSTi=np.ma.masked_array(nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape),mask=s.mask)
			phs[i+1]=plt.pcolormesh(LON,LAT,SSTi,cmap=cmap,vmin=vmin,vmax=vmax,shading='gouraud')
			plt.title(setup_names[i])
			ch=plt.colorbar(extend='both')
			SSTis.append(SSTi)
			plt.xticks([])
			plt.yticks([])
			Tmean[i+1]=np.hstack((Tmean[i+1],np.nanmean(SSTi[ivalid])))
		#fig.colorbar(ph1, ax=axes.ravel().tolist(),'horizontal')
		#fig.colorbar(ph1, ax=axes.ravel().tolist(),shrink=0.5,pad=0.5)
		#fig.colorbar(ph1, ax=axes.ravel().tolist(),pad=0.99)
		
		#for i,s in enumerate(setups.values()):
			plt.subplot(m,n,n+2+i)
			phs[i+1]=plt.pcolormesh(LON,LAT,SSTis[i]-SST,vmin=dmin,vmax=dmax,cmap=cmap2,shading='gouraud')
			plt.xticks([])
			plt.yticks([])

			#plt.title(setup_names[i] + ' - osisaf')
			plt.title('above - osisaf')
			ch=plt.colorbar(extend='both')
			#stats
			iuse2=iuse & ~np.isnan(SSTis[i])
			means[i][iuse2]+=SSTis[i][iuse2]
			count[i]+=iuse2
			biasi=(SSTis[i]-SST)[iuse2]
			bias[i][iuse2]+=biasi
			diffsqr[i][iuse2]+=(biasi)**2
			
		#t_ts.append(timei)
		day=(timei-t0)/np.timedelta64(1,'D')
		t_ts.append(day)	
		plt.subplot(3,1,3)
		plt.cla()
		for i in range(1+len(setups)):
			plt.plot(t_ts,Tmean[i])
			plt.grid()
		plt.ylabel('T mean [°C]')	
		plt.legend(names,ncol=4,loc='upper center',bbox_to_anchor=(0.5, 1.3))
		plt.tight_layout()
		plt.savefig(outdir+'{:04.0f}_ValidOsisaf'.format(day),dpi=dpivalue)

for i,s in enumerate(setups.values()):	
	bias[i]/=count[i]	
	rmse[i]=np.sqrt(diffsqr[i]/count[i])	
	means[i]/=count[i]	
mean/=count[0]	

for file in osisaffiles[:nt]:
	osisaf=xr.open_dataset(file)
	osisaf=osisaf.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
	#SST=osisaf['sea_surface_temperature'][0,:].values
	time=osisaf['time'].values

	for ti,timei in enumerate(time):
		SST=osisaf['sea_surface_temperature'][ti,:].values-273.15
		SST=np.ma.masked_array(SST,mask=s.mask | np.isnan(SST))
		iuse=SST.mask==False
		var[iuse]+=(SST[iuse]-mean[iuse])**2
		
		for i,s in enumerate(setups.values()):
			nci=s.nc.interp(time=str(timei)[:19])
			#SSTi=np.ma.masked_array(nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape),mask=s.mask)
			SSTi=nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape)
			iuse2=iuse & ~np.isnan(SSTi)
			vars[i][iuse2]+=(SSTi[iuse2]-means[i][iuse2])**2
			cov[i][iuse2]+=(SSTi[iuse2]-means[i][iuse2])*(SST[iuse2]-mean[iuse2])

for i,s in enumerate(setups.values()):
	cor[i]=	cov[i]/(np.sqrt(var*vars[i]))	

# bias cor rmse
plt.figure()
count[0]=count[1]
variables=[count,bias,rmse,cor]
varnames=['# valid','bias','rmse','cor']
clims=[(0,count[0].max()),(-2,2),(0,4),(0,1)]
#cmaps=[cmap,cmap2,cmap,cmap]
plt.clf()
i=0
for i,s in enumerate(setups.values()):
	for j,variable in enumerate(variables):
		plt.subplot(len(setups),4,4*i+1+j)
		plt.pcolormesh(LON,LAT,np.ma.masked_array(variable[i],mask=mask==False),vmin=clims[j][0],vmax=clims[j][1],shading='gouraud',cmap=cmap2)
		plt.colorbar(extend='both')
		#plt.clim(clims[j])
		plt.xticks([])
		plt.yticks([])
		if j==0:
			plt.ylabel(setup_names[i])
		if i==0:	
			plt.title(varnames[j])	
#plt.suptitle(str(time[0])[:10] + ' - ' + str(time[ti])[:10])		
#plt.suptitle(str(t0+np.timedelta64(np.int(t_ts[1]),'D'))[:10] + ' - ' + str(time[-1])[:10])				
plt.suptitle(str(t0+np.timedelta64(np.int(t_ts[1]),'D'))[:10] + ' - ' + str(timei)[:10])				
tstr=str(t0+np.timedelta64(np.int(t_ts[1]),'D'))[:10] + ' - ' + str(timei)[:10]
tstr=tstr.replace(' ','')
plt.savefig(outdir+'stats_ValidOsisaf'+tstr,dpi=dpivalue)

	
########
#timei=time[0]
#SST=osisaf['sea_surface_temperature'][0,:].values-273.15
#vmin=np.nanmin(SST)
#vmax=np.nanquantile(SST,0.99)
#
#phs={}
#plt.close('all')
#fig, axes = plt.subplots(nrows=m, ncols=n)
#plt.subplot(m,n,1)
#phs[0]=plt.pcolormesh(LON,LAT,np.ma.masked_array(SST,mask=s.mask),vmin=vmin,vmax=vmax)
#plt.title('osisaf')
#ch=plt.colorbar()
#SSTis=[]
#plt.suptitle(str(timei)[:19])
#
#Tmean={}
#Tmean[0]=np.nanmean(SST)
#for i,s in enumerate(setups.values()):
#	plt.subplot(m,n,2+i)
#	nci=s.nc.interp(time=timei)
#	SSTi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
#	phs[i+1]=plt.pcolormesh(LON,LAT,SSTi,cmap=cmap,vmin=vmin,vmax=vmax)
#	plt.title(setup_names[i])
#	ch=plt.colorbar(extend='both')
#	SSTis.append(SSTi)
##fig.colorbar(ph1, ax=axes.ravel().tolist(),'horizontal')
##fig.colorbar(ph1, ax=axes.ravel().tolist(),shrink=0.5,pad=0.5)
##fig.colorbar(ph1, ax=axes.ravel().tolist(),pad=0.99)
#
#for i,s in enumerate(setups.values()):
#	plt.subplot(m,n,n+2+i)
#	nci=s.nc.interp(time=timei)
#	SSTi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
#	phs[i+1]=plt.pcolormesh(LON,LAT,SSTis[i]-SST,vmin=dmin,vmax=dmax,cmap=cmap2)
#	plt.title(setup_names[i] + ' - osisaf')
#	ch=plt.colorbar(extend='both')
#	SSTis.append(SSTi)
#	Tmean[i+1]=np.nanmean(SSTi)
#	#
#plt.subplot(3,1,3)
#plt.tight_layout()
#
#
#
#
## plotting
#for ti,timei in enumerate(time[1:]):
#SST=osisaf['sea_surface_temperature'][ti+1,:].values-273.15
#SST=np.ma.masked_array(SST,mask=s.mask | np.isnan(SST))
#vmin=np.nanmin(SST)
#vmax=np.nanquantile(SST,0.99)
#
#
#
#phs[0].set_array(SST.flatten()[:-1])
#plt.suptitle(str(timei)[:19])
#plt.title('osisaf')
#ch=plt.colorbar()
#SSTis=[]
#
#
#Tmean={}
#Tmean[0]=np.nanmean(SST)
#for i,s in enumerate(setups.values()):
#	plt.subplot(m,n,2+i)
#	nci=s.nc.interp(time=timei)
#	SSTi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
#	phs[i+1]=plt.pcolormesh(LON,LAT,SSTi,cmap=cmap,vmin=vmin,vmax=vmax)
#	plt.title(setup_names[i])
#	ch=plt.colorbar(extend='both')
#	SSTis.append(SSTi)
##fig.colorbar(ph1, ax=axes.ravel().tolist(),'horizontal')
##fig.colorbar(ph1, ax=axes.ravel().tolist(),shrink=0.5,pad=0.5)
##fig.colorbar(ph1, ax=axes.ravel().tolist(),pad=0.99)
#
#for i,s in enumerate(setups.values()):
#	plt.subplot(m,n,n+2+i)
#	nci=s.nc.interp(time=timei)
#	SSTi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
#	phs[i+1]=plt.pcolormesh(LON,LAT,SSTis[i]-SST,vmin=dmin,vmax=dmax,cmap=cmap2)
#	plt.title(setup_names[i] + ' - osisaf')
#	ch=plt.colorbar(extend='both')
#	SSTis.append(SSTi)
#	Tmean[i+1]=np.nanmean(SSTi)
#	#
#plt.subplot(3,1,3)
#plt.tight_layout()
