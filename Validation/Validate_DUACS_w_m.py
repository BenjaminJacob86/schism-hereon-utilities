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

# Tested with Osisaf & Duacs

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
background=True
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
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/')
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/Lib/')
from schism import * # import schism functions
#from techit2 import * # import latex script
from TaylorDiagram import * # import taylordiagram
from data_and_model_classes import cmems
import pandas as pd
#import utide
from matplotlib.dates import date2num
import xarray as xr
#from numba import jit
import time
from matplotlib import path
plt.ion()

########## settings #################################

# GRIDED satelite products #############################
# satelite dir

#SateliteProdcut='OSISAF'
SateliteProdcut='DUACS'

# OSISAF:
if SateliteProdcut=='OSISAF':
	SATdir='/gpfs/work/ksddata/observation/remote/AVHRR/' #tgdir # directories (have to end with '/')
	#SATdir='/gpfs/work/ksddata/observation/remote/AVHRR_unzip/'
	varname='sea_surface_temperature'
	offset=-273.15
	anomaly=False										  # Satelite variable is anomalie and schism needs to
														  # be computed as such
	label='T mean [�C]'									  # LABEL FOR TIME SERIES PLOT

elif SateliteProdcut=='DUACS':
	# DUACS
	SATdir='/sciclone/data10/bjacob01/BlackSea/'  
	varname='sla'
	anomaly=True										  #	temporal mean will be removed from model		
	offset=0
	label='mean SLA [m]'
######################################################

####### SCHISM SETUP ################################

# schism setup(s) for cross validation TG stations are selected based on coverage of first setup
setupdir=['/gpfs/work/jacobb/data/SETUPS/NorthSea/NorthSea200mHelgoland/RERUN/']*2  # schism run directory, containts hgrid.* etc.
ncdir=[setupdir[0] + 'l_nudge_combine/'] 		  #   directory of schism nc output 
ncdir+=[setupdir[1] + 'combined_nudged2/'] 		  #   directory of schism nc output 

setupdir=['/sciclone/home10/yinglong/DISKS/vims20/BlackSea/RUN25b/']
ncdir=['/sciclone/home10/yinglong/DISKS/vims20/BlackSea/RUN25b/outputs/']

#setupdir=['/gpfs/work/jacobb/data/RUNS/BlackSea/2017/']
#ncdir=[setupdir[0] + 'combined/']



setup_names=['BlackSea2017']
varname_model='elev'
#varname_model='temp'
######################

################ CMEMS comparison #############
oceandir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/' #amm15 directory
use_amm=False						  # compare against amm15
pattern='TEM'						  # Pattern im AMM15 file names
varname_cmems='thetao'
###################

#year=2012							  # year to be analyesed 	 

######### image output ############################


vmin,vmax=0,16				# colorbar range
#dmin,dmax=-2.5,2.5          # colorbar range in difference plots
dmin,dmax=-0.2,0.2          # colorbar range in difference plots

outdir='/sciclone/data10/bjacob01/BlackSea'+SateliteProdcut+'_all_'+varname+'/'		  # output directory where images will be stored
pic_format='.png'
dpivalue=300

# colormaps
cmap=plt.cm.viridis			# Values
cmap2=plt.cm.jet			# Differences
cmap.set_over('m')
cmap.set_under('k')
cmap2.set_over('m')
cmap2.set_under('k')
if anomaly:				# use divergent colormap
	cmap=cmap2
######################################################

##### Latex Documentation ######################################
#put_pics_to_texdoc=True    										# images will be put in tex document
#latexname='OSISAF_valid.tex'									
#latextitle='SNS 2018'
#latextext='Tide Gauge validation of SNS ' + str(year) +'n against Emodnet data.'
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
	# Chunking along the dimension to be interpolated (1) is not yet supported. on W & M
	#s.nc=xr.open_mfdataset(schismfiles,concat_dim='time')
	s.nc=xr.open_mfdataset(schismfiles,concat_dim='time',chunks=None,)
	
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

# corresponding SAT	
print('loading satelite data')
#import dask
#dask.config.set({"array.slicing.split_large_chunks": True}) # problem on sciclone?


if SateliteProdcut=='DUACS':
	SATfiles=np.sort(np.concatenate([glob.glob(SATdir+'*.nc')]))
else:
	SATfiles=np.sort(np.concatenate([glob.glob(SATdir+'*{:d}*.nc'.format(year)) for year in years]))	
	
#SAT=xr.open_mfdataset(SATfiles) #memory laoding all
SAT=xr.open_dataset(SATfiles[0])
# reduce extend
try:
	SAT=SAT.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
	LON,LAT=np.meshgrid(SAT['lon'],SAT['lat'])
	sattype=1
except:	
	SAT=SAT.sel(longitude=slice(lonmin,lonmax),latitude=slice(latmin,latmax),time=slice(t0,t1))
	LON,LAT=np.meshgrid(SAT['longitude'],SAT['latitude'])
	sattype=2
time=SAT['time'].values
OBS=SAT[varname][0,:].values

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
	#s.nn=nnTree.query(list(zip(s.lon,s.lat)))[1]  # nn of schism to SAT
varnames=[varname_model]*len(setups)	


## compute anomalies for schism
if anomaly==True:
	print('computing temporal mean for anomalies')
	offsets=[ s.nc[varname_model].mean(axis=0).values  for s in setups.values()]
else:
	offsets=[ np.zeros(s.nnodes)  for s in setups.values()]

# corresponding cmems	
if use_amm:	
	cmemsfiles=np.sort(np.concatenate([glob.glob(oceandir+'*'+pattern+'*{:d}*.nc'.format(year)) for year in years]))
	extend=[lonmin,lonmax,latmin,latmax]
	varnames+=[varname_cmems]
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
	
	## compute anomalies for schism
	if anomaly==True:
		print('computing temporal mean for anomalies for cmems')
		offsets+=[amm15.nc[varname_model].mean(axis=0).flatten()]
	else:
		offsets+=[ np.zeros(amm15.LON.shape).flatten()]

		
		
################################ Start Computations and Plotting ################################		
ti=0
timei=time[ti]
Tmean={}
OBS=SAT[varname][ti,:].values+offset
OBS=np.ma.masked_array(OBS,mask=s.mask | np.isnan(OBS))
ivalid=OBS.mask==False
Tmean[0]=np.nanmean(OBS)


mean=np.zeros(OBS.shape)
var=np.zeros(OBS.shape)
means={i:np.zeros(OBS.shape) for i in range(len(setups))}
vars={i:np.zeros(OBS.shape) for i in range(len(setups))}
count={i:np.zeros(OBS.shape) for i in range(len(setups))}
bias={i:np.zeros(OBS.shape) for i in range(len(setups))}
rmse={i:np.zeros(OBS.shape) for i in range(len(setups))}
diffsqr={i:np.zeros(OBS.shape) for i in range(len(setups))}
cor={i:np.zeros(OBS.shape) for i in range(len(setups))}
stds={i:np.zeros(OBS.shape) for i in range(len(setups))}
cov={i:np.zeros(OBS.shape) for i in range(len(setups))}

names=['SAT']+setup_names	
for i,s in enumerate(setups.values()):
	#nci=s.nc.interp(time=str(timei)[:19]) # error w & m
	nci=s.nc.sel(time=str(timei)[:19])
	
	#OBSi=np.ma.masked_array(nci['temp'].values[s.nn].reshape(LON.shape),mask=s.mask)
	OBSi=np.ma.masked_array(nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape),mask=s.mask)
	Tmean[i+1]=np.nanmean(OBSi[ivalid])
	
day=(timei-t0)/np.timedelta64(1,'D')
t_ts=[day+(-np.timedelta64(1,'s'))/np.timedelta64(1,'D')	]	

# plotting
if not os.path.exists(outdir): os.mkdir(outdir)
n=1+len(setups)
m=3
fig, axes = plt.subplots(nrows=m, ncols=n)
nt=len(time)

nt=len(SATfiles)
#nt=10
fig1=plt.figure()
fig1.set_size_inches(12,9)
for file in SATfiles[:nt]:  # osisaf files  need to be iterated

	# reduce extend
	SAT=xr.open_dataset(file)
	if sattype==1:
		SAT=SAT.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
		LON,LAT=np.meshgrid(SAT['lon'],SAT['lat'])
		SAT=SAT.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
	else:	
		SAT=SAT.sel(longitude=slice(lonmin,lonmax),latitude=slice(latmin,latmax),time=slice(t0,t1))
		LON,LAT=np.meshgrid(SAT['longitude'],SAT['latitude'])
		SAT=SAT.sel(longitude=slice(lonmin,lonmax),latitude=slice(latmin,latmax),time=slice(t0,t1))
	time=SAT['time'].values
	#OBS=SAT[varname][0,:].values


	for ti,timei in enumerate(time):
		#if ti%10 == 0:
		#	print('{:d}/{:d}'.format(ti,nt))
		OBS=SAT[varname][ti,:].values+offset
		OBS=np.ma.masked_array(OBS,mask=s.mask | np.isnan(OBS))
		iuse=OBS.mask==False
		mean[iuse]+=OBS[iuse]
		ivalid=OBS.mask==False

		vmin=np.nanmin(OBS[ivalid])
		vmax=np.nanquantile(OBS[ivalid],0.99)
		
		# compute	
		Tmean[0]=np.hstack((Tmean[0],OBS.mean()))
		OBSis=[]
		for i,s in enumerate(setups.values()):
			#nci=s.nc.interp(time=str(timei)[:19])
			nci=s.nc.sel(time=str(timei)[:19])
			OBSi=(nci[varnames[i]].values.flatten()+offsets[i])[s.nn]
			OBSi=np.ma.masked_array(OBSi.reshape(LON.shape),mask=s.mask)
			OBSis.append(OBSi)
			Tmean[i+1]=np.hstack((Tmean[i+1],np.nanmean(OBSi[ivalid])))

			vmin=np.min((vmin,np.nanmin(OBSi[ivalid])))
			vmax=np.max((vmax,np.nanquantile(OBSi[ivalid],0.99)))
		
		# plot
		plt.clf()
		phs={}
		plt.subplot(m,n,1)
		phs[0]=plt.pcolormesh(LON,LAT,OBS,vmin=vmin,vmax=vmax,shading='gouraud',cmap=cmap)
		plt.xticks([])
		plt.yticks([])
		plt.title('SAT')
		ch=plt.colorbar()
		plt.suptitle(str(timei)[:19])
		
		for i,s in enumerate(setups.values()):
			plt.subplot(m,n,2+i)
			#nci=s.nc.interp(time=str(timei)[:19])
			nci=s.nc.sel(time=str(timei)[:19])
			OBSi=OBSis[i]
			phs[i+1]=plt.pcolormesh(LON,LAT,OBSi,cmap=cmap,vmin=vmin,vmax=vmax,shading='gouraud')
			plt.title(setup_names[i])
			ch=plt.colorbar(extend='both')
			plt.xticks([])
			plt.yticks([])
		
			plt.subplot(m,n,n+2+i)
			phs[i+1]=plt.pcolormesh(LON,LAT,OBSis[i]-OBS,vmin=dmin,vmax=dmax,cmap=cmap2,shading='gouraud')
			plt.xticks([])
			plt.yticks([])
			plt.title('above - SAT')
			ch=plt.colorbar(extend='both')
			#stats
			iuse2=iuse & ~np.isnan(OBSis[i])
			means[i][iuse2]+=OBSis[i][iuse2]
			count[i]+=iuse2
			biasi=(OBSis[i]-OBS)[iuse2]
			bias[i][iuse2]+=biasi
			diffsqr[i][iuse2]+=(biasi)**2
			
		#t_ts.append(timei)
		day=(timei-t0)/np.timedelta64(1,'D')
		t_ts.append(day)	
		plt.subplot(3,1,3)
		plt.cla()
		for i in range(1+len(setups)):
			plt.plot(t_ts,Tmean[i],linewidth=0.5)
			plt.grid()
		plt.ylabel(label)	
		plt.legend(names,ncol=4,loc='upper center',bbox_to_anchor=(0.5, 1.3),frameon=False)
		plt.tight_layout()
		plt.savefig(outdir+'{:04.0f}_Valid_SAT'.format(day),dpi=dpivalue)

## calculate statistics
for i,s in enumerate(setups.values()):	
	bias[i]/=count[i]	
	rmse[i]=np.sqrt(diffsqr[i]/count[i])	
	means[i]/=count[i]	
mean/=count[0]	

for file in SATfiles[:nt]:

	# reduce extend
	SAT=xr.open_dataset(file)
	if sattype==1:
		SAT=SAT.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
		LON,LAT=np.meshgrid(SAT['lon'],SAT['lat'])
		SAT=SAT.sel(lon=slice(lonmin,lonmax),lat=slice(latmin,latmax),time=slice(t0,t1))
	else:	
		SAT=SAT.sel(longitude=slice(lonmin,lonmax),latitude=slice(latmin,latmax),time=slice(t0,t1))
		LON,LAT=np.meshgrid(SAT['longitude'],SAT['latitude'])
		SAT=SAT.sel(longitude=slice(lonmin,lonmax),latitude=slice(latmin,latmax),time=slice(t0,t1))
	time=SAT['time'].values
	#OBS=SAT[varname][0,:].values

	for ti,timei in enumerate(time):
		OBS=SAT[varname][ti,:].values+offset
		OBS=np.ma.masked_array(OBS,mask=s.mask | np.isnan(OBS))
		iuse=OBS.mask==False
		var[iuse]+=(OBS[iuse]-mean[iuse])**2
		for i,s in enumerate(setups.values()):
			#nci=s.nc.interp(time=str(timei)[:19])
			nci=s.nc.sel(time=str(timei)[:19])
			#OBSi=np.ma.masked_array(nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape),mask=s.mask)
			#OBSi=nci[varnames[i]].values.flatten()[s.nn].reshape(LON.shape)
			OBSi=(nci[varnames[i]].values.flatten()+offsets[i])[s.nn]
			OBSi=np.ma.masked_array(OBSi.reshape(LON.shape),mask=s.mask)
			
			iuse2=iuse & ~np.isnan(OBSi)
			vars[i][iuse2]+=(OBSi[iuse2]-means[i][iuse2])**2
			cov[i][iuse2]+=(OBSi[iuse2]-means[i][iuse2])*(OBS[iuse2]-mean[iuse2])

for i,s in enumerate(setups.values()):
	cor[i]=	cov[i]/(np.sqrt(var*vars[i]))	

# bias cor rmse
fig1=plt.figure()
fig1.set_size_inches(12,9) 
#count[0]=count[1]
variables=[count,bias,rmse,cor]
varnames=['# valid','bias','rmse','cor']
#clims=[(0,count[0].max()),(-0.2,0.2),(0,0.3),(0,1)]
clims=[(0,count[0].max()),(-1.5,1.5),(0,2),(0,1)]
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
			plt.title('{:s} ({:.2f})'.format(varnames[j],np.nanmean(variable[i])))	
		else:	
			plt.title('avg= {:.2f}'.format(np.nanmean(variable[i])))	

#plt.suptitle(str(time[0])[:10] + ' - ' + str(time[ti])[:10])		
#plt.suptitle(str(t0+np.timedelta64(np.int(t_ts[1]),'D'))[:10] + ' - ' + str(time[-1])[:10])				
plt.suptitle(str(t0+np.timedelta64(np.int(t_ts[1]),'D'))[:10] + ' - ' + str(t0+np.timedelta64(np.int(t_ts[-1]),'D'))[:10])				
plt.savefig(outdir+'stats_ValidOsisaf.png'.format(day),dpi=dpivalue)
