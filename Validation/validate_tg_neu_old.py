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
# own and 3d party libraries
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
plt.ion()

########## settings #################################
tgdir='/gpfs/work/ksddata/observation/insitu/TideGauge/MyOcean/' #tgdir # directories (have to end with '/')
oceandir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/' #amm15 directory

# schism setup(s) for cross validation TG stations are selected based on coverage of first setup
setupdir=['/gpfs/work/jacobb/data/SETUPS/NorthSea/NorthSea200mHelgoland/RERUN/']  # schism run directory, containts hgrid.* etc.
ncdir=[setupdir[0] + 'combined/'] 		  #   directory of schism nc output 
use_station_in=[False]					  # True: use station output False: extract from netcdf (slow)

setupdir+=['/gpfs/work/jacobb/data/SETUPS/NorthSea/NorthSea200mHelgoland/RERUN2/']  # schism run directory, containts hgrid.* etc.
ncdir+=[setupdir[1] + 'combined_no_nudge/'] 		  #   directory of schism nc output 
setup_names=['schism_dwd','schism_era']
use_station_in+=[True]
######################
 
outdir=setupdir[-1]+'/TGvalid3/'	  # output directory where images will be stored

year=2018							  # year to be analyesed 	 
dtol=0.05           				  # distance tolerance in degree lon/lat towards tg stations 


skipdays=5							  # skip days in beginnig															
remove_mean=True  					  # remove temporal mean from Data and Model to compare 
use_amm=True						  # compare against amm15

#--- what to plot True/False:			
overview_map=True												
satistic_maps=True
full_timeseries=False 				      # for flull year data hard to see anything
first_two_Weeks=True 				      # zomed period 			
running_avg_window=np.timedelta64(25,'h') # running mean filter in hours
monthly_subplots=True
taylor_diag=True
consts=['M2','S2','M4','M6']									# tidal constituents to compute
tidal_bar_plots=False
tidal_difference_maps=False


# FIlter Criteria Exclude bad statons or add known offsets
exclude_stations=[]            #'HusumTG',
offset_stations={}			   #{'WittduenTG':-5,'BuesumTG':-5} # add value to stations 

###################

# images ############################
pic_format='.eps'
dpivalue=300
cmap='jet'
Rcircle=150 							# radius of colorcoded cricles used in scatter maps
limit_to_data=True   					# limit map plots to bounding box of data.

### Latex output
put_pics_to_texdoc=True    										# images will be put in tex document
latexname='SNS_amm15_valid.tex'									
latextitle='SNS 2018'
latextext='Tide Gauge validation of SNS ' + str(year) +'n against Emodnet data.'
#############################################################

######### load SCHISM setup   ##################################
cwd=os.getcwd()
setups={}
output={}
access={}
for i,folder in enumerate(setupdir):
	os.chdir(folder)
	s=schism_setup()
	s.nntree = cKDTree(list(zip(s.lon,s.lat))) 
	setups[i]=s
	if use_station_in[i]:
		if os.path.exists('station.in')	& os.path.exists(ncdir[i]+'staout_1'):
			staout=schism_station_output()
			output[i]=staout
		else:
			print('station output does not exist, exiting program')
			exit
	else:
		schismfiles=[] 
		for iorder in range(6): # check for schout_nc files until 99999
			schismfiles+=glob.glob(ncdir[i]+'schout_'+'?'*iorder+'.nc')
		nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
		schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
		nrs=list(np.asarray(nrs)[np.argsort(nrs)])
		s.nc=xr.open_mfdataset(schismfiles)
		access[i]=s.nc
	
	
# use cmems TG files  # pre select in domain
pattern='SSH'							 							# pattern to check in amm15 files 				
folders=np.sort(glob.glob('/gpfs/work/ksddata/observation/insitu/CMEMS/NorthWestShelf/TG/{:d}*'.format(year)))
files=glob.glob(folders[0]+'/*.nc')

# select stations withing domain (tolerance distance) # and constructe file acess and nearest neighbours
stations={'coord':[],'names':[],'TG':{}}
sources=setup_names.copy()
if use_amm:
	sources.append('amm')
for key in sources:	
	stations[key]={'nc':[],'time':[],'nn':[],'coord':[],'zeta':[],'names':[]}
names=[]	
for i,setup_name in enumerate(setup_names):
	s=setups[i]
	print('find neighbours for ' + setup_name )
	if use_station_in[i]:
		staout=output[i]
	for file in files:
		a=xr.open_dataset(file)
		coord=np.float(a.geospatial_lon_min),np.float(a.geospatial_lat_min)
		if use_station_in[i]:
			D=np.sqrt((staout.coords[:,0]-coord[0])**2+(staout.coords[:,1]-coord[1])**2)
			nn=np.argmin(D)
			lldists=D[nn]
		else:
			lldists,nn=s.nntree.query(coord)		
		# in domain
		if lldists < dtol:
			name=file[file.rindex('/')+1:].split('_')[-2]
			print('using ' + name)
			
			if i==0:
				ncfiles=[file.replace('{:d}01'.format(year),'{:d}{:02d}'.format(year,mon+1)) for mon in range(12) if os.path.exists(file.replace('{:d}01'.format(year),'{:d}{:02d}'.format(year,mon+1)))]
				ncs=MFDataset(ncfiles,aggdim='TIME')
				t=(ncs['TIME'][:])
				try:
					igood=(ncs['SLEV_QC'][:] < 3)[:,0] # select good data (ncs['TIME_QC'][:]<3) and
				except:
					if np.ma.isMaskedArray(ncs['SLEV'][:,0]):
						igood=ncs['SLEV'][:,0].mask==False # np.ones(len(t),bool)
					else:	
						igood=np.ones(len(t),bool)
				t=t[igood]
				t0=dt.datetime.strptime(ncs['TIME'].units[11:30],'%Y-%m-%dT%H:%M:%S')
				timeunit=ncs['TIME'].units[:ncs['TIME'].units.index(' ')]
				stations['TG'][name]={'time':np.asarray([t0+eval('dt.timedelta({:s}=float(ti))'.format(timeunit)) for ti in t]),'zeta':(ncs['SLEV'][:])[igood]}
				names.append(name)
				stations['coord'].append(coord)
			stations[setup_name]['names'].append(name)	
			stations[setup_names[i]]['coord'].append((s.lon[nn],s.lat[nn]))
			stations[setup_names[i]]['nn'].append(nn)
	print('done selecting and loading tide gauges, found: '+str(len(stations['coord'])) + ' stations meeting criterea')
########################### CMEMS data acces #####################################

## remove not share stations from schism files
# shared stations between setups  needs to be done flexible
#s.plot_domain_boundaries()
#stations[setup_names[i]]['coord']

#for i in range(len(a)):
#	print('dwd: ' + a[i] + 'era: ' + b[i] )

a=np.asarray(list(stations[setup_names[0]]['names']))
b=np.asarray(list(stations[setup_names[1]]['names']))

shared=[]
notshared=[]
remove=[]
for name in a:
	if name in b:
		shared.append(name)
	else:	
		notshared.append(name)
		remove.append(np.where(a==name)[0][0])
for key in ['nn','coord','names']:
	delete=[stations[setup_names[0]][key][ind] for ind in remove]
	for deltei in delete:
		stations[setup_names[0]][key].remove(deltei)  


delete=[stations['coord'][ind] for ind in remove]
for deltei in delete:
	stations['coord'].remove(deltei)  
delete=[names[ind] for ind in remove]
for deltei in delete:
	names.remove(deltei)
for deltei in delete:
	del stations['TG'][deltei]	
################################################################

a=np.asarray(list(stations[setup_names[0]]['names']))
b=np.asarray(list(stations[setup_names[1]]['names']))
for i in range(len(a)):
	print('TG: ' + list(stations['TG'].keys())[i] + 'dwd: ' + a[i] + ' era: ' + b[i] )




coords=np.asarray(stations['coord'])
######### load AMM15 setup ######################################
if use_amm:
	print('load cmems Data ')
	files=np.sort(glob.glob(oceandir+'*'+pattern+'*'+str(year)+'*'))
	amm15=cmems(SSHfile=files[0])
	amm15.varnames={'ssh':'zos'}
	lldists2,nn2=amm15.tree.query(stations['coord'])		
	ii,jj=np.unravel_index(nn2,amm15.LON.shape)
	amm15.nc=xr.open_mfdataset(files)
	stations['amm']['time']=amm15.nc.sel(time=slice("{:d}-12-31".format(year-1), "{:d}-01-01".format(year+1))).indexes['time'].to_datetimeindex()
	stations['amm']['zeta']=np.asarray([amm15.nc[amm15.varnames['ssh']][i,:].values[ii,jj] for i in range(len(amm15.nc['time']))])
	print('done load cmems Data ')
####################################

# load data from schism next neighbours to TG statsions
modeldt=[]
dates=[]
for i,setup_name in enumerate(setup_names):
	print(setup_name)
	s=setups[i]
	name=names[i]
	#ind=stations[setup_name]['names'].index(name)
	if use_station_in[i]:
		staout=output[i]
		modeldt.append(np.timedelta64(staout.time[1]-staout.time[0]))
		dates.append(np.asarray(staout.time,np.datetime64))   # convert to np.datetime64
		stations[setup_name]['zeta']=staout.station_out['elev'][:,stations[setup_name]['nn']]
		modeldt[i]=np.timedelta64(modeldt[i])
		stations[setup_name]['time']=dates[i]#np.asarray(dates,np.datetime64)
		lons,lats=staout.coords[:,0],staout.coords[:,1]
	else:
		s.nc=access[i]
		lons,lats=np.asarray(s.lon),np.asarray(s.lat)
		print('load SCHISM sea level time series at TG next neighbouring nodes')
		use_elev='elev' in s.nc.variables.keys()
		s.nc2=s.nc.sel(nSCHISM_hgrid_node=stations[setup_name]['nn'],nSCHISM_vgrid_layers=-1)
		n=len(schismfiles)	
		print(str(n) + 'files')
		if use_elev:	
			stations[setup_name]['zeta']=s.nc2['elev'].values
		else: # use zcor
			stations[setup_name]['zeta']=s.nc2['zcor'].values
		stations[setup_name]['time']=s.nc2['time'].values
		modeldt.append(stations[setup_name]['time'][1]-stations[setup_name]['time'][0])
		dates.append(stations[setup_name]['time'])
	print('load model Tide Gauge Data')


#if remove_mean:
#		latextext+=' Time Series were mean removed and so is rmse then. In Taylor it is always mean removed'
# plt.ion() has to be off to work in background mode
exclude_by_names=['dummy1','dummy2'] # list of Station names to remove from analysis
				     # to exclude bad station in second iterarion (names have to match those in the overviewmap)

# offset
for key in offset_stations.keys():	
	if key in names:
		stations[key][name]['zeta']+=offset_stations[key]
		latextext+=' station ' + name + ' was offseted by ' + str(offset_stations[key])
stations['amm']['time']=np.asarray(stations['amm']['time'],np.datetime64)
### remove mean
if remove_mean:
	print('removing mean')
	for key in list(stations.keys())[2:]:
		if key=='TG':
			stations['TG_mean']={}
			for name in names:
				stations['TG_mean'][name]=stations['TG'][name]['zeta'].mean()
				stations['TG'][name]['zeta']-=stations['TG_mean'][name]
		else:	
			stations[key]['zeta_mean']=np.nanmean(stations[key]['zeta'],axis=0)
			stations[key]['zeta']-=stations[key]['zeta_mean']
		latextext+=' tide gauges and models where mean removed for every step excluding bias computation'

		
####### temporal interpolation of model to data  and error statistics ##########
print('interploating to common time steps of data and calculating error statistics')
sources=['TG']+sources
for i,key in enumerate(sources):
	stations[key]['interp']={name:{'time':0,'zeta':0} for name in names}
	stations[key]['stats']={'std':[],'bias':[],'rmse':[],'cor':[],'std_rel':[]}
	for subkey in stations[key]['stats'].keys():
		stations[key]['stats'][subkey]={name:0 for name in names}

n=len(names)
time_schism=stations[setup_names[0]]['time'] # dirret acces seems faster
for i,name in enumerate(names):

	start = time.time()
	print(str(i)+'/'+str(n)+' '+name)
	stations['TG'][name]['time']=np.asarray(stations['TG'][name]['time'],np.datetime64)	
	Todates=stations['TG'][name]['time']
	ilast=(Todates<=time_schism[-1]).sum()
	ifirst=(Todates<=time_schism[0]).sum()#-1?
	ifirst=max(ifirst,0)
	Todates=Todates[ifirst:ilast]

	# interpolate to model time step but at maximum the temporla resolution of the refrence model	
	Todates2=[Todates[0]]
	for date in Todates:
		if (date-Todates2[-1]) >= modeldt[0]:
			Todates2.append(date)
	Todates2=np.asarray(Todates2)
	stations['TG']['interp'][name]['time']=Todates2
	a,ainb,bina=np.intersect1d(Todates, Todates2, assume_unique=False, return_indices=True)
	stations['TG']['interp'][name]['zeta']=stations['TG'][name]['zeta'][ifirst:ilast][ainb]
	
	std0=np.std(stations['TG']['interp'][name]['zeta'])
	stations['TG']['stats']['std'][name]=std0
	end = time.time()
	print("prepate TG Elapsed = %s" % (end - start))
	
	for key in sources[1:]:	
		#start = time.time()
		tin=(stations[key]['time']-time_schism[0])/(Todates2[1]-Todates2[0])
		tout=(Todates2-time_schism[0])/(Todates2[1]-Todates2[0])
		stations[key]['interp'][name]['time']=Todates2
		fintp=interp1d(tin, stations[key]['zeta'][:,i])
		stations[key]['interp'][name]['zeta']=fintp(tout)
		#end = time.time()
		#print("interp Elapsed = %s" % (end - start))

		#start = time.time()
			
		# calculate error  stats ###########################  SLOW	
		Zeta_tg=stations['TG']['interp'][name]['zeta'][:,0]
		Zetamod=stations[key]['interp'][name]['zeta']
		stations[key]['stats']['bias'][name]=stations[key]['zeta_mean'][i]-stations['TG_mean'][name]
		
		diff=Zeta_tg-Zetamod
		rmse=np.sqrt((diff**2).mean())
		std=np.std(stations[key]['interp'][name]['zeta']) # mean removed
		if np.isnan(np .nanmax(stations[key]['interp'][name]['zeta']) ):
			cor=np.nan
		else:	
			cor=np.corrcoef(Zetamod,Zeta_tg)[0,1]
		
		stations[key]['stats']['rmse'][name]=rmse # mean removed
		stations[key]['stats']['std'][name]=std
		stations[key]['stats']['cor'][name]=cor
		stations[key]['stats']['std_rel'][name]=std/std0
		#stations=calc_error_stats(stations,key)
		#end = time.time()
		#print("calc stats Elapsed = %s" % (end - start))


	
########## Plotting #####################################
print('plotting results')
if not os.path.exists(outdir): os.mkdir(outdir) 
##### Plot selected stations
coords=np.asarray(stations[setup_names[0]]['coord'])
xmin,ymin=np.min(coords,axis=0)
xmax,ymax=np.max(coords,axis=0)
x,y=coords[:,0],coords[:,1]
figures=[]
captions=[]
fname='0_TideGuageStationLocations.eps'
caption='Overview of Tide Gauge Stations.'
figures.append(fname)
lons,lats=np.asarray(setups[0].lon),np.asarray(setups[0].lat)
captions.append(caption)
if overview_map:
	plt.clf()
	s.plot_domain_boundaries()
	fig=plt.gcf()
	fig.set_size_inches(11,8,forward=True)
	plt.plot(lons[stations[setup_names[0]]['nn']],lats[stations[setup_names[0]]['nn']],'bo')
	#plt.plot( np.asarray(stations[setup_names[-1]]['coord'])[:,0],np.asarray(stations[setup_names[-1]]['coord'])[:,-1],'bo')
	for coord,name in zip(stations['coord'],names):
		lon,lat=coord
		plt.plot(lon,lat,'r+')
		plt.text(lon,lat,' '+name[:3],rotation=50, rotation_mode='anchor')
	plt.title('TG Stations')
	
	if limit_to_data:
		plt.xlim((xmin-1,xmax+1))
		plt.ylim((ymin-1,ymax+1))
	plt.tight_layout()	
	plt.savefig(outdir+fname,dpi=dpivalue)	
	plt.close()
###################################	


#name='HelgolandTG'
#ind=list(names).index(name)
#for key in setup_names:
#	plt.plot(stations[key]['time'],stations[key]['zeta'][:,ind],'.-')


####### Plot statistics #######################
names=np.asarray(names)
shortnames=np.asarray([name[:3] for name in names])
ivalids=[np.asarray(np.ones(len(coords)),bool)]*len(setup_names)

n=len(sources)+len(sources)%2
widths=np.asarray([ n*0.4/n/2 * (i+1) for i in range(np.int(n/2))])
widths=np.concatenate((-widths[::-1],widths))

colors=plt.cm.tab10(range(len(sources)))  #['b','r']
if use_amm:
	ivalid=~np.isnan(stations['amm']['zeta_mean'])
	ivalids.append(ivalid)	
plt.close()

Rcircles=[Rcircle*0.25**np.float(i) for i in range(len(sources)-1)]
phs={}
if satistic_maps:
	labels={'bias':'bias [m]','rmse':'rmse [m]','cor':'correlation [-]','std_rel':'relative std [-]'}
	for key in list(stations[setup_names[0]]['stats'].keys())[1:]:

		data=[np.asarray(list(stations[setup_names[0]]['stats'][key].values()))]
		ilarger=np.abs(data[0])>3 # station offset
		data[0][ilarger]-=np.floor(data[0][ilarger])
		vmin=data[0].min()
		vmax=data[0].max()
		if len(sources[1:]) > 1:
			for i,stp in enumerate(sources[1:]):
				k=i+1
				data.append(np.asarray(list(stations[stp]['stats'][key].values())))
				data[k]=np.ma.masked_array(data[1],mask=np.isnan(data[1]))
				ilarger=np.abs(data[k])>3 # station offset
				data[k][ilarger]-=np.floor(data[1][ilarger])
		vmin=np.nanquantile(np.concatenate(data),0.1)
		vmax=np.nanquantile(np.concatenate(data),0.99)
		if key=='cor':
			vmin=0
			vmax=1
			

		plt.clf()
		s.plot_domain_boundaries()
		fig=plt.gcf()
		fig.set_size_inches(11,8,forward=True)
		for nr,model in enumerate(sources[1:]):
			phs[nr]=plt.scatter(x[ivalids[nr]],y[ivalids[nr]],s=Rcircles[nr],c=data[nr][ivalids[nr]]+nr,vmin=vmin,vmax=vmax)			
		ph=plt.colorbar()
		ph.set_label(labels[key])
		plt.legend(phs.values(),sources[1:],loc='upper center',ncol=2)

		for coord,name in zip(stations['coord'],names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		plt.tight_layout()
		if limit_to_data:
			plt.xlim((xmin-.2,xmax+.2))
			plt.ylim((ymin-.2,ymax+.2))

		fname='1b_'+labels[key][:labels[key].index(' ')]+'_schism_amm'+pic_format
		plt.savefig(outdir+fname,dpi=dpivalue)	
		caption='Map of '+key+ ' between ' + str(Todates2[0]) + ' and ' + str(Todates2[-1])
		figures.append(fname)
		captions.append(caption)
		plt.close()

		# bar charts
		plt.clf()
		for nr,model in enumerate(sources[1:]):
			if nr==0:
				phs[nr]=plt.bar(np.where(ivalids[nr])[0],data[nr][ivalids[nr]],widths[nr],align='edge',label=model,color=colors[nr],tick_label=shortnames[ivalids[nr]])
				plt.grid()
			else:
				phs[nr]=plt.bar(np.where(ivalids[nr])[0],data[nr][ivalids[nr]],widths[nr],align='edge',label=model,color=colors[nr])
		plt.xticks(rotation=45)
		plt.ylabel(labels[key])
		plt.legend(phs.values(),sources[1:],loc='upper center',ncol=2)
		fname='1c_'+labels[key][:labels[key].index(' ')]+'_bar_plot.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)	
		plt.close()
		caption='Barplot of '+key+ ' between ' + str(Todates2[0]) + ' and ' + str(Todates2[-1])
		figures.append(fname)
		captions.append(caption)
###################################	
	

	
###### Plot first two weeks of data #####
inds=[]
if first_two_Weeks: 
	print('making plot of first two weeks')
	tmax=stations[setup_names[0]]['time'][0]+np.timedelta64(14,'D')
	#inds=stations[setup_names[0]]['time']<tmax
	for key in sources[1:]:
		inds.append(stations[key]['time']<tmax)
	
	for i,name in enumerate(names):
		inds2=stations['TG'][name]['time']<tmax
		plt.clf()
		plt.plot(stations['TG'][name]['time'][inds2],stations['TG'][name]['zeta'][inds2],'k-')
		for nr,key in enumerate(sources[1:]):
			plt.plot(stations[key]['time'][inds[nr]],stations[key]['zeta'][:,i][inds[nr]],'-')
		plt.legend(sources,ncol=4,loc='upper center')					
		plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.title(name)
		plt.grid()
		fname='3_TG'+name+pic_format
		plt.savefig(outdir+fname,dpi=dpivalue)
		caption='Time Series at Station ' + name
		figures.append(fname)
		captions.append(caption)
################################################	

## filter longterm TS
for key in sources[1:]:
	dates=stations[key]['time']
	zeta=stations[key]['zeta']
	zeta_filter=np.zeros(zeta.shape)
	for i,date in enumerate(dates):
		ind=np.abs(dates-date)<= running_avg_window/2	
		zeta_filter[i,:]=zeta[ind,:].mean(axis=0)
	stations[key]['zeta_filter']=zeta_filter.copy()

###slow 
#for i, name in enumerate(names):
#	start = time.time()
#	print('running average for ' + name)	
#	dates=stations['TG'][name]['time']
#	zeta=stations['TG'][name]['zeta']
#	# running mean
#	stations['TG'][name]['zeta_filter']=np.asarray([zeta[((dates[ti]-running_avg_window/2) <= dates) & (dates <= dates[ti]+running_avg_window/2)].mean() for ti in range(len(dates))])
#	end = time.time()
#	print("prepate TG Elapsed = %s" % (end - start))
#
###slow 
#wndw_half=running_avg_window/2
#for i, name in enumerate(names):
#	start = time.time()
#	print('running average for ' + name)	
#	dates=stations['TG'][name]['time']
#	zeta=stations['TG'][name]['zeta']
#	# running mean
#	zeta_filter=np.asarray([zeta[((dates[ti]-wndw_half) <= dates) & (dates <= dates[ti]+wndw_half)].mean() for ti in range(len(dates))])
#	stations['TG'][name]['zeta_filter']=zeta_filter
#	end = time.time()
#	print("prepate TG Elapsed = %s" % (end - start))

wndw_half=running_avg_window/2	
@jit(nopython=True)
def runnin_mean(dates,zeta,wndw_half):
	zeta_filter=np.zeros(zeta.shape)
	for ti in range(len(dates)):
		zeta_filter[ti]=np.mean(zeta[((dates[ti]-wndw_half) <= dates) & (dates <= dates[ti]+wndw_half)])	
	return zeta_filter
#np.abs(dates-dates[ti])<= windw_half
	
##slow 
wndw_half=running_avg_window/2
for i, name in enumerate(names):
	start = time.time()
	print('running average for ' + name)	
	dates=stations['TG'][name]['time']
	zeta=stations['TG'][name]['zeta']
	stations['TG'][name]['zeta_filter']=runnin_mean(dates,zeta,wndw_half)
	end = time.time()
	print("prepate TG Elapsed = %s" % (end - start))	
	
	
	
###### compare data	non normalized 
if monthly_subplots:
	print('making monthly subplots of complete time series')
	
	for i, name in enumerate(names):
		plt.clf()
		
		data=[stations['TG'][name]['zeta_filter'][:,0]]
		vmin=data[0].min()
		vmax=data[0].max()
		if len(sources[1:]) > 1:
			k=0
			for i,stp in enumerate(sources[1:]):
				k=i+1
				data.append( stations[stp]['zeta_filter'][:,i])
				data[k]=np.ma.masked_array(data[k],mask=np.isnan(data[k]))
		#vmin=np.nanquantile(np.concatenate(data),0.1)
		#vmax=np.nanquantile(np.concatenate(data),0.99)	
		vmin=np.nanmin(np.concatenate(data))
		vmax=np.nanmax(np.concatenate(data))	


		Time=stations['TG'][name]['time']
		plt.clf()
		for month in range(1,13): 
			plt.subplot(4,3,month)
			inds=np.asarray([date.astype(object).month ==month for date in Time ])
			plt.plot(Time[inds],data[0][inds],'b',linewidth=1.5)
			for nr,stp in enumerate(sources[1:]):
					time_model=stations[stp]['time']
					inds=np.asarray([np.int(str(date)[5:7])==month for date in time_model ])
					plt.plot(time_model[inds],data[nr+1][inds],'--')
			if month != 2:
				plt.title('month: '+str(month))
			else:
				plt.title(name)
			plt.grid()
			plt.ylim((vmin,vmax))
			if month%3 !=1:
				plt.tick_params(axis='y',labelleft=False)  
			plt.tick_params(axis='x',labelbottom=False)  	
			if month==12:
				plt.legend(sources,loc='upper center',ncol=3)
		fig=plt.gcf()
		fig.set_size_inches(8,8,forward=True)
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()	
		fname='4_TG'+name+pic_format
		plt.savefig(outdir+fname,dpi=dpivalue)
		plt.close()	
		caption='Monthly timeseries at station '+str(name) +'. Time series where filtered with a running window of '+ str(running_avg_window)
		figures.append(fname)
		captions.append(caption)	

if taylor_diag:
	
	key=sources[1]
	samples=[[ [stations[key]['stats']['std_rel'][name],stations[key]['stats']['cor'][name],name[:3]] for name in names ]]
	for key in sources[2:]:
		samples.append([ [stations[key]['stats']['std_rel'][name],stations[key]['stats']['cor'][name],name[:3]] for name in names ])

	plt.clf()
	dia=plotTaylor(samples[0],stdref=1,extend=True) #negative
	for nr,key in enumerate(sources[2:]):
		#Add models to Taylor diagram
		for i,(stddev, corrcoef, name) in enumerate(samples[nr]):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)

	fname='5_taylorA.eps'	
	figures.append(fname)
	caption='Taylor diagram of full signals interpolated to Tide Gauge Timesteps. Black: Schism, Red Amm15.'
	captions.append(caption)
	plt.savefig(outdir+fname,dpi=dpivalue)
	plt.close()
	
	plt.clf()
	dia=plotTaylor(samples[0],stdref=1,extend=False) #negative
	for nr,key in enumerate(sources[2:3]):
		#Add models to Taylor diagram
		for i,(stddev, corrcoef, name) in enumerate(samples[nr]):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
	fname='5_taylorB.eps'	
	figures.append(fname)
	caption='Taylor diagram of full signals interpolated to Tide Gauge Timesteps. Black: Schism, Red Amm15.'
	captions.append(caption)
	plt.savefig(outdir+fname,dpi=dpivalue)
	plt.close()

	
	ivlaid=~np.isnan(list(stations[sources[-1]]['stats']['bias'].values()))
	bias=[]
	rmse=[]
	cor=[]
	
	for stp in sources[1:]:
		bias.append(np.mean(np.abs(np.asarray(list(stations[stp]['stats']['bias'].values()))[ivlaid])))
		rmse.append(np.mean((np.asarray(list(stations[stp]['stats']['rmse'].values()))[ivlaid])))
		cor.append(np.mean((np.asarray(list(stations[stp]['stats']['cor'].values()))[ivlaid])))
	k=0
	plt.clf()
	xpos=np.arange(3)
	for k,stp in enumerate(sources[1:]):
		plt.bar(xpos+k/n,[bias[k],rmse[k],cor[k]],1/n,align='edge',label=stp)
	plt.xticks(xpos+2/n,('bias','mae','cor'))	
	fname='5.2_sation_average.eps'	
	plt.ylabel('station averaged values')
	figures.append(fname)
	caption='Statiscal properties averaged over commonly covered stations'
	captions.append(caption)
	plt.savefig(outdir+fname,dpi=dpivalue)
	plt.close()


	
# harmonic analysis
if tidal_bar_plots:
	Amps={'tg':dict.fromkeys(consts),'schism':dict.fromkeys(consts),'amm':dict.fromkeys(consts)}
	Phas={'tg':dict.fromkeys(consts),'schism':dict.fromkeys(consts),'amm':dict.fromkeys(consts)}
	
	for const in consts:
		for const in consts:
			for key in 'tg','schism','amm':
				Amps[key][const]=np.zeros(len(names))
				Phas[key][const]=np.zeros(len(names))
	
	
	for i,name in enumerate(names):			
		print(name)
		utout=utide.solve(date2num(Time[name]), u=np.asarray(ZETA[name][:,0]), v=None, lat=coords[i,1])
		utout_schism=utide.solve(date2num(dates), u=np.asarray(zeta_schism[:,i]), v=None, lat=coords[i,1])
		if True in np.isnan(zeta_amm[:,i]):
			utout_amm=utout_schism
			ammfactor=np.nan
		else:
			utout_amm=utide.solve(date2num(dates2), u=np.asarray(zeta_amm[:,i]), v=None, lat=coords[i,1])
			ammfactor=1.0
		for const in consts:
			ind=np.where(utout.name==const)
			Amps['tg'][const][i]=utout.A[ind]
			Phas['tg'][const][i]=utout.g[ind]
			ind=np.where(utout_schism.name==const)
			Amps['schism'][const][i]=utout_schism.A[ind]
			Phas['schism'][const][i]=utout_schism.g[ind]
			ind=np.where(utout_amm.name==const)
			Amps['amm'][const][i]=utout_amm.A[ind]*ammfactor
			Phas['amm'][const][i]=utout_amm.g[ind]*ammfactor
			
	for const in consts:
		plt.clf()
		plt.subplot(2,1,1)
		xbar=np.asarray(range(len(names)))
		ph1=plt.bar(xbar,Amps['tg'][const],width=0.25,align='edge',label='TideGauge',color='b',tick_label=shortnames)
		ph2=plt.bar(xbar+0.25,Amps['schism'][const],width=0.25,align='edge',label='SCHSIM',color='r') #,tick_label=valid_names
		ph3=plt.bar(xbar+0.5,Amps['amm'][const],width=0.25,align='edge',label='Amm15',color='k') #,tick_label=valid_names
		plt.xticks(rotation=45)
		plt.grid()
		plt.legend([ph1,ph2,ph3],('TG','SCHISM','Amm15'),ncol=3,loc='upper center')
		plt.ylabel(const +' Amplitude [m]')
		plt.tight_layout()
		#plt.xlabel('station')	

		plt.subplot(2,1,2)
		xbar=np.asarray(range(len(names)))
		ph1=plt.bar(xbar,Phas['tg'][const],width=0.25,align='edge',label='TideGauge',color='b',tick_label=shortnames)
		ph2=plt.bar(xbar+0.25,Phas['schism'][const],width=0.25,align='edge',label='SCHSIM',color='r') #,tick_label=valid_names
		ph3=plt.bar(xbar+0.5,Phas['amm'][const],width=0.25,align='edge',label='Amm15',color='k') #,tick_label=valid_names
		plt.xticks(rotation=45)
		plt.grid()
		#plt.legend([ph1,ph2,ph3],('TG','SCHISM','Amm15'),ncol=3)
		plt.ylabel(const +' Phase [deg]')
		plt.tight_layout()
		#plt.xlabel('station')	
		
		caption=const+' tidal amplitude (top) and phase (bottom)'
		fname=const+'amp_phase.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)
		plt.close()			
		figures.append(fname)
		captions.append(caption)

		
if tidal_difference_maps:
	for const in ['M2','M4']:		
		plt.clf()
		s.plot_domain_boundaries()
		damp=Amps['schism'][const]-Amps['tg'][const]
		damp2=Amps['amm'][const]-Amps['tg'][const]
		vmin=np.max((np.min((damp.min(),damp2.min())),-0.5))
		vmax=np.min((np.max((damp.max(),damp2.max())),0.5))
		ph1=plt.scatter(x,y,s=Rcircle,c=damp)
		plt.clim((vmin,vmax))
		ph2=plt.scatter(x[ivalid],y[ivalid],s=Rcircle*0.25,c=damp2[ivalid])
		plt.clim((vmin,vmax))

		plt.legend([ph1,ph2],('SCHISM','AMM15'),loc='upper center',ncol=2)
		ch=plt.colorbar()
		ch.set_label(const+ 'Amplitude difference')
		plt.set_cmap(cmap)
		for coord,name in zip(coords,names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		
		fname=const+'_amplitude_difference.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)
		caption='Map of '+const+ ' amplitude difference (model - data)'
		figures.append(fname)
		captions.append(caption)#


		plt.clf()
		s.plot_domain_boundaries()
		dphi=Phas['schism'][const]-Phas['tg'][const]
		dphi2=Phas['amm'][const]-Phas['tg'][const]
		dphi=-np.sign(dphi)*((np.abs(dphi)>180)*360-np.abs(dphi))
		dphi2=-np.sign(dphi2)*((np.abs(dphi2)>180)*360-np.abs(dphi2))
		vmin=np.min((dphi.min(),dphi2.min()))
		vmax=np.max((dphi.max(),dphi2.max()))
		ph1=plt.scatter(x,y,s=Rcircle,c=dphi)
		plt.clim((vmin,vmax))
		ph2=plt.scatter(x[ivalid],y[ivalid],s=Rcircle*0.25,c=dphi2[ivalid])
		plt.clim((vmin,vmax))

		plt.legend([ph1,ph2],('SCHISM','AMM15'),loc='upper center',ncol=2)
		ch=plt.colorbar()
		ch.set_label(const+ 'phase difference')
		plt.set_cmap(cmap)
		for coord,name in zip(coords,names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		fname=const+'phase_lag.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)	
		caption='Map of '+const+ ' phase lag (model -data)'
		figures.append(fname)
		captions.append(caption)
		
		
		
# bar charts
val=np.ma.masked_array(list(data2.values()),mask=mask_amm)

	
	
	
#rms=np.sqrt(1-R[name]**2)*stdmodel[name]
#rmse[name]
	
# create latex + pdf
if put_pics_to_texdoc:
		print('generating tex doc and pdf')
		os.chdir(outdir)
		#techit(latexname,latextitle,latextext)
		techit(latexname,latextitle,latextext,figures,captions)
		if os.path.isfile(latexname[:latexname.rindex('.')]+'.pdf'):
			 os.remove(latexname[:latexname.rindex('.')]+'.pdf')
		os.system('pdflatex '+latexname)
		print('done generating '+latexname[:latexname.rindex('.')]+'.pdf')


## write error stats to csv
M=[['station name','mean data','mean SCHISM','mean Amm','stddata','std schism','std Amm','bias schism','bias Amm','rmse schism','rmse Amm','correlation schism','correlation Amm']]
for name in names:
	M.append([name,ZETAmean[name],zeta_schismmean[name],zeta_ammmean[name],stddata[name],stdmodel[name],stdmodel2[name],bias[name],bias2[name],rmse[name],rmse2[name],R[name],R2[name]])
M=np.asarray(M)
with open(outdir+'errorstats.csv', 'w') as csvFile:
	writer = csv.writer(csvFile)
	for i in range(M.shape[0]):
		writer.writerow(M[i,:])
csvFile.close()
		
# hceck shift for coralation +-1 hour	
#output all data as pickle	
analysis={'stations':{'coords':coords,'names':names},
'mes':{'data':{'time':Time,'ssh':ZETA},'data_intp':{'time':timeintp,'ssh':zetaintp}},
'schism':{'data':{'time':dates,'ssh':zeta_schism},'data_intp':{'time':timeintp,'ssh':zetaintp_schism}},
'amm':{'data':{'time':dates2,'ssh':zeta_amm},'data_intp':{'time':timeintp,'ssh':zetaintp_amm}},
'mes_stats':{'mean':ZETAmean,'std':stddata},
'schism_stats':{'mean':zeta_schismmean,'std':stdmodel},
'amm_stats':{'mean':zeta_ammmean,'std':stdmodel2},
'schism_errors_stats':{'bias':bias,'rmse':rmse,'cor':R},
'amm_errors_stats':{'bias':bias2,'rmse':rmse2,'cor':R2}}

pickle.dump(analysis,open(outdir+"analysisdata","wb"))


obd_compare=False
if obd_compare:
	## compasrions at open boundary
	ibd=np.asarray(s.bdy_segments[0][::50])
	bdcoords=[(s.lon[i-1],s.lat[i-1]) for i in ibd]
	bdnn=amm15.tree.query(bdcoords)[1]
	ii2,jj2=np.unravel_index(bdnn,amm15.LON.shape)

	date=dates2[0]
	amm15.update(date)
	tt=amm15.nc['time'][:]
	ssh=amm15.nc[amm15.varnames['ssh']][:][:,ii2,jj2]

	for date in dates2[nt2:nt2*10+1:nt2]:
		print('loading cmems ' + str(date))
		amm15.update(date)
		tt=np.concatenate((tt,amm15.nc['time'][:]))
		ssh=np.concatenate((ssh,amm15.nc[amm15.varnames['ssh']][:][:,ii2,jj2]))
		
	ncs=MFDataset(schismfiles[:10])	
	ssh2=ncs['elev'][:,ibd]
	tt2=ncs['time'][:]


	date=dates2[0]
	amm15.update(date)
	tt=amm15.t0+dt2*np.arange(len(tt))
	tt2=reftime+dt.timedelta(seconds=modeldt)*np.arange(1,len(tt2)+1) 



	nsub=int(np.ceil(np.sqrt(len(bdnn))))
	plt.close()
	for i in range(nsub):
		plt.subplot(nsub,1,i+1)
		ph1=plt.plot(tt,ssh[:,i],'b')
		ph2=plt.plot(tt2,np.asarray(ssh2[:,i]),'r--')
		plt.legend(('Amm','SCHISM'),frameon=False)
		plt.title('bd node ' + str(ibd[i]))
		plt.ylabel('ssh')
		plt.xlabel('time')
	#plt.tight_layout()
	plt.gcf().autofmt_xdate()
	plt.savefig('bnd_comparisons',dpi=dpivalue)

