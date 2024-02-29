"""
Read multiple scenario (subfoldes) schism output and compute monthly statistics (min,mean,max quantiles)
for selected variables and save as analysis results
"""
from glob import glob
import os
import sys
#sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
from schism import* # import schism class to read grid structure

# parallel stuff
#%matplotlib inline
import xarray as xr
import matplotlib.pyplot as plt

from dask.diagnostics import ProgressBar
from dask_jobqueue import SLURMCluster
from distributed import Client, progress
from distributed.utils import tmpfile

import dask
import distributed
dask.__version__, distributed.__version__


##### setting ###################################


# take statistics over period:
date0='20170701'
date1='20170801'

# Folders of involded scenarios
basedir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/'#

# subfolder names

# labels for the experiments in the netcdf output
experiments=['Veg_REF','Veg_CNTRL','Veg_max','Veg_HE','Veg_LE']

## construct folders
dirs=['Veg_REF','Veg_CNTRL','Veg_max','Veg_HE','Veg_LE']

setupdir=[basedir+dir+'/' for dir in dirs]
setupdir[1]='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdirs=[setupdiri+'outputs_all/' for setupdiri in setupdir]#'outputs/'] # where the outputs are.
basedir=setupdir[0]

quantaties=['mean','min','max','std','quantile95'] # ['mean','min','max','std','quantile95'] potential quantaties (percentiles, eg.   'quantile95','quantile5', ...
###############################################

#varnames=['dryFlagNode','sigWaveHeight','dominantDirection','meanWavePeriod','zeroDowncrossPeriod','peakPeriod','frictionalVelocity','turbulentKineticEner','totalSuspendedLoad','meanDirSpreading','rmsOrbitalVelocity']
varnames=['dryFlagNode','elevation','sigWaveHeight','totalSuspendedLoad']

# reduce varlist selection to speed up loading
varlist=['out2d','turbulentKineticEner','totalSuspendedLoad','horizontalVelX','horizontalVelY']
vector_vars=['depthAverageVel','bottomStress'] # add variables with x and y component
#----------------------------

# time period - potenetial input argument
if len(sys.argv)>2:
	date0=int(sys.argv[1])
	date1=int(sys.argv[2])
	print(mon0)
	print(mon1)

print('extracting month {:d} - {:d}'.format(mon0,mon1))
print(varnames)


#output name
prefix='c_mon_Ana_mveg_jan_'
###################################################



################## code begins #################################

#####read grid_strtucture
os.chdir(basedir)
s=schism_setup()


# find the maximum commonly available stack between experiments ()
min_max_stack=278
key='out2d'
for ncdir in ncdirs:
	files=np.hstack([np.sort(glob.glob('{:s}{:s}_{:s}.nc'.format(ncdir,key,'?'*iorder))) for iorder in range(1,6)])
	print(files[-1])
	min_max_stack=np.minimum(int(files[-1][:-3].split('_')[-1]),min_max_stack)
#access=[schism_outputs_by_variable(ncdiri,max_stack=min_max_stack-1,varlist=varlist) for ncdiri in ncdirs]
access=[schism_outputs_by_variable(ncdiri,max_stack=-1,varlist=varlist) for ncdiri in ncdirs]
##############################################

#access=[]
#for ncdiri in ncdirs:
#	print(ncdiri)
#	access.append(schism_outputs_by_variable(ncdiri,max_stack=min_max_stack-1,varlist=varlist))

######### restict to commont time range for experiments ##########
varname='elevation'
nts=[len(accessi.get('elevation')) for accessi in access]
nt=np.min(nts)
print('selected folders have ' + str(nt) + ' time steps')
imin=np.argmin(nts)

#check latest file instead to be faster
tmin=access[imin].ds[access[imin].vardict[varname]]['time'][0].values #[0]
tmax=access[imin].ds[access[imin].vardict[varname]]['time'][nt-1].values #[0]
print( 'tmin: {:f} days'.format(tmin/86400))
print( 'tmax: {:f} days'.format(tmax/86400))
##################################################################


###################### identify vector variables #############################
vector_vars=[vari[:-1] for vari in access[0].vardict.keys() if vari[-1] =='Y'] # stack components for convenience
for variablefile in np.unique(list(access[0].vardict.values())):
	if variablefile not in ['hvel','wind']:
		for accessi in access:
			if variablefile in vector_vars:
				accessi.ds[variablefile][variablefile]=accessi.ds[variablefile][variablefile].sel(time=slice(tmin,tmax))
			else:
				accessi.ds[variablefile]=accessi.ds[variablefile].sel(time=slice(tmin,tmax))
#####################################################################	



############

###### Reults ##############################	
varnames+=vector_vars

# x and y compneints of vectr vars
for comp in ['X','Y']:
	varnames+=[tmp+comp for tmp in list(np.asarray(vector_vars)[np.asarray([0,1,3,4],int)])]
	#units+=['m/s','kg/m/s^2','m/s','g/l']
	#symbols+=['U','btress', 'Umean','WED']


# intialized output array	
results=dict.fromkeys(experiments)
for expi in results.keys():
	results[expi]=dict.fromkeys(varnames)
	for vari in results[expi].keys():
		results[expi][vari]=dict.fromkeys(quantaties)


# get time in datetime		
p=param('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_LE//')
reftime=np.asarray(dt.datetime(int(p.get_parameter('start_year')),
int(p.get_parameter('start_month')),
int(p.get_parameter('start_day')),
int(p.get_parameter('start_hour')),0,0),np.datetime64)		
a=access[0]
t=a.get('elevation').time.values
startdate=reftime+3600*np.timedelta64(1,'s')
dates=reftime+np.asarray(t,np.timedelta64(1,'s'))
years = dates.astype('datetime64[Y]').astype(int) + 1970
months = dates.astype('datetime64[M]').astype(int) % 12 + 1
days = (dates - dates.astype('datetime64[M]') + 1)/np.timedelta64(1,'D')


import pandas as pd
date0=np.asarray(pd.to_datetime(date0, format='%Y%m%d'),np.datetime64)
date1=np.asarray(pd.to_datetime(date1, format='%Y%m%d'),np.datetime64)


date0=dt.datetime(startdate)
dates



index=(dates>= date0) & (dates <= date1)
levels=[len(s.vgrid[1])-1,0]
label=['surface','bottom']
mode='w'
print('doing month ' + str(mon))

# schism time indices for selectio in seconds since starts (depends on schism version)
t0sec=t[index][0]
t1sec=t[index][-10]

ndry=0
# lGet dry entries - oop over scenarios
for i,expi in enumerate(results.keys()):
	ndry+=access[i].get('dryFlagNode').sel(time=slice(t0sec,t1sec)).values
ndry=ndry>0



filenameout=prefix+str(years[0])+'_{:s}_{:s}.nc'.format(str(date0)[:10],str(date1)[:10]) # outputfile
for ivar,vari in enumerate(list(results[expi].keys())):   # variavle loop
	
	for i,expi in enumerate(results.keys()): 			  # scenario loop

	# 2D  or surface + bottom  fields
		if 'nSCHISM_vgrid_layers' in  access[i].get(vari).dims:
			nrun=2
			results[expi][vari+'_surface']=dict.fromkeys(results[expi][vari].keys())
			results[expi][vari+'_bottom']=dict.fromkeys(results[expi][vari].keys())
		else: 
			nrun=1	
		for irun in range(nrun):
			
			# make subselection of time steps where all are wet for statistics
			
			
			if vari in vector_vars: # vector variables -calc absolute value
				isvec=True
				if nrun==1:						
					data=access[i].get(vari).sel(time=slice(t0sec,t1sec)).values
				else:
					data=access[i].get(vari).sel(time=slice(t[index][0],t[index][-1]),nSCHISM_vgrid_layers=levels[irun]).values
				u=data[0,:]
				v=data[1,:]
				data=np.sqrt(u*u+v*v)
				u=np.ma.masked_array(u,mask=ndry)
				v=np.ma.masked_array(v,mask=ndry)
			else:					# scalar variables
				isvec=False
				
				if nrun==1:						
					data=access[i].get(vari).sel(time=slice(t0sec,t1sec)).values
				else:
					data=results[expi][vari][quanti]=access[i].get(vari).sel(nSCHISM_vgrid_layers=levels[irun],time=slice(t0sec,t1sec)).values
					
			data=np.ma.masked_array(data,mask=ndry) #mask
		
			for iquant,quanti in enumerate(quantaties):	# statistical quantaty loop
				print('exporting '+ expi  +' '+vari+ ' ' + quanti)
				
				name=expi+'_'+vari+'_'+quanti	
				if 'quantile' in quanti:
					value=float(quanti.split('quantile')[1])/100					 #get quantile nominaö
					results[expi][vari][quanti]=tmp=np.quantile(data,axis=0,q=value)
					
					if nrun==1:
						exec('tmp=xr.Dataset(data_vars=dict({:s}=(["nSCHISM_hgrid_node"], tmp)))'.format(name))
					else:	
						exec('tmp=xr.Dataset(data_vars=dict({:s}_{:s}=(["nSCHISM_hgrid_node"], tmp)))'.format(name,label[irun]))
					if isvec: # macht das sinn
						utmp=np.quantile(u,axis=0,q=value)	
						vtmp=np.quantile(v,axis=0,q=value)	
						if nrun==1:
							exec('utmp=xr.Dataset(data_vars=dict({:s}=(["nSCHISM_hgrid_node"], utmp)))'.format(vari+'X_'+quanti))
							exec('vtmp=xr.Dataset(data_vars=dict({:s}=(["nSCHISM_hgrid_node"], vtmp)))'.format(vari+'Y_'+quanti))							
						else:
							exec('utmp=xr.Dataset(data_vars=dict({:s}_{:s}=(["nSCHISM_hgrid_node"], utmp)))'.format(vari+'Y_'+quanti,label[irun]))
							exec('vtmp=xr.Dataset(data_vars=dict({:s}_{:s}=(["nSCHISM_hgrid_node"], vtmp)))'.format(vari+'Y_'+quanti,label[irun]))		
							
				else:	# not quantile	
					if nrun==1:						
						exec('results[expi][vari][quanti]=tmp=xr.Dataset(data_vars=dict({:s}=(["nSCHISM_hgrid_node"], np.{:s}(data,axis=0))))'.format(name,quanti))
					else:		
						exec('results[expi]["{:s}_{:s}"][quanti]=tmp=xr.Dataset(data_vars=dict({:s}_{:s}_{:s}=(["nSCHISM_hgrid_node"], np.{:s}(data,axis=0))))'.format(vari,label[irun],quanti,expi+'_'+vari,label[irun],quanti))						
					if isvec: # macht das sinn fuer quantile?
						if nrun==1:
							exec('utmp=xr.Dataset(data_vars=dict({:s}=(["nSCHISM_hgrid_node"], np.{:s}(u,axis=0))))'.format(quanti+expi+'_'+vari+'X_',quanti))
							exec('vtmp=xr.Dataset(data_vars=dict({:s}=(["nSCHISM_hgrid_node"], np.{:s}(v,axis=0))))'.format(quanti+expi+'_'+vari+'Y_',quanti))							
						else:
							exec('utmp=xr.Dataset(data_vars=dict({:s}_{:s}=(["nSCHISM_hgrid_node"], np.{:s}(u,axis=0))))'.format(quanti+expi+'_'+vari+'Y_',label[irun],quanti))
							exec('vtmp=xr.Dataset(data_vars=dict({:s}_{:s}=(["nSCHISM_hgrid_node"], np.{:s}(v,axis=0))))'.format(quanti+expi+'_'+vari+'Y_',label[irun],quanti))	
						
				# write output		
				tmp.to_netcdf(filenameout,mode=mode)
				mode='a'	
				if isvec:
					utmp.to_netcdf(filenameout,mode=mode)
					vtmp.to_netcdf(filenameout,mode=mode)
					utmp.close()
					vtmp.close()
				tmp.close()

