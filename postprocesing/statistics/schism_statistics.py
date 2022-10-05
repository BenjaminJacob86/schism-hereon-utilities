
def calc_stats(rundir='./',ncdir='./outputs',varnames=['hvel','wind_speed','zcor','salt','temp'],lvl=-1,periods=['year']):
	""" Compute variable statistics min,max,mean, and variance grouped by temporal periods for variables indicated by list: varname. For variables with depth dimension staistics are calculate for sigma level lvl (-1:=surface)
	returns values:
	- stats: statistics dictionary by variable example:
		{'elev': {'min': {'year': None},
		'max': {'year': None},
		'mean': array([[ 0.07768171,  0.07824935,  0.07704046, ..., -0.13147822,
			  -0.14230849, -0.13129728]], dtype=float32),
		'var': array([[0.03823537, 0.03888669, 0.0381129 , ..., 0.03570636, 0.03540066,
			  0.03722837]], dtype=float32)},
		'years': array([2018])}
	- countdry: numpy array counting timesteps an element is dry for masking
	"""

	import os
	import netCDF4
	import sys
	from netCDF4 import Dataset
	import numpy as np
	from matplotlib import path
	from matplotlib import pyplot as plt
	import datetime as dt
	import glob
	from scipy.spatial import cKDTree
	import xarray as xr
	from matplotlib.path import Path
	from matplotlib.colors import LinearSegmentedColormap
	import pickle
	import datetime as dt
	sys.path.insert(0,'/gpfs/home/jacobb/code/python/')   #strand
	sys.path.insert(0,'/mnt/lustre01/pf/g/g260114/Programs/python/scripts/')  #mistral
	sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/')   #vims
	from schism import *



	os.chdir(rundir)
	os.chdir(rundir)
	s=schism_setup()
	
	## check for schout_nc files
	schismfiles=[]
	for iorder in range(6):
		schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
	nrs=[]  
	for file in schismfiles:
		nr=int(file[file.rfind('_')+1:file.index('.nc')])
		nrs.append(nr)
	schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
	nrs=list(np.asarray(nrs)[np.argsort(nrs)])
	ds=xr.open_mfdataset(schismfiles)


	t=ds['time'][:]
	basedate=t[0].base_date.split(' ')
	while '' in basedate:
		basedate.remove('')
	refdate=[]	
	for item in basedate:
			item=np.int(np.float(item))
			refdate.append(item)
	reftime=dt.datetime(refdate[0],refdate[1],refdate[2],refdate[3],refdate[4])
	dates=reftime+dt.timedelta(hours=1)+np.arange(len(t))*dt.timedelta(seconds= np.float((t[1]-t[0]))/1e9)

	
	yearly=ds.groupby('time.'+periods[0])
	means=yearly.mean('time')
	varname='dahv'
	monthly=ds.groupby('time.'+'month')
	mmeans=monthly.mean('time')
	
	
	nt,nnodes,nvert=ds['zcor'].shape
	
	stats=dict.fromkeys(varnames)
	quantity=['min','max','mean','var']
	for key in stats.keys():
		stats[key]=dict.fromkeys(quantity)
	for key2 in quantity:
		stats[key][key2]=dict.fromkeys(periods)

	yearly=ds.groupby('time.'+periods[0])
	means=yearly.mean('time')
	vars=yearly.var('time')
	mins=yearly.min('time')
	maxs=yearly.max('time')
	#medians=yearly.median('time')
	stats['years']=means.year.values
	
	for varname in varnames:
		if means[varname].shape==(len(stats['years']),nnodes,nvert):
			stats[varname]['mean']=means[varname][:,:,lvl].values
			#stats[varname]['median']=np.asarray(medians[varname][:,:,lvl])
			stats[varname]['var']=np.asarray(vars[varname][:,:,lvl])
			stats[varname]['min']=np.asarray(mins[varname][:,:,lvl])
			stats[varname]['max']=np.asarray(maxs[varname][:,:,lvl])
		elif means[varname].shape==(len(stats['years']),nnodes,nvert,2):		
			stats[varname]['mean']=np.asarray(means[varname][:,:,lvl,:])
			#stats[varname]['median']=np.asarray(medians[varname][:,:,lvl,:])
			stats[varname]['var']=np.asarray(vars[varname][:,:,lvl,:])
			stats[varname]['min']=np.asarray(mins[varname][:,:,lvl,:])
			stats[varname]['max']=np.asarray(maxs[varname][:,:,lvl,:])
			# from magnitude
			abs=np.sqrt(ds[varname][:,:,lvl,0]**2+ds[varname][:,:,lvl,1]**2)
			abs.groupby('time.'+periods[0])
			stats[varname]['magnitude_mean']=abs.mean('time').values
			#stats[varname]['magnitude_median']=abs.median('time').values
			stats[varname]['magnitude_var']=abs.var('time').values
			stats[varname]['magnitude_min']=abs.min('time').values
			stats[varname]['magnitude_max']=abs.max('time').values
			
		elif means[varname].shape==(len(stats['years']),nnodes,2):		
			stats[varname]['mean']=np.asarray(means[varname][:])
			stats[varname]['var']=np.asarray(vars[varname][:])
			# from magnitude
			abs=np.sqrt(ds[varname][:,:,0]**2+ds[varname][:,:,1]**2)
			abs=abs.groupby('time.'+periods[0])
			stats[varname]['magnitude_mean']=abs.mean('time').values
			#stats[varname]['magnitude_median']=abs.median('time').values
			stats[varname]['magnitude_var']=abs.var('time').values
			stats[varname]['magnitude_min']=abs.min('time').values
			stats[varname]['magnitude_max']=abs.max('time').values

		elif means[varname].shape==(len(stats['years']),nnodes):		
			stats[varname]['mean']=np.asarray(means[varname][:])
			stats[varname]['var']=np.asarray(vars[varname][:])


	countdry=ds['wetdry_elem'].sum(axis=0).values
	#s.plotAtnodes(np.sqrt(stats['elev']['var'][0,:]),mask=countdry)
	stats['setup']=s	 
	return stats,countdry