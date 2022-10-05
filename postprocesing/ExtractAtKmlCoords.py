"""
Extract time Series from SCHISM netCDF4 output
for locations specified within google earth kml files.
Put Time Series output to netcdf file
"""

#area trends
import os
import netCDF4
import sys
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
import matplotlib
#matplotlib.use('AGG')
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
from schism import *
import datetime as dt
import numpy as np	
import datetime as dt
from netCDF4 import Dataset
from matplotlib import pyplot as plt	
from netCDF4 import Dataset
from cftime import utime

# calculate spatial mean as function of time for defined polygons (kml from google earth)

########## settings #################################
# export OMP_NUM_THREADS=4 # before call


setupdir='/work/gg0028/g260114/RUNS/Europe/c1/flow_tweak/'
ncdir=setupdir+'/combined/'


outname='EuropeBasinExtracts.nc' # oztputbame fir betcdfukes

## location files

# local points as google Earth placemark kml
placemarkfile='MontioringCoordinates.kml'
#placemarkfile='' to turn off

# regions for basin average specified as google earth polygon kml
kmlfiles=['NorthAgeanSeaBSdomain.kml']
#kmlfiles=[] to turn off

# Eintire Basin as by name of official ocean areas within OceanDB.kml
ocnfile='/work/gg0028/g260114/RUNS/postproc/water_mass_balance/OceanDB.kml'  #ocean regions kml 
basins=['Mediterranean Sea','Sea of Marmara','Black Sea','Baltic Sea']
#basins=['Sea of Marmara','Black Sea']
# basins=[] # to turn off

# time otput
nspool = 216 #960 
ihfskip = 4320#9600 
nstep=np.int(ihfskip/nspool)

varnames=['elev','wind_speed','air_pressure','salt','temp']
#####################################################


def load_placemarks(file):
	kmlcoords=[]
	kmlnames=[]
	append=False
	with open(file) as f:
		for line in f.readlines():
			if '</Placemark>' in line:
				append=False
			elif '<Placemark>' in line:
				append=True
			if append: 
				if 'name' in line:
					kmlnames.append(line[line.index('>')+1:line.index('</')])
				elif 'coordinates' in line:
					kmlcoords.append(line[line.index('>')+1:line.index('</')].split(',')[:2])
	return {key:coord   for key,coord in zip(np.asarray(kmlnames),np.asarray(kmlcoords))}

def google_kml_coords(kmlfiles):
	"""
		extract coordinates from google earth polygon kml file
		as array of x,y tuples
	"""
	
	seas={}
	for kmlfile in kmlfiles:
		name=[kmlfile[:kmlfile.index('.')]][0]
		with open(kmlfile,'r') as infile:
			prevline=[]
			for line in infile:
				if '<coordinates>' in prevline:
					temp=line.replace('\t','').replace(' \n','').split(' ')
					coords=[]
					for coord in temp:
						coords.append(tuple(np.double(coord.split(',')[:2])))
					seas[name]=np.asarray(coords)
				prevline=line
	return seas

def kml_coords(kmlfile):
	"""
        extract coordinates from google earth polygon kml file
        as array of x,y tuples
	"""
	seas={}
	start = '<name>'
	end = '<\name>'
	string=''
	read_cords=False
	name=''
	
	with open(kmlfile,'r') as infile:
		prevline=[]
		for line in infile:
			
			if ('name' in line) and ('Events' not in line) and ('<Folder>' in prevline):
				name=line[line.find(start)+len(start):line.rfind(end)-len(end)-1]
			if ('<coordinates>' in prevline) and ('<Placemark>' not in prevline): #
				read_cords=True
				coords=[]
			if  ('</coordinates>' in line) and ('<Placemark>' not in line):
				seas[name]=np.asarray(coords)
				read_cords=False
			if read_cords:
				coords.append( (float(line.split(',')[0]),float(line.split(',')[1]) ) )
			prevline=line
	
	return seas
	
def extrac_basin_averaged_timeseries(hxarr,faces,A,eleminds,regname='Baltic Sea',varnames='zcor',nt=np.inf):
	""" extract spatial mean time series for define region of SCHISM grid using xarray acces """
	trinodes=faces[eleminds[regname]]
	nodes=np.unique(trinodes)
	trinodes_reduced=np.zeros(trinodes.shape,int)
	for i,node in enumerate(nodes):
		node
		trinodes_reduced[np.where(trinodes==node)]=i
	selection=hxarr.sel(nSCHISM_hgrid_node=nodes)
	
	if nt==np.inf:
		nt=len(selection['time'])
	t=selection['time'][:nt]
	
	if type(varnames) != list:
		varnames=[varnames,]
	basin_ts=dict.fromkeys(varnames)	
	for varname in varnames:
		if len(hxarr[varname].shape)==2:
			basin_ts[varname]=(np.asarray(selection[varname][:nt,:])[:,trinodes_reduced].mean(axis=2)*A[eleminds[regname]]).sum(axis=1)/A[eleminds[regname]].sum()
		else:
			basin_ts[varname]=(np.asarray(selection[varname][:nt,:,-1])[:,trinodes_reduced].mean(axis=2)*A[eleminds[regname]]).sum(axis=1)/A[eleminds[regname]].sum()
	
	return t,basin_ts

def datetime64_to_datetime(t):
	if len(t)==1:
		t=[t]
	return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00Z'))/ np.timedelta64(1, 's')) for ti in t])

#cwd=os.getcwd()
# modelgrid	
s=schism_setup()	
faces=s.nvplt
x=np.asarray(s.lon)
y=np.asarray(s.lat) 
# element centers
cx=np.mean(x[faces],axis=1)
cy=np.mean(y[faces],axis=1)
elcoords=[[cx[i],cy[i]] for i in range(len(cx))] # pooint pairs of nodes
elem_nn_tree = cKDTree(elcoords) # next neighbour search tree	     

# get element indices for basins
seas=kml_coords(ocnfile)	
if len(kmlfiles)>0:
	regions=google_kml_coords(kmlfiles)	
	key=list(regions.keys())
	seas={**seas,**regions}

################   net cdf access   ########################
## check for schout_nc files
schismfiles=[]
for iorder in range(6):
    iorder
    schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
# regular schismfiles for acces
nrs=[]  
for file in schismfiles:
    nr=int(file[file.rfind('_')+1:file.index('.nc')])
    nrs.append(nr)
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
nrs=list(np.asarray(nrs)[np.argsort(nrs)])

# initialize file acces
ds={}
ds['schism']=xr.open_mfdataset(schismfiles)
s.hgridnodes=np.arange(s.nnodes)

# put to nc
ab=datetime64_to_datetime([ds['schism']['time'][0].values])[0]
basedate=(ab.year,ab.month,ab.day,0,0,0)

nt=len(ds['schism']['time'])	
## load placemarks:
if len(placemarkfile)>0:
	placemarks=load_placemarks(placemarkfile)

	basin_ts=dict.fromkeys(placemarks.keys())
	nn=dict.fromkeys(placemarks.keys())
	for key in basin_ts.keys():
		basin_ts[key]={varname:[] for varname in varnames}
		nn[key]=s.find_nearest_node(placemarks[key][0],placemarks[key][1])-1

	ti=0
	while(ti<nt):
		for varname in varnames:
			if len(ds['schism'][varname].shape)==2:
				data_exert=ds['schism'][varname][ti:ti+nstep,s.hgridnodes].values
			elif ds['schism'][varname].shape == ds['schism']['zcor'].shape:	
				data_exert=ds['schism'][varname][ti:ti+nstep,s.hgridnodes,-1].values
			elif (len(ds['schism'][varname].shape)==3) & (ds['schism'][varname].shape[-1]==2):		#wind
				data_exert=ds['schism'][varname][ti:ti+nstep,s.hgridnodes,:].values
				
			for regname in basin_ts.keys():		
				trinodes=nn[regname]#faces[eleminds[regname]]
				basin_ts[regname][varname].append(data_exert[:,trinodes])
		if ti%(nstep*100)==0:
			print(ti/nstep*100)
		ti+=nstep	
		
	#basin_ts[regname][varname]=np.concatenate(basin_ts[regname][varname])		
	for key in basin_ts.keys():
			for varname in varnames:
				try:
					basin_ts[key][varname]=np.concatenate(basin_ts[key][varname])	
				except:
					pass

	nc = Dataset(outname,'w', format='NETCDF4') #'w' stands for write
	nc.createDimension('time',None)
	nc.createDimension('component',2)
	tv = nc.createVariable('time','f8',('time'))
	tv.long_name = 'Time'
	tv.standard_name = 'time'
	tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
	tv.base_date =  list(basedate)
	tout=datetime64_to_datetime(ds['schism']['time'][:])
	tv[:]=np.asarray([(ti-tout[0]).total_seconds()/86400 for ti in tout])
	for key in basin_ts.keys():
		for varname in varnames:
			outvarname=key.replace(' ','')+'_local_'+varname
			if basin_ts[key][varname].shape[-1]==2:
				tv = nc.createVariable(outvarname,'f8',('time','component'))
			else:
				tv = nc.createVariable(outvarname,'f8',('time'))
			tv[:]=basin_ts[key][varname]
	
### load polygons for basin average #quad
if len(basins)>0:
	A=[]
	for i in range(s.nvplt.shape[0]):
		nodes=s.nvplt[i,:]+1
		A.append(s.proj_area(nodes))
	A=np.asarray(A)
	Warea={regname:np.tile(A[eleminds[regname]],(nstep,1))/A[eleminds[regname]].sum() for regname in basins}	

	# extract basin averaged time series
	eleminds={}
	for ifile,tag in enumerate(basins):
		areaPoly=path.Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
		eleminds[tag]=areaPoly.contains_points(elcoords)
		
	# extract mean for region
	basin_ts2=dict.fromkeys(basins)
	for key in basin_ts2.keys():
		basin_ts2[key]={varname:[] for varname in varnames}
	ti=0
	while(ti<nt):
		for varname in varnames:
			if len(ds['schism'][varname].shape)==2:
				data_exert=ds['schism'][varname][ti:ti+nstep,:].values
			else:	
				data_exert=ds['schism'][varname][ti:ti+nstep,:,-1].values
				
			for regname in basin_ts2.keys():		
				trinodes=faces[eleminds[regname]]
				basin_ts2[regname][varname].append((data_exert[:,trinodes].mean(axis=2)*Warea[regname]).sum(axis=1))
		if ti%(nstep*10)==0:
			print((ti/nstep)*100)
		ti+=nstep

	for key in basin_ts2.keys():
			for varname in varnames:
				try:
					basin_ts2[key][varname]=np.concatenate(basin_ts2[key][varname])	
				except:
					pass		
	# convert schism time	
	t=ds['schism']['time'][:].values
	dates=datetime64_to_datetime(t)

	# put to nc
	if len(placemarkfile)==0:
		nc = Dataset(outname,'w', format='NETCDF4') #'w' stands for write
		nc.createDimension('time',None)
		tv = nc.createVariable('time','f8',('time'))
		tv.long_name = 'Time'
		tv.standard_name = 'time'
		tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
		tv.base_date =  list(basedate)
		tout=datetime64_to_datetime(ds['schism']['time'][:])
		tv[:]=np.asarray([(ti-tout[0]).total_seconds()/86400 for ti in tout])

	for key in basin_ts2.keys():
		for varname in varnames:
			outvarname=key.replace(' ','')+'_mean_'+varname
			tv = nc.createVariable(outvarname,'f8',('time'))
			tv[:]=basin_ts2[key][varname]
nc.close()	