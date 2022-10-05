import os
import netCDF4
import sys
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities/')
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/')
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/Lib/')
import matplotlib
matplotlib.use('AGG')
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
import xlwt
#import cftime
#import netcdftime
#from netcdftime import utime
from cftime import utime
########## settings #################################
# export OMP_NUM_THREADS=1 #

run_name='BlackSea24d_strand'            # Name of Excelfile Carrying bilances

# directories foor RUN and atmospheric forcing
setupdir='/gpfs/work/jacobb/data/RUNS/BlackSea/RUN24d' #'/sciclone/home10/yinglong/DISKS/vims20/BlackSea/RUN25b/'
sfluxdir=setupdir+'/sflux/'
sflux_in=setupdir+'/sflux/'

									    # https://www.eea.europa.eu/data-and-maps/data/europe-seas
cwd=os.getcwd()
ncdir=setupdir+'combined/'

#fluxfile='/sciclone/home10/yinglong/DISKS/vims20/BlackSea/RUN25b/outputs/flux.out' # new schism
fluxfile=ncdir+'flux.out'

outfolder='/gpfs/work/jacobb/data/RUNS/BlackSea/RUN24d/analysis'+run_name+'/'
kmlfile='/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/postprocesing/balances//OceanDB.kml'  #ocean regions kml file

# Areas to extract from boundary database
#areas=['Baltic Sea','Kattegat','Skagerrak','North Sea','English Channel',
#'Mediterranean Sea','Sea of Marmara',
#'Black Sea', 'Sea of Azov',] 
areas=['Black Sea', 'Sea of Azov',] 
########
years=range(2008,2008+1)     # years considered
reftime=dt.datetime(2008,6,1,0,0,0)#+dt.timedelta(days=214) # fake shift to end to have a fullyear
prec_years=np.arange(2008,2008+1)
evap_years=np.arange(2008,2008+1)

# analysis to be done for full years
# if not to full years start in summer build artifical year


use_atmo_from_schism=True   #True: get prec / evaporation from schism output

							#False: get from atmo
do_plot=True
# name of volume transport calculation sections in order of fluxflag.prop

locnames=['Asov','Bosporus N','Bosporus C','Bosporus S','Daradanelle N','Daradanelle C','Daradanelle S']
#locnames=['Asov','Bosporus N','Bosporus C','Bosporus S','Daradanelle N','Daradanelle C','Daradanelle S','Bosporus','Gibraltar','little Belt','Great Belt','Sound']

locdir=np.zeros(len(locnames))
locdir[:]=1
BasinStraits=dict.fromkeys(areas)
#BasinStraits['Mediterranean Sea']=[('Gibraltar',1)]   #+1 :flows ind, -1: flows pout (higher to lower number in fluxflag.prop)
BasinStraits['Black Sea']=[('Bosporus',-1)]
BasinStraits['Sea of Azov']=[('Asov',-1)]
BasinStraits['Sea of Marmara']=[('Bosporus',-1),('Daradanelle C',-1)]
#BasinStraits['Baltic Sea']=[('little Belt',1),('Great Belt',1),('Sound',1)]
#BasinStraits['Kattegat']=[('little Belt',-1),('Great Belt',-1),('Sound',-1)]
#BasinStraits['English CHannel']=[('NorthSea',-1)]
BasinStraits['ModelDomain']=None

# conversion factors
m3tokm3=1/1e9
rho=1000.0 # kg/m^3/
in_m_s=1.0/(rho)
in_mm_day=in_m_s*1000*86400

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
				name=line[line.find(start)+len(start):line.rfind(end)-len(end)-1]# voerher minus 1 warum geht nicht mehr? zaehlt jetzt eol vielleicht
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

###### load rivers 
def get_annual_mean_river_runoff(s,areas,years,startdate=dt.datetime(2007,6,1)):
	""" 	Extract river Data from flux.th and source files 	using schism class and domain polygons 	returns combined river inflow per basin averaged by year in m3/s
	uses source class.
	"""
	
	river_inflow={}
	######### load flux th	###################################
	if os.path.exists('flux.th'):
		q=np.loadtxt('flux.th')
		qt=q[:,0]
		qdates=startdate+np.arange(len(qt))*dt.timedelta(seconds= np.diff(q[:2,0])[0])
		Ryears=np.asarray([date.year for date in qdates])

		# associate river boundary node	with ocean region
		lon=np.asarray(s.lon)
		lat=np.asarray(s.lat)
		bdcenters=[]
		for bdy in s.bdy_segments:
				bdcenters.append( (np.mean(lon[np.asarray(bdy)-1]),np.mean(lat[np.asarray(bdy)-1])) )
		bdcenters=np.asarray(bdcenters)


		obdtype=[]
		bd_rnames=[]
		line='dummy'
		with open('bctides.in','r') as f:
			while len(line)>0:
				line=f.readline()
				if '4 4 4 4' in line:
					continue
				if 'river' in line:
					obdtype.append('river')
					line2=line.replace(' ','').replace('\n','')
					bd_rnames.append(line2[line2.rindex(':')+1:])
				elif 'ocean' in line:
					obdtype.append('ocean')
					bd_rnames.append('ocean')
				elif ('!' in line) and ('relax' not in line) :	
					line2=line.replace('\n','')
					bd_rnames.append(line2[line2.index('!')+1:])
					
		# in which basin are openbds ?
		basinbdrivers={}
		flux_indices={}
		for ifile,basin in enumerate(areas):
				areaPoly=Path(list(zip(seas[basin][:,0],seas[basin][:,1])))
				basinbdrivers[basin]=np.asarray(np.where(areaPoly.contains_points(bdcenters)))
				# bnd
				if len(basinbdrivers[basin][0]) > 0:
					flux_indices[basin]=basinbdrivers[basin]-obdtype[:basinbdrivers[basin][0][0]].count('ocean')+1 # one for time
				else:
					flux_indices[basin]=np.nan
	else:
		bd_rnames=[]
		Ryears=[]
					
	if os.path.exists('source_sink.in'):					
		src=sources(s=s,fname='source_sink.in',comments='!')	
		src.dates=startdate+np.arange(len(src.vtime))*dt.timedelta(seconds=np.float(src.vtime[1]-src.vtime[0]))
		src.years=np.asarray([date.year for date in src.dates])
		
		if len(Ryears)==0:
			Ryears=src.years

		# in which basin are source elements?
		src_indices={}
		for ifile,basin in enumerate(areas):
				areaPoly=Path(list(zip(seas[basin][:,0],seas[basin][:,1])))
				inds=np.asarray(np.where(areaPoly.contains_points(list(zip(src.x,src.y)))))
				if len(inds[0]) > 0:
					src_indices[basin]=inds[0]
				else:
					src_indices[basin]=[]
					
	else:
		src_indices={}	
		for ifile,basin in enumerate(areas):
			src_indices[basin]=[]
			src=None
			
	for basin in areas:
			if os.path.exists('flux.th'):
				from_bdriver=-np.mean(q[Ryears==years[0],:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
			else: 	
				from_bdriver=[]
			if len(src_indices[basin])> 0: # add from src
				from_srcriver=np.mean(src.vsource[src.years==years[0],:][:,src_indices[basin]] ,axis=0) #dunno bout bracket when not interactive
				river_inflow[basin]=np.hstack((from_bdriver, from_srcriver))
			else:
				river_inflow[basin]=from_bdriver
				
			for year in years[1:]:
				if os.path.exists('flux.th'):
					from_bdriver=-np.mean(q[Ryears==year,:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
				else:
					from_bdriver=[]				
				if len(src_indices[basin])> 0: # add from src
					from_srcriver=np.mean(src.vsource[src.years==year,:][:,src_indices[basin]] ,axis=0) #dunno bout bracket when not interactive				#from_srcriver=np.mean(src.vsource[src.years==year,:][src_indices[basin]] ,axis=1)[0]
					river_inflow_temp=np.hstack((from_bdriver, from_srcriver))
				else:
					river_inflow_temp=from_bdriver
				river_inflow[basin]=np.vstack((river_inflow[basin],river_inflow_temp))
			river_inflow[basin]=np.asarray(river_inflow[basin])
	print('done reading annual mean river runoff by basin')
	return river_inflow,Ryears,bd_rnames,src,src_indices	

	"""
	Extract river Data from flux.th and source files
	using schism class and domain polygons
	returns combined river inflow per basin averaged by year
	uses source class
	"""
	river_inflow={}
	######### load flux th	###################################
	q=np.loadtxt('flux.th')
	qt=q[:,0]
	qdates=startdate+np.arange(len(qt))*dt.timedelta(seconds= np.diff(q[:2,0])[0])
	Ryears=np.asarray([date.year for date in qdates])

	# associate river boundary node	with ocean region
	lon=np.asarray(s.x)
	lat=np.asarray(s.y)
	bdcenters=[]
	for bdy in s.bdy_segments:
			bdcenters.append( (np.mean(lon[np.asarray(bdy)-1]),np.mean(lat[np.asarray(bdy)-1])) )
	bdcenters=np.asarray(bdcenters)


	obdtype=[]
	bd_rnames=[]
	line='dummy'
	with open('bctides.in','r') as f:
		while len(line)>0:
			#try:
			line=f.readline()
			#except:
			#	pass
			#for line in f.readlines():
			#print(line)
			if 'river' in line:
				obdtype.append('river')
				line2=line.replace(' ','').replace('\n','')
				bd_rnames.append(line2[line2.rindex(':')+1:])
			elif 'ocean' in line:
				obdtype.append('ocean')
				bd_rnames.append('ocean')

			
	src=sources(s=s,fname='source_sink.in',comments='!')	
	src.dates=startdate+np.arange(len(src.vtime))*dt.timedelta(seconds=np.float(src.vtime[1]-src.vtime[0]))
	src.years=np.asarray([date.year for date in src.dates])

	# in which basin are openbds and source elements?
	basinbdrivers={}
	flux_indices={}
	src_indices={}
	for ifile,basin in enumerate(areas):
			areaPoly=Path(list(zip(seas[basin][:,0],seas[basin][:,1])))
			basinbdrivers[basin]=np.asarray(np.where(areaPoly.contains_points(bdcenters)))
			# bnd
			if len(basinbdrivers[basin][0]) > 0:
				flux_indices[basin]=basinbdrivers[basin]-obdtype[:basinbdrivers[basin][0][0]].count('ocean')+1 # one for time
			else:
				flux_indices[basin]=np.nan
			#src
			inds=np.asarray(np.where(areaPoly.contains_points(list(zip(src.x,src.y)))))
			if len(inds[0]) > 0:
				src_indices[basin]=inds
			else:
				src_indices[basin]=[]

	for basin in areas:
			from_bdriver=-np.mean(q[Ryears==years[0],:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
			if len(src_indices[basin])> 0: # add from src
				from_srcriver=np.mean(src.vsource[src.years==years[0],:][src_indices[basin]] ,axis=1)[0]
				river_inflow[basin]=np.hstack((from_bdriver, from_srcriver))
			else:
				river_inflow[basin]=from_bdriver
				
			for year in years[1:]:
				from_bdriver=-np.mean(q[Ryears==year,:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
				if len(src_indices[basin])> 0: # add from src
					from_srcriver=np.mean(src.vsource[src.years==year,:][src_indices[basin]] ,axis=1)[0]
					river_inflow_temp=np.hstack((from_bdriver, from_srcriver))
				else:
					river_inflow_temp=from_bdriver
				river_inflow[basin]=np.vstack((river_inflow[basin],river_inflow_temp))
			river_inflow[basin]=np.asarray(river_inflow[basin])
	print('done reading annual mean river runoff by basin')
	return river_inflow,Ryears,bd_rnames,src		
############################################################


if not (os.path.exists(outfolder)):
        os.mkdir(outfolder)


############################# load ocean areas ############################################
seas=kml_coords(kmlfile)
for key in areas:
	cx,cy=np.mean(seas[key],axis=0)
##################################	


####### load era data  ############################
if not use_atmo_from_schism:
	#files=[]
	#for year in range(20012,2013+1):
	#		year
	#		files+=glob.glob(sflux_in+'*'+str(year)+'*')
	#files=np.sort(files)	
	#evap_files=list(np.sort(glob.glob('/gpfs/work/jacobb/data/DATA/DOWNL_ERA5/*evap*nc')))
	#ls_files=list(np.sort(glob.glob('/gpfs/work/jacobb/data/DATA/DOWNL_ERA5/*lsheat*nc')))
	#netfiles=list(np.sort(glob.glob(sfluxdir+'*.nc')))

	#nc=Dataset(files[0])
	#era_lon=nc['longitude'][:]
	#era_lat=nc['latitude'][:]
	#era_lon,era_lat=np.meshgrid(nc['longitude'][:],nc['latitude'][:])
	#isocean=nc['sst'][0,:].mask==False	
	#t=nc['time'][:]
	#timestep=np.double(np.diff(t[0:2]))*3600.0
	#era_nn_tree = cKDTree(list(zip(era_lon.flatten(),era_lat.flatten())))
	##########################################################

	######## get element indices for basins and file acces fpr era data  ########################
	airfiles=list(np.sort(glob.glob(sfluxdir+'*air*.nc')))
	nc_era=Dataset(airfiles[0])
	era_lon=nc_era['lon'][:]
	era_lat=nc_era['lat'][:]
	era_nn_tree = cKDTree(list(zip(era_lon.flatten(),era_lat.flatten())))
	#isocean=nc_era['stmp'][0,:].mask==False	
	isocean=np.ones(nc_era['stmp'][0,:].shape,bool)
	#plt.pcolormesh(era_lon,era_lat,isocean)

	 ######## get element indices for basins and file acces fpr era data  ########################
	prcfiles=list(np.sort(glob.glob(sfluxdir+'*prc*.nc')))
	#nc_era=Dataset(prcfiles[0])
	#era_lon=nc_era['lon'][:]
	#era_lat=nc_era['lat'][:]
	#isocean=nc_era['stmp'][0,:].mask==False	
	#plt.pcolormesh(era_lon,era_lat,isocean)
	
	
	era_basins={}
	era_basin_ii_jj={}
	for ifile,tag in enumerate(areas):
			areaPoly=path.Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
			era_basins[tag]=areaPoly.contains_points(list(zip(era_lon.flatten(),era_lat.flatten()))) & isocean.flatten()
			ii,jj=np.unravel_index(np.where(era_basins[tag]),era_lon.shape)
			era_basin_ii_jj[tag]=(ii,jj)
#########################################################################

os.chdir(setupdir)
######### load setup #############################
s=schism_setup()
if do_plot:
	s.plot_domain_boundaries(append=1)

# extent_total_ares
xmin,ymin=np.min(s.lon),np.min(s.lat)
xmax,ymax=np.max(s.lon),np.max(s.lat)
seas['ModelDomain']=np.asarray(([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin])).T
areas+=['ModelDomain']
	
#quad
A=[]
for i in range(s.nvplt.shape[0]):
	nodes=s.nvplt[i,:]+1
	A.append(s.proj_area(nodes))
A=np.asarray(A)

#faces= np.asarray(s.nv)-1
faces=s.nvplt
x=np.asarray(s.lon)
y=np.asarray(s.lat) 

# element centers
cx=np.mean(x[faces],axis=1)
cy=np.mean(y[faces],axis=1)
elcoords=[[cx[i],cy[i]] for i in range(len(cx))] # pooint pairs of nodes
elem_nn_tree = cKDTree(elcoords) # next neighbour search tree	     

# get element indices for basins
eleminds={}
for ifile,tag in enumerate(areas):
	areaPoly=path.Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
	eleminds[tag]=areaPoly.contains_points(elcoords)
	

scale2year=86400*365  #m3/a


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
prcfiles=list(np.sort(glob.glob(sfluxdir+'*prc*.nc')))
if not use_atmo_from_schism:
	#ds['pure']=xr.open_mfdataset(files[5:])
	ds['pure']=xr.open_mfdataset(prcfiles)
	#ds['evap']=xr.open_mfdataset(evap_files)
	#ds['ls']=xr.open_mfdataset(ls_files)
	#ds['net']=xr.open_mfdataset(netfiles)
#########################################
volumina={}
basinareas={}
dmean={}
for area in areas:
	basinareas[area]=A[eleminds[area]].sum()
	dmean[area]=(A[eleminds[area]]/A[eleminds[area]].sum()*np.asarray(ds['schism']['depth'][0,:])[faces[eleminds[area]]].mean(axis=1)).sum()
	volumina[area]=dmean[area]*basinareas[area]


# check die konvertierung nicht workaround
tschism=ds['schism']['time'][:]
nc_temporary=netCDF4.Dataset(schismfiles[0])
ut=utime(nc_temporary['time'].units)
#ut=netcdftime.utime(nc_temporary['time'].units)
date0=ut.num2date(nc_temporary['time'][0])
nc_temporary.close()
# now works with offset but not gaps
schismdates=date0+np.arange(len(tschism))*dt.timedelta(seconds= np.float((tschism[1]-tschism[0]))/1e9)
schismyears=np.asarray([date.year for date in schismdates])

rho_water=1000 # kg/m3	
if not use_atmo_from_schism:
	# radiation precipitation
	year_means=ds['pure'].groupby('time.year').mean('time')
	means={}
	varname='prate'
	means[varname]=np.asarray(year_means[varname])/rho_water # / rhow_water when from air files#/timestep  #m/s
	#for varname in 'tp',:
	#	means[varname]=np.asarray(year_means[varname])/timestep  #m/s
	# evaporation
	#year_means=ds['evap'].groupby('time.year').mean('time')
	#for varname in 'e','es':
	#	varname
	#	means[varname]=np.asarray(year_means[varname]/timestep)  #m/s
	#schism_vars=['zcor']	
	#pickle.dump( means, open( run_name+"means_Pickle", "wb" ) )	
else:
	prec_years=np.unique(schismyears)
	evap_years=np.unique(schismyears)
	schism_vars=['zcor','precipitation','evaporation']	

schism_vars+=['solar_radiation','upward_longwave','downward_longwave','total_heat_flux']	
schism_means={}
year_means=ds['schism'].groupby('time.year').mean('time')
for varname in schism_vars[1:]:
	print(' doing ' + varname)
	if varname in  ['precipitation','evaporation']:
		factor=1.0/rho
	else:
		factor=1
	if len(year_means[varname].shape)==3:
		schism_means[varname]=np.asarray(year_means[varname][:,:,-1]*factor)
	else:		
		schism_means[varname]=np.asarray(year_means[varname][:]*factor)
pickle.dump( schism_means, open( run_name+"_schismmeans_Pickle", "wb" ) )
print('done extracting schism annual means')	

schism_means['total_flux']=year_means['solar_radiation']+year_means['total_heat_flux']  # including solar radiation

# plots
years=np.unique(schismyears)
for iyear,year in enumerate(years):
	for varname in ['solar_radiation','upward_longwave','downward_longwave','total_flux']:
		plt.clf()
		ph,ch=s.plotAtnodes(schism_means[varname][iyear])
		ch.set_label(varname + ' [W/m^2]')
		plt.title('Annual Mean '+ str(year))
		plt.tight_layout()
		plt.savefig(outfolder+'Annual_Mean_'+varname+'_'+str(year)+'.png',dpi=400)	

# load river data
years=np.unique(schismyears)
river_inflow,Ryears,bd_rnames,src,src_indices=get_annual_mean_river_runoff(s,areas,years,startdate=reftime)	

heat_by_basin=dict.fromkeys(['solar_radiation','upward_longwave','downward_longwave','total_flux'])
for key in heat_by_basin.keys():
	heat_by_basin[key]=dict.fromkeys(areas)

# spatially mean weighted precipitation
water_from_atmo={} # m3/s
water_to_atmo={} # m3/s
if not use_atmo_from_schism:
	nn2elem_center={}
	for basin in areas:
		nni=era_nn_tree.query(list(zip(cx[eleminds[basin]],cy[eleminds[basin]]) ))[1]
		nn2elem_center[basin]=nni
		water_from_atmo[basin]=[(means['prate'][i,:].flatten()[nn2elem_center[basin]]*A[eleminds[basin]]).sum() for i in range(means['prate'].shape[0]) ] # m^3/a
	#[(means['tp'][i,:].flatten()[nn2elem_center[basin]]*A[eleminds[basin]]/(A[eleminds[basin]].sum())).sum() for i in range(means['tp'].shape[0]) ]
		# water_from_atmo[basin]=[(means['tp'][i,:].flatten()[nn2elem_center[basin]]*A[eleminds[basin]]).sum() for i in range(means['tp'].shape[0]) ] # m^3/a
		#water_to_atmo[basin]=[(means['e'][i,:].flatten()[nn2elem_center[basin]]*A[eleminds[basin]]).sum() for i in range(means['e'].shape[0]) ] # m^3/a
		
		
	pickle.dump( water_from_atmo, open( run_name+"water_from_atmo", "wb" ) )	
	pickle.dump( water_to_atmo, open( run_name+"water_to_atmo", "wb" ) )	

	load=0
	if load:
		water_from_atmo=pickle.load(open( run_name+"water_from_atmo", "rb" ) )	
		water_to_atmo=pickle.load( open( run_name+"water_to_atmo", "rb" ) )	

else:
	
	for basin in areas:
		trimeansP=schism_means['precipitation'][:,faces[eleminds[basin],:]].mean(axis=2) # jahresmittel
		trimeansE=schism_means['evaporation'][:,faces[eleminds[basin],:]].mean(axis=2)
		water_from_atmo[basin]=[(trimeansP[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] # m^3/a
		water_to_atmo[basin]=[(trimeansE[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] # m^3/a

		# heatflux by basin		
		for varname in heat_by_basin.keys():
			if varname=='total_flux':
				trimeans=schism_means[varname].values[:,faces[eleminds[basin],:]].mean(axis=2)
			else:
				trimeans=schism_means[varname][:,faces[eleminds[basin],:]].mean(axis=2)
			heat_by_basin[varname][basin]=[(trimeans[i,:]*A[eleminds[basin]]/A[eleminds[basin]].sum()).sum() for i in range(trimeansP.shape[0]) ] # W/mm/a

				
#######   strait fluxes  ###############################
#
# if compiled without analysis module only total outflow
#
m=np.loadtxt(fluxfile)
nregs=len(m[0,:])-1
ntypes=np.where(np.diff(m[:100,0])>0)[0][0]+1
t1=m[::ntypes,0]
if ntypes==1:
	print('flux out without analysis module contains only total flux')
	flux1={'vol':{'tot':m[::ntypes,1:],'pos':m[1::ntypes,1:]*np.nan,'neg':m[2::ntypes,1:]*np.nan},'salt':{'pos':m[3::ntypes,1:]*np.nan,'neg':m[4::ntypes,1:]*np.nan}}
else:	
	flux1={'vol':{'tot':m[::ntypes,1:],'pos':m[1::ntypes,1:],'neg':m[2::ntypes,1:]},'salt':{'pos':m[3::ntypes,1:],'neg':m[4::ntypes,1:]}}
fluxdates=dt.datetime(2012,6,1)+np.arange(len(t1))*dt.timedelta(days=np.diff(t1[:2])[0])
flux_years=np.asarray([date.year for date in fluxdates])
basin_strais=[]
leg=['tot','pos','neg']

plot_trans=False
if plot_trans & ntypes > 1:
	plt.figure()
	for i in range(nregs):
		plt.subplot(nrows,ncols,i+1)
		plt.plot(t1, flux['vol']['tot'][:,i],'k')              
		plt.plot(t1, flux['vol']['pos'][:,i],'--')              
		plt.plot(t1, flux['vol']['neg'][:,i],'--')              
		plt.ylabel('vol [m^3/s]')
		plt.title(locnames[i])
		plt.legend(leg)
	plt.tight_layout()
	plt.savefig(outfolder+'flux_ts.png',dpi=300)
	plt.close()
########################################################




####### balancing ###############################
#for basin in areas:
#	plt.plot(years,river_inflow[basin].sum(axis=1)*scale2year*m3tokm3)
#plt.legend(areas)

# output excel	
row=0
col0=1
#labell=['Quantity [km3/a] \ Basin:','area [km^2]','precipitation','evaporation','River','Strait inflow','Strait outflow','Strait Flux tot','Netsum','volume 0','volume 1','Dv ']

labell=['Quantity [km3/a] \ Basin:','area [km^2]','precipitation','evaporation','River','Strait inflow','Strait outflow','Strait Flux tot','Netsum','','annual mean heatfluxes [W/m^2]','solar_radiation', 'upward_longwave', 'downward_longwave', 'total_flux']

#years=np.arange(2007,2010+1)
years=np.unique(schismyears)
row_flux_pos=labell.index('Strait inflow')
row_flux_neg=labell.index('Strait outflow')
row_flux_tot=labell.index('Strait Flux tot')


D=np.asarray(ds['schism']['depth'][0,:])

os.chdir(cwd)
wb = xlwt.Workbook()  # Workbook is created 
#years=np.asarray([2012,2013])
schism_unq_years=np.unique(schismyears)
for year in years:
	print('doing year ' + str(year))
	sheet1 = wb.add_sheet('Budget {:d}'.format(year))  # add_sheet is used to create 

	#i0,i1=np.where(schismyears==year)[0][[0,-1]]
	#labell[-3]='volume '+str(schismdates[i0])	
	#labell[-2]='volume '+str(schismdates[i1])	
	#mat=np.zeros((len(labell),len(basin)))

	# write labels	
	for irow,label in enumerate(labell):
		sheet1.write(irow, 0, label)
	for icol,basin in enumerate(areas):
		sheet1.write(row, col0+icol, basin)
	for icol,strait in enumerate(locnames):
		sheet1.write(0, len(areas)+1+icol, strait)	
	irow=labell.index('area [km^2]')		
	for icol,basin in enumerate(areas):
		sheet1.write(irow, col0+icol,basinareas[basin]/1000000)
		
		
	irow=labell.index('River')	
	for icol,basin in enumerate(areas):
		sheet1.write(irow, col0+icol, river_inflow[basin][np.where(years==year),:].sum()*scale2year*m3tokm3)
		#mat[irow,icol]=river_inflow[basin][np.where(years==year),:].sum()*scale2year*m3tokm3
		
	irow=labell.index('precipitation')	
	for icol,basin in enumerate(areas):
		sheet1.write(irow, col0+icol,water_from_atmo[basin][np.where(prec_years==year)[0][0]]*scale2year*m3tokm3) #
		#mat[irow,icol]=water_from_atmo[basin][np.where(prec_years==year)[0][0]]*scale2year*m3tokm3
	irow=labell.index('evaporation')	
	for icol,basin in enumerate(areas):
		sheet1.write(irow, col0+icol,-np.abs(water_to_atmo[basin][np.where(evap_years==year)[0][0]])*scale2year*m3tokm3) #*scale2year

	# add strait labels # this does not accpint direction prescriped
	for icol,strait in enumerate(locnames):
		#sheet1.write(0, len(areas)+1+icol, strait)
		if 	locdir[icol]==-1:
			sheet1.write(row_flux_neg, len(areas)+1+icol, -flux1['vol']['pos'][flux_years==year,icol].mean()*scale2year*m3tokm3)	
			sheet1.write(row_flux_pos, len(areas)+1+icol, -flux1['vol']['neg'][flux_years==year,icol].mean()*scale2year*m3tokm3)	
			sheet1.write(row_flux_tot, len(areas)+1+icol, -flux1['vol']['tot'][flux_years==year,icol].mean()*scale2year*m3tokm3)	
		else:
			sheet1.write(row_flux_pos, len(areas)+1+icol, flux1['vol']['pos'][flux_years==year,icol].mean()*scale2year*m3tokm3)	
			sheet1.write(row_flux_neg, len(areas)+1+icol, flux1['vol']['neg'][flux_years==year,icol].mean()*scale2year*m3tokm3)	
			sheet1.write(row_flux_tot, len(areas)+1+icol, flux1['vol']['tot'][flux_years==year,icol].mean()*scale2year*m3tokm3)	
	
	# add straits for basins
	irow=labell.index('Strait Flux tot')	
	for icol,basin in enumerate(areas):
		posflux=0
		negflux=0
		totflux=0
		if not BasinStraits[basin]==None:
			for item in BasinStraits[basin]:
				coli=locnames.index(item[0])
				if item[1]==-1:
					posflux+=-flux1['vol']['neg'][flux_years==year,coli].mean()
					negflux+=-flux1['vol']['pos'][flux_years==year,coli].mean()
					totflux+=-flux1['vol']['tot'][flux_years==year,coli].mean()
				else:
					posflux+=flux1['vol']['pos'][flux_years==year,coli].mean()
					negflux+=flux1['vol']['neg'][flux_years==year,coli].mean()
					totflux+=flux1['vol']['tot'][flux_years==year,coli].mean()
		sheet1.write(row_flux_pos, 1+icol, posflux*scale2year*m3tokm3)	
		sheet1.write(row_flux_neg, 1+icol, negflux*scale2year*m3tokm3)	
		sheet1.write(row_flux_tot, 1+icol, totflux*scale2year*m3tokm3)


		
#if year in schismyears:
#	
#	Z0=np.asarray(ds['schism']['zcor'][i0:(i0+4),:,-1]).mean(axis=0)
#	Z1=np.asarray(ds['schism']['zcor'][i1-3:i1+1,:,-1]).mean(axis=0)
#	V0={}
#	V1={}
#	for area in areas:
#		basinareas[area]=A[eleminds[area]].sum() # tages mittel
#		V0[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z0[faces[eleminds[area]]].mean(axis=1) ).sum()
#		V1[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z1[faces[eleminds[area]]].mean(axis=1) ).sum()
		
		#for icol,basin in enumerate(areas):
		#	sheet1.write(len(labell)-3, col0+icol,V0[basin]*m3tokm3)
		#	sheet1.write(len(labell)-2, col0+icol,V1[basin]*m3tokm3)
		#	sheet1.write(len(labell)-1, col0+icol,(V1[basin]-V0[basin])*m3tokm3)
			
	irow=labell.index('Netsum')
	Net=dict.fromkeys(areas)	
	for icol,basin in enumerate(areas):
		Net[basin]=0
		Net[basin]+=(water_from_atmo[basin][np.where(prec_years==year)[0][0]]-np.abs(water_to_atmo[basin][np.where(prec_years==year)[0][0]]))#*scale2year*m3tokm3)
		Net[basin]+=river_inflow[basin][np.where(years==year),:].sum()#*scale2year*m3tokm3
		totflux=0
		if not BasinStraits[basin]==None:
			for item in BasinStraits[basin]:
				coli=locnames.index(item[0])
				if item[1]==-1:
					totflux+=-flux1['vol']['tot'][flux_years==year,coli].mean()
				else:
					totflux+=flux1['vol']['tot'][flux_years==year,coli].mean()
		Net[basin]+=totflux
		sheet1.write(irow, 1+icol,Net[basin]*scale2year*m3tokm3)

	# heatflux vars
	for varname in heat_by_basin.keys():
		irow=labell.index(varname)		
		for icol,basin in enumerate(areas):
			sheet1.write(irow, col0+icol,heat_by_basin[varname][basin][np.where(schism_unq_years==year)[0][0]])



		
unqyears=np.unique(schismyears)
useyears=np.asarray([(schismyears==year).sum()/4 for year in unqyears])>360  # only accound years having more than 360 days

#useyear_atmo
years_used=unqyears[useyears]

sheet1 = wb.add_sheet('Mean Budget {:d} to {:d}'.format(unqyears[useyears][0],unqyears[useyears][-1]))  # add_sheet is used to create sheet. 

i0,i1=np.where(schismyears==year)[0][[0,-1]]
labell[-3]='volume '+str(schismdates[0])	
labell[-2]='volume '+str(schismdates[-1])	

for irow,label in enumerate(labell):
	sheet1.write(irow, 0, label)

for icol,basin in enumerate(areas):
	sheet1.write(row, col0+icol, basin)

irow=labell.index('River')	
for icol,basin in enumerate(areas):
	sheet1.write(irow, col0+icol, river_inflow[basin][useyears,:].mean(axis=0).sum()*scale2year*m3tokm3)

irow=labell.index('precipitation')	
for icol,basin in enumerate(areas):
	sheet1.write(irow, col0+icol,np.mean(np.asarray(water_from_atmo[basin])[useyears],axis=0)*scale2year*m3tokm3) #
	
irow=labell.index('evaporation')	
for icol,basin in enumerate(areas):
	sheet1.write(irow, col0+icol,np.mean(np.asarray(water_to_atmo[basin])[useyears])*scale2year*m3tokm3) #*scale2year

# add strait labels
iuse=[yeari in years_used for yeari in flux_years]
for icol,strait in enumerate(locnames):
	sheet1.write(0, len(areas)+1+icol, strait)
	
	if 	locdir[icol]==-1:
		sheet1.write(row_flux_neg, len(areas)+1+icol, -flux1['vol']['pos'][iuse,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_pos, len(areas)+1+icol, -flux1['vol']['neg'][iuse,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_tot, len(areas)+1+icol, -flux1['vol']['tot'][iuse,icol].mean()*scale2year*m3tokm3)	
	else:
		sheet1.write(row_flux_pos, len(areas)+1+icol, flux1['vol']['pos'][iuse,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_neg, len(areas)+1+icol, flux1['vol']['neg'][iuse,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_tot, len(areas)+1+icol, flux1['vol']['tot'][iuse,icol].mean()*scale2year*m3tokm3)	

# add straits for basins
irow=labell.index('Strait Flux tot')	
for icol,basin in enumerate(areas):
	posflux=0
	negflux=0
	totflux=0
	if not BasinStraits[basin]==None:
		for item in BasinStraits[basin]:
			coli=locnames.index(item[0])
			if item[1]==-1:
				posflux+=-flux1['vol']['neg'][iuse,coli].mean()
				negflux+=-flux1['vol']['pos'][iuse,coli].mean()
				totflux+=-flux1['vol']['tot'][iuse,coli].mean()
			else:
				posflux+=flux1['vol']['pos'][iuse,coli].mean()
				negflux+=flux1['vol']['neg'][iuse,coli].mean()
				totflux+=flux1['vol']['tot'][iuse,coli].mean()
	sheet1.write(row_flux_pos, 1+icol, posflux*scale2year*m3tokm3)	
	sheet1.write(row_flux_neg, 1+icol, negflux*scale2year*m3tokm3)	
	sheet1.write(row_flux_tot, 1+icol, totflux*scale2year*m3tokm3)

	

	

Z0=np.asarray(ds['schism']['zcor'][0:4,:,-1]).mean(axis=0)
Z1=np.asarray(ds['schism']['zcor'][-4:,:,-1]).mean(axis=0)
V0={}
V1={}
for area in areas:
	basinareas[area]=A[eleminds[area]].sum() # tages mittel
	V0[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z0[faces[eleminds[area]]].mean(axis=1) ).sum()
	V1[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z1[faces[eleminds[area]]].mean(axis=1) ).sum()

for icol,basin in enumerate(areas):
	sheet1.write(len(labell)-3, col0+icol,V0[basin]*m3tokm3)
	sheet1.write(len(labell)-2, col0+icol,V1[basin]*m3tokm3)
	sheet1.write(len(labell)-1, col0+icol,(V1[basin]-V0[basin])*m3tokm3)
		
irow=labell.index('Netsum')
Net=dict.fromkeys(areas)	
for icol,basin in enumerate(areas):
	Net[basin]=0
	Net[basin]+=np.mean(np.asarray(water_from_atmo[basin])[useyears])+np.mean(np.asarray(water_to_atmo[basin])[useyears])#*scale2year*m3tokm3)
	Net[basin]+=np.mean(river_inflow[basin][useyears,:],axis=0).sum()#*scale2year*m3tokm3
	totflux=0
	if not BasinStraits[basin]==None:
		for item in BasinStraits[basin]:
			coli=locnames.index(item[0])
			if item[1]==-1:
				totflux+=-flux1['vol']['tot'][iuse,coli].mean()
			else:
				totflux+=flux1['vol']['tot'][iuse,coli].mean()
	Net[basin]+=totflux
	sheet1.write(irow, 1+icol,Net[basin]*scale2year*m3tokm3)		
		
wb.save(outfolder+'OceanBudgets'+run_name+'.xls')	

print('done computing balance')

#for iyear,year in enumerate(years):
#	sheet1.write(iyear, 0, str(year))
	

