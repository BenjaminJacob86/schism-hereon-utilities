# script to calculate basin balances

import os
import netCDF4
import sys
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/')
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/Lib/')
import matplotlib



# bctides in need to be modified with comments to work:
# ocean  and river: comments ae need for the function calculating averAGE RIVER INFLOWS
#09/01/2008 00:00:00 GMT
#0 40. ntip
#0  nbfr
#9 nope
#204 0 2 0 0 !ocean Aegen Sea
#1.0528e4 !mean from sum of all rvrs
#2 0 1 1 2 !river: Kizil_Irmak,
#1. !relax T
#0.
#1. relax S
#2 0 1 1 2 !river: Coruh,
#1. !relax T
#0.
#1. relax S
#2 0 1 1 2 !river: Rioni
#1. !relax T
#0.
#1. relax S
#2 0 2 1 2 !river: Kubon
#-425.
#
#
bgmode=False
if bgmode:
	matplotlib.use('AGG')

from matplotlib import path
from matplotlib import pyplot as plt
if not bgmode:
	plt.ion()
	
from netCDF4 import Dataset
import numpy as np
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
#from cftime import utime
########## settings #################################
# export OMP_NUM_THREADS=1 #


##### SETTINGS ########################################
run_name='LockExchange'            # Name of Excelfile Carrying bilances

# directories foor RUN and atmospheric forcing
#setupdir='/work/gg0028/g260114/RUNS/BLACKSEA/RUN24f/' #'/work/gg0028/g260114/RUNS/BLACKSEA/RUN24d/RUN24d/new_code/constant_bd/' 
#setupdir='/work/gg0028/g260114/RUNS/BLACKSEA/RUN24f/' 
setupdir='/work/gg0028/g260114/RUNS/BLACKSEA/RUN24h/' 
setupdir='/work/gg0028/g260114/RUNS/BLACKSEA/LockExchange/'
sfluxdir=setupdir+'/sflux/'
sflux_in=setupdir+'/sflux/'

									    # https://www.eea.europa.eu/data-and-maps/data/europe-seas
cwd=os.getcwd()
#ncdir=setupdir+'outputs_all2/'
ncdir=setupdir+'outputs_all/'

#fluxfile='/sciclone/home10/yinglong/DISKS/vims20/BlackSea/RUN25b/outputs/flux.out' # new schism
fluxfile=ncdir+'flux.out'

outfolder=setupdir+'/analysis'+run_name+'/'
kmlfile='/work/gg0028/SCHISM/schism-hzg-utilities/postprocesing/balances/OceanDB.kml'  #ocean regions kml file

# Areas to extract from boundary database
#areas=['Baltic Sea','Kattegat','Skagerrak','North Sea','English Channel',
#'Mediterranean Sea','Sea of Marmara',
#'Black Sea', 'Sea of Azov',] 
areas=['Black Sea', 'Sea of Azov','Sea of Marmara','Aegean Sea'] 
########
years=range(2008,2008+1)     # years considered
reftime=dt.datetime(2008,9,1,0,0,0)#+dt.timedelta(days=214) # fake shift to end to have a fullyear
prec_years=np.arange(2008,2008+1)
evap_years=np.arange(2008,2008+1)

startdate=reftime
endtime=dt.datetime(2009,9,1) # None   if = None use last model timestep as end time

# analysis to be done for full years
# if not to full years start in summer build artifical year


use_atmo_from_schism=True   #True: get prec / evaporationRate from schism output

							#False: get from atmo
do_plot=True
# name of volume transport calculation sections in order of fluxflag.prop

# locations of defined throughflux regions
locnames=['Asov','Bosporus N','Bosporus C','Bosporus S','Daradanelle N','Daradanelle C','Daradanelle S']
#locnames=['Asov','Bosporus N','Bosporus C','Bosporus S','Daradanelle N','Daradanelle C','Daradanelle S','Bosporus','Gibraltar','little Belt','Great Belt','Sound']

locdir=np.zeros(len(locnames))
locdir[:]=1
BasinStraits=dict.fromkeys(areas)
#BasinStraits['Mediterranean Sea']=[('Gibraltar',1)]   #+1 :flows ind, -1: flows pout (higher to lower number in fluxflag.prop)
#BasinStraits['Black Sea']=[('Bosporus C',-1)]
BasinStraits['Black Sea']=[('Bosporus C',-1),('Asov',+1)]
BasinStraits['Sea of Azov']=[('Asov',-1)]
BasinStraits['Sea of Marmara']=[('Bosporus C',+1),('Daradanelle C',-1)]
BasinStraits['Aegean Sea']=[('Daradanelle C',+1)]

#BasinStraits['Baltic Sea']=[('little Belt',1),('Great Belt',1),('Sound',1)]
#BasinStraits['Kattegat']=[('little Belt',-1),('Great Belt',-1),('Sound',-1)]
#BasinStraits['English CHannel']=[('NorthSea',-1)]
BasinStraits['ModelDomain']=None

# conversion factors
m3tokm3=1/1e9
rho=1000.0 # kg/m^3/
in_m_s=1.0/(rho)
in_mm_day=in_m_s*1000*86400
scale2year=86400*365  #m3/a   # per second to per year
###################################### END Settints   ################################################

################################## FUnctions ############################################
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
def get_annual_mean_river_runoff(s,areas,years,startdate=dt.datetime(2007,6,1),enddate=dt.datetime(2007,6,1)):
	""" 	Extract river Data from flux.th and source files 	using schism class and domain polygons 	returns combined river inflow per basin averaged by year in m3/s
	uses source class.
	this also accounst constant ocean bd outflow as negative river
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
		
		
		tmax=(enddate-startdate).total_seconds()  # upper limiot for period mean - stats in addition to individual years, epsecially sim doe not start alwas at 01.01.yyyy

		obdtype=[]
		bd_rnames=[]
		line='dummy'
		qconst=[]
		with open('bctides.in','r') as f:
			while len(line)>0:
				line=f.readline()
				if '4 4 4 4' in line:
					continue
				if ('river' in line) or ('ocean' in line):  # work in progress
					if line.split()[2] == '2': #column vor vel forcing is constant volume 
						print(line)
						obdtype.append('const_river')		
						line2=line.replace('\n','')
						try:
							bd_rnames.append(line2[line2.rindex(':')+1:])
						except:
							bd_rnames.append(line2[line2.index('!')+1:])
						line=f.readline()
						if '!' in line:
							line=line.split('!')[0]
						qconst.append(np.float(line)) #read constant runuff value
					elif line.split()[2] == '1': #column vor vel forcing is flux.ts
						obdtype.append('river')
						line2=line.replace(' ','').replace('\n','')
						bd_rnames.append(line2[line2.rindex(':')+1:])
				elif 'ocean' in line:
					obdtype.append('ocean')
					bd_rnames.append('ocean')
				elif ('!' in line) and ('relax' not in line) :	
					line2=line.replace('\n','')
					bd_rnames.append(line2[line2.index('!')+1:])
		
		
		#to follow the logic of the code q read in from flux.th is modified, inserting the constant values as timeseries, in the order of where they would occur counting the order of boundaries
		q2=np.zeros((q.shape[0],q.shape[1]+obdtype.count('const_river')))
		q2[:,0]=-q[:,0] # positiver irver inflow
		
		iconst=0 # index within qconst
		ifluxth=1 # start 1 one becuase 0 is time
		for i,obd in enumerate(obdtype):
		
			if obd=='const_river':
				q2[:,i+1]=qconst[iconst]
				iconst+=1
			else:
				q2[:,i+1]=q[:,ifluxth]
				ifluxth+=1
		q=q2
		
		
		# in which basin are openbds ?
		basinbdrivers={} # column index for flux.th  indices of of openboundaries belonging to this basin
		flux_indices={}  # column index for flux.th
		for ifile,basin in enumerate(areas):
				areaPoly=Path(list(zip(seas[basin][:,0],seas[basin][:,1])))
				basinbdrivers[basin]=np.asarray(np.where(areaPoly.contains_points(bdcenters)))
				# bnd
				if len(basinbdrivers[basin][0]) > 0:
					#flux_indices[basin]=basinbdrivers[basin]-obdtype[:basinbdrivers[basin][0][0]].count('ocean')+1 # one for time
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
	
	years=np.unique(Ryears)	
	
	years=years[years<=endtime.year]
	for basin in areas:
			if os.path.exists('flux.th'):
				#from_bdriver=-np.mean(q[Ryears==years[0],:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
				from_bdriver=-np.mean(q[Ryears==years[0],:][:,basinbdrivers[basin]+1],axis=0)[0]
			else: 	
				from_bdriver=[]
			if len(src_indices[basin])> 0: # add from src
				from_srcriver=np.mean(src.vsource[src.years==years[0],:][:,src_indices[basin]] ,axis=0) #dunno bout bracket when not interactive
				river_inflow[basin]=np.hstack((from_bdriver, from_srcriver))
			else:
				river_inflow[basin]=from_bdriver
				
			for year in years[1:]:
			
				#if year > endtime.year:
				#	print( str(year) + 'contained in runn off is left out becuase it is past prescribed end date' + str(endtime.year))
				#	continue
					
			
				if os.path.exists('flux.th'):
					#from_bdriver=-np.mean(q[Ryears==year,:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
					from_bdriver=-np.mean(q[Ryears==year,:][:,basinbdrivers[basin]+1],axis=0)[0]
				else:
					from_bdriver=[]				
				if len(src_indices[basin])> 0: # add from src
					from_srcriver=np.mean(src.vsource[src.years==year,:][:,src_indices[basin]] ,axis=0) #dunno bout bracket when not interactive				#from_srcriver=np.mean(src.vsource[src.years==year,:][src_indices[basin]] ,axis=1)[0]
					river_inflow_temp=np.hstack((from_bdriver, from_srcriver))
				else:
					river_inflow_temp=from_bdriver
				river_inflow[basin]=np.vstack((river_inflow[basin],river_inflow_temp))

				
			if os.path.exists('flux.th'):
				#from_bdriver=-np.mean(q[Ryears==years[0],:][:,np.asarray(basinbdrivers[basin])-4+1],axis=0)[0]
				from_bdriver=-np.mean(q[qt<=tmax,:][:,basinbdrivers[basin]+1],axis=0)[0]
			else: 	
				print('this part is not in place yet')
				from_bdriver=[]
			if len(src_indices[basin])> 0: # add from src
				from_srcriver=np.mean(src.vsource[src.years==years[0],:][:,src_indices[basin]] ,axis=0) #dunno bout bracket when not interactive
				river_inflow_tmp=np.hstack((from_bdriver, from_srcriver))
			else:
				river_inflow_tmp=from_bdriver
			river_inflow[basin]=np.vstack((river_inflow[basin],river_inflow_tmp))  #np.asarray(river_inflow[basin])
	# referenced periods =	
	period=[]
	for year in years:
		date0=qdates[Ryears==year][0]
		date1=qdates[Ryears==year][-1]
		period.append('-'.join((date0.strftime('%Y%m%d'),date1.strftime('%Y%m%d'))))
	period.append('-'.join((startdate.strftime('%Y%m%d'),enddate.strftime('%Y%m%d'))))
	river_inflow['periods']=period
	
	
	
	Ryears=period
	print('done reading annual mean river runoff by basin')
	return river_inflow,Ryears,bd_rnames,src,src_indices,basinbdrivers,q2	

############################################################


if not (os.path.exists(outfolder)):
        os.mkdir(outfolder)
		

#fluxflag
#ff=np.loadtxt('fluxflag.prop')
#tcv=np.loadtxt('tvd.prop')
#plt.figure()
#s.plotAtelems(ff[:,1])
#s.plotAtelems(ff[:,1][s.nvplt2nvp])
		
		
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
	
#quad #calculate aereas
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

plt.figure()
for i,key in enumerate(eleminds.keys()):
	plt.subplot(3,2,i+1)	
	s.plotAtelems(np.asarray(eleminds[key],np.float)*1)
	plt.title(key)
plt.savefig(outfolder+'domain_parts.png',dpi=300)	




################   net cdf access   ########################
## check for schout_nc files

s.ds=schism_outputs_by_variable(ncdir).ds  #load access to xarray for variable files
#ds={}
#ds['schism']=xr.open_mfdataset(schismfiles)
if not use_atmo_from_schism:
	prcfiles=list(np.sort(glob.glob(sfluxdir+'*prc*.nc')))
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
	dmean[area]=(A[eleminds[area]]/A[eleminds[area]].sum()*np.asarray(s.depths)[faces[eleminds[area]]].mean(axis=1)).sum()
	volumina[area]=dmean[area]*basinareas[area]


	
# check die konvertierung nicht workaround
#time not in new variables
p=param()
reftime=dt.datetime(int(p.get_parameter('start_year')),
int(p.get_parameter('start_month')),
int(p.get_parameter('start_day')),
int(p.get_parameter('start_hour')),0,0)	


schismdates=reftime+s.ds['out2d'].time.values*dt.timedelta(seconds=1)
#schismdates=date0+np.arange(len(tschism))*dt.timedelta(seconds= np.float((tschism[1]-tschism[0]))/1e9)
schismyears=np.asarray([date.year for date in schismdates])


# end date of simulation
# add 	 dates to dataset for easier sorting within xarray
dates = xr.DataArray(name='dates',data=schismdates,dims=["time"])
s.ds['out2d']['dates']=dates
datestr=str(s.ds['out2d'].dates[-1].values)
datestr=datestr[:datestr.rindex('.')]
if endtime==None:
	endtime=dt.datetime.strptime(datestr,'%Y-%m-%dT%H:%M:%S')

river_inflow,Ryears,bd_rnames,src,src_indices,basin_bdrivers,q2=get_annual_mean_river_runoff(s,areas,years,reftime,endtime)	


use_prec_evap=p.get_parameter('isconsv')>0

rho_water=1000 # kg/m3	
if not use_atmo_from_schism:
	# radiation precipitationRate
	year_means=ds['pure'].groupby('time.year').mean('time')
	means={}
	varname='prate'
	means[varname]=np.asarray(year_means[varname])/rho_water # / rhow_water when from air files#/timestep  #m/s
	#for varname in 'tp',:
	#	means[varname]=np.asarray(year_means[varname])/timestep  #m/s
	# evaporationRate
	#year_means=ds['evap'].groupby('time.year').mean('time')
	#for varname in 'e','es':
	#	varname
	#	means[varname]=np.asarray(year_means[varname]/timestep)  #m/s
	#schism_vars=['zcor']	
	#pickle.dump( means, open( run_name+"means_Pickle", "wb" ) )	
else:
	prec_years=np.unique(schismyears)
	evap_years=np.unique(schismyears)
	if use_prec_evap:
		schism_vars=['elevation','precipitationRate','evaporationRate']	
	else:	
		schism_vars=['elevation']	
		
		
schism_vars+=[ var for var in ['solarRadiation','totalHeat'] if var in  list(s.ds['out2d'].keys())]	 #,'upward_longwave','downward_longwave'
schism_means={}
#year_means=s.ds['out2d'].groupby('time.year').mean('time') # not woriking because no explicit time in schism
   
#schismyears=   
#time_inds=dict.fromkeys(prec_years)  #to reference years
#for year in prec_years:
#time_inds[year]=np.where(schismyears==year)


# add 	 dates to dataset for easier sorting within xarray
dates = xr.DataArray(name='dates',data=schismdates,dims=["time"])
s.ds['out2d']['dates']=dates
	
#year_means=s.ds['out2d'].groupby('time.year').mean('time') # not woriking because no explicit time in schism	
year_means=s.ds['out2d'].groupby('dates.year').mean('time')	


# cannot slice select using new introduced variable dates 
# convert to time indeices
tmax=(endtime-startdate).total_seconds()
#period_means=s.ds['out2d'].sel(time=slice(np.asarray(reftime,np.datetime64),np.asarray(endtime,np.datetime64))).mean('time')		
period_means=s.ds['out2d'].sel(time=slice(0,tmax)).mean('time')		
	
	
for varname in schism_vars[1:]:
	print(' doing ' + varname)
	if varname in  ['precipitationRate','evaporationRate']:
		factor=1.0/rho  # [kg/m/m/s]  /  (kg/m^3)/  -> m/s
	else:
		factor=1
	if len(year_means[varname].shape)==3:
		#schism_means[varname]=np.asarray(year_means[varname][:,:,-1]*factor)
		schism_means[varname]=np.asarray(year_means[varname][:,:,-1].values*factor)
		# append average for whole simulation period
		schism_means[varname]=np.vstack((schism_means[varname],period_means[varname][:,:,-1].values*factor))
	else:		
		#schism_means[varname]=np.asarray(year_means[varname][:]*factor)
		schism_means[varname]=np.asarray(year_means[varname][:].values*factor)
		schism_means[varname]=np.vstack((schism_means[varname],period_means[varname][:].values*factor))
#pickle.dump( schism_means, open( run_name+"_schismmeans_Pickle", "wb" ) )
print('done extracting schism annual means')	


if 'totalHeat' in  list(year_means.keys()):
	schism_means['total_flux']=year_means['solarRadiation']+year_means['totalHeat']  # including solar radiation
	schism_means['total_flux']=np.vstack((schism_means['total_flux'],period_means['solarRadiation']+period_means['totalHeat'])) #append perid mean

# plots
years=np.hstack((np.unique(schismyears),9999))
for iyear,year in enumerate(years):
	for varname in ['solarRadiation','upward_longwave','downward_longwave','total_flux']:
		if (varname in schism_vars) or (  varname in list(schism_means.keys())):
			plt.clf()
			ph,ch,ax=s.plotAtnodes(schism_means[varname][iyear])
			ch.set_label(varname + ' [W/m^2]')
			#plt.title('Annual Mean '+ str(year))
			plt.title('Annual Mean '+ Ryears[iyear])
			plt.tight_layout()
			plt.savefig(outfolder+'Annual_Mean_'+varname+'_'+str(year)+'.png',dpi=400)	

			plt.savefig(outfolder+'Annual_Mean_'+varname+'_'+str(year)+'.png',dpi=400)				
			
			
# load river data
years=np.unique(schismyears)


if 'total_flux' in schism_means.keys():
	heat_by_basin=dict.fromkeys(['solarRadiation','total_flux'])
else:	
	heat_by_basin=dict.fromkeys(['solarRadiation'])
for key in heat_by_basin.keys():
	heat_by_basin[key]=dict.fromkeys(areas)

# spatially mean weighted precipitationRate
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

elif use_prec_evap:
	
	for basin in areas:
		trimeansP=schism_means['precipitationRate'][:,faces[eleminds[basin],:]].mean(axis=2) # jahresmittel
		trimeansE=schism_means['evaporationRate'][:,faces[eleminds[basin],:]].mean(axis=2)
		water_from_atmo[basin]=[(trimeansP[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] ## m/s *m^2  -> m^3/s    m^3/a
		water_to_atmo[basin]=[(trimeansE[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] # m^3/a

		# heatflux by basin		
		for varname in heat_by_basin.keys():
			if varname=='total_flux':
				#trimeans=schism_means[varname].values[:,faces[eleminds[basin],:]].mean(axis=2) # value call was done earlier
				trimeans=schism_means[varname][:,faces[eleminds[basin],:]].mean(axis=2)
			else:
				trimeans=schism_means[varname][:,faces[eleminds[basin],:]].mean(axis=2)
			heat_by_basin[varname][basin]=[(trimeans[i,:]*A[eleminds[basin]]/A[eleminds[basin]].sum()).sum() for i in range(trimeansP.shape[0]) ] # W/mm/a
else:		# just filly Prec Evap wioth dummy value 0 as the run was without	

	dummy=np.zeros((len(years)+1,s.nnodes))
	for basin in areas:
		trimeansP=dummy[:,faces[eleminds[basin],:]].mean(axis=2) # jahresmittel
		trimeansE=dummy[:,faces[eleminds[basin],:]].mean(axis=2)
		water_from_atmo[basin]=[(trimeansP[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] ## m/s *m^2  -> m^3/s    m^3/a
		water_to_atmo[basin]=[(trimeansE[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] # m^3/a

		# heatflux by basin		
		for varname in heat_by_basin.keys():
			if varname=='total_flux':
				#trimeans=schism_means[varname].values[:,faces[eleminds[basin],:]].mean(axis=2) # value call was done earlier
				trimeans=schism_means[varname][:,faces[eleminds[basin],:]].mean(axis=2)
			else:
				trimeans=schism_means[varname][:,faces[eleminds[basin],:]].mean(axis=2)
			heat_by_basin[varname][basin]=[(trimeans[i,:]*A[eleminds[basin]]/A[eleminds[basin]].sum()).sum() for i in range(trimeansP.shape[0]) ] # W/mm/a

			
#for iyear,year in enumerate(years):
#	for varname in ['precipitationRate',]:
#		if (varname in schism_vars) or (  varname in list(schism_means.keys())):
#			plt.clf()
#			ph,ch,ax=s.plotAtnodes(schism_means[varname][iyear])
#			ch.set_label(varname + ' [W/m^2]')
#			#plt.title('Annual Mean '+ str(year))
#			plt.title('Annual Mean '+ Ryears[iyear])
#			plt.tight_layout()
#			
#			
			
				
#######   strait fluxes  ###############################
#
# if compiled without analysis module only total outflow
#
m=np.loadtxt(fluxfile)

if m.shape==(0,):
	print('there are no outputs for crossection volume fluxes setting all to zero')
	
	nregs=len(locnames)-1
	ntypes=5 #  np.where(np.diff(m[:100,0])>0)[0][0]+1
	#t1=m[::ntypes,0]
	t1=np.loadtxt('flux.th')[:,0]
	dummy=np.ones((len(t1),nregs+1))*np.nan
	if ntypes==1:
		print('flux out without analysis module contains only total flux')
		flux1={'vol':{'tot':dummy,'pos':dummy,'neg':dummy},'salt':{'pos':dummy,'neg':dummy}}
	else:	
		flux1={'vol':{'tot':dummy,'pos':dummy,'neg':dummy},'salt':{'pos':dummy,'neg':dummy}}
	fluxdates=reftime+np.arange(len(t1))*dt.timedelta(seconds=np.diff(t1[:2])[0])
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
else:	
	nregs=len(m[0,:])-1
	ntypes=np.where(np.diff(m[:100,0])>0)[0][0]+1
	t1=m[::ntypes,0]
	if ntypes==1:
		print('flux out without analysis module contains only total flux')
		flux1={'vol':{'tot':m[::ntypes,1:],'pos':m[1::ntypes,1:]*np.nan,'neg':m[2::ntypes,1:]*np.nan},'salt':{'pos':m[3::ntypes,1:]*np.nan,'neg':m[4::ntypes,1:]*np.nan}}
	else:	
		flux1={'vol':{'tot':m[::ntypes,1:],'pos':m[1::ntypes,1:],'neg':m[2::ntypes,1:]},'salt':{'pos':m[3::ntypes,1:],'neg':m[4::ntypes,1:]}}
	fluxdates=reftime+np.arange(len(t1))*dt.timedelta(days=np.diff(t1[:2])[0])
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


# total domain total period balance

water_from_atmo[basin]=[(trimeansP[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ]



#water_from_atmo[basin]=[(trimeansP[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] ## m/s *m^2  -> m^3/s    m^3/a
#water_to_atmo[basin]=[(trimeansE[i,:]*A[eleminds[basin]]).sum() for i in range(trimeansP.shape[0]) ] 

#-1 total period
basin='ModelDomain'
years=np.hstack((np.unique(schismyears),9999))	
for iyear,year in enumerate(years):
	Pbasin=trimeansP[iyear,:]*A[eleminds[basin]]*scale2year*m3tokm3 ## m/s *m^2  -> m^3/s    m^3/a
	Ebasin=trimeansE[iyear,:]*A[eleminds[basin]]*scale2year*m3tokm3  # km3/a
	for varname,var in zip(['precipitationRate','EvaporationRate'],[Pbasin,Ebasin]):
			plt.clf()
			ph,ch=s.plotAtelems(var)
			ch.set_label(varname + ' [km3/a] at triangle')
			#plt.title('Annual Mean '+ str(year))
			plt.title('Annual Mean '+ Ryears[iyear])
			plt.tight_layout()
			plt.savefig(outfolder+varname+'_'+str(years[-1]),dpi=300)	

#-1 total period
m2mm=1000
basin='ModelDomain'
Pbasin=trimeansP[-1,:]*scale2year*m2mm   #m3tokm3 ## m/s *m^2  -> m^3/s    m^3/a
Ebasin=trimeansE[-1,:]*scale2year*m2mm   #m3tokm3  # km3/a
years=np.hstack((np.unique(schismyears),9999))	
for iyear,year in enumerate(years):
	Pbasin=trimeansP[iyear,:]*scale2year*m2mm   #m3tokm3 ## m/s *m^2  -> m^3/s    m^3/a
	Ebasin=trimeansE[iyear,:]*scale2year*m2mm   #m3tokm3  # km3/a
	for varname,var in zip(['precipitationRate','EvaporationRate','Prec-Evap  Rate'],[Pbasin,Ebasin,Pbasin-Ebasin]):
			plt.clf()
			ph,ch,ax=s.plotAtnodes(var)
			ch.set_label(varname + ' [mm/a]')
			#plt.title('Annual Mean '+ str(year))
			plt.title('Annual Mean '+ Ryears[iyear])
			plt.tight_layout()
			plt.savefig(outfolder+varname+'mm_'+str(years[-1]),dpi=300)	
			

			
			

basin='ModelDomain'
P=water_from_atmo[basin][-1]*scale2year*m3tokm3
E=water_to_atmo[basin][-1]*scale2year*m3tokm3
R=river_inflow[basin][-1]*scale2year*m3tokm3
Rin=R[R>0].sum()
Bdout=R[R<0].sum()			
			
plt.figure()
plt.clf()
ph1=plt.bar([0,2],[P,Rin])
ph2=plt.bar([1,3],[E,-Bdout],color='r')
plt.xticks([0,1,2,3,4])
plt.gca().set_xticklabels(['Prec','Evap','River','bdout','sum'],rotation=45)
S=P-E+Rin+Bdout
fac=(-1)**(1-S>0)
plt.bar([4],[S*fac],color=[['b','r'][0+(fac<0)]])
plt.legend((ph1,ph2),['positive','negative'])
plt.ylabel('flux into model domain [km^3/a]')
plt.grid()
plt.title('averaging period ' + reftime.strftime('%Y%m%d')+ '-' + endtime.strftime('%Y%m%d'))
plt.tight_layout()
for i,val in enumerate([P,E,Rin,Bdout,S]):
	plt.text(i-0.1,15,str(round(val,2)),rotation=90,fontsize=14)
plt.savefig(outfolder+varname+'mm_'+str(years[-1]),dpi=300)		



# output excel	
row=0
col0=1
#labell=['Quantity [km3/a] \ Basin:','area [km^2]','precipitationRate','evaporationRate','River','Strait inflow','Strait outflow','Strait Flux tot','Netsum','volume 0','volume 1','Dv ']

labell=['Quantity [km3/a] \ Basin:','area [km^2]','precipitationRate','evaporationRate','River','Strait inflow','Strait outflow','Strait Flux tot','Netsum','','annual mean heatfluxes [W/m^2]','solarRadiation', 'upward_longwave', 'downward_longwave', 'total_flux']


#years=np.arange(2007,2010+1)
years=np.unique(schismyears)
row_flux_pos=labell.index('Strait inflow')
row_flux_neg=labell.index('Strait outflow')
row_flux_tot=labell.index('Strait Flux tot')


#D=np.asarray(ds['schism']['depth'][0,:])
D=np.asarray(s.depths)

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
		
	irow=labell.index('precipitationRate')	
	for icol,basin in enumerate(areas):
		sheet1.write(irow, col0+icol,water_from_atmo[basin][np.where(prec_years==year)[0][0]]*scale2year*m3tokm3) #
		#mat[irow,icol]=water_from_atmo[basin][np.where(prec_years==year)[0][0]]*scale2year*m3tokm3
	irow=labell.index('evaporationRate')	
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



# for the entire period:

year=9999
print('doing year ' + str(year))
sheet1 = wb.add_sheet('Budget {:s}-{:s}'.format( reftime.strftime('%Y%m%d'),endtime.strftime('%Y%m%d')))  # add_sheet is used to create 

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
	sheet1.write(irow, col0+icol, river_inflow[basin][-1,:].sum()*scale2year*m3tokm3)
	#mat[irow,icol]=river_inflow[basin][np.where(years==year),:].sum()*scale2year*m3tokm3
	
irow=labell.index('precipitationRate')	
for icol,basin in enumerate(areas):
	sheet1.write(irow, col0+icol,water_from_atmo[basin][-1]*scale2year*m3tokm3) #
	#mat[irow,icol]=water_from_atmo[basin][np.where(prec_years==year)[0][0]]*scale2year*m3tokm3
irow=labell.index('evaporationRate')	
for icol,basin in enumerate(areas):
	sheet1.write(irow, col0+icol,-np.abs(water_to_atmo[basin][-1])*scale2year*m3tokm3) #*scale2year

# add strait labels # this does not accpint direction prescriped


iperiod=(fluxdates>=reftime) & (fluxdates<=endtime) # determined period
for icol,strait in enumerate(locnames):
	#sheet1.write(0, len(areas)+1+icol, strait)
	if 	locdir[icol]==-1:
		sheet1.write(row_flux_neg, len(areas)+1+icol, -flux1['vol']['pos'][iperiod,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_pos, len(areas)+1+icol, -flux1['vol']['neg'][iperiod,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_tot, len(areas)+1+icol, -flux1['vol']['tot'][iperiod,icol].mean()*scale2year*m3tokm3)	
	else:
		sheet1.write(row_flux_pos, len(areas)+1+icol, flux1['vol']['pos'][iperiod,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_neg, len(areas)+1+icol, flux1['vol']['neg'][iperiod,icol].mean()*scale2year*m3tokm3)	
		sheet1.write(row_flux_tot, len(areas)+1+icol, flux1['vol']['tot'][iperiod,icol].mean()*scale2year*m3tokm3)	

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
				posflux+=-flux1['vol']['neg'][iperiod,coli].mean()
				negflux+=-flux1['vol']['pos'][iperiod,coli].mean()
				totflux+=-flux1['vol']['tot'][iperiod,coli].mean()
			else:
				posflux+=flux1['vol']['pos'][iperiod,coli].mean()
				negflux+=flux1['vol']['neg'][iperiod,coli].mean()
				totflux+=flux1['vol']['tot'][iperiod,coli].mean()
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
	Net[basin]+=(water_from_atmo[basin][-1]-np.abs(water_to_atmo[basin][-1]))#*scale2year*m3tokm3)
	Net[basin]+=river_inflow[basin][-1,:].sum()#*scale2year*m3tokm3
	totflux=0
	if not BasinStraits[basin]==None:
		for item in BasinStraits[basin]:
			coli=locnames.index(item[0])
			if item[1]==-1:
				totflux+=-flux1['vol']['tot'][iperiod,coli].mean()
			else:
				totflux+=flux1['vol']['tot'][iperiod,coli].mean()
	Net[basin]+=totflux
	sheet1.write(irow, 1+icol,Net[basin]*scale2year*m3tokm3)

# heatflux vars
for varname in heat_by_basin.keys():
	irow=labell.index(varname)		
	for icol,basin in enumerate(areas):
		sheet1.write(irow, col0+icol,heat_by_basin[varname][basin][-1])
			

Z0=s.ds['out2d']['elevation'][0:4,:].mean(axis=0).values   # 6 hourly outputs, daily mean for volume at start and end
Z1=s.ds['out2d']['elevation'][-4:,:].mean(axis=0).values  

plt.subplot(3,1,1)
s.plotAtnodes(Z0)
plt.subplot(3,1,2)
s.plotAtnodes(Z1)
plt.subplot(3,1,3)
s.plotAtnodes(Z1-Z0)

#daily mean elv

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
		
	sheet1.write(irow, 1+icol,Net[basin]*scale2year*m3tokm3)					
			
wb.save(outfolder+'OceanBudgets'+run_name+'2xls')	
print('done computing balance')
				




				
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
		
		

# ignore belod for now


		
		
unqyears=np.hstack((np.unique(schismyears)))

#useyears=np.asarray([(schismyears==year).sum()/4 for year in unqyears])>360  # only accound years having more than 360 days
useyears=np.ones((len(unqyears)),np.bool)  #np.hstack((unqyears>0,True)) # last true tells including the total period balance



#useyear_atmo
years_used=unqyears[useyears] #+[True,] # last true tells including the total period balance

sheet1 = wb.add_sheet('Mean Budget {:d} to {:d}'.format(unqyears[useyears][0],unqyears[useyears][-1]))  # add_sheet is used to create sheet. 

year=years[-1]
i0,i1=np.where(schismyears==year)[0][[0,-1]]
labell[-3]='volume '+str(schismdates[0])	
labell[-2]='volume '+str(schismdates[-1])	

for irow,label in enumerate(labell):
	sheet1.write(irow, 0, label)

for icol,basin in enumerate(areas):
	sheet1.write(row, col0+icol, basin)

irow=labell.index('River')	
for icol,basin in enumerate(areas):
	#sheet1.write(irow, col0+icol, river_inflow[basin][useyears,:].mean(axis=0).sum()*scale2year*m3tokm3)
	sheet1.write(irow, col0+icol, river_inflow[basin][useyears,:].mean(axis=0).sum()*scale2year*m3tokm3)

irow=labell.index('precipitationRate')	
for icol,basin in enumerate(areas):
	sheet1.write(irow, col0+icol,np.mean(np.asarray(water_from_atmo[basin])[useyears],axis=0)*scale2year*m3tokm3) #
	
irow=labell.index('evaporationRate')	
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

	






#Z0=np.asarray(ds['schism']['zcor'][0:4,:,-1]).mean(axis=0)   # group by year Iguess
#Z1=np.asarray(ds['schism']['zcor'][-4:,:,-1]).mean(axis=0)

#Z0=np.asarray(ds['schism']['zcor'][0:4,:,-1]).mean(axis=0)   # group by year Iguess
#Z1=np.asarray(ds['schism']['zcor'][-4:,:,-1]).mean(axis=0)

Z0=s.ds['out2d']['elevation'][0:4,:].mean(axis=0).values   # 6 hourly outputs, daily mean for volume at start and end
Z1=s.ds['out2d']['elevation'][-4:,:].mean(axis=0).values  

plt.subplot(3,1,1)
s.plotAtnodes(Z0)
plt.subplot(3,1,2)
s.plotAtnodes(Z1)
plt.subplot(3,1,3)
s.plotAtnodes(Z1-Z0)

#daily mean elv

V0={}
V1={}
for area in areas:
	basinareas[area]=A[eleminds[area]].sum() # tages mittel
	V0[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z0[faces[eleminds[area]]].mean(axis=1) ).sum()
	V1[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z1[faces[eleminds[area]]].mean(axis=1) ).sum()

# for plooint conv km3
#	
#V0={}
#V1={}
#for area in areas:
#	basinareas[area]=A[eleminds[area]].sum() # tages mittel
#	V0[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z0[faces[eleminds[area]]].mean(axis=1) ).sum()*m3tokm3
#	V1[area]=(A[eleminds[area]]*D[faces[eleminds[area]]].mean(axis=1)+Z1[faces[eleminds[area]]].mean(axis=1) ).sum()*m3tokm3
#km3	

#dV=dict.fromkeys(V0)
#for key in dV.keys():
#		dV[key]=(V1[key]-V0[key])
	
#dV=dict.fromkeys(V0)
#for key in dV.keys():
#		dV[key]=(V1[key]-V0[key])/m3tokm3


	
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
	

	
#volume loss over time
# volum curves for basin and sea level


fluxdt=90
for basin in areas:
	plt.figure(1)
	plt.clf()
	print(basin)
	# surface fluxes and elevation
	Bs=np.zeros(s.ds['out2d']['elevation'].shape[0])
	BsPE=np.zeros(s.ds['out2d']['elevation'].shape[0])
	#BP=np.zeros(s.ds['out2d']['elevation'].shape[0])
	#BE=np.zeros(s.ds['out2d']['elevation'].shape[0])
	weights=A[eleminds[basin]]/A[eleminds[basin]].sum()

	# Net evapotranspiration
	if use_prec_evap:
		PmE=s.ds['out2d']['precipitationRate']-s.ds['out2d']['evaporationRate']
	else:	
		PmE=s.ds['out2d']['elevation']*0			
	for ti in range(s.ds['out2d']['elevation'].shape[0]):
		PmEi=PmE[ti,:].values	
		BsPE[ti]=(PmEi[faces[eleminds[basin]]].mean(axis=1)*A[eleminds[basin]]).sum()
		#BP[ti]=(s.ds['out2d']['precipitationRate'][ti,:].values[faces[eleminds[basin]]].mean(axis=1)*A[eleminds[basin]]).sum()
		#BE[ti]=(s.ds['out2d']['evaporationRate'][ti,:].values[faces[eleminds[basin]]].mean(axis=1)*A[eleminds[basin]]).sum()
		ci=s.ds['out2d']['elevation'][ti,:].values	
		Bs[ti]=(ci[faces[eleminds[basin]]].mean(axis=1)*weights ).sum()
	BsPE*=in_m_s  ##    m^3/s  # convert mass in volume
	#BP*=in_m_s 
	#BE*=in_m_s 
	Vbs=basinareas[basin]*Bs   # basin from sea level
	
	# river inflow
	Dt=21600
	Qr=-np.squeeze(q2[:,1+basin_bdrivers[basin]])
	
	m=np.loadtxt('flux.th')
	riverdates=reftime+m[:,0]*dt.timedelta(seconds=1)
	Dt_river=(np.diff(riverdates[:2])/dt.timedelta(seconds=1))[0]
	
	
	# interpolate river time to time
	ToDates=s.ds['out2d']['time'].values
	indates=m[:,0] 
	
	if  len(Qr.shape) > 1:
		RiverCumsum=np.cumsum(Qr.sum(axis=1))*Dt_river*m3tokm3
	elif len(Qr.shape) == 1:	
		RiverCumsum=np.cumsum(Qr)*Dt_river*m3tokm3
	RiverCumsumIntp=np.interp(ToDates,indates,RiverCumsum)
	PEcumsum=np.cumsum(BsPE)*Dt*m3tokm3	
	
	# determine straitfluxes into river
	if BasinStraits[basin]==None:
		Total=0
	else:
		keys,facs=[],[]
		for item in BasinStraits[basin]:
			keyi,faci=item
			keys.append(keyi)
			facs.append(faci)	
		
		Rtemp=dict.fromkeys(keys)	
		Total=0
		for key,fac in zip(keys,facs):
			Rtemp[key]=np.cumsum(flux1['vol']['tot'][:,locnames.index(key)])*fluxdt*fac*m3tokm3
			Total+=Rtemp[key]

	plt.plot(dates,Vbs*m3tokm3,'k',linewidth=2)
	plt.plot(dates,PEcumsum)
	plt.plot(riverdates,RiverCumsum)
	plt.ylabel('cummulated dV [km^3]')
	leg=['VolCeta','P-E','River']
	# rivers	
	if BasinStraits[basin]!=None: # ad stratfluxes
		for key in Rtemp.keys():
			leg+=['Strait: '+ key,]
			plt.plot(fluxdates,Rtemp[key])
		indates=np.asarray((fluxdates-reftime)/dt.timedelta(seconds=1),np.float)
		TotalIntp=np.interp(ToDates,indates,Total)
		plt.plot(dates,RiverCumsumIntp+PEcumsum+TotalIntp)	
		leg+=['(P-E)+River+Straits',]
	else:	
		leg+=['(P-E)+River',]
		plt.plot(dates,RiverCumsumIntp+PEcumsum)	
		
	plt.legend(leg)
	plt.title(basin)
	plt.xlim(reftime,endtime)
	plt.tight_layout()
	plt.savefig(outfolder+'CumSumCurve'+basin.replace(' ','_')+'.png',dpi=300)
	
	
	
	plt.figure(2)
	plt.clf()
	P=water_from_atmo[basin][-1]*scale2year*m3tokm3
	E=water_to_atmo[basin][-1]*scale2year*m3tokm3
	Str=TotalIntp[-1]/ToDates[-1]*scale2year
	if basin=='ModelDomain':
		Str=Bdout.copy()
	Rin=RiverCumsumIntp[-1]/ToDates[-1]*scale2year
	ph1=plt.bar([0,2],[P,Rin])
	ph2=plt.bar([1],[E],color='r')
	fac=(-1)**(1-Str>0)
	plt.bar([3],[Str*fac],color=[['b','r'][0+(fac<0)]])
	plt.xticks([0,1,2,3,4])
	plt.gca().set_xticklabels(['Prec','Evap','River','StraitNET','sum'],rotation=45)
	S=P-E+Rin+Str
	fac=(-1)**(1-S>0)
	plt.bar([4],[S*fac],color=[['b','r'][0+(fac<0)]])
	plt.legend((ph1,ph2),['positive','negative'])
	plt.ylabel('flux into model domain [km^3/a]')
	plt.grid()
	plt.title(['averaging period ' + reftime.strftime('%Y%m%d')+ '-' + endtime.strftime('%Y%m%d'), basin])
	plt.tight_layout()
	for i,val in enumerate([P,E,Rin,Bdout,S]):
		plt.text(i-0.1,15,str(round(val,2)),rotation=90,fontsize=14)
	plt.savefig(outfolder+'Bars'+basin.replace(' ','_')+'.png',dpi=300)
	
	
###### flux plot###############
add_fluxplots=True
if add_fluxplots:
	MS=3 # markersize

	# Black Sea 5 x 10^4
	# Baltic Sea 2.5 10^5
	#endtime=endtime


	t=t1
	nrows=int(np.ceil(np.sqrt(nregs)))
	ncols=int(np.round(np.sqrt(nregs)))

	n=len(t)
	plt.close('all')
	fig=plt.figure()
	for i in range(nregs):
		plt.subplot(nrows,ncols,i+1)
		tot=flux1['salt']['pos'][:,i]+flux1['salt']['neg'][:,i]
		plt.plot(t, np.cumsum(tot[:n]),'k')              
		plt.plot(t, np.cumsum(flux1['salt']['pos'][:,i][:n]),'--')              
		plt.plot(t, np.cumsum(flux1['salt']['neg'][:,i][:n]),'--')              
		plt.ylabel('vol salt [kg m^3/s]')
		plt.title(locnames[i])
		plt.legend(['tot','pos','neg'],frameon=False)
	plt.tight_layout()
	plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
	fig.subplots_adjust(top=0.9)
	plt.savefig(outfolder+'flux1_salt_cumts.png',dpi=300)
	plt.close()


	fig.clf()
	for i in range(nregs):
		plt.subplot(nrows,ncols,i+1)
		plt.plot(t, np.cumsum(flux1['vol']['tot'][:,i][:n]),'k')              
		plt.plot(t, np.cumsum(flux1['vol']['pos'][:,i][:n]),'--')              
		plt.plot(t, np.cumsum(flux1['vol']['neg'][:,i][:n]),'--')              
		plt.ylabel('vol [m^3/s]')
		plt.title(locnames[i])
		plt.legend(['tot','pos','neg'],frameon=False)
	plt.tight_layout()
	plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
	fig.subplots_adjust(top=0.9)
	plt.savefig(outfolder+'flux1_cumts.png',dpi=300)
	plt.close()


	# plt.figure()
	# for i in range(nregs):
		# plt.subplot(nrows,ncols,i+1)
		# plt.plot(t, np.cumsum(flux1['salt']['pos'][:,i]),'--')              
		# plt.plot(t, np.cumsum(flux1['salt']['neg'][:,i]),'--')              
		# plt.plot(t, np.cumsum(flux1['salt']['pos'][:,i]-np.cumsum(flux1['vol']['neg'][:,i])),'k')              
		# plt.ylabel('vol salt [kg m^3/s]')
		# plt.title(locnames[i])
		# plt.legend(['pos','neg'],frameon=False)
	# plt.tight_layout()
	# plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
	# plt.savefig('flux1_saltts.png',dpi=300)
	# plt.close()

	fig.clf()
	for i in range(nregs):
		plt.subplot(nrows,ncols,i+1)
		plt.plot(flux1['vol']['tot'][:,i],flux1['vol']['pos'][:,i][:n],',',markersize=MS)              
		plt.plot(flux1['vol']['tot'][:,i],-flux1['vol']['neg'][:,i][:n],',',markersize=MS)              
		plt.xlabel('total flux1 [m^3/s]' )
		plt.ylabel('partial flux1 [m^3/s]')
		plt.legend(['pos','neg'],frameon=False)
		plt.title(locnames[i])
	plt.tight_layout()
	plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
	fig.subplots_adjust(top=0.9)
	plt.tight_layout()
	plt.savefig(outfolder+'flux1_diag.png',dpi=300)
	plt.close()

	fig.clf()
	for i in range(nregs):
		plt.subplot(nrows,ncols,i+1)
		tot=flux1['salt']['pos'][:,i]-flux1['salt']['neg'][:,i]
		plt.plot(tot,flux1['salt']['pos'][:,i][:n],',',markersize=MS)              
		plt.plot(tot,-flux1['salt']['neg'][:,i][:n],',',markersize=MS)              
		plt.xlabel('total salt \n flux1 [kg m^3/s]' )
		plt.ylabel('partial salt \n flux1 [kg m^3/s]')
		plt.legend(['pos','neg'],frameon=False)
		plt.title(locnames[i])
	plt.tight_layout()
	plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
	fig.subplots_adjust(top=0.9)
	plt.savefig(outfolder+'flux1_diag_salt.png',dpi=300)
	plt.close()	
		
		
	####


	i=locnames.index('Bosporus C')

	i0d=dt.datetime(2008,9,1)
	i1d=dt.datetime(2008,11,1)


	i0d=dt.datetime(2009,7,1)
	i1d=dt.datetime(2009,9,1)


	i0=np.where(fluxdates>=i0d)[0][0]
	i1=np.where(fluxdates>=i1d)[0][0]

	n=i1

	#plt.clf()
	#plt.plot(fluxdates[:n],flux1['vol']['pos'][:,i][:n]*scale2year*m3tokm3,',',markersize=MS)              
	#plt.plot(fluxdates[:n],-flux1['vol']['neg'][:,i][:n]*scale2year*m3tokm3,',',markersize=MS)             
	#plt.ylabel('<Q> [km³/a]')

	plt.clf()
	plt.plot(fluxdates[:n],flux1['vol']['pos'][:,i][:n]*scale2year*m3tokm3)              
	plt.plot(fluxdates[:n],-flux1['vol']['neg'][:,i][:n]*scale2year*m3tokm3)             
	plt.ylabel('Q km³/a]')
	plt.legend(['outflow','inflow'])



	fout=flux1['vol']['pos'][:,i][i0:i1].mean()*scale2year*m3tokm3
	fin=-flux1['vol']['neg'][:,i][i0:i1].mean()*scale2year*m3tokm3
	fnet=flux1['vol']['tot'][:,i][i0:i1].mean()*scale2year*m3tokm3
	plt.clf()
	plt.bar([0,1,2],[fin,fout,fnet])
	plt.xticks([0,1,2])
	plt.ylabel('<Q> [km³/a]')
	for i,val in enumerate([fin,fout,fnet]):
		plt.text(i-0.1,15,str(round(val,2)),rotation=90,fontsize=14)
	plt.gca().set_xticklabels(['inflow','outflow','Netout'],rotation=45)
	plt.title(locnames[i]+ ': averaging period '+'-'.join((fluxdates[i0].strftime('%Y%m%d'),fluxdates[i1].strftime('%Y%m%d'))))
	plt.savefig(outfolder+'inoutFlorBars.png',dpi=300)




	i0d=dt.datetime(2008,9,1)
	i1d=dt.datetime(2008,11,1)


	i0d=dt.datetime(2009,7,1)
	i1d=dt.datetime(2009,9,1)


	i0=np.where(fluxdates>=i0d)[0][0]
	i1=np.where(fluxdates>=i1d)[0][0]

	n=i1
	MS=8
	plt.clf()
	for i in range(nregs):
		plt.subplot(nrows,ncols,i+1)
		plt.plot(flux1['vol']['tot'][:,i][:n],flux1['vol']['pos'][:,i][:n],',',markersize=MS)              
		plt.plot(flux1['vol']['tot'][:,i][:n],-flux1['vol']['neg'][:,i][:n],',',markersize=MS)              
		plt.xlabel('total flux1 [m^3/s]' )
		plt.ylabel('partial flux1 [m^3/s]')
		plt.legend(['pos','neg'],frameon=False)
		plt.title(locnames[i])
	plt.tight_layout()
	#plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
	plt.suptitle('-'.join((fluxdates[i0].strftime('%Y%m%d'),fluxdates[i1].strftime('%Y%m%d'))))
	#fig.subplots_adjust(top=0.9)
	plt.tight_layout()
	plt.savefig(outfolder+'flux1_diag.png',dpi=300)
	plt.close()
		
		
		
		

	#straits	
		