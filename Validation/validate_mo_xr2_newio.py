#export OMP_NUM_THREADS=1 #
import os
import netCDF4
import sys
import matplotlib
#matplotlib.use('Agg') # backend
background=False
from matplotlib import pyplot as plt
if background:
	matplotlib.use('Agg') # backend
else:
	plt.ion()	
from matplotlib import path
import datetime as dt
import glob
from scipy.spatial import cKDTree

import dask
import distributed
dask.__version__, distributed.__version__

import xarray as xr
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
import pickle
import datetime as dt
from netCDF4 import Dataset,MFDataset
from matplotlib import gridspec  # unequal subplots
from scipy.interpolate import interp1d
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/Lib/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')

from matplotlib.dates import DateFormatter
myFmt = DateFormatter("%m-%d")

from techit import * # import latex script
#from techit2 import * # import latex script
from TaylorDiagram import * # import taylordiagram
from schism import * # import schism functions
#from data_and_model_classes import cmems

# validate agains marnet data as stored as netcdf by sebastian

# Marnet Data: netcd Variables in files
#<station_name>_year.nc
# WTxm Wasser temperature at depth x m
# WSxm Wasser Salnity  at depth x m    

# try parallel?

#### I have added 1 day to compensate error in falsely set param.in start date
## updated hard coded timesteps with automatically determined ones

########## settings #################################
# directories
mardir='/gpfs/work/ksddata/observation/insitu/Marnet_Sebastian/' 					  # myocean
mardir='/work/gg0028/ksddata/insitu/CMEMS/NorthWestShelf/MO/'

# list of model setups to validate first setup is reference for station selections
#setupdir=['/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/',]
setupdir=['/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/','/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/HeatTest/']

#ncdir=setupdir+'/outputs01/'
year=2017		    # currently data between 2011 and 2018				
dtol=0.05           # distance tolerance in degree lon/lat for station selection 

ncdir=[setupdir[0] + 'outputs_all/'] 		  #   directory of schism nc output or 
ncdir+=[setupdir[1] + 'outputs_all/']

setup_names=['Veg_CNTRL','Veg_CNTRL_Heatcheck'] #,'GNU','ParamReset']


#oceandir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/'# amm15 dir

outdir=setupdir[0]+'mo_valid4/'	   # output directory where images will be stored
if not os.path.exists(outdir): os.mkdir(outdir) 

Tmax=45
Smax=49
add_amm15=False  											# add plot comparing with Amm15 | Amm15 profiles 
use_station_in=[False,False]					  # True: use station output False: extract from netcdf (slow)
put_pics_to_texdoc=True    										# images will be put in tex document
latexname='HRBallje.tex' #'GB1_5Ems_marnet.tex'										# in alphanumerical order 
latextitle='HRBallje 2017 vs Marnet' #'GB 2011 + FW 1.5EMS  in EFWS vs Marnet'
latextext= str(year) + ' validation against Marnet'

# reduce color scales to qunatile rnage , bypass outliers affectiong scale
quantile_range=True
quantiles=[0.01, 0.99]

## image
def set_fontsize(fac=1): 
	SMALL_SIZE = 10 *fac
	MEDIUM_SIZE = 12 *fac 
	BIGGER_SIZE = 14 *fac 
	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure

######### load SCHISM setup   ##################################
cwd=os.getcwd()
setups={}
output={}
access={}

setups=[]
newio=[]
for i,folder in enumerate(setupdir):
	os.chdir(folder)
	#s=schism_setup()
	setups.append(schism_setup())
	setups[i].nntree = cKDTree(list(zip(setups[i].lon,setups[i].lat))) 

	#lon,lat,D=np.asarray(s.lon),np.asarray(s.lat),np.asarray(s.depths)
	#lon[D>dmin]=9999
	#s.nntree = cKDTree(list(zip(lon,lat))) 
	#setups[i]=s
	if use_station_in[i]:
		if os.path.exists('station.in'): #& os.path.exists(setupdir[i]+'staout_1'):
			staout=schism_station_output()
			output[i]=staout
		else:
			print('station output does not exist, exiting program')
			exit()
	else:
		schismfiles=[] 
		if len(glob.glob(ncdir[i]+'schout_*.nc'))==0: #
			#new io
			#s.nc=schism_outputs_by_variable()
			#s.nc=schism_outputs_by_variable(ncdir)
			access[i]=schism_outputs_by_variable(ncdir[i],varlist=['out2d','salinity','temperature','zCoordinates']) #
			newio.append(True)
		else:	#old io schout_
			#for iorder in range(6): # check for schout_nc files until 99999
			#	schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
			#nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
			#schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
			#nrs=list(np.asarray(nrs)[np.argsort(nrs)])
			setups[i].nc=schism_output2(ncdir[i]).nc
			#setups[i].nc=xr.open_mfdataset(schismfiles)
			access[i]=setups[i].nc
			newio.append(False)
	
# use cmems TG files  # pre select in domain
# pattern to check in amm15 files 		
# access[0].ds[access[0].vardict['temperature']]

if use_station_in[0]:
	T=staout.time
	T=np.asarray([np.datetime64(ti) for ti in T])
else:
	if newio[0]:
		os.chdir(setupdir[0]) #[0]
		p=param('')
		reftime=np.asarray(dt.datetime(np.int(p.get_parameter('start_year')),\
	np.int(p.get_parameter('start_month')),\
	np.int(p.get_parameter('start_day')),\
	np.int(p.get_parameter('start_hour')),0,0),np.datetime64)
		T=reftime+access[0].ds['out2d']['time'].values *np.timedelta64(1,'s')
	else:
		#T=s.nc['time'].values
		T=access[0]['time'].values


######### initialize marnet Data located in SCHISM Domain #######

if len(glob.glob((mardir+'20??*'))) > 0: #ear folders
	foldertype='time'
	folders=np.sort(glob.glob(mardir+'{:d}*'.format(year)))
	files=glob.glob(folders[0]+'/*.nc')
	nfolders=1
else:
	foldertype='station'
	folders=glob.glob(mardir+'/*/')
	nfolders=len(folders)
	files=glob.glob(folders[0]+'/*.nc')


# select stations withing domain (tolerance distance) # and constructe file acess and nearest neighbours
names=[]	
stations={'coord':[],'names':[],'MO':{}}	
sources=setup_names.copy()
for key in sources:	
	stations[key]={'nc':[],'time':[],'nn':[],'coord':[],'temp':[],'salt':[],'D':[],'names':[]}

for i,setup_name in enumerate(setup_names):
	stations[setup_name]['MO']=dict.fromkeys(names) # place holder for interpolated to stationd data
	s=setups[i]
	print('find neighbours for ' + setup_name )
	if use_station_in[i]:
		staout=output[i]
	for folder in folders[:nfolders]:
		nfiles=1 if foldertype == 'station' else len(files) # files are eother timer or station
		if foldertype=='station':  # read new stations otherwise stations are in date folder
			files=glob.glob(folder+'/*NO*.nc')
		for file in files[:nfiles]: 
			#if 'Helgo' in file:
			#	break
			a=xr.open_dataset(file,decode_times=False)
			coord=np.float(a.geospatial_lon_min),np.float(a.geospatial_lat_min)
			if use_station_in[0]:
				D=np.sqrt((staout.coords[:,0]-coord[0])**2+(staout.coords[:,1]-coord[1])**2)
				nn=np.argmin(D)
				lldists=D[nn]
			else:
				lldists,nn=s.nntree.query(coord)		
			# in domain
			if lldists < dtol:
				#break
				if foldertype=='station':
					name=folder[folder[:-1].rindex('/')+1:-1]
				else:
					name=file[file.rindex('/')+1:].split('_')[-2]
				
				if True: # (i==0):
					if nfiles > 1:
						ncfiles=[file.replace('{:d}01'.format(year),'{:d}{:02d}'.format(year,mon+1)) for mon in range(12) if os.path.exists(file.replace('{:d}01'.format(year),'{:d}{:02d}'.format(year,mon+1)))]
					else:	
						ncfiles=np.sort(files)
					try:
						ncs=MFDataset(ncfiles,aggdim='TIME')
						#if np.max([ 'SAL' in  key for key in  ncs.variables.keys() ]):
						#	print('has salt')
						#	from IPython import embed; embed()
					except:
						continue
					t=(ncs['TIME'][:])

					#platform_code
					try:
						if np.ma.isMaskedArray(ncs['TIME_QC'][:]):
							igood=(ncs['TIME_QC'][:] < 3)[:] & ncs['TIME_QC'][:].mask==False
						else:
							igood=(ncs['TIME_QC'][:] < 3)[:] # select good data (ncs['TIME_QC'][:]<3) and
					except:
						try:
							if np.ma.isMaskedArray(ncs['SLEV'][:,0]):
								igood=ncs['SLEV'][:,0].mask==False # np.ones(len(t),bool)
							else:	
								igood=np.ones(len(t),bool)
						except:
							#pass
							continue
					
					t0=dt.datetime.strptime(ncs['TIME'].units[11:30],'%Y-%m-%dT%H:%M:%S')
					timeunit=ncs['TIME'].units[:ncs['TIME'].units.index(' ')]
					t=t[igood]
					dates=np.asarray([t0+eval('dt.timedelta({:s}=float(ti))'.format(timeunit)) for ti in t],np.datetime64)
					
					# temporal overlap
					#shaere interval
					keep=(dates>T[0]-np.timedelta64(2,'h')) & (dates<T[-1]+np.timedelta64(2,'h'))
					if keep.max(): #and not (name in exclude_stations):
						print('adding station '+name)
						#keep=(dates>T[0]-np.timedelta64(2,'h')) & (dates<T[-1]+np.timedelta64(2,'h'))
						t=t[keep]
						dates=dates[keep]
						
						#stations={'coord':[],'names':[],'MO':{}}	
						
						#dates,uinds=np.unique(dates,return_index=True)
						if 'TEMP' in ncs.variables.keys():
							temp=ncs['TEMP'][:][igood][keep]#[:,0]
							#temp=temp[uinds]
						else:
							temp=np.nan*np.ones(len(ncs['DEPH'][0,:]))
						if 'PSAL' in ncs.variables.keys():
							salt=ncs['PSAL'][:][igood][keep]#[:,0]
							#salt=salt[uinds]
						else:
							salt=np.nan*np.ones(len(ncs['DEPH'][0,:]))	
						
						stations['MO'][name]={'time':dates,'T':temp,'S':salt,'D':ncs['DEPH'][0,:]}
						
						if i==0:
							names.append(name)
							stations['coord'].append(coord)
						stations[setup_name]['names'].append(name)	
						if use_station_in[i]:
							stations[setup_names[i]]['coord'].append((staout.coords[nn,0],staout.coords[nn,1]))
						else:
							stations[setup_names[i]]['coord'].append((s.lon[nn],s.lat[nn]))
						stations[setup_names[i]]['nn'].append(nn)
						# select nearest neighbour
					else:
						print('not adding station '+name)
	print('done selecting and loading tide gauges, found: '+str(len(stations['coord'])) + ' stations meeting criterea')

stations['names']=list(stations['MO'].keys())
	
###### SCHISM read next neighbours profiles #############################
for i,setup_name in enumerate(setup_names):
	os.chdir(ncdir[i])
	if use_station_in[i]:
		if os.path.exists('station.in')	& os.path.exists(ncdir+'staout_1'):
			staout=schism_station_output()
			if 0: # based on coordiantes
				coords2=coords.copy()	
				coords=[]
				nn2=[]
				lldists2=[]
				for coord in coords2:
					D=np.sqrt((staout.coords[:,0]-coord[0])**2+(staout.coords[:,1]-coord[1])**2)
					nn2.append(np.argmin(D))
					lldists2.append(D[nn2[-1]])
				#np.asarray(lldists2)<=dtol]
			# based on name
			keep=[]
			for i,stat in enumerate(staout.stations):
				for name in names:
					if name in stat:
						keep.append(i)
			keep=np.asarray(keep)
			staout.coords=staout.coords[keep,:]
			staout.stations=np.asarray(staout.stations)[keep]
			for key in list(staout.station_out.keys())[:-1]:
				if	type(staout.station_out[key])==np.ndarray:
					staout.station_out[key]=staout.station_out[key][keep,:]
			s.plot_domain_boundaries()
			plt.plot(staout.coords[:,0],staout.coords[:,1],'ko')
			
		else:
			print('station output does not exist, exiting program')
			exit
	#from dask.diagnostics import ProgressBar
	#with ProgressBar():
	#	s.nc=xr.open_mfdataset(schismfiles[:5],combine='by_coords')
	else:
			# load profiles at nn
		nns=stations[setup_names[i]]['nn']
		stations[setup_names[i]]['profiles']=dict.fromkeys(('salinity','temperature','zCoordinates','time'))
		if newio[i]:
			for varname in 'salinity','temperature','zCoordinates':	
				accessi=access[i].ds[access[i].vardict[varname]][varname]	
				accessi=accessi.sel(nSCHISM_hgrid_node=nns)
				print('load {:s} profiles for setup {:s}'.format(varname,setup_names[i]))
				# use parallel writing
				stations[setup_names[i]]['profiles'][varname]=accessi.values #.values
			os.chdir(setupdir[i]) #[0]
			p=param('')
			reftime=np.asarray(dt.datetime(np.int(p.get_parameter('start_year')),\
			np.int(p.get_parameter('start_month')),\
			np.int(p.get_parameter('start_day')),\
			np.int(p.get_parameter('start_hour')),0,0),np.datetime64)
			T=reftime+access[0].ds['out2d']['time'].values *np.timedelta64(1,'s')
			stations[setup_names[i]]['time']=T
	
		else:
			setups[i].nc=setups[i].nc.sel(nSCHISM_hgrid_node=nns)
			oldnames={'salinity':'salt','temperature':'temp','zCoordinates':'zcor','time':'time'}
			for varname in 'salinity','temperature','zCoordinates','time':
				stations[setup_names[i]]['profiles'][varname]=setups[i].nc[oldnames[varname]].values
################################################################################################				



#modeldt=[np.diff(stations[setup_name]['time'][:2])/np.timedelta64(1,'s') for setup_name in setup_names]

## load data from schism next neighbours to TG statsions
#modeldt=[]
#dates=[]
#for i,setup_name in enumerate(setup_names):
#	print(setup_name)
#	s=setups[i]
#	if use_station_in[i]:
#		staout=output[i]
#		modeldt.append(np.timedelta64(staout.time[1]-staout.time[0]))
#		dates.append(np.asarray(staout.time,np.datetime64))   # convert to np.datetime64
#		stations[setup_name]['zeta']=staout.station_out['elev'][:,stations[setup_name]['nn']]
#		modeldt[i]=np.timedelta64(modeldt[i])
#		stations[setup_name]['time']=dates[i]#np.asarray(dates,np.datetime64)
#		lons,lats=staout.coords[:,0],staout.coords[:,1]
#	else:
#	
#		if newio[i]:
#			s.nc=access[i].get('elevation')
#			#s.nc.time[:]=T
#			s.nc2=s.nc.sel(nSCHISM_hgrid_node=stations[setup_name]['nn'])		
#			#,nSCHISM_vgrid_layers=-1
#			stations[setup_name]['time']=T
#			#stations[setup_name]['zeta']=s.nc2['elevation'].values
#			stations[setup_name]['zeta']=s.nc2[:].values
#		else:
#			#break
#			s.nc=access[i]
#			use_elev='elev' in s.nc.variables.keys()
#			s.nc2=s.nc.sel(nSCHISM_hgrid_node=stations[setup_name]['nn'],nSCHISM_vgrid_layers=-1)
#			n=len(schismfiles)	
#			print(str(n) + 'files')
#			if use_elev:	
#				stations[setup_name]['zeta']=s.nc2['elev'].values
#			else: # use zcor
#				stations[setup_name]['zeta']=s.nc2['zcor'].values
#			stations[setup_name]['time']=s.nc2['time'].values
#		modeldt.append(stations[setup_name]['time'][1]-stations[setup_name]['time'][0])
#		dates.append(stations[setup_name]['time'])
#	print('load model Tide Gauge Data')
#
	
# initiate file acces

#ncs={tag:xr.open_dataset(file) for tag,file in zip(names,marnet)}       # create netcdf acces
##############################################################################

for coord,name in zip(stations['coord'],stations['names']):
	xi,yi=coord
	plt.plot(xi,yi,'ro')
	namei=name
	if len(name)>5:
		namei=name[0]+name[1:-2][::2]+name[-2:]
	plt.text(xi,yi,' '+namei)
plt.title('Marnet Stations')
s.plot_domain_boundaries(append=True,latlon=True)
plt.savefig(outdir+'0_MarnetStations.png',dpi=300)


#ncks -v salt,temp,zcor schout_1.nc out.nc
# ncks -v time,salt,temp,zcor -d nSCHISM_hgrid_node,1,2,3 schout_1.nc out.nc

# escape ncks nccat

# loop in windows
#for ti in range(wndw,nt,wndw):
#	print(ti)
#	salt=np.concatenate((salt,s.nc['salt'][ti:ti+wndw,:,:][:,nn,:]),axis=0)
#	temp=np.concatenate((temp,s.nc['temp'][ti:ti+wndw,:,:][:,nn,:]),axis=0)
#	zcor=np.concatenate((zcor,s.nc['zcor'][ti:ti+wndw,:,:][:,nn,:]),axis=0)

# prevent error of un compiled data
#salt=np.ma.masked_array(salt,mask=[salt>1e36])
#temp=np.ma.masked_array(temp,mask=[temp>1e36])
#zcor=np.ma.masked_array(zcor,mask=[zcor>1e36])

#TT=np.tile(time[:nhours],[s.nc['zcor'].shape[2],1]).swapaxes(0,1)
#TT=np.tile(time,[s.nc['zcor'].shape[2],1]).swapaxes(0,1)
#sampmar=np.int(np.ceil(Tres/TresMarnet))  # subsample marnet as necessary for schism
#sampmod=np.int(np.ceil(Tres/TresMarnet))

# higher resolution marnet or model

# Plot Mar net stations
# monthly   - salt



def plot_profile_comp(t_data,data,t_mod,data_mod,setup_name,tag,label='Temperature',unit='[°C]'):
	# plot
	vmin=np.nanmin([data_mod.flatten().min(),data.flatten().min()]) 
	vmax=np.nanmax([data_mod.flatten().min(),data.flatten().max()]) 
	plt.clf()
	gs = gridspec.GridSpec(10, 1)
	plt.subplot(gs[0:2])
	plt.plot(t_data[0,:],data[0,:])#sampmar:nhours:sampmar
	plt.plot(t_mod,data_mod[:,-1])
	plt.legend((tag,setup_name),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5),frameon=False)	
	plt.ylabel(' '.join((label[:9]+'surf',unit)))
	plt.tick_params(labelbottom=False)
	plt.xlim((t_mod[0,0],t_mod[-1,0]))
	plt.grid()
	plt.subplot(gs[3:10])
	# plot schism profile
	plt.pcolor(t_mod,z_mod,data_mod,vmin=vmin,vmax=vmax)
	plt.xlim((t_mod[0,0],t_mod[-1,0]))
	plt.ylabel('H [m]')
	plt.title(tag)
	
	#if data.mask.min()==False:
	ph=plt.scatter(t_data,z_data,30,c=data,vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
	plt.gca().xaxis.set_major_formatter(myFmt)
	plt.gcf().autofmt_xdate()
	ch=plt.colorbar(orientation='horizontal',pad=0.2)
	ch.set_label(' '.join((label,unit)))
	plt.tight_layout()


def vert_interp(zs,z_mod,data):	
	nz=z_mod.shape[-1]
	intp_data=np.ones((z_mod.shape[0],len(zs)))*np.nan
	depths=np.ones((z_mod.shape[0],len(zs)))*np.nan
	# convert depths positive downard
	
	if (zs[:]>0).max():
		print('warning expecting data z values < 0 applying minus')
		zs=-zs
	for iz,zi in enumerate(zs):
		if (not np.ma.is_masked(zi)) and (zi <=0):
			# positive dpeths values 
			#diffs=np.abs(zi)+z_mod
			#ibelow=(diffs<0).sum(axis=1)-1
			#iabove=ibelow+1-(ibelow+1==nz)  # use 

			diffs=zi-z_mod
			ibelow=(diffs>0).sum(axis=1)-1
			iabove=ibelow+1-(ibelow+1==nz)  # use 
			
			#diffs[range(len(ibelow)),ibelow]
			#d1=np.abs(z_mod[range(len(iabove)),iabove]+zi) # oppsing signs
			#d2=np.abs(-z_mod[range(len(iabove)),ibelow]-zi)
			d1=np.abs(z_mod[range(len(iabove)),iabove]-zi) # oppsing signs
			d2=np.abs(-z_mod[range(len(iabove)),ibelow]+zi)
			#w1a=1/d1/(1/d1+1/d2)
			#w2a=1/d2/(1/d1+1/d2)
			w1=1-d1/(d1+d2)
			w2=1-w1
			
			depths[:,iz]=w1*z_mod[range(len(iabove)),iabove]+w2*z_mod[range(len(iabove)),ibelow]
			intp_data[:,iz]=w1*data[range(len(iabove)),iabove]+w2*data[range(len(iabove)),ibelow]
	
	return depths,intp_data			
	

	
doplots=1
if doplots:
	plt.figure(figsize=(8,6))
	for setup_name in setup_names:
		nz=stations[setup_name]['profiles']['temperature'].shape[-1]
		for i,tag in enumerate(names):
			print(tag)

			#zs=-stations['MO'][tag]['D']
			zs=stations['MO'][tag]['D']
			#salinity data available  at mooring
			obs_time=stations['MO'][tag]['time'] 
			mod_time=stations[setup_name]['time']
			data_inds=((mod_time[0] <=obs_time) & (mod_time[-1] >=obs_time ))
			#data_inds=((mod_time[0] <= obs_time) & (mod_time[-1] >= obs_time )) 
			obs_time=obs_time[data_inds]
			
			depth=stations['MO'][tag]['D']
			z_data,t_data=np.meshgrid(-depth,obs_time)
			
			#time interp
			#inds=((obs_time[0] <=time) & (obs_time[-1] >=time )) # overlap
			# select model data one before and after observations for interpolation
			#inds=np.arange(np.where(mod_time < obs_time[0])[-1][0],np.where(mod_time > obs_time[-1])[0][0]+1)
			inds=np.arange(np.where(mod_time <= obs_time[0])[-1][0],np.where(mod_time >= obs_time[-1])[0][0]+1)
			
			t_mod=np.tile(mod_time,(nz,1)).swapaxes(0,1)[inds,:]
			z_mod=stations[setup_name]['profiles']['zCoordinates'][inds,i,:]
			
			# temporal interpolation times
			ttemp=mod_time[inds]
			tin=(ttemp-obs_time[0])/np.timedelta64(1,'s') #(ttemp[1]-ttemp[0])
			tout=(obs_time-obs_time[0])/np.timedelta64(1,'s')#(time[1]-time[0])# /(Todates[1]-Todates[0])
			stations['MO'][tag]['time']=stations['MO'][tag]['time'][data_inds]

			stations[setup_name]['MO'][tag]={'time':obs_time,'T':[],'S':[],'D':depth}	
			# salinity
			if  np.isnan(stations['MO'][tag]['S']).min()==False:
				#plt.pause(0.1)
				data=stations['MO'][tag]['S'][data_inds,:]
				#schism	
				stations['MO'][tag]['S']=data
				data_mod=stations[setup_name]['profiles']['salinity'][inds,i,:]
				#z_mod=stations[setup_name]['profiles']['zCoordinates'][inds,i,:]
				plot_profile_comp(t_data,data,t_mod,data_mod,setup_name,tag,label='Salinity',unit='[psu]')
				plt.savefig(outdir+'1_Profile_salt_'+setup_name+tag+'.png',dpi=300)			
			
				# vertical interpolation
				depths,intp_data = vert_interp(zs,z_mod,data_mod)
				# temporal interpolation
				fintp=interp1d(tin, intp_data, axis=0)
								#stations[setup_name]['MO'][tag]['S']=fintp(tout)
				interpolated=fintp(tout)
				if len(stations['MO'][tag]['time']) != len(interpolated):
					from IPython import embed; embed()	      
				stations[setup_name]['MO'][tag]['S']=np.ma.masked_array(interpolated,mask=np.isnan(interpolated))
			
							
			# temperature
			#type(data)==np.ndarray or 
			
			if np.isnan(stations['MO'][tag]['T']).min()==False:
				#plt.pause(0.1)
				data=stations['MO'][tag]['T'][data_inds,:]
				stations['MO'][tag]['T']=data
				data_mod=stations[setup_name]['profiles']['temperature'][inds,i,:]
				#z_mod=stations[setup_name]['profiles']['zCoordinates'][inds,i,:]
				plot_profile_comp(t_data,data,t_mod,data_mod,setup_name,tag,label='Temperature',unit='[°C]')
				plt.savefig(outdir+'1_Profile_temp_'+setup_name+tag+'.png',dpi=300)

				# vertical interpolation
				depths,intp_data = vert_interp(zs,z_mod,data_mod)
				# temporal interpolation
				fintp=interp1d(tin, intp_data, axis=0)
				interpolated=fintp(tout)
				if len(stations['MO'][tag]['time']) != len(interpolated):
					from IPython import embed; embed()					
				stations[setup_name]['MO'][tag]['T']=np.ma.masked_array(interpolated,mask=np.isnan(interpolated))

for tag in stations[setup_name]['MO'].keys():				
	T=stations[setup_name]['MO'][tag]['T']
	S=stations[setup_name]['MO'][tag]['S']
	stations[setup_name]['MO'][tag]['T']=np.ma.masked_array(T,mask=np.isnan(T))
	stations[setup_name]['MO'][tag]['S']=np.ma.masked_array(S,mask=np.isnan(S))

			
# statistics ---
# interp model to data
T_schism={name:0 for name in names}
S_schism={name:0 for name in names}
Ztemp_schism={name:0 for name in names}
Zsalt_schism={name:0 for name in names}
timeout={name:0 for name in names}

	
# Same for amm15
class Amm15:
	def __init__(self, Tfile=None, Sfile=None):
		self.dir=Tfile[:Tfile.rfind('/')+1]
		self.tnc=Dataset(Tfile)
		self.snc=Dataset(Sfile)
		self.t0=dt.datetime.strptime(self.tnc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.tnc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.tnc['time'][:2]))[0]))*np.arange(len(self.tnc['time']))
		self.lon,self.lat=self.tnc['lon'][:],self.tnc['lat'][:]
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		coords=[]
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	

	def get_slab(self,time):		
		nn=np.argsort(np.abs(self.ts-time))[:2] # interpolated linearly
		self.deltaT=np.abs(self.ts[nn]-time)
		dts=np.asarray([dti.total_seconds() for dti in self.deltaT])
		w=1/dts/(1/dts).sum()
		self.T=self.tnc['thetao'][nn[0],:]*w[0]+self.tnc['thetao'][nn[1],:]*w[1]
		self.t=time #self.ts[nn]
		self.S=self.snc['so'][nn[0],:]*w[0]+self.snc['so'][nn[1],:]*w[1]
		
	def get_hov(self,coords):		
		s.nn_moc=self.tree.query(coords)[1]
		ii,jj=np.unravel_index(nn_moc,moc.LON.shape)
		
		self.Sprofile=moc.snc['so'][:][:,:,ii,jj]
		self.Tprofile=moc.tnc['thetao'][:][:,:,ii,jj]
		
	# my ocean
	def gname(self,date):	
		date2=date+dt.timedelta(days=1)	
		return 'metoffice_foam1_amm15_NWS_TEM_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day), 'metoffice_foam1_amm15_NWS_SAL_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day)

	def close(self):
		self.snc.close()
	def update(self,date):
		tname,sname=self.gname(date)
		Tfile=self.dir+tname
		Sfile=self.dir+sname
		self.tnc.close()
		self.snc.close()
		self.tnc=Dataset(Tfile)
		self.snc=Dataset(Sfile)
		self.t0=dt.datetime.strptime(self.tnc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.tnc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.tnc['time'][:2]))[0]))*np.arange(len(self.tnc['time']))

		
if add_amm15:		
	# my ocean
	tfiles=np.sort(glob.glob(oceandir+'*TEM*nc'))
	sfiles=np.sort(glob.glob(oceandir+'*SAL*nc'))
	moc=cmems(Sfile=sfiles[0],Tfile=tfiles[0])
	#nn_moc=moc.tree.query(coords)[1]
	#ii,jj=np.unravel_index(nn_moc,moc.LON.shape)
	moc.tnc=xr.open_mfdataset(tfiles)
	moc.snc=xr.open_mfdataset(sfiles)
	
	# Load Profiles at coordinates
	# Define target latitude and longitude (where weather stations might be)
	target_lon = xr.DataArray(np.asarray(coords)[:,0], dims="points")
	target_lat = xr.DataArray(np.asarray(coords)[:,1], dims="points")
	# Retrieve data at the grid cells nearest to the target latitudes and longitudes
	moc.tnc = moc.tnc.sel(lon=target_lon, lat=target_lat, method="nearest")
	moc.snc = moc.snc.sel(lon=target_lon, lat=target_lat, method="nearest")
	#load TOdata
	Z15=moc.tnc['depth'].values
	T15=moc.tnc['thetao'].values
	S15=moc.snc['so'].values
	time15=moc.tnc['time'].values

	
# temp
plt.figure()#(figsize=(25,33)) # full dina 4 widht in inch
#compare={tag:{'temp':np.nan,'salt':np.nan} for tag in names}
compare={setup_name: {tag:{'T':np.nan,'S':np.nan} for tag in names} for setup_name in setup_names}
for setup_name in setup_names:
	for key in compare[setup_name].keys():
		for key2 in compare[setup_name][key].keys():
			compare[setup_name][key][key2]={tag:np.nan for tag in  ['bias','rmse','cor','std1','std2']}
			
		
	plt.clf()
	for i,tag in enumerate(names):

		obs_time=stations['MO'][tag]['time']
		D=stations['MO'][tag]['D']
		labels=[' T [deg C]',' S [psu]']
		for var in ['T','S']:
			stations[setup_name]['MO'][tag]['time']
			
			Tobs=stations['MO'][tag][var]
			#Sobs=stations['MO'][tag]['S']
			
			if (~np.isnan(Tobs)).sum()==0:
				continue
			Tmod=stations[setup_name]['MO'][tag][var]
			#Smod=stations[setup_name]['MO'][tag]['S']

			#Tmod.shape == Tobs.shape
			#zmat=np.tile(stations['MO'][tag]['D'],[len(obs_time),1]).T
			
			plt.clf()
			plt.pause(0.001)
			plt.subplot(2,2,1)
			if quantile_range:
				#vmin=np.quantile(np.concatenate((T_schism[tag].flatten(),T[tag].flatten())),quantiles[0])
				#vmax=np.quantile(np.concatenate((T_schism[tag].flatten(),T[tag].flatten())),quantiles[1])
				Tmasked=np.ma.concatenate((Tmod.flatten(),Tobs.flatten()))
				validvalues=Tmasked[Tmasked.mask==False]
				vmin=np.quantile(validvalues,quantiles[0])
				vmax=np.quantile(validvalues,quantiles[1])

				#plt.pcolor(obs_time[tag],zmat,T[tag],vmin=vmin,vmax=vmax)
				if (Tobs.shape)[1]==1:
					plt.scatter(obs_time,stations['MO'][tag]['D']*np.ones(Tobs.shape),c=Tobs.T,vmin=vmin,vmax=vmax)
				else:	
					plt.pcolor(obs_time,stations['MO'][tag]['D'],Tobs.T,vmin=vmin,vmax=vmax)
			else:
				#plt.pcolor(obs_time[tag],zmat,T[tag])
				plt.pcolor(obs_time,stations['MO'][tag]['D'],Tobs.T)
				vmin, vmax = plt.gci().get_clim()
			plt.colorbar()
			plt.title(tag + ' T [deg C]')
			plt.tick_params(labelbottom=False)
			plt.subplot(2,2,2)
			if (Tobs.shape)[1]==1:
				plt.scatter(obs_time,stations['MO'][tag]['D']*np.ones(Tobs.shape),c=Tmod.T,vmin=vmin,vmax=vmax)		
			else:
				plt.pcolor(obs_time,stations['MO'][tag]['D'],Tmod.T,vmin=vmin,vmax=vmax)		
			#plt.pcolor(obs_time[tag],zmat,T_schism[tag],vmin=vmin,vmax=vmax)
			plt.colorbar()
			plt.title('schism' + labels[var=='S']) #' T [deg C]'
			plt.subplot(2,2,3)
			if quantile_range:
				#vmax=np.quantile(np.abs(T[tag]-T_schism[tag]),quantiles[1])
				delta=Tmod-Tobs
				vmax=np.quantile(delta[delta.mask==False],quantiles[1])
				#np.nanquantile(np.abs(delta).flatten(),quantiles[1])
				vmin=-vmax
				#plt.pcolor(obs_time[tag],zmat,T[tag]-T_schism[tag],vmin=vmin,vmax=vmax,cmap='jet')
				if (Tobs.shape)[1]==1:
					plt.scatter(obs_time,stations['MO'][tag]['D']*np.ones(Tobs.shape),c=delta,vmin=vmin,vmax=vmax,cmap='jet')		
				else:
					plt.pcolor(obs_time,stations['MO'][tag]['D'],(Tmod-Tobs).T,vmin=vmin,vmax=vmax,cmap='jet')
			else:
				#plt.pcolor(obs_time[tag],zmat,T[tag]-T_schism[tag],cmap='jet')
				plt.pcolor(obs_time,stations['MO'][tag]['D'],(Tmod-Tobs).T,cmap='jet')
			plt.colorbar()
			plt.title(tag + ' - schism')
			plt.gcf().autofmt_xdate()
			mng = plt.get_current_fig_manager()
			plt.tight_layout()
			plt.subplot(2,2,4)
			#bias=(T_schism[tag]-T[tag]).mean(axis=1)
			bias=(delta).mean(axis=0)
			#rmse=np.sqrt(((T[tag]-T_schism[tag])**2).mean(axis=1))
			rmse=np.sqrt(((delta)**2).mean(axis=0))
			#R=np.asarray([np.corrcoef(T[tag][j,:],T_schism[tag][j,:])[0,1] for j in range(T[tag].shape[0])])
			R=np.asarray([np.corrcoef(Tobs[:,j],Tmod[:,j])[0,1] for j in range(Tobs.shape[1])])
			#if len(T[tag][T[tag].mask==False])<1:
			if (type(Tobs) == np.ma.masked_array) and (len(Tobs[Tobs.mask==False])<1):
				bias*=np.nan
				rmse*=np.nan
				R*=np.nan
			
		
			ind=range(len(bias))
			if (type(Tobs) == np.ma.masked_array): 
				ind=bias.mask==False
			
			plt.bar(D[ind]-0.33,bias[ind],color='g', label='bias',width=0.33)
			plt.bar(D[ind], rmse[ind],color='b', label='RMSE',width=0.33)
			plt.bar(D[ind]+0.33, R[ind],color='k', label='COR',width=0.33)


			#plt.title('total Bias/rmse {:.2f}/{:.2f}'.format((Tmod-Tobs).mean(),np.sqrt( ((T[tag]-T_schism[tag])**2).quantile95()))) 
			
			compare[setup_name][tag][var]['bias']=bias		
			compare[setup_name][tag][var]['rmse']=rmse
			#compare[tag]['temp']['std1']=np.std(T[tag],axis=1)
			compare[setup_name][tag][var]['std1']=np.std(Tobs,axis=0)
			#compare[tag]['temp']['std2']=np.std(T_schism[tag],axis=1)	
			compare[setup_name][tag][var]['std2']=np.std(Tmod,axis=0)	
			compare[setup_name][tag][var]['cor']=R
			for zi,biasi,rmsei,cori in zip(D,bias,rmse,R):
				#plt.text(zi-0.6,biasi,'{:.2f}'.format(biasi),rotation=90 )
				#plt.text(zi+0.6,rmsei,'{:.2f}'.format(rmsei),rotation=90 )
				plt.text(zi-0.33,biasi,'{:.2f}'.format(biasi),rotation=90 )
				plt.text(zi,rmsei,'{:.2f}'.format(rmsei),rotation=90 )
				plt.text(zi+0.33,cori,'{:.2f}'.format(cori),rotation=90 )
			plt.xlabel('depth')
			plt.legend(('bias','rmse','cor'),frameon=False)
			plt.grid()
			plt.tight_layout()
			plt.savefig(outdir+setup_name+'3_'+tag+var+'_statistics.png',dpi=300)		


import pandas as pd	
def write_table(varname,stations,compare,names,setup_name):
	d = {'station': names}	
	f=open(setup_name+varname+'.txt','w')
	f.write('quant\stat')
	for tag in names:
		f.write('\t {:8s}'.format(tag))
	f.write('\n')	
	#zs=np.unique(np.hstack([np.asarray(Zsalt[tag]) for tag in Zsalt.keys() ]))
	zs=np.unique(np.hstack([np.asarray(stations['MO'][tag]['D']) for tag in stations['MO'].keys() ]))
	for key in compare[tag][varname].keys():
		for zi in zs:
			vals=[]
			f.write(key + ' {:8.2f}'.format(zi))
			for tag in names:
				zstat=stations['MO'][tag]['D']
				if zi in zstat: #Zsalt[tag]:
					#val=np.asarray(compare[tag][varname][key])[Zsalt[tag]==zi][0]
					try:
						val=np.asarray(compare[tag][varname][key])[zstat==zi][0]
					except:
						val=np.nan
						#from IPython import embed; embed()
				else:
					val=np.nan
				f.write('\t {:8.4f}'.format(val))
				vals.append(val)
			vals=np.asarray(vals)
			vals[vals==0.0]=np.nan
			d[key+str(np.int(zi))]=vals
			f.write('\n')		
	f.close()
	return pd.DataFrame(data=d)
	
	

#dfsalt=write_table('salt',Zsalt,compare,names)
#dftemp=write_table('temp',Ztemp,compare,names)
dfsalt=dict.fromkeys(setup_names)
dftemp=dict.fromkeys(setup_names)
for setup_name in setup_names:
	dfsalt[setup_name]=write_table('S',stations,compare[setup_name],names,setup_name)
	dftemp[setup_name]=write_table('T',stations,compare[setup_name],names,setup_name)

# not working currently
##zs=np.unique(np.hstack([np.asarray(Zsalt[tag]) for tag in Zsalt.keys() ]))
#zs=np.unique(np.hstack([np.asarray(stations['MO'][tag]['D']) for tag in
## is this used?
#for key in compare[tag]['salt'].keys():
#	for zi in zs:
#		vals=[]
#		for tag in names:
#			zsstat=stations['MO'][tag]['D']
#			if zi in zsstat: #Zsalt[tag]:
#				#val=np.asarray(compare[tag]['salt'][key])[Zsalt[tag]==zi][0]
#				val=np.asarray(compare[tag]['salt'][key])[zsstat==zi][0]
#			else:
#				val=np.nan
#			vals.append(val)
#		d[key+str(np.int(zi))]=vals
#

for setup_name in setup_names:
	set_fontsize(fac=1.8)		
	zs=np.asarray(np.unique(np.hstack([np.asarray(stations['MO'][tag]['D']) for tag in stations['MO'].keys() ])),int)
	xbar=range(len(names))
	plt.figure()
	for zi in zs:
		phs=[]
		
		plt.clf()
		plt.pause(0.001)
		phs.append(plt.bar(xbar,dfsalt[setup_name]['bias'+str(np.int(zi))].values,width=-0.25,align='edge',label='bias',tick_label=names))
		phs.append(plt.bar(xbar,dfsalt[setup_name]['rmse'+str(np.int(zi))].values,width=0.25,align='edge',label='bias',tick_label=names))
		plt.grid()
		plt.legend(('bias','rmse'))
		plt.ylabel('salinity [-]')
		plt.tight_layout()
		plt.savefig('salt_error_bars_{:d}.png'.format((zi)),dpi=300)
	

	#zs=np.asarray(np.unique(np.hstack([np.asarray(Ztemp[tag]) for tag in Ztemp.keys() ])),int)
	#xbar=range(len(names))
	for zi in zs:
		phs=[]
		plt.clf()
		plt.pause(0.001)
		phs.append(plt.bar(xbar,dftemp[setup_name]['bias'+str(np.int(zi))].values,width=-0.25,align='edge',label='bias',tick_label=names))
		phs.append(plt.bar(xbar,dftemp[setup_name]['rmse'+str(np.int(zi))].values,width=0.25,align='edge',label='bias',tick_label=names))
		plt.grid()
		plt.legend(('bias','rmse'))
		plt.ylabel('temperature [°C]')
		plt.tight_layout()
		plt.savefig('temp_error_bars_{:d}.png'.format((zi)),dpi=300)
	
set_fontsize(fac=1)		

colors=plt.cm.tab10(range(len(setup_names)))  

# top botoopm
label=['top','bottom']
for nr,i in enumerate([0,-1]):
	compare0=compare[setup_names[0]]
	samples=[[ [compare0[tag]['T']['std2'][i]/compare0[tag]['T']['std1'][i],compare0[tag]['T']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'] for tag in names if ~np.isnan(compare0[tag]['T']['bias']).max()]]
	dia=plotTaylor(samples=samples[0])

	#add data
	for isetup,key in enumerate(setup_names[1:]):
		compare1=compare[setup_names[isetup+1]]
		samples.append([ [compare1[tag]['T']['std2'][i]/compare1[tag]['T']['std1'][i],compare1[tag]['T']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'] for tag in names if ~np.isnan(compare1[tag]['T']['bias']).max()])	
	
	
	for isetup,key in enumerate(setup_names[1:]):
		#Add models to Taylor diagram
		for i,(stddev, corrcoef, name) in enumerate(samples[isetup+1]):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)	
	
	plt.suptitle('Temperature ' +label[nr])
	plt.savefig(outdir+'5_taylor_temp'+label[nr]+'.png',dpi=300)
	plt.close()

# top botoopm
label=['top','bottom']
for nr,i in enumerate([0,-1]):

	compare0=compare[setup_names[0]]
	samples=[[ [compare0[tag]['S']['std2'][i]/compare0[tag]['S']['std1'][i],compare0[tag]['S']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'] for tag in names if ~np.isnan(compare0[tag]['S']['bias']).max()]]
	dia=plotTaylor(samples=samples[0])

	#add data
	for isetup,key in enumerate(setup_names[1:]):
		compare1=compare[setup_names[isetup+1]]
		samples.append([ [compare1[tag]['S']['std2'][i]/compare1[tag]['S']['std1'][i],compare1[tag]['S']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'] for tag in names if ~np.isnan(compare1[tag]['S']['bias']).max()])	
	
	
	for isetup,key in enumerate(setup_names[1:]):
		#Add models to Taylor diagram
		for i,(stddev, corrcoef, name) in enumerate(samples[isetup+1]):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)	

	plt.suptitle('Salinity ' +label[nr])
	plt.savefig(outdir+'5_taylor_salt'+label[nr]+'.png',dpi=300)
	plt.close()


	
# comparison with amm15	
if add_amm15: 
	 
	T15={name:0 for name in names}
	S15={name:0 for name in names}
	obs_time15={name:0 for name in names}
	Z15=moc.snc['depth'][:].values	
	
	
	for i,name in enumerate(names):
		S15[tag]=moc.snc['so'][:,:,i]
		T15[tag]=moc.snc['thetao'][:,:,i]

	############### plot itner comparison
	plt.figure()
	for i,tag in enumerate(names):
		print(tag)
		
		# salt grid for hovmoeller N
		XX=-np.tile(Zsalt[tag].T,[nhours,1])
		YY=np.tile(obs_time[tag][23:23+nhours],[len(Zsalt[tag]),1]).swapaxes(0,1)
		
		# salt
		# range
		data=S[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([salt[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([salt[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		#gs = gridspec.GridSpec(10, 1)
		#plt.subplot(gs[0:2])
		#plt.plot(obs_time[tag][23:23+nhours],data[0,23:23+nhours])
		#plt.plot(time[:nhours],salt[:,i,-1])
		#plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		#plt.ylabel('Surf Salt [psu]')
		#plt.tick_params(labelbottom=False)

		#plt.subplot(gs[3:10])
		
		plt.subplot(2,1,1)
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],salt[:,i,:],vmin=vmin,vmax=vmax)
		#ch=plt.colorbar(orientation='horizontal')
		#ch.set_label('Salinity [psu]')
		plt.ylabel('H [m]')
		plt.title(tag + 'SCHISM')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.subplot(2,1,2)
		# plot amm15 profile
		iend=(T15[tag][:,0].mask==False).sum()
		plt.pcolormesh(obs_time15[tag],-Z15[:iend],S15[tag][:iend,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Salinity [psu]')
		plt.ylabel('H [m]')
		plt.title(tag + ' Amm15')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		
		plt.savefig(outdir+'Profile_salt2_'+tag+'.png',dpi=300)
		
		
		# temp
		# salt grid for hovmoeller 
		XX=-np.tile(Ztemp[tag].T,[nhours,1])
		YY=np.tile(obs_time[tag][23:23+nhours],[len(Ztemp[tag]),1]).swapaxes(0,1)

		
		data=T[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([temp[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([temp[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		#gs = gridspec.GridSpec(10, 1)
		#plt.subplot(gs[0:2])
		#plt.plot(obs_time[tag][23:23+nhours],data[0,23:23+nhours])
		#plt.plot(time[:nhours],temp[:,i,-1])
		#plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		#plt.ylabel('Tsurf [deg C]')
		#plt.tick_params(labelbottom=False)

		#plt.subplot(gs[3:10])
		plt.subplot(2,1,1)
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],temp[:,i,:],vmin=vmin,vmax=vmax)
		#ch=plt.colorbar(orientation='horizontal')
		#ch.set_label('Temperature [deg C]')
		plt.ylabel('H [m]')
		plt.title(tag + 'SCHISM')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		
		plt.subplot(2,1,2)
		# plot amm15 profile
		plt.pcolormesh(obs_time15[tag],-Z15[:iend],T15[tag][:iend,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Temperature [deg C]')
		plt.ylabel('H [m]')
		plt.title(tag + ' Amm5')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.savefig(outdir+'Profile_temp2_'+tag+'.png',dpi=300)

#
## top botoopm
#label=['top','bottom']
#for nr,i in enumerate(range(10)):#[0,-1]
#	samples=[ [compare[tag]['T']['std2'][i]/compare[tag]['T']['std1'][i],compare[tag]['T']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'] for tag in names if (i<=len(stations['MO'][tag]['D'])) & (~np.isnan(compare[tag]['T']['bias']).max())]
#	plotTaylor(samples=samples)
#	plt.suptitle('Temperature ' +label[nr])
#	plt.savefig(outdir+'6_all_taylor_temp'+label[nr]+'.png',dpi=300)
#	plt.close()		
#		
#
## top botoopm
#label=['top','bottom']
#samples=[]
#for nr,i in enumerate(range(10)):#[0,-1]
#	for tag in names:
#		try:
#			samples.append( [compare[tag]['T']['std2'][i]/compare[tag]['T']['std1'][i],compare[tag]['T']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'])
#		except:
#			pass
#plotTaylor(samples=samples)		
##compare[tag]['T']['bias'].max()]
#plotTaylor(samples=samples)
#plt.suptitle('Temperature ' +label[nr])		
#plt.savefig(outdir+'6_all_taylor_temp'+label[nr]+'.png',dpi=300)		
#
## top botoopm
#label=['top','bottom']
#samples=[]
#for nr,i in enumerate(range(10)):#[0,-1]
#	for tag in names:
#		try:
#			samples.append( [compare[tag]['S']['std2'][i]/compare[tag]['S']['std1'][i],compare[tag]['S']['cor'][i] ,tag+' '+str(stations['MO'][tag]['D'][i])+'m'])
#		except:
#			pass
#plotTaylor(samples=samples)		
##plt.suptitle('Salinity ' +label[nr])		
#plt.savefig(outdir+'6_all_taylor_salt'+label[nr]+'.png',dpi=300)			
		
# create latex + pdf
if put_pics_to_texdoc:
		print('generating tex doc and pdf')
		os.chdir(outdir)
		techit(latexname,latextitle,latextext)
		if os.path.isfile(latexname[:latexname.rindex('.')]+'.pdf'):
			 os.remove(latexname[:latexname.rindex('.')]+'.pdf')
		os.system('pdflatex '+latexname)
		print('done generating '+latexname[:latexname.rindex('.')]+'.pdf')
