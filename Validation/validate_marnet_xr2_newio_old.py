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
import xarray as xr
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
import pickle
import datetime as dt
from netCDF4 import Dataset,MFDataset
from matplotlib import gridspec  # unequal subplots
from scipy.interpolate import interp1d
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities//Lib/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')


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


#### I have added 1 day to compensate error in falsely set param.in start date
## updated hard coded timesteps with automatically determined ones

########## settings #################################
# directories
mardir='/gpfs/work/ksddata/observation/insitu/Marnet_Sebastian/' 					  # myocean
mardir='/work/gg0028/ksddata/insitu/CMEMS/NorthWestShelf/MO/'
setupdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir=setupdir+'/outputs_all/'
year=2017	    # currently data between 2011 and 2018				
dtol=0.05           										   # distance tolerance in degree lon/lat 


#oceandir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/'# amm15 dir

outdir=setupdir+'/marnet_valid/'	   # output directory where images will be stored
if not os.path.exists(outdir): os.mkdir(outdir) 

Tmax=45
Smax=49
add_amm15=False  											# add plot comparing with Amm15 | Amm15 profiles 
use_station_in=False					  # True: use station output False: extract from netcdf (slow)
put_pics_to_texdoc=True    										# images will be put in tex document
latexname='HRBallje.tex' #'GB1_5Ems_marnet.tex'										# in alphanumerical order 
latextitle='HRBallje 2011 vs Marnet' #'GB 2011 + FW 1.5EMS  in EFWS vs Marnet'
latextext= str(year) + ' validation against Marnet'

# reduce color scales to qunatile rnage , bypass outliers affectiong scale
quantile_range=True
quantiles=[0.1, 0.99]

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



#plt.ion() 
######### load SCHISM setup   ##################################
cwd=os.getcwd()
os.chdir(setupdir)
s=schism_setup()
s.nntree = cKDTree(list(zip(s.lon,s.lat))) 


######### initialize marnet Data located in SCHISM Domain #######
marnet=glob.glob(mardir+'*_{:d}.nc'.format(year))	
#marnet=glob.glob(mardir+'*{:d}*.nc'.format(year))	
print('searching coordinates for available data in year {:d}'.format(year))
coords=[]
stations=[]
for station in marnet:
	nc=Dataset(station)
	coords.append((nc['ppos'][1],nc['ppos'][0]))
	nc.close()
print('done selecting tide gauges, found: '+str(len(coords)) + ' stations meeting criterea')
#----  reduce to data in domain as given by dtol -----------
lldists,nn=s.nntree.query(coords)		
nn=nn[lldists<=dtol]
marnet=np.asarray(marnet)[lldists<=dtol]
coords=np.asarray(coords)[lldists<=dtol]	
names=[ file[file.rindex('/')+1:file.rindex('_')]  for file in marnet]
###### SCHISM read next neighbours profiles #############################
ref=dt.datetime.now()
if use_station_in:
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
#		
else:
	schismfiles=[] 
	for iorder in range(6): # check for schout_nc files until 99999
		schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
	nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
	# default SLOW - why dask problems
	schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
	#s.nc=xr.open_mfdataset(schismfiles,combine='by_coords')
	#alt
	all_dsets = [xr.open_dataset(fname).chunk() for fname in schismfiles]
	s.nc=xr.concat(all_dsets, dim='time')
	#ds_concat = 
	#%%timeit ts = ds_concat.data[:, 0, 0].load()
	
	#nrs=list(np.asarray(nrs)[np.argsort(nrs)])
	s.nc2=s.nc.sel(nSCHISM_hgrid_node=nn)
	time=s.nc['time'][:].values
################################################################################################



	
# initiate file acces

ncs={tag:xr.open_dataset(file) for tag,file in zip(names,marnet)}       # create netcdf acces
##############################################################################

coords=[]
for item in ncs.items():
	tag,nc=item
	yi,xi=nc['ppos'][:]
	coords.append((xi,yi))
	plt.plot(xi,yi,'ro')
	plt.text(xi,yi,' '+tag)
plt.title('Marnet Stations')
s.plot_domain_boundaries(append=True)
plt.savefig(outdir+'0_MarnetStations.png',dpi=300)

# Data put Data in TIme Depth Matrix
T={name:0 for name in names}
S={name:0 for name in names}
Ztemp={name:0 for name in names}
Zsalt={name:0 for name in names}
Time={name:0 for name in names}
for item in ncs.items():
	tag,nc=item
	ncv=nc.variables
	Ti=[]
	Zi_temp=[]
	Zi_salt=[]
	Si=[]
	Time[tag]=nc['time'][:]
	for key in ncv.keys():
		if 'WT' in key:
			print(key)
			Zi_temp.append(float(key[2:key.rindex('m')]))
			Ti+=[ncv[key][:]]
		if 'WS' in key:
			print(key)
			Si+=[ncv[key][:]]
			Zi_salt.append(float(key[2:key.rindex('m')]))
	Ztemp[tag]=np.asarray(Zi_temp)			
	Zsalt[tag]=np.asarray(Zi_salt)			
	#T[tag]=np.asarray(Ti)			
	#S[tag]=np.asarray(Si)			
	T[tag]=np.ma.masked_array(Ti,(np.asarray(Ti)==-10) | (np.asarray(Ti)>Tmax))			
	S[tag]=np.ma.masked_array(Si,(np.asarray(Si)<0) | (np.asarray(Si)>Smax))			

	
	
#S[tag]
#T[tag]



nt=len(s.nc['time'])
wndw=np.int(nt/len(schismfiles))
Tres=np.diff(s.nc['time'][:2])[0]    #seconds schism temporal resolution 
nhours=np.int((s.nc['time'][-1]-s.nc['time'][0])/np.timedelta64(1,'h'))+1 # acutally  steps ' length hours schism
TresMarnet=(Time[names[0]][1]-Time[names[0]][0]) #.total_seconds() #seconds resolution marnet hours

samplingMarnet=np.float(Tres/TresMarnet)
if samplingMarnet < 1:
	samplingSchism=np.int(1/samplingMarnet)
	samplingMarnet=1
else:
	samplingSchism=1
	samplingMarnet=np.int(samplingMarnet)
# resample schism to marnet / marnet to schisn
#if Tres < TresMarnet:
#	s.nc=s.nc.interp(time=Time)
#else



# this is slow
# load schism data
salt=s.nc2['salt'][:].values 
temp=s.nc2['temp'][:].values  
zcor=s.nc2['zcor'][:].values

# Fastest method to load ts ?temp
#import time
#tstart = time.time()
#for i in range(5):
#	temp1=s.nc['temp'][:wndw,:,:][:,nn,:].values	#[:wndw,:,:][:,nn,:]  
#	tend = time.time()
#print(tend - tstart)
#tstart = time.time()  #fstest
#for i in range(5):
#	temp2=s.nc2['temp'][:wndw,:,:].values	#[:wndw,:,:][:,nn,:] 
#	tend = time.time()
#print(tend - tstart)
#for i in range(5):
#	temp2=s.nc3['temp'][:wndw,:,:][:,nn,:]	#[:wndw,:,:][:,nn,:] 
#	tend = time.time()
#print(tend - tstart)

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
salt=np.ma.masked_array(salt,mask=[salt>1e36])
temp=np.ma.masked_array(temp,mask=[temp>1e36])
zcor=np.ma.masked_array(zcor,mask=[zcor>1e36])

#TT=np.tile(time[:nhours],[s.nc['zcor'].shape[2],1]).swapaxes(0,1)
TT=np.tile(time,[s.nc['zcor'].shape[2],1]).swapaxes(0,1)
sampmar=np.int(np.ceil(Tres/TresMarnet))  # subsample marnet as necessary for schism
#sampmod=np.int(np.ceil(Tres/TresMarnet))

# higher resolution marnet or model

#check time overlap otherwise remove station
for i,tag in enumerate(names):
	print('checking temporal overlap for station ' + tag)
	inds=((time[0] <=Time[tag] ) & (time[-1] >=Time[tag] )).values
	if inds.sum()< 1:
		print('no overlap for station removing ' + tag)
		names.remove(tag)
# Plot Mar net stations
# monthly   - salt
doplots=1
if doplots:
	#plt.figure()
	plt.figure(figsize=(8,6))
	for i,tag in enumerate(names):
		print(tag)
		#YY=np.tile(Time[tag][sampmar:nhours:sampmar],[len(Zsalt[tag]),1]).swapaxes(0,1)
		YY=np.tile(Time[tag][::samplingMarnet],[len(Zsalt[tag]),1]).swapaxes(0,1)
		XX=-np.tile(Zsalt[tag].T,[YY.shape[0],1])		
		
		# salt
		# range
		data=S[tag]
		#data=np.ma.masked_array(data,mask=(data<0) | (data > 49)) # mask  Salt
		#vmin=np.nanmin([salt[:,i,:].flatten().min(),data[:,sampmar:nhours:sampmar].swapaxes(0,1).flatten().min()])
		#vmax=np.nanmax([salt[:,i,:].flatten().max(),data[:,sampmar:nhours:sampmar].swapaxes(0,1).flatten().max()])
		#inds=(time[0] <=Time[tag][::samplingMarnet] ) & (time[-1] >=Time[tag][::samplingMarnet] )
		inds=((time[0] <=Time[tag] ) & (time[-1] >=Time[tag] )).values
		
		# reduce to model time overlap
		data=data[:,inds]
		YY2=YY[inds,:]
		XX2=XX[inds,:]
		Time2=Time[tag][inds].copy()
		
		vmin=np.nanmin([salt[:,i,:].flatten().min(),data[:,::samplingMarnet].swapaxes(0,1).flatten().min()])
		vmax=np.nanmax([salt[:,i,:].flatten().max(),data[:,::samplingMarnet].swapaxes(0,1).flatten().max()])
		
		# range
		#plt.close()
		plt.pause(0.001)
		plt.clf()
		gs = gridspec.GridSpec(10, 1)
		plt.subplot(gs[0:2])
		#plt.plot(Time[tag][sampmar:nhours:sampmar],data[0,sampmar:nhours:sampmar])
		#plt.plot(Time[tag][::samplingMarnet],data[0,::samplingMarnet])
		plt.plot(Time2[::samplingMarnet],data[0,::samplingMarnet])
		
		#plt.plot(time[:nhours],salt[:,i,-1])#
		#plt.plot(time[::samplingSchism],salt[::samplingSchism,i,-1])
		plt.plot(time,salt[:,i,-1])
		plt.xlim((time[0],time[-1]))
		
		plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		plt.ylabel('Surf Salt [psu]')
		plt.tick_params(labelbottom=False)

		plt.subplot(gs[3:10])
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],salt[:,i,:],vmin=vmin,vmax=vmax)
		#plt.pcolor(TT,zcor[:,i,:],salt[:,i,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Salinity [psu]')
		plt.ylabel('H [m]')
		plt.title(tag)
		
		# add marnet data as scatter
		sampling=6
		#colors=data[:,sampmar:nhours:sampmar].swapaxes(0,1)[::sampling,:].flatten()
		colors=data[:,::samplingMarnet].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			#ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
			ph=plt.scatter(YY2[::sampling,:].flatten()[ivalid],XX2[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.savefig(outdir+'1_Profile_salt_'+tag+'.png',dpi=300)
		
		# temperature
		#YY=np.tile(Time[tag][sampmar:nhours:sampmar],[len(Ztemp[tag]),1]).swapaxes(0,1)
		YY=np.tile(Time[tag][::samplingMarnet],[len(Ztemp[tag]),1]).swapaxes(0,1)
		XX=-np.tile(Ztemp[tag].T,[YY.shape[0],1])		
		
		data=T[tag]
		
		data=data[:,inds]
		YY2=YY[inds,:]
		XX2=XX[inds,:]
		Time2=Time[tag][inds].copy()
		
		
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.nanmin([temp[:,i,:].flatten().min(),data[:,::samplingMarnet].swapaxes(0,1).flatten().min()]) #sampmar:nhours:sampmar
		vmax=np.nanmax([temp[:,i,:].flatten().max(),data[:,::samplingMarnet].swapaxes(0,1).flatten().max()])#sampmar:nhours:sampmar
		
		# range
		
		plt.pause(0.001)
		plt.clf()
		gs = gridspec.GridSpec(10, 1)
		plt.subplot(gs[0:2])
		plt.plot(Time2[::samplingMarnet],data[0,::samplingMarnet])#sampmar:nhours:sampmar
		plt.plot(time[:],temp[:,i,-1])
		plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		plt.ylabel('Tsurf [deg C]')
		plt.tick_params(labelbottom=False)

		plt.subplot(gs[3:10])
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],temp[:,i,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Temperature [deg C]')
		plt.ylabel('H [m]')
		plt.title(tag)
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,::samplingMarnet].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY2[::sampling,:].flatten()[ivalid],XX2[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.savefig(outdir+'1_Profile_temp_'+tag+'.png',dpi=300)

	
# statistics ---
# interp model to data
T_schism={name:0 for name in names}
S_schism={name:0 for name in names}
Ztemp_schism={name:0 for name in names}
Zsalt_schism={name:0 for name in names}
timeout={name:0 for name in names}
for i,tag in enumerate(names):
	
	# interpolate to layers
	
	# salt
	Si=[]
	Zi=[]
	for zi in Zsalt[tag]:
		diffs=zi+zcor[:,i,:]
		ibelow=(diffs<0).sum(axis=1)-1
		iabove=ibelow+1
		diffs[range(len(ibelow)),ibelow]
		d1=zcor[range(len(iabove)),i,iabove]+zi
		d2=-zcor[range(len(iabove)),i,ibelow]-zi
		w1=1/d1/(1/d1+1/d2)
		w2=1/d2/(1/d1+1/d2)
		D=w1*zcor[range(len(iabove)),i,iabove]+w2*zcor[range(len(iabove)),i,ibelow]
		Zi.append(D)
		Si.append(w1*salt[range(len(iabove)),i,iabove]+w2*salt[range(len(iabove)),i,ibelow])
	Zsalt_schism[tag]=np.asarray(Zi)
	S_schism[tag]=np.asarray(Si)
	
	# temp
	Ti=[]
	Zi=[]
	for zi in Ztemp[tag]:
		diffs=zi+zcor[:,i,:]
		ibelow=(diffs<0).sum(axis=1)-1
		iabove=ibelow+1
		diffs[range(len(ibelow)),ibelow]
		d1=zcor[range(len(iabove)),i,iabove]+zi
		d2=-zcor[range(len(iabove)),i,ibelow]-zi
		w1=1/d1/(1/d1+1/d2)
		w2=1/d2/(1/d1+1/d2)
		D=w1*zcor[range(len(iabove)),i,iabove]+w2*zcor[range(len(iabove)),i,ibelow]
		Zi.append(D)
		Ti.append(w1*temp[range(len(iabove)),i,iabove]+w2*temp[range(len(iabove)),i,ibelow])
	Ztemp_schism[tag]=np.asarray(Zi)
	T_schism[tag]=np.asarray(Ti)

	# temporal interpolation to data
	Todates=np.asarray(Time[tag])
	ilast=(Todates<=time[nt-1]).sum()
	ifirst=(Todates<=time[0]).sum()
	Todates=Todates[ifirst:ilast]


	tin=(time-time[0])/(time[1]-time[0])
	tout=(Todates-time[0])/(time[1]-time[0])# /(Todates[1]-Todates[0])
	#tin=[(timei-time[0]).total_seconds() for timei in time]
	#tout=[(timei-time[0]).total_seconds() for timei in Todates]

	timeout[tag]=Todates	
	fintp=interp1d(tin[:nt], S_schism[tag], axis=1)
	S_schism[tag]=fintp(tout)
	fintp=interp1d(tin[:nt], T_schism[tag], axis=1)
	T_schism[tag]=fintp(tout)
	
	# limit data to smae range
	# Data put Data in TIme Depth Matrix
	T[tag]=T[tag][:,ifirst:ilast]
	S[tag]=S[tag][:,ifirst:ilast]
	Time[tag]=Todates
	

	
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
compare={tag:{'temp':np.nan,'salt':np.nan} for tag in names}
for key in compare.keys():
	for key2 in compare[key].keys():
		compare[key][key2]={tag:np.nan for tag in  ['bias','rmse','cor','std1','std2']}
plt.clf()
for i,tag in enumerate(names):
	zmat=np.tile(Ztemp[tag],[len(Time[tag]),1]).T
	plt.clf()
	plt.pause(0.001)
	plt.subplot(2,2,1)
	if quantile_range:
		vmin=np.quantile(np.concatenate((T_schism[tag].flatten(),T[tag].flatten())),quantiles[0])
		vmax=np.quantile(np.concatenate((T_schism[tag].flatten(),T[tag].flatten())),quantiles[1])
		plt.pcolor(Time[tag],zmat,T[tag],vmin=vmin,vmax=vmax)
	else:
		plt.pcolor(Time[tag],zmat,T[tag])
		vmin, vmax = plt.gci().get_clim()
	plt.colorbar()
	plt.title(tag + ' T [deg C]')
	plt.tick_params(labelbottom=False)
	plt.subplot(2,2,2)
	plt.pcolor(Time[tag],zmat,T_schism[tag],vmin=vmin,vmax=vmax)
	plt.colorbar()
	plt.title('schism' + ' T [deg C]')
	plt.subplot(2,2,3)
	if quantile_range:
		vmax=np.quantile(np.abs(T[tag]-T_schism[tag]),quantiles[1])
		vmin=-vmax
		plt.pcolor(Time[tag],zmat,T[tag]-T_schism[tag],vmin=vmin,vmax=vmax,cmap='jet')
	else:
		plt.pcolor(Time[tag],zmat,T[tag]-T_schism[tag],cmap='jet')
	plt.colorbar()
	plt.title(tag + ' - schism')
	plt.gcf().autofmt_xdate()
	mng = plt.get_current_fig_manager()
	plt.tight_layout()
	plt.subplot(2,2,4)
	bias=(T_schism[tag]-T[tag]).mean(axis=1)
	rmse=np.sqrt(((T[tag]-T_schism[tag])**2).mean(axis=1))
	R=np.asarray([np.corrcoef(T[tag][j,:],T_schism[tag][j,:])[0,1] for j in range(T[tag].shape[0])])
	if len(T[tag][T[tag].mask==False])<1:
		bias*=np.nan
		rmse*=np.nan
		R*=np.nan

		
	plt.bar(Ztemp[tag]-0.5,bias, label='bias')
	plt.bar(Ztemp[tag]+0.5, rmse,color='g', label='RMST')
	plt.title('total Bias/rmse {:.2f}/{:.2f}'.format((T[tag]-T_schism[tag]).mean(),np.sqrt( ((T[tag]-T_schism[tag])**2).mean()))) 
	compare[tag]['temp']['bias']=bias		
	compare[tag]['temp']['rmse']=rmse
	compare[tag]['temp']['std1']=np.std(T[tag],axis=1)
	compare[tag]['temp']['std2']=np.std(T_schism[tag],axis=1)	
	compare[tag]['temp']['cor']=R
	for zi,biasi,rmsei in zip(Ztemp[tag],bias,rmse):
		plt.text(zi-0.6,biasi,'{:.2f}'.format(biasi),rotation=90 )
		plt.text(zi+0.6,rmsei,'{:.2f}'.format(rmsei),rotation=90 )
	plt.xlabel('depth')
	plt.legend(('bias','rmse'),frameon=False)
	plt.grid()
	plt.tight_layout()
	plt.savefig(outdir+'3_'+tag+'Temp_statistics.png',dpi=300)


# salt
plt.clf()
for i,tag in enumerate(names):
	#zmat=np.tile(Zsalt[tag],[nhours-1,1]).T
	zmat=np.tile(Zsalt[tag],[len(Time[tag]),1]).T
	plt.pause(0.001)
	plt.clf()
	plt.subplot(2,2,1)
	if quantile_range:
		vmin=np.quantile(np.concatenate((S_schism[tag].flatten(),S[tag].flatten())),quantiles[0])
		vmax=np.quantile(np.concatenate((S_schism[tag].flatten(),S[tag].flatten())),quantiles[1])
		plt.pcolor(Time[tag],zmat,S[tag],vmin=vmin,vmax=vmax)
	else:
		plt.pcolor(Time[tag],zmat,S[tag])
		vmin, vmax = plt.gci().get_clim()
	plt.colorbar()
	plt.title(tag+ ' S [psu]')
	
	plt.tick_params(labelbottom=False)
	plt.subplot(2,2,2)
	plt.pcolor(Time[tag],zmat,S_schism[tag],vmin=vmin,vmax=vmax)
	plt.colorbar()
	plt.title('schism'+ ' S [psu]')
	plt.subplot(2,2,3)
	if quantile_range:
		vmax=np.quantile(np.abs(S[tag]-S_schism[tag]),quantiles[1])
		vmin=-vmax
		plt.pcolor(Time[tag],zmat,S[tag]-S_schism[tag],vmin=vmin,vmax=vmax,cmap='jet')
	else:
		plt.pcolor(Time[tag],zmat,S[tag]-S_schism[tag],cmap='jet')
	plt.colorbar()
	plt.title(tag + ' - schism')
	plt.gcf().autofmt_xdate()
	mng = plt.get_current_fig_manager()
	plt.tight_layout()
	plt.subplot(2,2,4)
	bias=(S_schism[tag]-S[tag]).mean(axis=1)
	rmse=np.sqrt(((S[tag]-S_schism[tag])**2).mean(axis=1))
	plt.bar(Zsalt[tag]-0.5,bias, label='bias')
	plt.bar(Zsalt[tag]+0.5, rmse,color='g', label='RMST')
	plt.title('total Bias/rmse {:.2f}/{:.2f}'.format((S[tag]-S_schism[tag]).mean(),np.sqrt( ((S[tag]-S_schism[tag])**2).mean()))) 

	for zi,biasi,rmsei in zip(Zsalt[tag],bias,rmse):
		plt.text(zi-0.6,biasi,'{:.2f}'.format(biasi),rotation=90 )
		plt.text(zi+0.6,rmsei,'{:.2f}'.format(rmsei),rotation=90 )
	plt.xlabel('depth')
	plt.legend(('bias','rmse'),frameon=False)
	plt.grid()
	plt.tight_layout()
	plt.savefig(outdir+'3_'+tag+'Salt_statistics.png',dpi=300)

	#bias=(S[tag]-S_schism[tag]).mean(axis=1)
	#rmse=np.sqrt(((S[tag]-S_schism[tag])**2).mean(axis=1))
	R=np.asarray([np.corrcoef(S[tag][j,:],S_schism[tag][j,:])[0,1] for j in range(S[tag].shape[0])])

	if len(S[tag][S[tag].mask==False])<1:
		bias*=np.nan
		rmse*=np.nan
		R*=np.nan
	
	
	compare[tag]['salt']['bias']=bias		
	compare[tag]['salt']['rmse']=rmse		
	compare[tag]['salt']['std1']=np.std(S[tag],axis=1)
	compare[tag]['salt']['std2']=np.std(S_schism[tag],axis=1)	
	compare[tag]['salt']['cor']=R




varname='salt'
import pandas as pd	
def write_table(varname,Zsalt,compare,names):
	d = {'station': names}	
	f=open(varname+'.txt','w')
	f.write('quant\stat')
	for tag in names:
		f.write('\t {:8s}'.format(tag))
	f.write('\n')	
	zs=np.unique(np.hstack([np.asarray(Zsalt[tag]) for tag in Zsalt.keys() ]))
	for key in compare[tag]['salt'].keys():
		for zi in zs:
			vals=[]
			f.write(key + ' {:8.2f}'.format(zi))
			for tag in names:
				if zi in Zsalt[tag]:
					val=np.asarray(compare[tag][varname][key])[Zsalt[tag]==zi][0]
					
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
	
	

dfsalt=write_table('salt',Zsalt,compare,names)
dftemp=write_table('temp',Ztemp,compare,names)

for key in compare[tag]['salt'].keys():
	for zi in zs:
		vals=[]
		for tag in names:
			if zi in Zsalt[tag]:
				val=np.asarray(compare[tag]['salt'][key])[Zsalt[tag]==zi][0]
			else:
				val=np.nan
			vals.append(val)
		d[key+str(np.int(zi))]=vals


set_fontsize(fac=1.8)		
zs=np.asarray(np.unique(np.hstack([np.asarray(Zsalt[tag]) for tag in Zsalt.keys() ])),int)
xbar=range(len(names))
plt.figure()
for zi in zs:
	phs=[]
	
	plt.clf()
	plt.pause(0.001)
	phs.append(plt.bar(xbar,dfsalt['bias'+str(np.int(zi))].values,width=-0.25,align='edge',label='bias',tick_label=names))
	phs.append(plt.bar(xbar,dfsalt['rmse'+str(np.int(zi))].values,width=0.25,align='edge',label='bias',tick_label=names))
	plt.grid()
	plt.legend(('bias','rmse'))
	plt.ylabel('salinity [-]')
	plt.tight_layout()
	plt.savefig('salt_error_bars_{:d}.png'.format((zi)),dpi=300)
	

zs=np.asarray(np.unique(np.hstack([np.asarray(Ztemp[tag]) for tag in Ztemp.keys() ])),int)
xbar=range(len(names))
for zi in zs:
	phs=[]
	plt.clf()
	plt.pause(0.001)
	phs.append(plt.bar(xbar,dftemp['bias'+str(np.int(zi))].values,width=-0.25,align='edge',label='bias',tick_label=names))
	phs.append(plt.bar(xbar,dftemp['rmse'+str(np.int(zi))].values,width=0.25,align='edge',label='bias',tick_label=names))
	plt.grid()
	plt.legend(('bias','rmse'))
	plt.ylabel('temperature [°C]')
	plt.tight_layout()
	plt.savefig('temp_error_bars_{:d}.png'.format((zi)),dpi=300)
	
set_fontsize(fac=1)		

# top botoopm
label=['top','bottom']
for nr,i in enumerate([0,-1]):
	samples=[ [compare[tag]['temp']['std2'][i]/compare[tag]['temp']['std1'][i],compare[tag]['temp']['cor'][i] ,tag+' '+str(Zsalt[tag][i])+'m'] for tag in names ]
	plotTaylor(samples=samples)
	plt.suptitle('Temperature ' +label[nr])
	plt.savefig(outdir+'5_taylor_temp'+label[nr]+'.png',dpi=300)
	plt.close()

# top botoopm
label=['top','bottom']
for nr,i in enumerate([0,-1]):
	samples=[ [compare[tag]['salt']['std2'][i]/compare[tag]['salt']['std1'][i],compare[tag]['salt']['cor'][i] ,tag+' '+str(Zsalt[tag][i])+'m'] for tag in names ]
	plotTaylor(samples=samples)
	plt.suptitle('Salinity ' +label[nr])
	plt.savefig(outdir+'5_taylor_salt'+label[nr]+'.png',dpi=300)	
	plt.close()	
	
# comparison with amm15	
if add_amm15: 
	 
	T15={name:0 for name in names}
	S15={name:0 for name in names}
	Time15={name:0 for name in names}
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
		YY=np.tile(Time[tag][23:23+nhours],[len(Zsalt[tag]),1]).swapaxes(0,1)
		
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
		#plt.plot(Time[tag][23:23+nhours],data[0,23:23+nhours])
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
		plt.pcolormesh(Time15[tag],-Z15[:iend],S15[tag][:iend,:],vmin=vmin,vmax=vmax)
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
		YY=np.tile(Time[tag][23:23+nhours],[len(Ztemp[tag]),1]).swapaxes(0,1)

		
		data=T[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([temp[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([temp[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		#gs = gridspec.GridSpec(10, 1)
		#plt.subplot(gs[0:2])
		#plt.plot(Time[tag][23:23+nhours],data[0,23:23+nhours])
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
		plt.pcolormesh(Time15[tag],-Z15[:iend],T15[tag][:iend,:],vmin=vmin,vmax=vmax)
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


# create latex + pdf
if put_pics_to_texdoc:
		print('generating tex doc and pdf')
		os.chdir(outdir)
		techit(latexname,latextitle,latextext)
		if os.path.isfile(latexname[:latexname.rindex('.')]+'.pdf'):
			 os.remove(latexname[:latexname.rindex('.')]+'.pdf')
		os.system('pdflatex '+latexname)
		print('done generating '+latexname[:latexname.rindex('.')]+'.pdf')
