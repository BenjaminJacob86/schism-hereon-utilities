"""
Compare and interpolate model setups to each other
"""
#export OMP_NUM_TtimHREADS=4
import os
import netCDF4
import sys
import csv
import matplotlib
#matplotlib.use('Agg') # backend
import datetime as dt
import glob
from scipy.spatial import cKDTree
import datetime as dt
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
from matplotlib import path
# own and 3d party libraries
#sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
#sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
# strand
#sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
#sys.path.insert(0,'/gpfs/work/jacobb/data/shared/SCHISM/validation/lib')
#sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
#sys.path.insert(0,'/gpfs/work/jacobb/data/shared/SCHISM/validation/lib')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hereon-utilities/')#validation/lib')

#from techit import * # import latex script
#from TaylorDiagram import * # import taylordiagram
#from data_and_model_classes import Amm15,cmems
import xarray as xr
from schism import *
plt.ion()
# subplot surface slab of SCHISM and Amm15 and their difference 
# iterating over timesteps outputing images for video
# including time series at given location coords
		
		

########## settings ##################################################
# directories (have to end with '/')
# SCHIMS grid files
setupdir=['/work/gg0877/g260200/GB_sth_HydroSed_Feb2020_NWF/','/work/gg0877/g260200/GB_sth_HydroSed_Feb2020_EWF/','/work/gg0877/g260200/GB_sth_HydroSed_Feb2020_FIT/']

names=setup_names=['NWF','EWF','FIT']


oceandir=['/work/gg0028/g260114/SETUPS/NWBS/setup/Forcing/data/nrt.cmems-du.eu/Core/BLKSEA_ANALYSISFORECAST_PHY_007_001/all/']

#ncdir=[setupdir[0]+'/outputs_all/']
ncdir=[setupdir[i] + 'outputs_all/' for i in range(len(setup_names))]  


outdir='/work/gg0028/g260114/PROJECTS/OLAMUR/comps/' #setupdir[0]+'comps/'
if not os.path.exists(outdir): os.mkdir(outdir) 
#plt.ion()
#names=['CMEMS','NWBS']

#varnames=['temp','ssh','salt']						#varnames=['ssh','salt','temp'] ['ssh',] if only one has to have ,
varnames=['ssh',]
#varnames=['temp','salt']
min_max={'ssh':(-1.5,1.5),'salt':(0,35),'temp':(2,24)}	# axis range for variables
difflims={'ssh':1,'salt':4,'temp':4}                # # axis limits +- in difference plot
dthours={'ssh':3,'salt':12,'temp':12}				# make plot each dthours hour
ndays_ts={'ssh':2,'salt':30,'temp':30}             # nr of days depicted in time series subplot

# considrede time periods and steps for vriables
#vartimes={
#'ssh':{'startdate':dt.datetime(2020,1,1,1,0),'enddate':dt.datetime(2020,3,1,1,0),'step[hours]':1},
#'salt':{'startdate':dt.datetime(2020,1,1,12,0),'enddate':dt.datetime(2020,12,31,12,0),'step[hours]':24},
#'temp':{'startdate':dt.datetime(2020,1,1,12,0),'enddate':dt.datetime(2020,12,31,12,0),'step[hours]':24},
#}

vartimes={
'ssh':{'startdate':dt.datetime(2020,2,1,1,0),'enddate':dt.datetime(2020,3,1,1,0),'step[hours]':1},
'salt':{'startdate':dt.datetime(2020,2,1,12,0),'enddate':dt.datetime(2020,3,21,12,0),'step[hours]':24*3},
'temp':{'startdate':dt.datetime(2020,2,1,12,0),'enddate':dt.datetime(2020,3,21,12,0),'step[hours]':24*3},
}



# coords for time series
#lon=7.6548217184545075
#lat=53.75554953775793
# helgoland
#lat,lon = 54.181453837423874, 7.8974252215073095
limx=((5.11,10.41))	#((-1.14,9.84))
limy=((53.038,55.63))	#((49.7,56.21))


#lat,lon = 54.181453837423874, 7.8974252215073095
#OWF center
lonmlat=7.7135470958629035, 54.36756716567857
limx=((5.11,10.41))	#((-1.14,9.84))
limy=((53.038,55.63))	#((49.7,56.21))


dthour=1
ndays=365 # Run durations

# apperence
cm=plt.cm.jet # colormap
ax_limitEQ_commonDataRange=True    # minmize ax limits to common Extend which is not Land mask in any of the models
cm.set_over(color='m', alpha=None)
cm.set_under(color='k', alpha=None)
cm.set_bad(color='gray', alpha=None)
##############################################################################		

# open file acess schism
#def init_schism_file_access(ncdir):
#	schismfiles=glob.glob(ncdir+'*.nc')	
#	nrs=[]
#	for file in schismfiles:
#		nrs.append(int(file[file.rindex('_')+1:file.rindex('.')]))
#	schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
#	ds=xr.open_mfdataset(schismfiles[:3])
#	return ds
#dss={name:init_schism_file_access(ncdiri) for name,ncdiri in zip(names[1:],ncdir)}
#
# test


###### open file acess reference model potentially different files per vairables

#ls /gpfs/work/jacobb/data/SETUPS/GermanBight/GermanBight_HR_Ballje/gcoast2/201105*
#GCWAVE_COUPLED_02_1h_OUTPUT_20110531.nc
#GCWAVE_COUPLED_02_6h_OUTPUT_20110531.nc

# patterns 
#folders=np.sort(['/gpfs/work/jacobb/data/SETUPS/GermanBight/GermanBight_HR_Ballje/gcoast2/20101231'] + glob('/gpfs/work/jacobb/data/SETUPS/GermanBight/GermanBight_HR_Ballje/gcoast2/2011????') +
#['/gpfs/work/jacobb/data/SETUPS/GermanBight/GermanBight_HR_Ballje/gcoast2/20120101'])
#

## myocean
#ocfiles=np.sort(glob.glob(oceandir+'*.nc'))
#date=reftime #dt.datetime(2006,1,2,0,0,0)
#
#if 1: # amm 15
#	Tfile,Sfile,SSHfile=[oceandir + namei for namei in Amm15.gnameExtern(date)]
#	moc=Amm15(Tfile,Sfile,SSHfile)
#	#moc=cmems(Tfile,Sfile,SSHfile)
#	moc.gname(date)		
#elif False:	# myocean
#	ocfiles=np.sort(glob.glob(oceandir+'*.nc'))
#	date=reftime #dt.datetime(2006,1,2,0,0,0)
#	#Tfile,Sfile,SSHfile=[oceandir + namei for namei in Amm15.gnameExtern(date)]
#	#moc=Amm15(Tfile,Sfile,SSHfile)
#	#moc.gname(date)		
#	Tfile=Sfile=SSHfile=ocfiles[0] #manual
#	moc=cmems(Tfile)
#else :# gcoast
#	#datum folders
#	folders=np.sort(['/gpfs/work/jacobb/data/SETUPS/GermanBight/GermanBight_HR_Ballje/gcoast2/201105??'])
#	#['/gpfs/work/jacobb/data/SETUPS/GermanBight/GermanBight_HR_Ballje/gcoast2/20120101'])
#
#	# gather filelist
#	ssh_pattern='GCWAVE_COUPLED_02_1h_OUTPUT' # shh 
#	ssh_files=[glob.glob(folder+'/*'+ssh_pattern+'*') for folder in folders]
#	# temp / salt
#	temp_pattern='GCWAVE_COUPLED_02_6h_OUTPUT'
#	temp_files=[glob.glob(folder+'/*'+temp_pattern+'*') for folder in folders]
#	salt_pattern='GCWAVE_COUPLED_02_6h_OUTPUT'
#	salt_files=temp_files#[glob(folder+'/*'+salt_pattern+'*') for folder in folders]
#	#uv_pattern='GCWAVE_COUPLED_02_3h_OUTPUT' # uv
#	#uv_files=[glob(folder+'/*'+uv_pattern+'*') for folder in folders]
#	# establish file access
#	ds={'ssh':xr.open_mfdataset(ssh_files),
#	'salt':xr.open_mfdataset(salt_files)} #,'uv':xr.open_mfdataset(uv_files)
#	ds['temp']=ds['salt']
#

#############

def datetime64_to_datetime(t):
	if len(t)==1:
		t=[t]
	return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00Z'))/ np.timedelta64(1, 's')) for ti in t])
		

class cmems():

  def __init__(self,file):
    self.nc = netCDF4.Dataset(file)
    print(file)	
#    tv = tnc.variables
#    snc = netCDF4.Dataset('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/woa13_decav_s01_04v2.nc')
    sv = self.nc
#    latslice=slic
#    lonslice=slice(648,860)
    self.lon = sv['longitude'][:]#[lonslice]
    self.lat = sv['latitude'][:]#[latslice]
    self.d = -sv['depth'][:]
    self.time = sv['time'][:]
    self.timeunits = sv['time'].units
    self.tidx = 0
    self.s = sv['so'][:,:,:,:]#[:,:,latslice,lonslice]
    self.t = sv['thetao'][:,:,:,:]#[:,:,latslice,lonslice]
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    self.LON,self.LAT = np.meshgrid(self.lon,self.lat)
    self.u = sv['uo'][:,:,:,:]#[:,:,latslice,lonslice]
    self.v = sv['vo'][:,:,:,:]#[:,:,latslice,lonslice]
    self.ssh = sv['zos'][:,:,:]#[:,:,latslice,lonslice]
	
	
    self.ts=netCDF4.num2date(sv['time'][:],sv['time'].units)	
    self.t0=self.ts[0]
	
    #snc.close()
    self.tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.varnames={'salt':'so','temp':'thetao','ssh':'zos','u':'uo','v':'vo'}
	
  def update(self,file):
    self.nc.close() 
    self.nc  = netCDF4.Dataset(file)
    sv = self.nc.variables
    self.time = sv['time'][:]
    self.timeunits = sv['time'].units
    self.tidx = 0
    self.s = sv['so'][:,:,:,:]#[:,:,latslice,lonslice]
    self.t = sv['thetao'][:,:,:,:]#[:,:,latslice,lonslice]
    self.u = sv['uo'][:,:,:,:]#[:,:,latslice,lonslice]
    self.v = sv['vo'][:,:,:,:]#[:,:,latslice,lonslice]
    self.ssh = sv['zos'][:,:,:]#[:,:,latslice,lonslice]
	
    self.ts=netCDF4.num2date(sv['time'][:],sv['time'].units)	
    self.t0=self.ts[0]

    #snc.close()
	
  def get_slab(self,time,level=0,varname='salt'):		
    #nn=np.argsort(np.abs(self.ts-time))[:2] # interpolated linearly
    nn=[0]
    self.deltaT=np.abs(self.ts[nn]-time)
    dts=np.asarray([dti.total_seconds() for dti in self.deltaT])
    w=1/dts/(1/dts).sum()
    w[np.isnan(w)]=1.0
    self.t=time 
    if self.nc[self.varnames[varname]].shape==self.nc[self.varnames['salt']].shape:
        self.slab=self.nc[self.varnames[varname]][nn[0],level,:]#*w[0]+self.ncs[varname][self.varnames[varname]][nn[1],level,:]*w[1]
    else:
        self.slab=self.nc[self.varnames[varname]][nn[0],:]#*w[0]+self.ncs[varname][self.varnames[varname]][nn[1],:]*w[1]
	
  def get_hov(self,coords,varname):		
    """ Generate hovmoeller data for nextneighbour coordinates of coords """
    s.nn_self=self.tree.query(coords)[1]
    ii,jj=np.unravel_index(nn_self,self.LON.shape)
    if len(self.ncs[varname][self.varnames[varname]].shape)==4:
        self.profile[varname]=self.ncs[varname][self.varnames[varname]][:][::ii,jj]
	
class gcoast():

  def __init__(self,ds):
    self.nc = ds
#    tv = tnc.variables
    #snc = netCDF4.Dataset('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/woa13_decav_s01_04v2.nc')
    sv = self.nc
    self.lon = sv['salt']['lon'][:].values#[lonslice]
    self.lat = sv['salt']['lat'][:].values#[latslice]
    self.d = -sv['salt']['depth'][:].values
    self.time = sv['salt']['time'][:].values

    self.s = sv['salt']['vosaline'][0,:,:,:]#[:,:,latslice,lonslice]
    self.t = sv['temp']['votemper'][0,:,:,:]#[:,:,latslice,lonslice]
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    self.LON,self.LAT = np.meshgrid(self.lon,self.lat)
    if 'uv' in sv.keys():
        self.u = sv['uo'][:,:,:,:]#[:,:,latslice,lonslice]
        self.v = sv['vo'][:,:,:,:]#[:,:,latslice,lonslice]
    self.ssh = sv['ssh']['sossheig'][0,:,:]#[:,:,latslice,lonslice]
	
    #snc.close()
    self.tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.varnames={'salt':'vosaline','temp':'votemper','ssh':'sossheig','u':'uo','v':'vo'}
	
  def get_slab(self,time,level=0,varname='salt'):		
    self.varname=varname
    if varname=='ssh':
        self.slab=ds[varname].interp(time=np.datetime64(date))[self.varnames[varname]][:].values	
    else:
        self.slab=ds[varname].interp(time=np.datetime64(date))[self.varnames[varname]][level,:].values
    self.slab=np.ma.masked_invalid(self.slab)
    self.date=time	
  def get_hov(self,coords,varname):		
    """ Generate hovmoeller data for nextneighbour coordinates of coords """
    s.nn_self=self.tree.query(coords)[1]
    ii,jj=np.unravel_index(nn_self,self.LON.shape)
    if len(self.ncs[varname][self.varnames[varname]].shape)==4:
        self.profile[varname]=self.ncs[varname][self.varnames[varname]][:][::ii,jj]
  def plot_slab(self):		
    self.ph=plt.pcolormesh(self.lon,self.lat,self.slab)  
    self.ch=plt.colorbar()
    self.ch.set_label(self.varname)
    plt.title(str(self.date))		
	
  def plot_update(self):	
    slab=self.slab.flatten()
    slab=slab[~np.isnan(slab)][:-1]
    self.ph.set_array(slab)
    plt.gcf().canvas.draw()
    plt.gcf().canvas.flush_events() # flush gui events
	
cwd=os.getcwd()
######### load setups   #######################################
#def update(self,date):
#       try:
#            self.nc.close()
#        except:
#            pass
#        (date-self.reftime)
#        self.file=self.ncdir+'schout_'+str(int(np.ceil((date-self.reftime).total_seconds()/(self.dt*self.nt))))+'.nc'     
#        self.nc=Dataset(self.file)
#        self.date=self.reftime+dt.timedelta(seconds=self.nc['time'][0]-self.dt)
#        self.time=self.date+dt.timedelta(seconds=self.dt)*np.arange(1,self.nt+1) # here wase an error before        

# linear time interpolation
# get slab
def get_slab(self,time,varname,layer=-1):
		""" linear interpolate onto slab """
		# sel two closest time steps
		time=np.datetime64(time)
		nn=np.argsort(np.abs(self.nc.time-time))[:2].values
		deltaT=np.abs(np.asarray(self.nc.time[nn]-time,float)/1e9) #nano secinds
		w=1-deltaT/deltaT.sum()
		
		if len(self.nc[varname].shape)==2:
			self.slab=(w[0]*self.nc[varname][nn[0],:]+w[1]*self.nc[varname][nn[1],:]).values
		elif len(self.nc[varname].shape)==3:
			self.slab=(w[0]*self.nc[varname][nn[0],:,layer]+w[1]*self.nc[varname][nn[1],:,layer]).values
			#self.slab=np.ma.masked_array(self)
		#np.asarray([np.float(self.nc.time[nni]-time)/1e9 for nni in nn])
		self.wetdry=np.round(w[0]*self.nc['wetdry_elem'][nn[0],:]+w[1]*self.nc['wetdry_elem'][nn[1],:]).values
		self.ti=time
		
#
#def get_slab(self,date,varname,layer=-1):
#		#np.datetime64(time)  Out[811]: numpy.datetime64('2020-01-02T00:00:00.000000')
#		#model.nc['time'][23]
#		self.slab=self.nc[varname].interp(time=np.datetime64(date))[:,layer].values
#		self.wetdry=np.round(self.nc['wetdry_elem'].interp(time=np.datetime64(date)).values)
#		self.date=time	
#		self.ti=self.date

# schism setups
schism_setups=[]
for i_setup,diri in enumerate(setupdir):
	os.chdir(diri)
	exec('s{:d}=schism_setup()'.format(i_setup))
	exec('x=np.asarray(s{:d}.lon)'.format(i_setup))
	exec('y=np.asarray(s{:d}.lat)'.format(i_setup))

	# interpolatiion
	exec('s{:d}.nntree = cKDTree(list(zip(s{:d}.lon,s{:d}.lat))) '.format(i_setup,i_setup,i_setup))

	# initiate file acces
	schismfiles=glob.glob(ncdir[i_setup]+'*.nc')	
	nrs=[]
	for file in schismfiles:
		nrs.append(int(file[file.rindex('_')+1:file.rindex('.')]))
	schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
	#exec('s{:d}.nc=Dataset(schismfiles[0])'.format(i_setup))
	exec('s{:d}.nc=xr.open_mfdataset(schismfiles)'.format(i_setup))
	
	# set reference time
	#exec("reftime=dt.datetime.strptime(s{:d}.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')".format(i_setup))
	exec("Y,m,d,h,M=s{:d}.nc['time'].base_date.split()".format(i_setup))
	exec("reftime=dt.datetime(int(Y),int(m),int(d),int(float(h)),int(float(M)),0)".format(i_setup))
	
	#time=reftime+dt.timedelta(days=0.5)

	# enhance schism for conveinet load
	exec('s{:d}.ncdir=ncdir[i_setup]'.format(i_setup))		
	exec('s{:d}.reftime=reftime'.format(i_setup))
	exec('s{:d}.file="{:s}"'.format(i_setup,schismfiles[0]))
	#exec('s{:d}.nc=Dataset(s{:d}.file)'.format(i_setup,i_setup))
	exec('s{:d}.dates=datetime64_to_datetime(s{:d}.nc["time"])'.format(i_setup,i_setup))
	
	#exec('s{:d}.dt=np.diff(s{:d}.nc["time"][:2])[0]'.format(i_setup,i_setup))
	exec('s{:d}.dt=np.diff(s{:d}.dates[:2])[0].total_seconds()'.format(i_setup,i_setup))
	
	exec('s{:d}.nt=len(s{:d}.nc["time"])'.format(i_setup,i_setup))
	#exec('s{:d}.date=s{:d}.reftime+dt.timedelta(seconds=s{:d}.nc["time"][0]-s{:d}.dt)'.format(i_setup,i_setup,i_setup,i_setup))
	#exec('s{:d}.update=update'.format(i_setup))
	exec('s{:d}.get_slab=get_slab'.format(i_setup))
	exec('schism_setups.append(s{:d})'.format(i_setup))
##################################################################################################

	
################### Load schism ############################################
# myocean
#files=np.sort(glob.glob(oceandir+'*'+pattern+'*'+str(year)+'*'))
#files=np.sort(glob.glob(oceandir+'*.nc'))
#file=files[0]
#amm15=cmems(SSHfile=file)
#date=dt.datetime(2012,6,1,0,0,0)
date=s0.dates[0]	
if 1:
	Tfile,Sfile,SSHfile=[oceandir + namei for namei in Amm15.gnameExtern(date)]
	moc=Amm15(Tfile,Sfile,SSHfile)
	moc.gname(date)		
	moc.vardict={'ssh':'ssh','temp':'temp','salt':'salt'}
elif 0:
	# myocean
	ocfiles=np.sort(glob.glob(oceandir+'*.nc'))
	date=reftime #dt.datetime(2006,1,2,0,0,0)
	#Tfile,Sfile,SSHfile=[oceandir + namei for namei in Amm15.gnameExtern(date)]
	#moc=Amm15(Tfile,Sfile,SSHfile)
	#moc.gname(date)		
	Tfile=Sfile=SSHfile=ocfiles[0] #manual
	moc=cmems(Tfile)

	nc=Dataset(ocfiles[0])
	ocdates=netCDF4.num2date(nc['time'][0],nc['time'].units)+dt.timedelta(hours=24)*np.arange(len(ocfiles))
	nc.close()

if 0: #gcoast
	############################################################

	############# load data
	#date=reftime+dt.timedelta(seconds=3600)
	
	#file=ocfiles[np.argmin(np.abs(ocdates-time))]	
	#moc.update(date)
	#moc.update(file)
	moc=gcoast(ds)
	moc.vardict={'ssh':'ssh','temp':'temp','salt':'salt'}
	moc.get_slab(date,level=0,varname='ssh')
	#moc.plot_slab()


models=[moc]
struct=[1]
for s in schism_setups:
	s.vardict={'ssh':'elev','temp':'temp','salt':'salt'}
	s.get_slab(s,date,'elev')
	models+=[s]
	struct.append(0)

# cross interpolation
nmodels=len(models)
nx=nmodels
ny=nx
for j,model in enumerate(models):
	if 'schism' in str(model):
		bdnodes=[]
		x,y=np.asarray(model.lon),np.asarray(model.lat)
		for ocean, land in zip(model.bdy_segments,model.land_segments):
			bdnodes+=ocean #+land[1:]
			bdnodes+=land[1:]
		bdnodes=np.asarray(bdnodes)-1
		model.p=path.Path([(x[bdnodes][i],y[bdnodes][i]) for i in range(len(bdnodes))])
#oi=b.p		
#plt.plot(oi.vertices[:,0],oi.vertices[:,1])
		
####### nn between models ###############################################################
for model in models:
	model.isin={i:None for i in range(len(models))}
	model.isin2d={i:None for i in range(len(models))}
	model.masks={i:None for i in range(len(models))}
	model.nn={i:None for i in range(len(models))}
	model.delta={i:None for i in range(len(models))}

for j, model2 in enumerate(models):
	for i,model in enumerate(models):
		#if model2 != model:
		if j > i:
		#if j != i:
			if not 'schism' in str(type(model)):
				if not 'schism' in str(type(model2)):
					pass
				else: # model 1 struct model 2 unstruct
					model2.init_node_tree(latlon=True)
					coords=list(zip(model.LON.flatten(),model.LAT.flatten()))
					nn=model2.node_tree_latlon.query(coords)[1]# nn to all model1 nodes
					model2.nn[i]=nn.reshape(model.LON.shape)
					model2.masks[i]=~model2.p.contains_points(coords)  # 
			else:   # adapt to plot at model 1 below part not checked
				if not 'schism' in str(type(model2)):
					coords=list(zip(model2.LON.flatten(),model2.LAT.flatten())) # coords of model
					#
					#model.nn[j]=model.nntree.query(np.asarray(coords)[model.isin[j]])[1]
					#model.isin[j]=model.p.contains_points(coords) # check Amm15 lies within schism
					#model.isin2d[j]=np.reshape(model.isin[j],model2.LON.shape) #schism domain # area of model2
					#model.nn[j]=model.nntree.query(np.asarray(coords)[model.isin[j]])[1]
					
					
				else:
					#coords=list(zip(model2.lon,model2.lat))
					#model.isin[j]=model.p.contains_points(coords) # check Amm15 
					#model.nn[j]=model.nntree.query(np.asarray(coords)[model.isin[j]])[1]
					# to confusing with sub indices -> use nn and mask:
					coords=list(zip(model.lon,model.lat))
					model2.nn[i]=model2.nntree.query(np.asarray(coords))[1]
					model2.masks[i]=~model2.p.contains_points(coords)  # 
#############################################################################################


# test the nn trees
#a,b,c=models[0],models[1],models[2]
#
#plt.pcolormesh(a.LON,a.LAT,np.asarray(b.depths)[b.nn[0]])
#
#mplot=np.ma.masked_array(np.asarray(b.depths)[b.nn[0]],mask=b.masks[0])
#plt.clf()
#plt.pcolormesh(a.LON,a.LAT,mplot)
#
#mplot=np.ma.masked_array(np.asarray(c.depths)[c.nn[0]],mask=c.masks[0])
#plt.clf()
#plt.pcolormesh(a.LON,a.LAT,mplot)

#plt.clf()
#mplot=np.ma.masked_array(np.asarray(c.depths)[c.nn[1]],mask=c.masks[1])
#b.plotAtnodes(mplot,mask=c.masks[1][b.nvplt].max(axis=1))


#b.plotAtelems(mplot,mask=c.masks[1][b.nvplt].max(axis=1))
#b.plotAtelems(mplot,mask=c.masks[1])
#.max(axis=1)




# get cosest elements for dry check
for model in models[1:]:
	cx=np.asarray(model.x)[model.nvplt].mean(axis=1)
	cy=np.asarray(model.y)[model.nvplt].mean(axis=1)
	model.elemtree=cKDTree(list(zip(cx,cy)))
	model.nn_element=model.elemtree.query(coords)[1]
######## make plots ###########################################


# plot timeseries ######
moc.nntree=moc.tree  # adapt nn tree name as schisms
nns={name:model.nntree.query((lon,lat))[1] for name,model in zip(names,models)}


###########################################################
######### Variable and time iteration #####################
elems=np.asarray(list(s.nvdict))
for varname in varnames:

	#times=[s0.dates[0]+i*dt.timedelta(hours=dthour) for i in range(2*ndays*int(24/dthour))]
	ntimes=np.int((vartimes[varname]['enddate']-vartimes[varname]['startdate'])/dt.timedelta(hours=vartimes[varname]['step[hours]']))+1
	times=[vartimes[varname]['startdate']+ ti*dt.timedelta(hours=vartimes[varname]['step[hours]']) for ti in range(ntimes)]
	time=times[0]
	
	plt.close('all')
	vmin,vmax=min_max[varname]
	difflim=difflims[varname]
	print('ploting '+ varname)
	
	# data1 data2 run1 run2  run2-run1
	fig, axes = plt.subplots(nrows=nmodels,ncols=nmodels)
	fig.set_dpi(200)
	plt.tight_layout()
	#figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
	

	outdir2=outdir+varname+'/'
	if not os.path.exists(outdir2): os.mkdir(outdir2)
	
	# plot values
	phs=[]
	updata=[]
	for i,model in enumerate(models):
		i,model
		plt.subplot(ny,ny,i+1)
		plt.set_cmap('jet')
		if not 'schism' in str(type(model)):
			model.update(time)
			model.get_slab(time=time,varname=model.vardict[varname])
			exec('ph{:d}=plt.pcolormesh(model.LON,model.LAT,model.slab)'.format(i))
			#plt.colorbar()
			plt.colorbar(extend='both')
		else:
			model.get_slab(model,time=time,varname=model.vardict[varname])
			#exec('ph{:d},ch=s.plotAtelems(s.slab,cmap=cm,mask=None,extend="both")'.format(i))
			exec('ph{:d},ch=model.plotAtelems(model.slab,cmap=cm,mask=None)'.format(i))
		eval('phs.append(ph{:d})'.format(i))
		plt.title(names[i])
		plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
		plt.clim((vmin,vmax))
		if ax_limitEQ_commonDataRange:
			plt.xlim(limx)
			plt.ylim(limy)
		updata.append(model.slab)
		
	# plot differences
	count=len(models)-1
	pmodels=models.copy()  # models interpolated to and used for plotting
	for i,model in enumerate(models[:int(np.ceil(nx/2))]):
		for j, model2 in enumerate(models):
			#if model2 != model:
			if j > i:
				plt.subplot(ny,ny,(i+1)*nx+j+1)
				plt.set_cmap('jet')
				count+=1
				if not 'schism' in str(type(model)):
					if not 'schism' in str(type(model2)):
						#exec('ph{:d}=plt.pcolormesh(model.lon,model.lat,model.T[0,:,:])'.format(i))
						#plt.colorbar()
						pass
					else:
						mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
						exec('ph{:d}=plt.pcolormesh(model.LON,model.LAT,mplot)'.format(count))
						#eval('phs.append(ph{:d})'.format((i+1)*nx+j+1))
						eval('phs.append(ph{:d})'.format(count))
						plt.colorbar(extend='both')
						updata.append(mplot)
						struct.append(1)
						pmodels.append(model)
				elif 'schism' in str(type(model2)):
					mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
					#model.plotAtelems(mplot,mask=model.masks[1])
					exec('ph{:d},ch=model.plotAtelems(mplot,cmap=cm,mask=model2.masks[i])'.format(count))
					eval('phs.append(ph{:d})'.format(count))
					updata.append(mplot)
					struct.append(0)
					pmodels.append(model)

				plt.title(names[j]+' - '+names[i])
				plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
				plt.gca().tick_params(labelbottom=False,labelleft=False)
				plt.clim((-difflim,difflim))
				if ax_limitEQ_commonDataRange:
					plt.xlim(limx)
					plt.ylim(limy)

				
	# update labels         
	ttls=[[] for model in models]
	for i,model in enumerate(models):
		plt.subplot(ny,ny,i+1)
		if not 'schism' in str(type(model)):
			ttls[i]=plt.title(names[i]+ ' ' + str(model.t)[5:])
		else:
			ttls[i]=plt.title(names[i] + ' ' + str(model.ti)[10:])

	count=len(models)-1
	th1=[]
	for i,model in enumerate(models[:int(np.ceil(nx/2))]):
		for j, model2 in enumerate(models):
			#if model2 != model:
			if j > i:
				th1.append([])
				plt.subplot(ny,ny,(i+1)*nx+j+1)
				plt.set_cmap('jet')
				count+=1
				th1[count-len(models)]=plt.text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[count]),np.nanmean(np.abs(updata[count]))) )
	
	# time series data
	plt.subplot(nx,ny,1)
	plt.plot(lon,lat,'r+')
	 # update figures in loop
	 
	tdata={name:[] for name in names} 
	ydata={name:[] for name in names} 
	for i,model in enumerate(models):
		if not 'schism' in str(type(model)):
			tdata[names[i]].append(model.t)
			ydata[names[i]].append(model.slab.flatten()[nns[names[i]]])
		else:		
			tdata[names[i]].append(model.ti)
			ydata[names[i]].append(model.slab[nns[names[i]]])
			
	plt.subplot(nx,ny,ny+1)
	ax1 = plt.gca()
	plt.ylim((vmin,vmax))
	for i,name in enumerate(names):
		exec('tsph{:d}=plt.plot(np.asarray(tdata[name]),np.asarray(ydata[name]),linewidth=2)'.format(i))
	plt.legend(names,loc='lower center',ncol=2,frameon=False)
	plt.grid()
	plt.gcf().autofmt_xdate()
	plt.ylim((-difflim,difflim))
	
	#ax2 = ax1.twinx()
	#g3,=ax2.plot(np.asarray(tdata2),np.asarray(ydata2)-np.asarray(ydata),'k')
	#plt.ylim((-difflim,difflim))
	#plt.ylabel('difference')
	
	#lns = [g1,g2,g3]
	#ax.legend(lns, labs, loc=0)
	#plt.legend(lns,['Amm15','SCHISM','diff'],loc='lower center',ncol=2,frameon=False)
	#xlim=ax1.get_xlim()
	
	#################
	for ti,time in enumerate(times):
		t0=dt.datetime.now()
		print('doing plot'+str(ti)+': '+str(time))
		#update data
		for inr,model in enumerate(models):
			if struct[inr]:
				model.update(time)
				model.get_slab(time,0, model.vardict[varname])
			else:			
				model.get_slab(model,time,varname=model.vardict[varname],layer=-1)
				model.wetdry=np.asarray(model.wetdry,bool)
				model.slab[model.nvplt[model.wetdry[model.nvplt2nvp],:]]=np.nan # set dry elemts nodes tp nan
			updata[inr]=model.slab
		
		#update data - differences
		for i,model in enumerate(models[:int(np.ceil(nx/2))]):
			for j, model2 in enumerate(models):
				if model2 != model:
					if not 'schism' in str(type(model)):
						if not 'schism' in str(type(model2)):
							pass
						else:
							#mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
							inr+=1
							updata[inr]=(np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab)
					elif 'schism' in str(type(model2)):
						#mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
						inr+=1
						updata[inr]=(np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab)

		#plt.clf()					
		#ab=plt.pcolormesh(moc.LON,moc.LAT,updata[0])
		#ab.set_array(updata[3])
		# wet dry masking
		#t_nn=np.argmin(np.abs(s.time-time))
		#idry=s.nc['wetdry_elem'][t_nn,:][s.nvplt2nvp]==1
		#dryElem=(s.nc['wetdry_elem'][t_nn,:][s.nvplt2nvp][nn_element])==1
		#slab=s.slab[nn]  # mask if node not member of wet elment
		#slab[dryElem]=np.nan
		#moc.Mdiff[isin2d]=slab
		
		
		#update plot
		l=0
		
		
		for phi,data,structi,model in zip(phs,updata,struct,pmodels):
			l=+1
			if structi==1:
				phi.set_array(data[:-1,:-1].flatten())
			else:
				#phi.set_array(data[s.nvplt[:,:3]].mean(axis=1))
				phi.set_array(np.ma.masked_array(data[model.nvplt[:,:3]].mean(axis=1),mask=model.wetdry)) #ignore mask
				#phi.set_array(np.ma.masked_array(data[model.nvplt[:,:3]].mean(axis=1),mask=None))
				
		for ii,thi in enumerate(th1):
			thi.set_text('bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[nx+ii]),np.nanmean(np.abs(updata[nx+ii]))))
		
		# update labels         
		
		# update labels       
		alldata=np.concatenate([updata[0][~updata[0].mask].flatten()]+[updata[i].flatten() for i in range(1,len(models))]) 
		vmin,vmax=np.floor(np.nanquantile(alldata,0.1)),np.ceil(np.nanquantile(alldata,0.99))
		for i,model in enumerate(models):
			plt.subplot(ny,ny,i+1)
			if not 'schism' in str(type(model)):
				ttls[i]=plt.title(names[i]+ ' ' + str(model.t)[5:])
			else:
				ttls[i]=plt.title(names[i] + ' ' + str(model.ti)[10:])
			# update clims
			eval('ph{:d}.set_clim((vmin,vmax))'.format(i))
			
#		deltay=np.asarray(ydata2)-np.asarray(ydata)
#		g3.set_ydata(deltay)

		for i,model in enumerate(models):
			if not 'schism' in str(type(model)):
				tdata[names[i]].append(model.t)
				ydata[names[i]].append(model.slab.flatten()[nns[names[i]]])
			else:		
				tdata[names[i]].append(model.ti)
				ydata[names[i]].append(model.slab[nns[names[i]]])
		for i,name in enumerate(names):
			exec('tsph{:d}[0].set_xdata(np.asarray(tdata[name]))'.format(i))
			exec('tsph{:d}[0].set_ydata(np.asarray(ydata[name]))'.format(i))
			
		#plt.subplot(nx,ny,ny+1)
		#ax1 = plt.gca()
		ax1.set_xlim(time-dt.timedelta(days=ndays_ts[varname]),time+dt.timedelta(days=0.5))
		#ax2.set_xlim(time-dt.timedelta(days=ndays),time+dt.timedelta(days=0.5))
		ax1.set_ylim( np.floor(np.min([np.min(ydata[name]) for name in names])), np.ceil(np.max([np.max(ydata[name]) for name in names])))
		#ax1.set_ylim( np.floor(np.min((np.min(ydata),np.min(ydata2)))), np.ceil(np.max((np.max(ydata),np.max(ydata2)))))
		#dval=np.ceil(np.max(np.abs(deltay)))
		#ax2.set_ylim((-dval,dval))
		
		#ax1.set_ylim( np.floor(np.min((np.min(ydata),np.min(ydata2)))), np.ceil(np.max((np.max(ydata),np.max(ydata2)))))
		
		#fig.canvas.draw()
		plt.savefig(outdir2+'{:04d}_intercomp'.format(ti)+'_'+varname,dpi=400)
		print('took '+str((dt.datetime.now()-t0).total_seconds())+' s')

