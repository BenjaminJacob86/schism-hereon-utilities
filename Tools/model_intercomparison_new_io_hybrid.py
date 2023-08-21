"""
# subplot surface slab of SCHISM and Amm15 and their difference 
# iterating over timesteps outputing images for video
# including time series at given location coords
"""

#export OMP_NUM_TtimHREADS=4
import os, sys
import numpy as np
import netCDF4
import csv
interactive=True
import matplotlib
if interactive:
	from matplotlib import pyplot as plt
	plt.ion()
else:	
	matplotlib.use('Agg')
	#matplotlib.use('Cairo')
	
import datetime as dt
import xarray as xr
import glob
from scipy.spatial import cKDTree
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
from matplotlib import path

# @strand 
#sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
#sys.path.insert(0,'/gpfs/work/jacobb/data/shared/SCHISM/validation/lib')

# @levante
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/g260114/SETUPS/NWBS/setup/newcode/NWBS/schism-hzg-utilities/') 
#add schism python class to path
from schism import *

import warnings
warnings.filterwarnings("ignore")	

#sys.path.insert(0,'/work/gg0028/g260192/SEARECAP/validationScripts/')
#from validationFunctions import *



########## settings ##################################################
# directories (have to end with '/')
# SCHIMS grid files
setupdir=['/work/gg0028/g260114/SETUPS/NWBS2b',]
#setupdir=['/work/gg0028/g260114/SETUPS/NWBS/setup/newcode/NWBS_danube_repair/']
#setupdir+=['/work/gg0028/g260114/SETUPS/NWBS/setup/newcode/NWBS_danube_repair/NWBS_vgrid_re/Gauss_filter_grid/']


oceandir=['/work/gg0028/g260114/SETUPS/NWBS/setup/Forcing/data/nrt.cmems-du.eu/Core/BLKSEA_ANALYSISFORECAST_PHY_007_001/all/']

ncdir=[setupdir[0]+'/outputs_all/'] # where the outputs are.
#ncdir+=[setupdir[1]+'outputs/'] # where the outputs are.
#ncdir=[setupdir[0]+'outputs_hotstart_not_freshened/']

#outdir='/work/gg0028/g260192/SEARECAP/validationTest'

outdir=setupdir[0]+'imageComp/' #'/work/gg0028/g260114/postproc/modelcomp/comp4/'
if not os.path.exists(outdir): os.mkdir(outdir) 

names=['CMEMS','NWBS2b']#,'NWBS_vgrid_re','NWBS_GaussFilter',]

#varnames=['temp','ssh','salt']	#varname ['ssh',] if only one has to have ,
#varnames=['temp','salt']
varnames=['ssh',]
sdictvnames = {'temp':'temperature','ssh':'zCoordinates','salt':'salinity'}
min_max={'ssh':(-1,1),'salt':(0,25),'temp':(5,25)}	# axis range for variables
difflims={'ssh':1,'salt':4,'temp':4}     # # axis limits +- in difference plot
dthours={'ssh':3,'salt':12,'temp':12}	# make plot each dthours hour
ndays_ts={'ssh':2,'salt':30,'temp':30}             # nr of days depicted in time series subplot

# considrede time periods and steps for vriables
year0=2021
year1=2022
vartimes={'ssh':{'startdate':dt.datetime(year0,11,2,1,0),'enddate':dt.datetime(year1,11,2,1,0),'step[hours]':1},\
'salt':{'startdate':dt.datetime(year0,11,2,12,0),'enddate':dt.datetime(year1,2,1,12,0),'step[hours]':24},\
'temp':{'startdate':dt.datetime(year0,11,2,12,0),'enddate':dt.datetime(year1,2,1,12,0),'step[hours]':24},\
} # limit later by s0.dates


# coords for time series
# helgoland
lat,lon = 43, 29
limx=((27.3,31.8))	#((-1.14,9.84))
limy=((41,47))	#((49.7,56.21))

dthour=1
ndays=25 # Run durations


# apperence
cm=plt.cm.turbo # colormap
cmDiff=plt.cm.bwr
ax_limitEQ_commonDataRange=True    # minmize ax limits to common Extend which is not Land mask in any of the models
cm.set_over(color='m', alpha=None)
cm.set_under(color='k', alpha=None)
cm.set_bad(color='gray', alpha=None)


strdate=str(dt.datetime.now())[:19]
#f=open(''.join(('model_comp',strdate.replace(' ','_'),'.log')),'w')


#################################################################################################










############# Functions / Classes ###############################

def datetime64_to_datetime(t):
  if len(t)==1:
    t=[t]
  return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')) for ti in t])

def schism_time_to_datetime(t, baseTime):
  if len(t)==1:
    t=[t]
  return np.asarray([ dt.timedelta(seconds=ti) + baseTime for ti in t])		

class cmems():

	def __init__(self, files, varname):
		self.vardicts={'salt':'so','temp':'thetao','ssh':'zos','u':'uo','v':'vo'}
		varn = self.vardicts[varname]
		#self.nc = netCDF4.Dataset(file)
		sv = xr.open_mfdataset(files)
		ndims = len(sv.dims)
		self.lon = sv['lon'][:].data#[lonslice]
		self.lat = sv['lat'][:].data#[latslice]
		if ndims == 4:
			self.d = -sv['depth'][:].data
			exec("self.{:s} = sv['{:s}'][:,:,:,:]".format(varname, self.vardicts[varname]))
		elif ndims == 3: 
			exec("self.{:s} = sv['{:s}'][:,:,:]".format(varname, self.vardicts[varname]))
		self.ts = datetime64_to_datetime(sv['time'][:].data)
		self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
		self.LON,self.LAT = np.meshgrid(self.lon,self.lat)
		self.t0 = sv['time'][:].data.astype(str)[0:19]
		self.tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
		#exec("self.{:s}.ndims = ndims".format(varname))
		sv.close()
	
	def update_var(self, files, varname):
		varn = self.vardicts[varname]
		sv = xr.open_mfdataset(files)
		ndims = len(sv.dims)
		nt = sv['time'][:].shape[0]
		if self.ts.shape[0] != nt:
			print("updating cmems with a different time coordinate size")
		if ndims == 4:
			self.d = -sv['depth'][:].data
			exec("self.{:s} = sv['{:s}'][:,:,:,:].data".format(varname, self.vardicts[varname]))
		elif ndims == 3: 
			exec("self.{:s} = sv['{:s}'][:,:,:].data".format(varname, self.vardicts[varname]))
		#exec("self.{:s}.ndims = ndims".format(varname))
		sv.close()
	
	def get_slab(self, time, varname, ndims=4, layer=0):
		if (varname == 'ssh'):
			ndims=3
		nn=np.argsort(np.abs(self.ts-time))[:2]
		self.t = time
		#(1) in case of time == model time
		if (self.ts[nn[0]] - time).seconds == 0:
			if ndims == 3:
				exec("self.slab = self.{:s}[nn[0], :, :]".format(varname))
			elif ndims == 4:
				exec("self.slab = self.{:s}[nn[0], layer, :, :]".format(varname))
		#(2) interpolate in time
		else: 
			deltaT = np.asarray([ np.abs(ti - time).total_seconds() for ti in self.ts[nn]]) #seconds
			w = (1 - deltaT/deltaT.sum())	
			if ndims == 3:
				exec("self.slab = (w[0] * self.{:s}[nn[0], :, :] + w[1] * self.{:s}[nn[1], :, :])".format(varname,varname))
			elif ndims == 4:
				exec("self.slab = (w[0] * self.{:s}[nn[0], layer, :, :] + w[1] * self.{:s}[nn[1], layer, :, :])".format(varname,varname))
	
#  def get_hov(self,coords,varname):		
#    """ Generate hovmoeller data for nextneighbour coordinates of coords """
#    s.nn_self=self.tree.query(coords)[1]
#    ii,jj=np.unravel_index(nn_self,self.LON.shape)
#    if len(self.ncs[varname][self.varnames[varname]].shape)==4:
#        self.profile[varname]=self.ncs[varname][self.varnames[varname]][:][::ii,jj]
	

######## load setups   #######################################
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
def get_slab(self, time, varname, layer=-1):
		""" linear interpolate onto slab """
		# sel two closest time steps
		#time=np.datetime64(time)
		nn=np.argsort(np.abs(self.dates-time))[:2]
		
		#(1) in case of time == model time
		if (self.dates[nn[0]] - time).seconds == 0:
			self.slab = self.nc[varname][varname][nn[0], :, layer].values
			self.wetdry = self.nc['out2d']['dryFlagElement'][nn[0],:].values
			
		#(2) interpolate in time
		else: 
			#deltaT = np.asarray([ np.abs(ti - time).total_seconds() for ti in self.dates]) #seconds
			#w = (1 - deltaT/deltaT.sum())[nn]
			deltaT = np.asarray([ np.abs(ti - time).total_seconds() for ti in self.dates])[nn] #seconds
			w = (1 - deltaT/deltaT.sum())#

			
			self.slab= (w[0]*self.nc[varname][varname][nn[0], :, layer] \
			              + w[1]*self.nc[varname][varname][nn[1], :, layer]).values
			self.wetdry = np.round( w[0]*self.nc['out2d']['dryFlagElement'][nn[0],:] \
			              + w[1]*self.nc['out2d']['dryFlagElement'][nn[1],:] ).values
		self.ti = time
#####################################################################



		
########################### Program start #################################################################	

# schism setups
schism_setups=[]
for i_setup,diri in enumerate(setupdir):
	os.chdir(diri)
	exec("s{:d}=schism_setup(vgrid_file='vgrid.in.old')".format(i_setup)) #due to format change necessary for new i/o can only load the old vgrid.in
	exec('x=np.asarray(s{:d}.lon)'.format(i_setup))
	exec('y=np.asarray(s{:d}.lat)'.format(i_setup))

	# interpolatiion
	exec('s{:d}.nntree = cKDTree(list(zip(s{:d}.lon,s{:d}.lat))) '.format(i_setup,i_setup,i_setup))

	# initiate file access
	# by variables (new output format since 2022)
	vars=glob.glob('{}*_1.nc'.format(ncdir[i_setup])) # which variable files exist
	vars=[var[var.rindex('/')+1:var.rindex('_')] for var in vars]

	files=dict.fromkeys(vars)
	ds=dict.fromkeys(vars)
	for key in files.keys():
		files[key]=np.hstack([np.sort(glob.glob('./outputs/{:s}_{:s}.nc'.format(key,'?'*iorder))) for iorder in range(1,6)])
		files[key]=files[key][:-4]
		ds[key]=xr.open_mfdataset(files[key][:-1])

	#get the base time (t0) --> maybe this will change 
	p=param()
	reftime=dt.datetime(np.int(p.get_parameter('start_year')),\
	np.int(p.get_parameter('start_month')),\
	np.int(p.get_parameter('start_day')),\
	np.int(p.get_parameter('start_hour')),0,0)
		
	#time=reftime+dt.timedelta(days=0.5)

	# enhance schism for conveinet load
	exec('s{:d}.nc=ds'.format(i_setup))
	exec('s{:d}.ncdir=ncdir[i_setup]'.format(i_setup))		
	exec('s{:d}.reftime=reftime'.format(i_setup))
	exec('s{:d}.file_by_var=files'.format(i_setup))
	#exec('s{:d}.dates=datetime64_to_datetime(s{:d}.nc["time"])'.format(i_setup,i_setup))
	exec("s{:d}.dates=schism_time_to_datetime(ds['out2d'].time.data,reftime)".format(i_setup))
	exec("s{:d}.time=ds['out2d'].time.data".format(i_setup))
	exec("s{:d}.time_units='seconds since {}'".format(i_setup,reftime.strftime('%Y-%m-%d %H:%M:%S')))
	exec('s{:d}.dt=np.diff(s{:d}.dates[:2])[0].total_seconds()'.format(i_setup,i_setup))
	exec('s{:d}.nt=len(s{:d}.dates)'.format(i_setup,i_setup))

	exec('s{:d}.get_slab=get_slab'.format(i_setup)) #add the function to s0
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
if 0:
	Tfile,Sfile,SSHfile=[oceandir + namei for namei in Amm15.gnameExtern(date)]
	moc=Amm15(Tfile,Sfile,SSHfile)
	moc.gname(date)		
	moc.vardict={'ssh':'ssh','temp':'temp','salt':'salt'}
elif 1:
	# myocean
	Tfiles = np.sort(glob.glob(oceandir[0]+'*TEM*.nc'))
	Sfiles = np.sort(glob.glob(oceandir[0]+'*PSAL*.nc'))
	SSHfiles = np.sort(glob.glob(oceandir[0]+'*ASLV*.nc'))
	#date=reftime
	moc=cmems(Tfiles, 'temp') # this dask handle
	moc.update_var(Tfiles, 'temp')  #  if loaded they get nd array - quickfix so load again temp to have variables at same type
	moc.update_var(Sfiles, 'salt')  
	moc.update_var(SSHfiles, 'ssh')

	#nc=Dataset(ocfiles[0])
	#ocdates=netCDF4.num2date(nc['time'][0],nc['time'].units)+dt.timedelta(hours=24)*np.arange(len(ocfiles))
	#nc.close()

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
	s.vardict = sdictvnames
	s.get_slab(s, date, s.vardict['temp'])
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
#a,b=models[0],models[1]
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


#b.plotAtelems(mplot,mask=a.masks[0][b.nvplt].max(axis=1))
#b.plotAtelems(mplot,mask=a.masks[1])
#.max(axis=1)


#ok here ---

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

#time_array = s0.dates[5::12]
#time_array = s0.dates[5::48]


 # limit date range later by s0.dates
for key in vartimes.keys():
	vartimes[key]['startdate']=np.minimum(vartimes[key]['startdate'],s0.dates[0])
	vartimes[key]['enddate']=np.minimum(vartimes[key]['enddate'],s0.dates[-1])

# replot appears faster now	
#for varname in varnames:
#
#	ntimes=np.int((vartimes[varname]['enddate']-vartimes[varname]['startdate'])/dt.timedelta(hours=vartimes[varname]['step[hours]']))+1
#	times=[vartimes[varname]['startdate']+ ti*dt.timedelta(hours=vartimes[varname]['step[hours]']) for ti in range(ntimes)]
#	time=times[0]
#	
#	plt.close('all')
#	vmin,vmax=min_max[varname]
#	difflim=difflims[varname]
#	print('ploting '+ varname)
#	#f.write('ploting '+ varname)
#	
#	######### INITIAL PLOT FOR VARIABLE #######################################
#	# data1 data2 run1 run2  run2-run1
#	fig, axes = plt.subplots(nrows=nmodels,ncols=nmodels)
#	fig.set_dpi(200)
#	plt.suptitle( str(s0.dates[0]))
#	plt.tight_layout()
#	#figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
#	
#
#	outdir2=outdir+varname+'/'
#	if not os.path.exists(outdir2): os.mkdir(outdir2)
#	
#	# plot values
#	phs=[]
#	updata=[]
#	for i,model in enumerate(models):
#		i,model
#		plt.subplot(ny,ny,i+1)
#		plt.set_cmap('jet')
#		if not 'schism' in str(type(model)):
#			#model.update(time)
#			#model.get_slab(time=time,varname=model.vardict[varname])
#			model.get_slab(time,varname,layer=0)
#			exec('ph{:d}=plt.pcolormesh(model.LON,model.LAT,model.slab)'.format(i))
#			#plt.colorbar()
#			plt.colorbar(extend='both')
#		else:
#			model.get_slab(model,time=time,varname=model.vardict[varname])
#			exec('ph{:d},ch=model.plotAtelems(model.slab,cmap=cm,mask=None)'.format(i))
#		eval('phs.append(ph{:d})'.format(i))
#		plt.title(names[i])
#		plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
#		plt.clim((vmin,vmax))
#		if ax_limitEQ_commonDataRange:
#			plt.xlim(limx)
#			plt.ylim(limy)
#		updata.append(model.slab)
#		
#	# plot differences
#	count=len(models)-1
#	pmodels=models.copy()  # models interpolated to and used for plotting
#	for i,model in enumerate(models[:int(np.ceil(nx/2))]):
#		for j, model2 in enumerate(models):
#			#if model2 != model:
#			if j > i:
#				plt.subplot(ny,ny,(i+1)*nx+j+1)
#				plt.set_cmap('jet')
#				count+=1
#				if not 'schism' in str(type(model)):
#					if not 'schism' in str(type(model2)):
#						#exec('ph{:d}=plt.pcolormesh(model.lon,model.lat,model.T[0,:,:])'.format(i))
#						#plt.colorbar()
#						pass
#					else:
#						mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
#						exec('ph{:d}=plt.pcolormesh(model.LON,model.LAT,mplot)'.format(count))
#						#eval('phs.append(ph{:d})'.format((i+1)*nx+j+1))
#						eval('phs.append(ph{:d})'.format(count))
#						plt.colorbar(extend='both')
#						updata.append(mplot)
#						struct.append(1)
#						pmodels.append(model)
#				elif 'schism' in str(type(model2)):
#					mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
#					#model.plotAtelems(mplot,mask=model.masks[1])
#					exec('ph{:d},ch=model.plotAtelems(mplot,cmap=cm,mask=model2.masks[i])'.format(count))
#					eval('phs.append(ph{:d})'.format(count))
#					updata.append(mplot)
#					struct.append(0)
#					pmodels.append(model)
#
#				plt.title(names[j]+' - '+names[i])
#				plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
#				plt.gca().tick_params(labelbottom=False,labelleft=False)
#				plt.clim((-difflim,difflim))
#				if ax_limitEQ_commonDataRange:
#					plt.xlim(limx)
#					plt.ylim(limy)
#
#				# axes is currently retured array of arrays
#				# rowwise so index axes[i][j]
#				
#				
#	# update labels        
#	count=0	
#	ttls=[[] for model in models]
#	for i,model in enumerate(models):
#		#plt.subplot(ny,ny,i+1) overwrites plot now  seemingly due to change in matplolib
#		if not 'schism' in str(type(model)):
#			#ttls[i]=plt.title(names[i]+ ' ' + str(model.t)[5:])
#			ttls[i]=axes[0][i].set_title(names[i]+ ' ' + str(model.t)[5:])
#		else:
#			#ttls[i]=plt.title(names[i] + ' ' + str(model.ti)[10:])
#			ttls[i]=axes[0][i].set_title(names[i] + ' ' + str(model.ti)[10:])
#		count+=1
#
#	# display bias and mae error at timestep	
#	count=len(models)-1
#	th1=[]
#	for i,model in enumerate(models[:int(np.ceil(nx/2))]):
#		for j, model2 in enumerate(models):
#			#if model2 != model:
#			if j > i:
#				th1.append([])
#				#plt.subplot(ny,ny,(i+1)*nx+j+1)
#				#plt.set_cmap('jet')
#				
#				count+=1
#				#th1[count-len(models)]=plt.text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[count]),np.nanmean(np.abs(updata[count]))) )
#				th1[count-len(models)]=axes[i+1][j].text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[count]),np.nanmean(np.abs(updata[count]))) )
#	
#	# Plot time series data
#	#plt.subplot(nx,ny,1)
#	#plt.plot(lon,lat,'r+')
#	axes[0][0].plot(lon,lat,'r+')
#	 # update figures in loop
#	 
#	tdata={name:[] for name in names} 
#	ydata={name:[] for name in names} 
#	for i,model in enumerate(models):
#		if not 'schism' in str(type(model)):
#			tdata[names[i]].append(model.t)
#			ydata[names[i]].append(model.slab.flatten()[nns[names[i]]]) # add values
#		else:		
#			tdata[names[i]].append(model.ti)
#			ydata[names[i]].append(model.slab[nns[names[i]]])
#			
#	#plt.subplot(nx,ny,ny+1)
#	#ax1 = plt.gca()
#	#ax1 = plt.gca()
#	#plt.ylim((vmin,vmax))
#	for i,name in enumerate(names):
#		#exec('tsph{:d}=plt.plot(np.asarray(tdata[name]),np.asarray(ydata[name]),linewidth=2)'.format(i))
#		exec('tsph{:d}=axes[1][0].plot(np.asarray(tdata[name]),np.asarray(ydata[name]),".-",linewidth=2)'.format(i))
#	axes[1][0].legend(names,loc='lower center',ncol=2,frameon=False) # replaced plt with axes handle
#	axes[1][0].grid()
#	#axes[1][0].gcf().autofmt_xdate() this not yet working
#	axes[1][0].set_ylim((-difflim,difflim))
#	ax1=axes[1][0]
#
#	########################## ######### ######### ######### ######### ######### 
#	
#	# TIME LOOP only update date arays in plots for efficency ##################
#	for ti,time in enumerate(times):
#		t0=dt.datetime.now()
#		print('doing plot'+str(ti)+': '+str(time))
#		#f.write('doing plot \n'+str(ti)+': '+str(time)+'\n')
#		#f.write('doing plot {:d} \n'.format(ti))
#		#f.flush()
#		#update data
#		for inr,model in enumerate(models):
#			if struct[inr]:
#				#model.update(time)
#				model.get_slab(time,varname,layer=0)
#				#model.get_slab(time,0, model.vardict[varname])
#			else:			
#				model.get_slab(model,time,varname=model.vardict[varname],layer=-1)
#				model.wetdry=np.asarray(model.wetdry,bool)
#				model.slab[model.nvplt[model.wetdry[model.nvplt2nvp],:]]=np.nan # set dry elemts nodes tp nan
#			updata[inr]=model.slab
#		
#		#update data - differences
#		for i,model in enumerate(models[:int(np.ceil(nx/2))]):
#			for j, model2 in enumerate(models):
#				if model2 != model:
#					if not 'schism' in str(type(model)):
#						if not 'schism' in str(type(model2)):
#							pass
#						else:
#							
#							inr+=1
#							updata[inr]=(np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab)
#					elif 'schism' in str(type(model2)):
#						
#						inr+=1
#						updata[inr]=(np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab)
#
#		l=0
#		
#		
#		for phi,data,structi,model in zip(phs,updata,struct,pmodels):
#			l=+1
#			if structi==1:
#				#phi.set_array(data[:-1,:-1].values.flatten())
#				phi.set_array(data.flatten()) #.values
#			else:
#				#phi.set_array(data[s.nvplt[:,:3]].mean(axis=1))
#				phi.set_array(np.ma.masked_array(data[model.nvplt[:,:3]].mean(axis=1),mask=model.wetdry)) #ignore mask
#				#phi.set_array(np.ma.masked_array(data[model.nvplt[:,:3]].mean(axis=1),mask=None))
#				
#		for ii,thi in enumerate(th1):
#			thi.set_text('bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[nx+ii]),np.nanmean(np.abs(updata[nx+ii]))))
#		
#		# update labels         
#		
#		# update labels       
#		updata[0]=np.ma.masked_array(updata[0],mask=np.isnan(updata[0]))
#		alldata=np.concatenate([updata[0][~updata[0].mask].flatten()]+[updata[i].flatten() for i in range(1,len(models))]) 
#		vmin,vmax=np.floor(np.nanquantile(alldata,0.1)),np.ceil(np.nanquantile(alldata,0.99))
#		for i,model in enumerate(models):
#			#plt.subplot(ny,ny,i+1)
#			if not 'schism' in str(type(model)):
#				#ttls[i]=plt.title(names[i]+ ' ' + str(model.t)[5:])
#				ttls[i]=axes[0][i].set_title(names[i]+ ' ' + str(model.t)[5:])
#			else:
#				#ttls[i]=plt.title(names[i] + ' ' + str(model.ti)[10:])
#				ttls[i]=axes[0][i].set_title(names[i] + ' ' + str(model.ti)[10:])
#			# update clims
#			eval('ph{:d}.set_clim((vmin,vmax))'.format(i))
#			
##		deltay=np.asarray(ydata2)-np.asarray(ydata)
##		g3.set_ydata(deltay)
#
#		# this mayve bottleneck
#		for i,model in enumerate(models):
#			if not 'schism' in str(type(model)):
#				tdata[names[i]].append(model.t)
#				ydata[names[i]].append(model.slab.flatten()[nns[names[i]]])
#			else:		
#				tdata[names[i]].append(model.ti)
#				ydata[names[i]].append(model.slab[nns[names[i]]])
#		for i,name in enumerate(names):
#			exec('tsph{:d}[0].set_xdata(np.asarray(tdata[name]))'.format(i))
#			exec('tsph{:d}[0].set_ydata(np.asarray(ydata[name]))'.format(i))
#			
#		#plt.subplot(nx,ny,ny+1)
#		#ax1 = plt.gca()
#		ax1.set_xlim(time-dt.timedelta(days=ndays_ts[varname]),time+dt.timedelta(days=0.5))
#		#ax2.set_xlim(time-dt.timedelta(days=ndays),time+dt.timedelta(days=0.5))
#		ax1.set_ylim( np.floor(np.min([np.min(ydata[name]) for name in names])), np.ceil(np.max([np.max(ydata[name]) for name in names])))
#		#ax1.set_ylim( np.floor(np.min((np.min(ydata),np.min(ydata2)))), np.ceil(np.max((np.max(ydata),np.max(ydata2)))))
#		#dval=np.ceil(np.max(np.abs(deltay)))
#		#ax2.set_ylim((-dval,dval))
#		
#		#ax1.set_ylim( np.floor(np.min((np.min(ydata),np.min(ydata2)))), np.ceil(np.max((np.max(ydata),np.max(ydata2)))))
#		plt.suptitle(str(time))
#		fig.canvas.draw()
#		print('took '+str((dt.datetime.now()-t0).total_seconds())+' s')
#		plt.savefig(outdir2+'{:04d}_intercomp'.format(ti)+'_'+varname,dpi=400)
#		print('inc write took '+str((dt.datetime.now()-t0).total_seconds())+' s')
#		#f.write=('took '+str((dt.datetime.now()-t0).total_seconds())+' s\n')
#		#f.flush()
#
########### re plot #########

for varname in varnames:
	break
	ntimes=np.int((vartimes[varname]['enddate']-vartimes[varname]['startdate'])/dt.timedelta(hours=vartimes[varname]['step[hours]']))+1
	times=[vartimes[varname]['startdate']+ ti*dt.timedelta(hours=vartimes[varname]['step[hours]']) for ti in range(ntimes)]
	time=times[0]

	outdir2=outdir+varname+'/'
	if not os.path.exists(outdir2): os.mkdir(outdir2)

	vmin,vmax=min_max[varname]
	difflim=difflims[varname]
	plt.close('all')
	fig, axes = plt.subplots(nrows=nmodels,ncols=nmodels)
	fig.set_dpi(200)
	for ti,time in enumerate(times):
		#plt.close('all')
		print('ploting '+ varname)
		#f.write('ploting '+ varname)
		t0=dt.datetime.now()
		######### INITIAL PLOT FOR VARIABLE #######################################
		# data1 data2 run1 run2  run2-run1
		#fig, axes = plt.subplots(nrows=nmodels,ncols=nmodels)
		#fig.set_dpi(200)
		plt.suptitle( str(s0.dates[0]))
		plt.tight_layout()
		#figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
		
		# plot values
		phs=[]
		updata=[]
		for i,model in enumerate(models):
			i,model
			plt.subplot(ny,ny,i+1)
			plt.set_cmap('jet')
			if not 'schism' in str(type(model)):
				model.get_slab(time,varname,layer=0)
				exec('ph{:d}=plt.pcolormesh(model.LON,model.LAT,model.slab)'.format(i))
				plt.colorbar(extend='both')
				plt.title(names[i]+ ' ' + str(model.t)[5:])
			else:
				model.get_slab(model,time=time,varname=model.vardict[varname])
				exec('ph{:d},ch=model.plotAtelems(model.slab,cmap=cm,mask=None)'.format(i))
				plt.title(names[i] + ' ' + str(model.ti)[10:])
				
			eval('phs.append(ph{:d})'.format(i))
			#plt.title(names[i])
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
							#updata.append(mplot)
							#struct.append(1)
							pmodels.append(model)
					elif 'schism' in str(type(model2)):
						mplot=np.ma.masked_array(model2.slab[model2.nn[i]],mask=model2.masks[i])-model.slab
						#model.plotAtelems(mplot,mask=model.masks[1])
						exec('ph{:d},ch=model.plotAtelems(mplot,cmap=cm,mask=model2.masks[i])'.format(count))
						eval('phs.append(ph{:d})'.format(count))
						#updata.append(mplot)
						#struct.append(0)
						pmodels.append(model)
					#axes[i+1][j].text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[count]),np.nanmean(np.abs(updata[count]))) )
					axes[i+1][j].text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(mplot),np.nanmean(np.abs(mplot))) )
					plt.title(names[j]+' - '+names[i])
					plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
					plt.gca().tick_params(labelbottom=False,labelleft=False)
					plt.clim((-difflim,difflim))
					if ax_limitEQ_commonDataRange:
						plt.xlim(limx)
						plt.ylim(limy)
		
					# axes is currently retured array of arrays
					# rowwise so index axes[i][j]
					
					
		## update labels        
		#count=0	
		#ttls=[[] for model in models]
		#for i,model in enumerate(models):
		#	#plt.subplot(ny,ny,i+1) overwrites plot now  seemingly due to change in matplolib
		#	if not 'schism' in str(type(model)):
		#		#ttls[i]=plt.title(names[i]+ ' ' + str(model.t)[5:])
		#		ttls[i]=axes[0][i].set_title(names[i]+ ' ' + str(model.t)[5:])
		#	else:
		#		#ttls[i]=plt.title(names[i] + ' ' + str(model.ti)[10:])
		#		ttls[i]=axes[0][i].set_title(names[i] + ' ' + str(model.ti)[10:])
		#	count+=1
		
		# display bias and mae error at timestep	
		#count=len(models)-1
		#th1=[]
		#for i,model in enumerate(models[:int(np.ceil(nx/2))]):
		#	for j, model2 in enumerate(models):
		#		#if model2 != model:
		#		if j > i:
		#			th1.append([])
		#			#plt.subplot(ny,ny,(i+1)*nx+j+1)
		#			#plt.set_cmap('jet')
		#			
		#			count+=1
		#			#th1[count-len(models)]=plt.text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[count]),np.nanmean(np.abs(updata[count]))) )
		#			th1[count-len(models)]=axes[i+1][j].text( limx[0],limy[0]-0.1,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[count]),np.nanmean(np.abs(updata[count]))) )
		
		# Plot time series data
		#plt.subplot(nx,ny,1)
		#plt.plot(lon,lat,'r+')
		axes[0][0].plot(lon,lat,'r+')
		# update figures in loop
		
		tdata={name:[] for name in names} 
		ydata={name:[] for name in names} 
		for i,model in enumerate(models):
			if not 'schism' in str(type(model)):
				tdata[names[i]].append(model.t)
				ydata[names[i]].append(model.slab.flatten()[nns[names[i]]]) # add values
			else:		
				tdata[names[i]].append(model.ti)
				ydata[names[i]].append(model.slab[nns[names[i]]])
				
		#plt.subplot(nx,ny,ny+1)
		#ax1 = plt.gca()
		#ax1 = plt.gca()
		#plt.ylim((vmin,vmax))
		for i,name in enumerate(names):
			#exec('tsph{:d}=plt.plot(np.asarray(tdata[name]),np.asarray(ydata[name]),linewidth=2)'.format(i))
			exec('tsph{:d}=axes[1][0].plot(np.asarray(tdata[name]),np.asarray(ydata[name]),".-",linewidth=2)'.format(i))
		axes[1][0].legend(names,loc='lower center',ncol=2,frameon=False) # replaced plt with axes handle
		axes[1][0].grid()
		#axes[1][0].gcf().autofmt_xdate() this not yet working
		axes[1][0].set_ylim((-difflim,difflim))
		ax1=axes[1][0]
		
		########################## ######### ######### ######### ######### ######### 
		#plt.suptitle(str(time))
		#fig.canvas.draw()
		print('took '+str((dt.datetime.now()-t0).total_seconds())+' s')
		plt.savefig(outdir2+'{:04d}_intercomp'.format(ti)+'_'+varname,dpi=400)
		for i,model in enumerate(models[:int(np.ceil(nx/2))]):
			try:
				axes[0][i].clear()	
			except:
				pass
			for j, model2 in enumerate(models):
				try:
					axes[i+1][j].clear()
				except:
					pass
		print('inc write took '+str((dt.datetime.now()-t0).total_seconds())+' s')