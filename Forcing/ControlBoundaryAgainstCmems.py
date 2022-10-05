#!/gpfs/home/jacobb/anaconda3/bin/python
"""
 Control forcing genertated with forcing script against cmems amm15 nearest neighbours
 for amplitude and phase shifts
"""
import os
import sys
import matplotlib
#matplotlib.use('Agg') # backend
import datetime as dt
import glob
from scipy.spatial import cKDTree
import datetime as dt
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')
from schism import * # import schism functions
import xarray as xr
plt.ion()

############### SETTINGS ####################################
rundir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/'
ncdir=rundir+'combined/'
amm15dir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/'

forcingbds=[0]         # nrs all of boundaries used with thnc
openboundary=[0]   # nr of openboundary to check start counting from zero
each_nodes=100     # make time series comparison for <each_nodes> nodes along the boundary
offset_days=2      # start days after start
ndays=1.5            # consider days since model start
##########################################################


######### load SCHISM setup   ##################################
cwd=os.getcwd()
os.chdir(rundir)
s=schism_setup()
s.nntree = cKDTree(list(zip(s.lon,s.lat))) 
schismfiles=[] 
for iorder in range(6): # check for schout_nc files until 99999
	schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
nrs=list(np.asarray(nrs)[np.argsort(nrs)])
s.nc=xr.open_mfdataset(schismfiles)


#In [14]: thnc_time[0]
#Out[14]: datetime.datetime(2018, 1, 3, 0, 0)


### load forcing from boundary forcing file
fileE='elev2D.th.nc'
nc=xr.open_dataset(fileE)
thnc_time=nc['time'][:].values
keep=(offset_days*86400 <= thnc_time)  & (thnc_time <= (offset_days + ndays)*86400)
thnc_time=thnc_time[keep]
thnc_all=np.squeeze(nc['time_series'][keep,:,:,:].values)
thnc_by_frcbd=[]
nc.close()	
istart,iend=0,0
thncnodes=[]   # indices for thnc for each forcing boundary
frcbdnodes=[]  # indices wirth respect to entire schism grid
     # indics for indiviudal boundaries within thnc 
for nr,i in enumerate(forcingbds): #range(len(stp.bdy_segments)):
	iend+=len(s.bdy_segments[i])
	inds=np.asarray(range(istart,iend))
	thncnodes.append(inds)
	thnc_by_frcbd.append(thnc_all[:,inds])
	frcbdnodes.append(s.bdy_segments[i])
	istart=iend

### laod data from schism
s.nc=xr.open_mfdataset(schismfiles)
schism_time_all=s.nc['time'][:].values
modeldt=s.nc['time'][1]-s.nc['time'][0]
# convert datetime formats  considered period (somehow the time in different netcdf not compatible)
a=[s.nc['time'][0]-modeldt+np.timedelta64(dt.timedelta(seconds=thnc_time[i])) for i in range(len(thnc_time))]
thnc_time=np.asarray([dt.datetime.utcfromtimestamp((a[i] - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')) for i in range(len(a))])
ts = (s.nc['time'][0]-modeldt - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
t0=dt.datetime.utcfromtimestamp(ts)+dt.timedelta(days=offset_days)
ts = (s.nc['time'][0]+np.timedelta64(dt.timedelta(days=ndays)) - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
t1=dt.datetime.utcfromtimestamp(ts)+dt.timedelta(days=offset_days)

s.nc=s.nc.sel(time=slice(np.datetime64(t0),np.datetime64(t1)),nSCHISM_hgrid_node=np.concatenate(np.asarray(frcbdnodes)-1))
model_at_bd_all=s.nc['elev'].values
model_by_frcbd=[model_at_bd_all[:,nodes] for nodes in thncnodes]
tmodel=s.nc['time'].values

###### load amm15 setup at nn nodes #########
files=np.sort(glob.glob(amm15dir+'*SSH*nc'))
amm15_nc=xr.open_mfdataset(files)
amm15_nc=amm15_nc.sel(time=slice(t0,t1))
# nearest neighbours
allnodes=np.concatenate(frcbdnodes)-1
lon,lat = np.asarray(s.lon)[allnodes], np.asarray(s.lat)[allnodes]
coords=list(zip(lon,lat))
LON,LAT=np.meshgrid(amm15_nc['lon'],amm15_nc['lat'])
sshtest=amm15_nc['zos'][0,:].values
LON2=LON
LON2[np.isnan(sshtest)]=9999
nntree = cKDTree(list(zip(LON2.flatten(),LAT.flatten()))) 
lldists2,nn2=nntree.query(coords)		
ii,jj=np.unravel_index(nn2,LON.shape)
zeta_amm=np.asarray([amm15_nc['zos'][i,:].values[ii,jj] for i in range(len(amm15_nc['time']))])
time_amm=amm15_nc.indexes['time'].to_datetimeindex().values
plon=[LON2.flatten()[nn2][nodes] for nodes in thncnodes ]
plat=[LAT.flatten()[nn2][nodes] for nodes in thncnodes ]
zeta_amm_by_bdy=[zeta_amm[:,nodes] for nodes in thncnodes ]


##### plot ##################
plt.close('all')
plt.figure(dpi=300)
for bd in range(len(forcingbds)):
	consider_nodes=np.arange(0,len(thncnodes[bd]),each_nodes) 
	npoints=len(consider_nodes)
	m=np.int(np.sqrt(npoints))
	n=m+1
	plt.clf()
	plt.subplot(m,n,1)
	s.plot_domain_boundaries(append=True)
	plt.tight_layout()
	for i,node in enumerate(consider_nodes):
		plt.subplot(m,n,1)
		# thnc
		plt.plot(s.lon[frcbdnodes[bd][node]],s.lat[frcbdnodes[bd][node]],'ko')
		plt.plot(plon[bd][node],plat[bd][node],'b+')
		#plt.text(s.lon[frcbdnodes[bd][node]],s.lat[frcbdnodes[bd][node]],'bdn {:d}'.format(node))
		#plt.text(plon[bd][node],plat[bd][node],'nn amm {:d}'.format(node))
		plt.subplot(m,n,2+i)
		plt.title('bdn {:d}'.format(node))
		plt.plot(time_amm,zeta_amm_by_bdy[bd][:,consider_nodes[i]],'.-')
		plt.plot(thnc_time,thnc_by_frcbd[bd][:,consider_nodes[i]],'r.--')
		plt.plot(tmodel,model_by_frcbd[bd][:,consider_nodes[i]],'k.:')
		plt.grid()
		if i ==	len(consider_nodes)-1:
			plt.legend(['nn amm15','thnc','schism output'],frameon=False)
	plt.suptitle('bd {:d}'.format(bd))
	plt.gcf().autofmt_xdate()
	plt.savefig('bd{:d}_cmems_comp'.format(bd),dpi=300)