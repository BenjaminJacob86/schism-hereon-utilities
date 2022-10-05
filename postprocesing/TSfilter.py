import os
import netCDF4
import sys
#sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities/')
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
from scipy import signal
from scipy.signal import butter, lfilter
########## settings #################################
# export OMP_NUM_THREADS=4 # call before python


#rundir='/gpfs/work/jacobb/data/RUNS/BlackSea/2017/'#'/gpfs/work/jacobb/data/SETUPS/BlackSea/'
rundir='/work/gg0028/g260114/RUNS/Europe/c1/flow_tweak/'
ncdir=rundir+'combined/'#'/gpfs/work/jacobb/data/SETUPS/BlackSea/combined/'



os.chdir(rundir)
s=schism_setup()

## check for schout_nc files
schismfiles=[]
for iorder in range(6):
    iorder
    schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
nrs=[]  
for file in schismfiles:
    nr=int(file[file.rfind('_')+1:file.index('.nc')])
    nrs.append(nr)
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
nrs=list(np.asarray(nrs)[np.argsort(nrs)])

# initialize file acces
ds=xr.open_mfdataset(schismfiles)

#get reftime
t=ds['time'][:]
basedate=t[0].base_date.split(' ')
while '' in basedate:
	basedate.remove('')
refdate=[]	
for item in basedate:
		item=np.int(np.float(item))
		refdate.append(item)
reftime=dt.datetime(refdate[0],refdate[1],refdate[2],refdate[3],refdate[4])
dates=reftime+np.arange(len(t))*dt.timedelta(seconds= np.float((t[1]-t[0]))/1e9)

def calc_stats(ds,varnames=['hvel','wind_speed','zcor','salt','temp'],lvl=-1,periods=['year']):
	""" Compute variable statistics min,max,mean, and variance grouped by temporal periods for variables indicated by list: varname. For variables with depth dimension staistics are calculate for sigma level lvl (-1:=surface)
	"""
	
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
			stats[varname]['mean']=np.asarray(means[varname][:,:,lvl])
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
	
	return stats

# schism 	
stats=calc_stats(ds,varnames=['elev'])

# duacs
#duacs=xr.open_dataset('/gpfs/work/jacobb/data/DATA/Download/blacksea/dataset-duacs-rep-blacksea-merged-allsat-phy-l4_1580991903572.nc')
duacs=xr.open_dataset('duacs2016to2018.nc')

duacs_yearly=duacs.groupby('time.'+'year')
duacs_means=duacs_yearly.mean('time')
duacs_vars=duacs_yearly.var('time')

def datetime64_to_datetime(t):
	if len(t)==1:
		t=[t]
	return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00Z'))/ np.timedelta64(1, 's')) for ti in t])
dates_duacs=datetime64_to_datetime(duacs['time'])



# anomalies
nt=len(duacs['time'])
#duacs_var=((duacs['sla']**2).sum(axis=0)/nt).values
lon=duacs['longitude'][:]
lat=duacs['latitude'][:]
cmap=plt.cm.jet
cmap.set_bad('w',0)
duacs_var=np.ma.masked_array(duacs_var,mask=duacs_var==0.0)

lon2d,lat2d=np.meshgrid(lon,lat)

s.init_node_tree() 

nn=s.node_tree_latlon.query(list(zip(lon2d.flatten(),lat2d.flatten())))[1]
nn2d=np.reshape(nn,lon2d.shape)

year=2012
for year in range(2012,2018):
	print(year)
	plt.clf()
	iyear=np.where(duacs_vars['sla'].year==year)[0][0]
	plt.subplot(2,2,1)
	plt.pcolormesh(lon,lat,np.sqrt(duacs_vars['sla'][iyear,:]),cmap=cmap)
	xlim=plt.xlim()
	ylim=plt.ylim()
	ch=plt.colorbar()
	ch.set_label('std SLA')
	plt.clim((0,0.25))
	vmin, vmax = plt.gci().get_clim()
	plt.title('duacs')
	plt.subplot(2,2,2)
	iyear=np.where(stats['years']==year)[0][0]
	s.plotAtnodes(np.sqrt(stats['elev']['var'][iyear]))
	#ch=plt.colorbar()
	ch.set_label('std SLA')
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.clim((vmin,vmax))
	plt.title('schism')
	
	slainterp=np.ma.masked_array(np.sqrt(stats['elev']['var'][iyear])[nn2d],mask=np.isnan(duacs_vars['sla'][iyear,:]))
	bias=slainterp-np.sqrt(duacs_vars['sla'][iyear,:])
	plt.subplot(2,2,3)
	iyear=np.where(stats['years']==year)[0][0]
	#plt.pcolormesh(lon,lat,slainterp,cmap=cmap)
	plt.pcolormesh(lon,lat,np.sqrt(duacs_vars['sla'][iyear,:])+bias.mean(),cmap=cmap)
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.clim((vmin,vmax))
	plt.colorbar()
	plt.title('ducas + mean bias')
	
	plt.subplot(2,2,4)
	plt.pcolormesh(lon,lat,bias,cmap=cmap)
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.suptitle(year)
	ch=plt.colorbar()
	ch.set_label('Delta std SLA')
	plt.tight_layout()
	plt.title('SCHISM - duacs')
	plt.savefig('schism_duacs_SLA_{:d}.png'.format(year))


year=2012
for year in range(2012,2018):
	print(year)
	plt.clf()
	iyear=np.where(duacs_means['sla'].year==year)[0][0]
	plt.subplot(2,2,1)
	plt.pcolormesh(lon,lat,(duacs_means['sla'][iyear,:]),cmap=cmap)
	xlim=plt.xlim()
	ylim=plt.ylim()
	ch=plt.colorbar()
	ch.set_label('mean SLA')
	#plt.clim((0,0.05))
	vmin, vmax = plt.gci().get_clim()
	plt.title('duacs')
	plt.subplot(2,2,2)
	iyear=np.where(stats['years']==year)[0][0]
	s.plotAtnodes((stats['elev']['mean'][iyear]))
	#ch=plt.colorbar()
	ch.set_label('mean SLA')
	plt.xlim(xlim)
	plt.ylim(ylim)
	#plt.clim((vmin,vmax))
	plt.title('schism')
	
	slainterp=np.ma.masked_array((stats['elev']['mean'][iyear])[nn2d],mask=np.isnan(duacs_means['sla'][iyear,:]))
	bias=slainterp-(duacs_means['sla'][iyear,:])
	plt.subplot(2,2,3)
	iyear=np.where(stats['years']==year)[0][0]
	#plt.pcolormesh(lon,lat,slainterp,cmap=cmap)
	plt.pcolormesh(lon,lat,(duacs_means['sla'][iyear,:])+bias.mean(),cmap=cmap)
	plt.xlim(xlim)
	plt.ylim(ylim)
	#plt.clim((vmin,vmax))
	plt.colorbar()
	plt.title('ducas + mean bias')
	
	plt.subplot(2,2,4)
	plt.pcolormesh(lon,lat,bias,cmap=cmap)
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.suptitle(year)
	ch=plt.colorbar()
	ch.set_label('Delta std SLA')
	plt.tight_layout()
	plt.title('SCHISM - duacs')
	plt.savefig('schism_duacs_mean_SLA_{:d}.png'.format(year))
	

# anomaly at instance	
for day in range(0,365*5,10):
	date=dt.datetime(2012,6,1)	+dt.timedelta(days=day)
	ischism=np.argmin(np.abs(dates-date))
	iduacs=np.argmin(np.abs(dates_duacs-date))
	plt.clf()
	plt.subplot(2,2,1)
	plt.pcolormesh(lon,lat,np.sqrt(duacs['sla'][iduacs,:]),cmap=cmap)
	xlim=plt.xlim()
	ylim=plt.ylim()
	ch=plt.colorbar()
	ch.set_label('SLA')
	plt.clim((-0.5,0.5))
	vmin, vmax = plt.gci().get_clim()
	plt.title('duacs')
	plt.subplot(2,2,2)
	s.plotAtnodes(ds['elev'][ischism,:]-ds['elev'][ischism,:].mean())
	#ch=plt.colorbar()
	ch.set_label('std SLA')
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.clim((vmin,vmax))
	plt.title('schism')
	plt.savefig('{:03d}_spat_anomaly'.format(day))
	
##

cmems_files=list(np.sort(glob.glob('/gpfs/work/jacobb/data/DATA/Download/blacksea/nrt.cmems-du.eu/Core/BLKSEA_ANALYSIS_FORECAST_PHYS_007_001/bs-cmcc-ssh-an-fc-d/2017/all/*.nc')))

cmems_files=list(np.sort(glob.glob('/work/gg0028/g260114/Data/cmems/my.cmems-du.eu/Core/BLKSEA_REANALYSIS_PHYS_007_004/sv04-bs-cmcc-ssh-rean-d/all/*.nc')))

cmems=xr.open_mfdataset(cmems_files)
periods=['year']

yearly=cmems.groupby('time.year')
means=yearly.mean('time')
#duacs_yearly=duacs.groupby('time.'+'year')
#duacs_means=duacs_yearly.mean('time')

means[varname][1]
duacs_means

vars=yearly.var('time')

varname=['sossheig']
cmems_stats=dict.fromkeys(varname)
quantity=['min','max','mean','var']
for key in cmems_stats.keys():
	cmems_stats[key]=dict.fromkeys(quantity)
cmems_stats[varname[0]]['mean']=means[varname[0]].values
cmems_stats[varname[0]]['var']=vars[varname[0]].values




np.asarray(vars[varname][:,:,lvl])
cmems_stats[varname]['var']=np.asarray(vars[varname])

#mins=yearly.min('time')



nt2=len(cmems['time'])
#cmems_var=(cmems['zos']-cmems['zos'].mean(axis=0))**2/nt2
cmems_var=(cmems['zos'].var(axis=0)).values
cmems_mean=(cmems['zos'].mean(axis=0)).values
lon2=cmems['lon']
lat2=cmems['lat']


# next neighbour interp
tree0=cKDTree(list(zip(s.lon,s.lat)))
LON2,LAT2=np.meshgrid(lon2,lat2)
tree=cKDTree(list(zip(LON2.flatten(),LAT2.flatten())))
dists0,coords0=tree0.query(list(zip(LON2.flatten(),LAT2.flatten())))
dists,coords=tree.query(list(zip(np.asarray(s.lon)[coords0],np.asarray(s.lat)[coords0])))



for year in range(2013,2017):
	plt.clf()
	plt.subplot(2,2,1)
	iyear0=np.where(means.year.values==year)[0]
	plt.pcolor(lon2,lat2,cmems_stats[varname[0]]['mean'][iyear0,:][0,:],cmap=cmap)
	ch=plt.colorbar()
	plt.clim((-0.4,0.1))
	vmin, vmax = plt.gci().get_clim()
	plt.title('cmems')
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch.set_label('<ssh [m]>')
	plt.subplot(2,2,2)
	#plt.pcolor(lon2,lat2,cmems_mean,cmap=cmap)
	iyear=np.where(stats['years']==year)[0][0]
	ph,ch3=s.plotAtnodes(stats['elev']['mean'][iyear,:])
	plt.title('schism')
	plt.suptitle('mean of '+str(year))
	plt.clim((vmin,vmax))
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch3.set_label('<ssh [m]>')
	schismNN=np.ma.masked_array(np.reshape(stats['elev']['mean'][iyear][coords0][coords],LON2.shape),mask=np.isnan(cmems_stats[varname[0]]['mean'][0,:]))
	bias=schismNN-cmems_stats[varname[0]]['mean'][iyear0,:][0,:]
	plt.subplot(2,2,3)
	iyear0=np.where(means.year.values==year)[0]
	plt.pcolor(lon2,lat2,cmems_stats[varname[0]]['mean'][iyear0,:][0,:]+bias.mean(),cmap=cmap)
	ch=plt.colorbar()
	plt.clim((-0.5,0.1))
	vmin, vmax = plt.gci().get_clim()
	plt.title('cmems + <bias>')
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch.set_label('<ssh [m]>')
	plt.subplot(2,2,4)
	plt.pcolor(lon2,lat2,bias,cmap=cmap)
	plt.title('schism - cmems')
	ch2=plt.colorbar()
	#plt.clim((-0.15,0.15))
	plt.clim((-0.5,0))
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch3.set_label('ssh [m]')
	plt.clim()
	plt.tight_layout()
	plt.savefig('schism_vs_cmems_{:d}.png'.format(year),dpi=300)

	
	
	
for year in range(2013,2017):
	plt.clf()
	plt.subplot(2,2,1)
	iyear0=np.where(vars.year.values==year)[0]
	plt.pcolor(lon2,lat2,np.sqrt(cmems_stats[varname[0]]['var'][iyear0,:][0,:]),cmap=cmap)
	ch=plt.colorbar()
	plt.clim((0,0.25))
	vmin, vmax = plt.gci().get_clim()
	plt.title('cmems')
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch.set_label('std ssh [m]')
	plt.subplot(2,2,2)
	#plt.pcolor(lon2,lat2,cmems_var,cmap=cmap)
	iyear=np.where(stats['years']==year)[0][0]
	ph,ch3=s.plotAtnodes(np.sqrt(stats['elev']['var'][iyear,:]))
	plt.title('schism')
	plt.suptitle('var of '+str(year))
	plt.clim((vmin,vmax))
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch3.set_label('var ssh [m]')
	schismNN=np.ma.masked_array(np.reshape(stats['elev']['var'][iyear][coords0][coords],LON2.shape),mask=np.isnan(cmems_stats[varname[0]]['var'][0,:]))
	bias=np.sqrt(schismNN)-np.sqrt(cmems_stats[varname[0]]['var'][iyear0,:][0,:])
	plt.subplot(2,2,3)
	iyear0=np.where(vars.year.values==year)[0]
	plt.pcolor(lon2,lat2,np.sqrt(cmems_stats[varname[0]]['var'][iyear0,:][0,:])+bias.mean(),cmap=cmap)
	ch=plt.colorbar()
	plt.clim((0,0.25))
	vmin, vmax = plt.gci().get_clim()
	plt.title('cmems + <bias>')
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch.set_label('std ssh [m]')
	plt.subplot(2,2,4)
	plt.pcolor(lon2,lat2,bias,cmap=cmap)
	plt.title('schism - cmems')
	ch2=plt.colorbar()
	#plt.clim((-0.15,0.15))
	plt.clim((-0.25,0.25))
	plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	ch3.set_label('std ssh [m]')
	plt.clim()
	plt.tight_layout()
	plt.savefig('schism_vs_cmems_std_{:d}.png'.format(year),dpi=300)	
	
	
	
	
xp=34.5
yp=43.2	
plt.clf()
ph,ch3=s.plotAtnodes(stats['elev']['mean'][iyear,:])
plt.plot(xp,yp,'ko')	
plt.clim((-0.5,0.2))
plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
	



	
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
butter_bandpass(lowcut, highcut, fs, order=5)
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    #y = lfilter(b, a, data) # not necissarily phase zerp
    y = signal.filtfilt(b, a, data)
    return y


def butter_pass(cut, fs, order=5,type='high'):
	nyq = 0.5 * fs
	cut = cut / nyq
	b, a = butter(order, cut, btype=type)
	return b, a
def butter_pass_filter(data, cut, fs, order=5,type='high'):
    b, a = butter_pass(cut, fs, order=order,type=type)
    #y = lfilter(b, a, data) # not necissarily phase zerp
    y = signal.filtfilt(b, a, data)
    return y

	
	

dUfilt=dUfilt.transpose()
äax.plot(tdata[indx],dUfilt[indx],'b',linewidth=1)

plt.plot(dates,dUfilt)
plt.plot(dates,dUfilt+dUfilt_low,'--')


#import xrscipy.signal as dsp  :-(


	

nnschism=tree0.query((xp,yp))[1] #tree0=cKDTree(list(zip(s.lon,s.lat)))
#tree=cKDTree(list(zip(LON2.flatten(),LAT2.flatten())))	
nncmems=tree.query((xp,yp))[1] #tree0=cKDTree(list(zip(s.lon,s.lat)))
ii,jj=np.unravel_index(nncmems,LON2.shape)
tscenter=cmems['sossheig'][:,ii,jj].values
tcmems=datetime64_to_datetime(cmems['time'].values)
tscenter_schism=ds['elev'][:,nncmems].values

fig=plt.figure()
plt.clf()
plt.plot(tcmems,tscenter)
plt.plot(dates,tscenter_schism)
plt.xlim(dates[0],dates[-1])
plt.legend(['cmems','schism'])
plt.ylabel('ssh [m]')
plt.grid()
fig.autofmt_xdate()
plt.tight_layout()
plt.savefig('basin_center_ts.png',dpi=300)


S=tscenter_schism
fs=1/np.diff(dates[:2])[0].total_seconds()
#fs=1/24/3600 #1 day in Hz (sampling frequency)
nyquist = fs / 2 # 0.5 times the sampling frequency
cuttofdays=180
lowcut=1/cuttofdays*nyquist*24*3600


cuttofdays=180
lowcut=1/cuttofdays*nyquist*24*3600


#cutoff=0.1 # fraction of nyquist frequency, here  it is 5 days
print('cutoff= ',1/cutoff*nyquist*24*3600,' days') #cutoff=  4.999999999999999  days
#b, a = signal.butter(5, cutoff, btype='lowpass') #low pass filter
b, a = signal.butter(5, cutoff, btype='highpass') #low pass filter
dUfilt = signal.filtfilt(b, a, S)
dUfilt=np.array(dUfilt)

c, d = signal.butter(5, cutoff, btype='lowpass') #low pass filter
dUfilt_low = signal.filtfilt(c, d, S)
dUfilt_low=np.array(dUfilt_low)
plt.plot(dates,dUfilt_low)

#cuttofdays=180
#lowcut=1/cuttofdays*nyquist*24*3600
#cuttofdays=1
#highcut=1/cuttofdays*nyquist*24*3600
days=180
lowcut=1/(days*86400)
days=1
highcut=1/(days*86400)


y=butter_bandpass_filter(S, lowcut, highcut, fs, order=5)
y2=butter_pass_filter(tscenter, lowcut, 1/86400, order=5,type='high')
y3=butter_pass_filter(tscenter, lowcut, 1/86400, order=5,type='low')

fig=plt.figure()
plt.clf()
plt.plot(tcmems,y2)
plt.plot(dates,y)
plt.legend(['cmems (180d )','schism (180d 1d)'])
plt.ylabel('ssh [m]')
plt.grid()
plt.tight_layout()
plt.xlim(dates[0],dates[-1])
fig.autofmt_xdate()
plt.savefig('basin_center_ts_filter.png',dpi=300)


from scipy.fft import fft

fft=scipy.fft
# Number of sample points
N = 600
# sample spacing
T = 1.0 / 800.0
x = np.linspace(0.0, N*T, N)
y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
yf = fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
import matplotlib.pyplot as plt
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
plt.show()



fig=plt.figure()
plt.clf()
plt.plot(tcmems,np.log(y2))
plt.plot(dates,np.log(y))
plt.legend(['cmems (180d )','schism (180d 1d)'])
plt.ylabel('log - ssh [m]')
plt.grid()
plt.tight_layout()
plt.xlim(dates[0],dates[-1])
fig.autofmt_xdate()
plt.savefig('basin_center_ts_filter_log.png',dpi=300)






# daily mean?

plt.contour()

ph,ch3=s.plotAtnodes(stats['elev']['mean'][1,:])
cs=plt.tricontour(s.lon,s.lat,stats['elev']['mean'][1,:], np.arange(-0.2,0.00,0.02,), linewidths=0.5, colors='k')
plt.clabel(cs)
s.plot_domain_boundaries(append=True)




plt.figure()
plt.pcolormesh(lon2,lat2,schismNN)



dists,coords=tree.query(list(zip(LON2.flatten(),LAT2.flatten())))


ii,jj=np.unravel_index(coords,LON2.shape)

plt.subplot(2,2,1)
plt.pcolor(lon,lat,duacs_var,cmap=cmap)
plt.title('duacs')
ch=plt.colorbar()
vmin, vmax = plt.gci().get_clim()
vmax+=0.005
plt.clim((vmin,vmax))
plt.subplot(2,2,2)
plt.pcolor(lon2,lat2,cmems_var,cmap=cmap)
plt.title('cmems')
ch2=plt.colorbar()
plt.clim((vmin,vmax))
plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
plt.subplot(2,2,3)
ph,ch3=s.plotAtnodes(stats['elev']['var'][1,:])
plt.title('schism')
plt.suptitle('variance of 2017 ssh(anomalies)')
plt.clim((vmin,vmax))
plt.subplot(2,2,3)
plt.axis((27.0625, 41.9375, 40.0625, 46.9375))
plt.tight_layout()
plt.savefig('variance_ssh.png',dpi=300)

stats0=calc_stats(ds2,varnames=['elev'])

stats_atmo=calc_stats(ds,varnames=['sensible_flux','latent_heat','upward_longwave','downward_longwave','total_heat_flux'])

plt.figure()
for i,key in enumerate(stats_atmo.keys()):
	plt.subplot(3,2,1+i)
	#s.plotAtnodesGeo(stats_atmo[key]['mean'][iyear])
	s.plotAtnodes(stats_atmo[key]['mean'][iyear])
	plt.title('Ann. mean ' +key + '[W/m^2] ')
plt.suptitle(stats['years'][iyear])
plt.tight_layout()
stats_atmo2=calc_stats(ds2,varnames=['sensible_flux','latent_heat','upward_longwave','downward_longwave','total_heat_flux'])

plt.figure()
for i,key in enumerate(stats_atmo.keys()):
	plt.subplot(3,2,1+i)
	#s.plotAtnodesGeo(stats_atmo[key]['mean'][iyear])
	s.plotAtnodes(stats_atmo2[key]['mean'][iyear])
	plt.title('Ann. mean ' +key + '[W/m^2] ')
plt.suptitle('2008')
plt.tight_layout()
stats_atmo2=calc_stats(ds2,varnames=['sensible_flux','latent_heat','upward_longwave','downward_longwave','total_heat_flux'])

plt.figure()
for i,key in enumerate(stats_atmo.keys()):
	plt.subplot(3,2,1+i)
	#s.plotAtnodesGeo(stats_atmo[key]['mean'][iyear])
	s.plotAtnodes(stats_atmo[key]['mean'][iyear]-stats_atmo2[key]['mean'][iyear])
	plt.title('Ann. mean ' +key + '[W/m^2] ')
	plt.gca().set_xticklabels([])
plt.suptitle('20017-2008:')
plt.tight_layout()
stats_atmo2=calc_stats(ds2,varnames=['sensible_flux','latent_heat','upward_longwave','downward_longwave','total_heat_flux'])

#sum=stats_atmo['downward_longwave']['mean'][iyear]-stats_atmo['upward_longwave']['mean'][iyear]-
#stats_atmo['latent']['mean'][iyear]-stats_atmo['latent']['mean'][iyear]
#radiation missing all 0 bcause

for i,key in enumerate(stats_atmo.keys()):
	plt.figure()
	plt.subplot(2,1,1)
	s.plotAtnodes(stats_atmo2[key]['mean'][iyear])
	vmin0, vmax0 = plt.gci().get_clim()
	plt.title(2008)
	plt.subplot(2,1,2)
	s.plotAtnodes(stats_atmo[key]['mean'][iyear])
	plt.title(2017)
	vmin1, vmax1 = plt.gci().get_clim()
	plt.suptitle('Ann. mean ' +key + '[W/m^2] ')
	plt.subplot(2,1,1)
	plt.clim((np.min((vmin0,vmin1)),np.max((vmax0,vmax1))))
	plt.subplot(2,1,2)
	plt.clim((np.min((vmin0,vmin1)),np.max((vmax0,vmax1))))
	plt.tight_layout()	

	
	
stats_atmo2=calc_stats(ds2,varnames=['sensible_flux','latent_heat','upward_longwave','downward_longwave','total_heat_flux'])


bstats=calc_stats(ds,varnames=['temp','salt'])
bstats2=calc_stats(ds2,varnames=['temp','salt'])

varnames=['hvel','zcor','salt','temp']# ,'hvel']
#stats[varname]['mean']=np.asarray(means[varname][:,:,:])
#stats[varname]['var']=np.asarray(vars[varname][:,:,:])

iyear=1
periods=['year']	
yearly=ds.groupby('time.'+periods[0])
spatial_variability=(ds['wind_speed'][:,:,0]**2+ds['wind_speed'][:,:,1]**2).mean(axis=0).values-(stats['wind_speed']['mean']**2).sum(axis=2)
u=stats['wind_speed']['mean'][iyear,:,0]
v=stats['wind_speed']['mean'][iyear,:,1]
lon=np.asarray(s.lon)
lat=np.asarray(s.lat)




ph,ch=s.plotAtnodes(spatial_variability[iyear,:])
ch.set_label('wind variability [m^2 s^-2]')
plt.clim((0,50))
n=40 # arrows shown along one axis
xlim=plt.xlim()
ylim=plt.ylim()
x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/n)
y=np.arange(ylim[0],ylim[1],(ylim[1]-ylim[0])/n)
X, Y = np.meshgrid(x,y)
positions = np.vstack([X.ravel(), Y.ravel()])
s.init_node_tree()
d,qloc=s.node_tree_latlon.query(positions.transpose())
xref,yref=np.asarray(plt.axis())[[1,2]] +  np.diff(plt.axis())[[0,2]]*[- 0.2, 0.1]
vmax=2
quiver=plt.quiver(np.concatenate((lon[qloc],(xref,))),np.concatenate((lat[qloc],(yref,))),np.concatenate((u[qloc],(vmax,))),np.concatenate((v[qloc],(0,))),scale=5*vmax,scale_units='inches') #self.maxfield.get()
arrowlabel=plt.text(xref,yref*1.01,'\n'*3+str(vmax)+' m/s')
plt.title(stats['years'][iyear])


ph,ch=s.plotAtnodes(stats['zcor']['mean'][iyear,:]+0.06)
ch.set_label('mean ssh [m')

yearly=ds.groupby('time.'+periods[0])	
	
plt.clf()	
iyear=1
varname='zcor'
synonym=dict.fromkeys(varnames)
synonym['zcor']='ssh'
synonym['temp']='sst'
synonym['salt']='sss'
synonym['hvel']='surf vel'

plt.subplot(2,2,1)
ph,ch=s.plotAtnodes(stats[varname]['mean'][iyear,:])
plt.title(stats['years'][iyear])
ch.set_label('mean '+synonym[varname])
plt.subplot(2,2,2)
ph,ch=s.plotAtnodes(np.sqrt(stats[varname]['var'][iyear,:]))
ch.set_label('std '+synonym[varname])
plt.subplot(2,2,3)
ph,ch=s.plotAtnodes(stats[varname]['mean'][iyear,:])
plt.title(stats['years'][iyear])
ch.set_label('mean '+synonym[varname])
plt.subplot(2,2,4)
ph,ch=s.plotAtnodes(np.sqrt(stats[varname]['var'][iyear,:]))
ch.set_label('std '+synonym[varname])

plt.suptitle(stats['years'][iyear])
plt.subplot(2,2,1)
plt.title('mean')
plt.subplot(2,2,2)
plt.title('variability')
plt.subplot(2,2,3)
plt.clim((-0.07,0.07))
plt.title('limits as Stanev and Ricker')
plt.subplot(2,2,4)
plt.clim((0.07,0.15))
plt.title('limits as Stanev and Ricker')
plt.savefig('blacksea_'+synonym[varname]+'_'+str(stats['years'][iyear])+'.png',dpi=300)


###
varname='wind_speed'
ph,ch=s.plotAtnodes((stats[varname]['magnitude_var'][iyear,:])
iyear=1
for varname in varnames[2:]:
	plt.clf()	
	plt.subplot(2,1,1)
	ph,ch=s.plotAtnodes(stats[varname]['mean'][iyear,:])
	plt.title(stats['years'][iyear])
	ch.set_label('mean '+synonym[varname])
	plt.subplot(2,1,2)
	ph,ch=s.plotAtnodes(np.sqrt(stats[varname]['var'][iyear,:]))
	ch.set_label('std '+synonym[varname])
	plt.suptitle(stats['years'][iyear])
	plt.savefig('blacksea_'+synonym[varname]+'_'+str(stats['years'][iyear])+'.png',dpi=300)