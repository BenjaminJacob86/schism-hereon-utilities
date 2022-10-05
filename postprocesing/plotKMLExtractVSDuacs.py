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
from cftime import utime
########## settings #################################
# export OMP_NUM_THREADS=4 # call before python

#rundir='/gpfs/work/jacobb/data/RUNS/BlackSea/2017/'#'/gpfs/work/jacobb/data/SETUPS/BlackSea/'
rundir='/work/gg0028/g260114/RUNS/Europe/c1/flow_tweak/'
ncdir=rundir+'combined/'#'/gpfs/work/jacobb/data/SETUPS/BlackSea/combined/'


##############
fname='EuropeBasinExtracts.nc'
plt.ion()
nc=Dataset(fname)   
t0=dt.datetime.strptime(nc['time'].units[11:],'%Y-%m-%d %H:%M:%S')
dates=t0+np.arange(1,len(nc['time'])+1)*dt.timedelta(hours=6)
year=np.asarray([ti.year for ti in dates])
#names=[key[:key.index('_')]  for key in nc.variables.keys() if key != 'time' ]
names=[key[:key.rindex('_')]  for key in nc.variables.keys() if key != 'time' ]




duacs=Dataset('duacs_2012_to_2018.nc') #xr.open_dataset('duacs_2012_to_2018.nc')
refLon=34+19/60+16.57/3600
refLat=43+11/60+54.71/3600

ilat=np.argmin(np.abs(duacs['latitude'][:]-refLat ))
ilon=np.argmin(np.abs(duacs['longitude'][:]-refLon ))
ut=utime(duacs['time'].units)
dates_duacs=ut.num2date(duacs['time'][:])
#plt.plot(dates_duacs,duacs['sla'][:,ilat,ilon],'m',linewidth=2)

plt.figure()
plt.pcolormesh(duacs['longitude'],duacs['latitude'],duacs['sla'][0,:])
plt.plot(duacs['longitude'][ilon],duacs['latitude'][ilat],'ro')
plt.colorbar()

plt.figure()
plt.plot_date(dates_duacs,np.nanmean(np.nanmean(duacs['sla'],axis=1),axis=1))

meanref=dates[3::4]
daymeans=dict.fromkeys(list(nc.variables.keys())[1:])
monmeans=dict.fromkeys(list(nc.variables.keys())[1:])
yearmeans=dict.fromkeys(list(nc.variables.keys())[1:])
years=np.asarray([date.year for date in dates])		
months=np.asarray([date.month for date in dates])		

#for varname in nc.variables.keys():
varname='BlackSea_local_elev'
demean=True
if 1:
	if varname != 'time':
		daymeans[varname]=(nc[varname][3:-3]+nc[varname][4:-2]+nc[varname][5:-1]+nc[varname][6:])[::4]/4
		ts=nc[varname][:]
		tsmean=[]
		if demean==True:
			daymeans[varname]-=ts.mean()
			ts-=ts.mean()
		if ts.shape[-1]==2:
			yearmeans[varname]=np.asarray([ts[years==year].mean(axis=0) for year in range(2012,2017)])
			for year in range(2012,2017):
				for mon in range(1,13):
					tsmean.append(ts[(years==year) & (months==mon), :].mean(axis=0))
		else:
			yearmeans[varname]=np.asarray([ts[years==year].mean(axis=0) for year in range(2012,2017)])
			for year in range(2012,2017):
				for mon in range(1,13):
					tsmean.append(ts[(years==year) & (months==mon) ].mean(axis=0))
		monmeans[varname]=np.asarray(tsmean)

monref=[]		
for year in range(2012,2017):
	for mon in range(1,13):
		monref.append(dt.datetime(year,mon,15))
monref=np.asarray(monref)

yearref=np.asarray([dt.datetime(year,6,15) for year in range(2012,2017)])
#ivalid=~np.isnan(monmeans['AlboranSea_mean_elev'])
ivalid=~np.isnan(monmeans[varname])
#ivalid=~np.isnan(monmeans[list(monmeans.keys())[0]])
fit=np.polyfit(np.arange(len(monref))[ivalid]+1,monmeans[varname][ivalid],1)


# daily means
plt.clf()
count=0
varname='BlackSea_local_elev'
#for i,varname in enumerate(list(nc.variables.keys())[1:]):
i=list(nc.variables.keys()).index(varname)
if True:
	#if 'elev' in varname:
	if len( daymeans[varname].shape)==1:
		plt.clf()
		plt.plot(meanref[:-1],daymeans[varname],linewidth=1)
		plt.plot(monref,monmeans[varname],linewidth=3)
		plt.plot(yearref,yearmeans[varname],'k',linewidth=2)
		fit=np.polyfit(np.arange(len(monref))[ivalid]+1,monmeans[varname][ivalid],1)
		plt.plot(monref[ivalid],np.polyval(fit,np.arange(len(monref))+1)[ivalid],'r')
		plt.text(monref[ivalid][-15],plt.ylim()[0]+np.diff(plt.ylim())*0.9,'y= {:f} x + {:f}'.format(fit[0],fit[1]),color='r')
		plt.title(names[i])
		plt.grid()
		#plt.ylim((-0.8,0))
		if 'elev' in varname: 
			plt.ylabel('ssh anomaly [m]')
		elif 'salt' in varname: 
			plt.ylabel('salinity [psu]')
		elif 'temp' in varname: 
			plt.ylabel('tempertature [degC]')
		elif 'pressure' in varname: 
			plt.ylabel('pressure [Pa]')
		#plt.legend(['daily mean','monthly mean','annual mean','monthly trend'],ncol=2,loc='upper left')
		#plt.tight_layout()
		plt.savefig(varname+'.png',dpi=300)

		
		
plt.plot_date(dates,ts,'b',linewidth=2)
plt.plot_date(dates_duacs,duacs['sla'][:,ilat,ilon],'m',linewidth=2)		
plt.grid()
plt.legend(['schism ','duacs'])
plt.title('SLA anomalie in central Black Sea ')

plt.legend(['daily mean','monthly mean','annual mean','monthly trend','DUACS'],ncol=3,loc='upper left')		
plt.savefig(varname+'_duacs_'+'.png',dpi=300)
!mkdir salt temp elev wind pressure
!mv *temp.png temp/
!mv *salt.png salt/
!mv *elev.png elev/
!mv *pressure.png pressure/
!mv *speed.png wind/

		
## daily means
#plt.clf()
#count=0
#for i,varname in enumerate(list(nc.variables.keys())[1:]):
#	if 'elev' in varname:
#		count+=1
#		plt.subplot(3,2,count)
#		plt.plot(meanref[:-1],daymeans[varname],linewidth=1)
#		plt.plot(monref,monmeans[varname],linewidth=3)
#		plt.plot(yearref,yearmeans[varname],'k',linewidth=2)
#		fit=np.polyfit(np.arange(len(monref))[ivalid]+1,monmeans[varname][ivalid],1)
#		plt.plot(monref[ivalid],np.polyval(fit,np.arange(len(monref))+1)[ivalid],'r')
#		plt.text(monref[ivalid][-15],plt.ylim()[1]*0.9,'y= {:f} x + {:f}'.format(fit[0],fit[1]),color='r')
#		plt.title(names[i])
#		plt.grid()
#		plt.ylim((-0.8,0))
#		if i%2==0:
#			plt.ylabel('ssh [m]')
#plt.legend(['daily mean','monthly mean','annual mean','monthly trend'],ncol=4)
#plt.tight_layout()