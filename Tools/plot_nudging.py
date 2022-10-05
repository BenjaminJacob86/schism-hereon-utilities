import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from matplotlib.path import Path
from schism import *
import xarray as xr
from glob import glob
plt.ion()
s=schism_setup()

include=np.loadtxt('include.gr3',skiprows=2,max_rows=s.nnodes)[:,-1]
include=np.asarray(include,bool)

S=xr.open_dataset('SAL_nu.nc')
T=xr.open_dataset('TEM_nu.nc')

S2=xr.open_mfdataset(['tsfiles/S_{:d}.nc'.format(ti) for ti in range(1,10)])
T2=xr.open_mfdataset(['tsfiles/T_{:d}.nc'.format(ti) for ti in range(1,10)])

lon0,lon1=np.min(s.lon),np.max(s.lon)
lat0,lat1=np.min(s.lat),np.max(s.lat)


S2=S2.sel(lon=slice(lon0-1,lon1+1),lat=slice(lat0-1,lat1+1))
T2=T2.sel(lon=slice(lon0-1,lon1+1),lat=slice(lat0-1,lat1+1))

ti=0
sv=np.zeros(s.nnodes)
tv=np.zeros(s.nnodes)

sv[include]=S['tracer_concentration'][ti,:,-1,0]
tv[include]=T['tracer_concentration'][ti,:,-1,0]

plt.figure()
plt.subplot(2,2,1)
s.plotAtnodes(tv)
plt.clim((0,12))
plt.subplot(2,2,2)
s.plotAtnodes(sv)
plt.clim((10,33))
plt.subplot(2,2,3)
T2['thetao'][0,0,:].plot(cmap=plt.cm.jet,vmin=0,vmax=12)
plt.clim((0,12))
plt.subplot(2,2,4)
S2['so'][0,0,:].plot(cmap=plt.cm.jet,vmin=10,vmax=30)
plt.clim((10,33))
plt.tight_layout()
plt.savefig('controlPlotSurface',dpi=300)



sv[include]=S['tracer_concentration'][ti,:,0,0]
tv[include]=T['tracer_concentration'][ti,:,0,0]
plt.figure()
plt.subplot(1,2,1)
s.plotAtnodes(tv)
plt.clim((0,12))
plt.subplot(1,2,2)
s.plotAtnodes(sv)


plt.figure()
plt.subplot(2,2,1)
s.plotAtnodes(tv)
plt.clim((0,12))
plt.subplot(2,2,2)
s.plotAtnodes(sv)
plt.clim((10,33))
plt.subplot(2,2,3)
T2['thetao'][0,0,:].plot(cmap=plt.cm.jet,vmin=0,vmax=12)
plt.clim((0,12))
plt.subplot(2,2,4)
S2['so'][0,0,:].plot(cmap=plt.cm.jet,vmin=10,vmax=30)
plt.clim((10,33))
plt.tight_layout()
plt.savefig('controlPlotBottom',dpi=300)