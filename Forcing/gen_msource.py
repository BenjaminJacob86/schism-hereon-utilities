# generate temperature / salinity input for river 
# from next neighbour atmosphere temperatures from
#sflux
import os
import sys
sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib')
import matplotlib
#matplotlib.use('AGG')
from schism import *
from netCDF4 import Dataset
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from matplotlib import pyplot as plt
import xarray
import netCDF4
import datetime as dt
import glob 

setupdir='./'
atmodir=setupdir+'sflux/'
# load schism setup
os.chdir(setupdir)
s=schism_setup()

# get elem_ids
elemids=[]
with open('source_sink.in') as f:
	nr=np.int(f.readline().split('!')[0])
	for i in range(nr):
		elemids.append(np.int(f.readline().split('!')[0]))
elemids=np.asarray(elemids)

# tri only
#nv=np.asarray(s.nv)
lon=np.asarray(s.lon)
lat=np.asarray(s.lat)
#xcoords=lon[nv-1].mean(axis=1)[elemids]
#ycoords=lat[nv-1].mean(axis=1)[elemids]

cx=[]
cy=[]
nv=np.asarray(s.nv)
for elem in nv:
	cx.append(lon[np.asarray(elem)-1].mean())
	cy.append(lat[np.asarray(elem)-1].mean())
cx=np.asarray(cx)[elemids-1]
cy=np.asarray(cy)[elemids-1]
xcoords=cx
ycoords=cy

s.plot_domain_boundaries()
plt.plot(cx,cy,'ko')

m=np.loadtxt('vsource.th')
dt_flux=np.diff(m[:2,0])

# atmosphere grid next neighbours		
print('determine next neighbours on atmosphere grid')
nc=Dataset(atmodir+'sflux_air_1.0001.nc')
ncv=nc.variables			
tree=cKDTree(list(zip(ncv['lon'][:].flatten(),ncv['lat'][:].flatten()))) # nextneighbour look up tree
d,inds=tree.query(list(zip(xcoords,ycoords)))
ii,jj=np.unravel_index(inds,ncv['lat'].shape)

# atmosphere timestep	
dt_atmo=np.round(np.diff(ncv['time'][:2])*86400)
step=np.int(dt_flux/dt_atmo) # day means
try:
	startdate = netCDF4.num2date(ncv['time'][0],ncv['time'].units)
except:
	startdate = dt.datetime(ncv['time'].base_date[0],ncv['time'].base_date[1],ncv['time'].base_date[2],ncv['time'].base_date[1],ncv['time'].base_date[3])
nc.close()

# load temp data from atmo
print('load surface temperature from atmosphere at river nn nodes')
files=np.sort(glob.glob(atmodir + '*air*.nc'))
xr=xarray.open_mfdataset(list(files),concat_dim='time')
nt,nx,ny=xr['stmp'].shape	
ti=np.arange(step/2,nt,step,int) # half to get to noon in units of sflux
blubb=np.asarray(xr['stmp'][ti-1,ii,jj])[:,range(len(ii)),range(len(jj))]
xr.close()


ts=np.column_stack((m[:,0],blubb[:m.shape[0],:],np.zeros((m.shape[0],m.shape[1]-1))))


ts[:,1:np.int((ts.shape[1]-1)/2+1)]-=273.15
# set limits
ts[:,1:m.shape[1]]=np.minimum(np.maximum(ts[:,1:m.shape[1]],-1),25)


print('write  output')

plt.clf()
s.plot_domain_boundaries()
plt.scatter(cx,cy,c=ts[0,1:np.int((ts.shape[1]-1)/2+1)],s=20)
ch=plt.colorbar()
ch.set_label('Temp [degC]')
np.savetxt('vsource_tempFromSflux.th',ts)
#np.savetxt( ('TEM_1.th_{}_{}'.format(str(startdate),str(startdate+dt.timedelta(ti[-1]/86400)))).replace(' ','_'),np.column_stack((ti,blubb)) )
dates=startdate+ ti * dt.timedelta(seconds=1)
fig=plt.figure()
plt.plot(dates,ts[:,1:])
plt.grid()
plt.tight_layout()
fig.autofmt_xdate()
plt.savefig('river_temperatures_from_air.eps',format='eps',dpi=600)
#a=np.asarray(xr['stmp'][:,ii[0],jj[0]])
#plt.plot(np.arange(nt)/24,a)
print('write done')
