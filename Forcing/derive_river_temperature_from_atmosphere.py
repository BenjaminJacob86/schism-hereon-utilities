# generate temperature input for river runoff 
# from next neighbour atmosphere temperature

""" Overwrite temperature in in sources (msource.th) with values derived
from air temperature
"""

import os
import sys
sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib')
import matplotlib
matplotlib.use('AGG')
from schism import *
from netCDF4 import Dataset
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from matplotlib import pyplot as plt
import xarray
import netCDF4
import datetime as dt
import glob 



setupdir='/mnt/lustre01/work/gg0028/g260114/RUNS/Europe/Europe6/'
atmodir=setupdir+'sflux/'

class sources:		
	"""	load source fluxes for schism setup """
	
	def __init__(self,s,fname='source_sink.in',comments='!'):
		self.i_src=np.asarray(np.loadtxt(fname,comments=comments)[1:-1],int)-1
		#self.elem=np.asarray(s.nv)[self.i_src,:]-1	
		#self.x=np.asarray(s.x)[self.elem].mean(axis=1)
		#self.y=np.asarray(s.y)[self.elem].mean(axis=1)
		# also for quads:
		self.x=[]
		self.y=[]
		x=np.asarray(s.x)
		y=np.asarray(s.y)
		for ielem in self.i_src:
			elem=np.asarray(s.nv[ielem])
			self.x.append(x[elem-1].mean())
			self.y.append(y[elem-1].mean())
		self.x=np.asarray(self.x)	
		self.y=np.asarray(self.y)	
		
		
		with open(fname) as f:
			n=np.int(f.readline().split(comments)[0])
			#self.names=[f.readline().split(comments)[1].split(',')[0].split('river')[1] for i in range(n)]
			self.names=[f.readline().split(comments)[1] for i in range(n)]
			
		# volume
		try:
			self.vsource=np.asarray(np.loadtxt('vsource.th',comments=comments),int)
			self.vtime=self.vsource[:,0]
			self.vsource=self.vsource[:,1:]
			
			# concentration i.e temperature
			self.msource=np.asarray(np.loadtxt('msource.th',comments=comments),int)
			self.mtime=self.msource[:,0]
			self.msource=self.msource[:,1:]	
		except:
			pass


os.chdir(setupdir)

# load schism setup
s=schism_setup()
srcs=sources(s)
dt_flux=np.diff(srcs.vtime[:2])


# atmosphere grid next neighbours		
print('determine next neighbours on atmosphere grid')
nc=Dataset(atmodir+'sflux_air_1.0001.nc')
ncv=nc.variables			
tree=cKDTree(list(zip(ncv['lon'][:].flatten(),ncv['lat'][:].flatten()))) # nextneighbour look up tree
d,inds=tree.query(list(zip(srcs.x,srcs.y)))
ii,jj=np.unravel_index(inds,ncv['lat'].shape)

# atmosphere timestep	
dt_atmo=np.round(np.diff(ncv['time'][:2])*86400)
step=np.int(dt_flux/dt_atmo) # day means
startdate = netCDF4.num2date(ncv['time'][0],ncv['time'].units)
nc.close()


# load temp data from atmo
print('load surface temperature from atmosphere at river nn nodes')
files=np.sort(glob.glob(atmodir + '*air*.nc'))
xr=xarray.open_mfdataset(list(files))
nt,nx,ny=xr['stmp'].shape	
ti=np.arange(step/2,nt,step,int) # half to get to noon in units of sflux
#ti=np.arange(step/2,365*24,step)
blubb=np.asarray(xr['stmp'][ti-1,ii,jj])[:,range(len(ii)),range(len(jj))]
xr.close()
#ti=np.arange(step/2,365*24,step)
ti=(ti-step/2)*3600  # convert to seconds starting at 0
blubb=np.maximum(blubb-273.15,0) #convert to celcius

msource_out=np.zeros((msource.shape[0],msource.shape[1]+1])
msource_out[:,1:]=msource
msource_out[:,1:len(srcs.x)+1]=blubb # update temperature values
msource_out[:,0]=ti

print('write  output')
np.savetxt('msource_new.th',msource_out)

#np.savetxt( ('TEM_1.th_{}_{}'.format(str(startdate),str(startdate+dt.timedelta(ti[-1]/86400)))).replace(' ','_'),np.column_stack((ti,blubb)) )
#dates=startdate+ ti * dt.timedelta(seconds=1)
#fig=plt.figure()
#plt.plot(dates,blubb)
#plt.grid()
#plt.tight_layout()
#fig.autofmt_xdate()
#plt.savefig('river_temperatures_from_air.eps',format='eps',dpi=600)
#a=np.asarray(xr['stmp'][:,ii[0],jj[0]])
#plt.plot(np.arange(nt)/24,a)

print('write done')
