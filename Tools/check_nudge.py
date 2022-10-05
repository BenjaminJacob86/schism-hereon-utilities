import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
from schism import *
import xarray as xr
plt.ion()



s=schism_setup(vgrid_file='vgrid.in.old')

include=np.loadtxt('include.gr3',skiprows=2,max_rows=s.nnodes)



#include=np.loadtxt('TEM_nudge.gr3',skiprows=2,max_rows=s.nnodes)
#s.plotAtnodes(include[:,-1]>0.0)
s.plotAtnodes(include[:,-1])

dst=xr.open_dataset('TEM_nu.nc')
dss=xr.open_dataset('SAL_nu.nc')
Tnu=dst['tracer_concentration'][0,:,-1,0].values
Snu=dss['tracer_concentration'][0,:,-1,0].values
T=np.zeros(s.nnodes)
S=np.zeros(s.nnodes)
T[include[:,-1]>0]=Tnu
S[include[:,-1]>0]=Snu

plt.clf()
plt.subplot(2,1,1)
s.plotAtnodes(T)
plt.subplot(2,1,2)
s.plotAtnodes(S)
plt.savafig('nudgeTest.png',dpi=300)