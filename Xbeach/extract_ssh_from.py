#!/usr/bin/env python
"""
Extract bd forcing for xbeach from SCHISM
using nn interpolation
"""

__author__ = "Benjamin Jacob"
__copyright__   = "Copyright 04\2022 - Helmholtz-Zentrum Hereon GmbH"
__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"

import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
import dask
dask.delayed()
import xarray as xr
plt.ion()
import time


################# settings #############################################################################

rundir='/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/'
ncdir=rundir+'combined/'


#I only need four points for the periods: 2021.06.01~2021.09.30, with time interval 1200 s
#Time    (5968200, 379270)    (5968200, 384510)    (5954200, 384510)    (5954200, 379270)

# coordinates # here utm as schism grid
# x,y pairs
#coords=[(5968200, 379270),    (5968200, 384510),    (5954200, 384510),    (5954200, 379270)] #corner coordinates from xbeach

#pairs of x,y coordinates
coords=[(5968200, 379270)[::-1],    (5968200, 384510)[::-1],    (5954200, 384510)[::-1],    (5954200, 379270)[::-1]] #corner coordinates from

########################################################################################

# load schism grid and netcdf acess
os.chdir(rundir)
s=schism_setup()
s.ds=schism_output2(ncdir).nc # load handle

# select neatest negibouts
s.init_node_tree(latlon=False) # initialie ckd tree 
Pnn=s.node_tree_xy.query(coords)[1] 
s.dsnn=s.ds.sel(nSCHISM_hgrid_node=Pnn) ## sele nearest neighbours within data set

#control location
plt.close('all')	
s.plot_domain_boundaries(latlon=False,append=True)
for i,coord in enumerate(coords):
	plt.plot(coord[0],coord[1],'r+')
	plt.text(coord[0],coord[1],'P'+str(i))
xmin,ymin=np.min(np.asarray(coords),axis=0)
xmax,ymax=np.max(np.asarray(coords),axis=0)
plt.xlim((xmin-5000,xmax+5000))
plt.ylim((ymin-5000,ymax+5000))
plt.savefig('timeseries_location_export4xbeach.png',dpi=300)

#inerpolate time

t0=s.ds['time'][0].values
t1=s.ds['time'][-1].values
dtout=np.timedelta64(1200,'s') # output time step (linear interpolation)
timeout=np.arange(t0,t1,dtout)
dsout=s.dsnn['elev'].interp(time=timeout)

elevout=dsout.values
tout= (timeout-timeout[0])/np.timedelta64(1,'s') #/dtout

headertext='time('+str(t0)[:10]+'-'+str(t1)[:10]+')'+''.join([' P{:d}{:s}'.format(i,str(coords[i])) for i in range(len(coords))])
np.savetxt('SSH_out'+str(t0)[:10]+'-'+str(t1)[:10]+'.ascii',dataout,header=headertext)

plt.close('all')	
plt.clf()
plt.subplot(2,1,1)
plt.plot(timeout,elevout)
plt.subplot(2,1,2)
plt.plot(timeout[-72:],elevout[-72:])
plt.legend(('P0','P1','P2','P3'))
plt.tight_layout()
plt.savefig('timeseries_export4xbeach.png',dpi=300)