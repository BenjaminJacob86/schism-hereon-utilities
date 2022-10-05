# libraries
import os
import netCDF4
import sys
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import netCDF4
import pickle 
plt.ion()

#own libraries
sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/gpfs/work/jacobb/data/PreProc/Setup2Other/')
from schism import *


os.chdir('/gpfs/work/jacobb/data/RUNS/BlackSea/RUN24d/')
s0=schism_setup()

os.chdir('/gpfs/work/jacobb/data/SETUPS/NWBlackSea/setup2')
s1=schism_setup()

nvrt0=np.asarray([s0.vgrid[key].count()    for key in s0.vgrid.keys()])
nvrt1=np.asarray([s1.vgrid[key].count()    for key in s1.vgrid.keys()])


s1.plotAtnodes(nvrt1)
ax=plt.axis()
plt.figure()
s0.plotAtnodes(nvrt0)
plt.axis(ax)

bdnodes=np.asarray(s1.bdy_segments[0])-1



d=np.asarray(s1.depths)[bdnodes]

lon=np.asarray(s1.lon)[bdnodes]
lat=np.asarray(s1.lat)[bdnodes]
i=0
dm=np.vstack([s1.vgrid[bdnodes[i]+1].filled(-1)*d[i] for i in range(len(bdnodes))])
lat2=np.tile(lat,(dm.shape[1],1)).T

# new grid
plt.clf()
plt.plot(lat2,dm)

s0.init_node_tree(latlon=True)
d0=np.asarray(s0.depths)[nn0]

nn0=s0.node_tree_latlon.query(list(zip(lon,lat)))[1]
dm0=np.vstack([s0.vgrid[nn0[i]+1].filled(-1)*d0[i] for i in range(len(nn0))])
lat02=np.tile(lat,(dm0.shape[1],1)).T


ax,phs=plt.subplots(2,1,sharex=True,sharey=True)
phs[0].plot(lat2,dm0)
phs[1].plot(lat2,dm)

plt.subplot(2,1,1)
plt.plot(lat2,dm0)
plt.title('OG BS')
plt.subplot(2,1,2)
plt.plot(lat2,dm)
plt.title('NW BS')


