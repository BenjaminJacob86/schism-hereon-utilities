"""
 Hotstart: overwrite salinty values in defined polygin
 with zero to correct interpolation from coarser productts
"""
import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *
from netCDF4 import Dataset
from matplotlib import path


ncfilein='GB_hot_jun202301.nc'
ncfileout='GB_hot_jun202301_river.nc'


s=schism_setup()
if not (os.path.exists('rivers.reg')):
	exit()

riverPoly=np.loadtxt('river.reg',skiprows=3)
riverPoly=path.Path(list(zip(riverPoly[:,0],riverPoly[:,1])))
setsalt0=riverPoly.contains_points(list(zip(s.x,s.y)))


# read file
nc=Dataset(ncfilein)
tr_nd=nc['tr_nd'][:]
tr_nd[setsalt0,:,1]=0
nc.close()
os.rename(ncfilein,ncfilein+'bck_up')
s.write_hotstart(tr_nd,filename=ncfilein)

# plot
suffix=['surface', 'bottom' ]
inds=[-1, 0]
# plot surf bottmom layer
for i, lab in zip(inds,suffix): 
	ph,ch=s.plotAtnodes(s[:,i])
	ch.set_label('salinity g/kg')
	plt.title(lab)
	plt.savefig('salinity_correct'+lab,dpi=600)
	plt.close()