"""
set boundary connected nodes to a minimum value to prevent dry bound error message
"""
import glob
import sys 
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities')
from schism import *
import datetime as dt	
import os
plt.ion()

/work/gg0028/g260114/SETUPS/NWBS/setup/newcode/Forcing/data/nrt.cmems-du.eu/Core/BLKSEA_ANALYSISFORECAST_PHY_007_001/all

dmin=3
outname='hgrid_minbd_depth{:.2f}.gr3'.format(dmin)
s=schism_setup()

elemTree=s.nvplt+1 # pure trianlge (quads splitted) element tree counting nodes from 1

nodes2deepen=[]
for bdy in s.bdy_segments:
	for bdnde in np.asarray(bdy):
		ind=np.where((bdnde == elemTree).sum(axis=1))[0]
		nodes2deepen=np.hstack((nodes2deepen,np.unique(elemTree[ind,:])))
nodes2deepen=np.unique(nodes2deepen)	
nodes2deepen-=1 # make 0 based again for depth counting

s.depths=np.asarray(s.depths)
nodes2deepen=np.asarray(nodes2deepen,int)
s.depths[nodes2deepen]=np.maximum(s.depths[nodes2deepen],dmin)
s.dump_hgridgr3(filename=outname)
os.rename('hgrid.gr3','hgrid.gr3_opbd_not_deepened')
os.rename(outname,'hgrid.gr3')
