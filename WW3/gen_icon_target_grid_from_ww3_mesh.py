""""
Generate Icon target mesh for global 13km grid based on ww3 mesh.
ww3mesh extend is rounded and target mesh in 1/8th degree 

target mesh is of type:
# CDO grid description file for regular grid Europe source ICON_GLOBAL
gridtype = lonlat
xsize = 1201
ysize = 601
xfirst = -75.0
xinc = 0.125
yfirst = 5.0
yinc = 0.125

"""

import os
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')
# own libraries 
from data_and_model_classes import WW3_mesh

ww3mesh='/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd_2017/longrun2/check_mergeVerena/NBSext_bl.msh'
ww3=WW3_mesh(ww3mesh)


xinc = 0.125
yinc = 0.125

# resolution
icon_target_mesh=ww3mesh[ww3mesh.rindex('/')+1:ww3mesh.rindex('.')]+'_target.txt'

xfirst=np.floor(np.min(ww3.x))
yfirst=np.floor(np.min(ww3.y))
xlast=np.ceil(np.max(ww3.x))
ylast=np.ceil(np.max(ww3.y))
xsize=np.int(np.ceil(xlast-xfirst)/xinc+1)
ysize=np.int(np.ceil(ylast-yfirst)/yinc+1)

with open(icon_target_mesh,'w') as f:
	f.write('# CDO grid description file for regular grid Europe source ICON_GLOBAL to {:s}\n'.format(icon_target_mesh))
	f.write('gridtype = lonlat\n')
	f.write('xsize = {:d}\n'.format(xsize))
	f.write('ysize = {:d}\n'.format(ysize))
	f.write('xfirst = {:f}\n'.format(xfirst))
	f.write('xinc = {:f}\n'.format(xinc))
	f.write('yfirst = {:f}\n'.format(yfirst))
	f.write('yinc = {:f}\n'.format(yinc))
	
	
