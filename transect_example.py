#export OMP_NUM_THREADS=1 # call before python
import os
import sys
import matplotlib
from matplotlib import pyplot as plt
background=False
if background:
	matplotlib.use('Agg') # backend
else:
	plt.ion()	
import datetime as dt
import glob
from scipy.spatial import cKDTree
import datetime as dt
from scipy.interpolate import interp1d
# own and 3d party libraries
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')
from schism import * # import schism functions
import pandas as pd
import utide
from matplotlib.dates import date2num
import xarray as xr
from numba import jit
import time
plt.ion()


setupdir='/gpfs/work/jacobb/data/RUNS/BlackSea/RUN24d/evap/' # schism run directory, containts hgrid.* etc.
ncdir=setupdir + 'combined/'


ncdir2=setupdir + 'combined_galberin_test_no_evap/'
ncdir='/gpfs/work/jacobb/data/RUNS/BlackSea/RUN24d/no_evap/combined_compiled_with_prec_evap/'

# images ############################
pic_format='.png'
dpivalue=300
cmap='jet'
Rcircle=150 							# radius of colorcoded cricles used in scatter maps
limit_to_data=True   					# limit map plots to bounding box of data.

# Font sizes
Fincrease=1.5
#SMALL_SIZE = 10*Fincrease
SMALL_SIZE = 8*Fincrease
MEDIUM_SIZE = 10*Fincrease
BIGGER_SIZE = 12*Fincrease
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure
#############################################################

######### load SCHISM setup   ##################################
s=schism_setup()

s.nc=schism_output2(ncdir)
s.nc2=schism_output2(ncdir2)

#bpfile='/gpfs/work/jacobb/data/LOCATIONS/Blacksea/BP/alonbosporus_ll.bp'
bpfile='/gpfs/work/jacobb/data/LOCATIONS/Blacksea/BP/bosporus.bp'

s.nc.nc
s.trans=bp_file(bpfile)
s.trans=bp_transect(s,s.nc.nc,s.trans.x,s.trans.y)
s.trans2=bp_transect(s,s.nc2.nc,s.trans.x,s.trans.y)


s.plot_domain_boundaries(latlon=False)
plt.plot(s.trans.x,s.trans.y,'r-')

ti=239
z=s.trans.ds['zcor'][ti,:,:].values
t=s.trans.ds['temp'][ti,:,:].values
salt=s.trans.ds['salt'][ti,:,:].values
time=s.trans.ds['time'][ti].values
plt.pcolor(s.trans.l,z,t)
#plt.title(time)

ualong,uacross=s.trans.proj_hvel_along_accros(s.trans.ds,s.trans,s.trans.ds['hvel'][ti,:].values)
ualong2,uacross2=s.trans.proj_hvel_along_accros(s.trans2.ds,s.trans,s.trans2.ds['hvel'][ti,:].values)

1.0528e4*1.5

plt.figure()
plt.pcolor(s.trans.l,z,salt)
ch=plt.colorbar(extend='both')
plt.title('salt')


leg=['Galperin -0.56', 'Galperin -0.14']
plt.figure()
plt.clf()
for i,val in enumerate([ualong,ualong2]):
	plt.subplot(2,2,i+1)
	plt.pcolor(s.trans.l,z,val)
	plt.set_cmap('RdBu_r')
	plt.clim((-1,1))
	ch=plt.colorbar(extend='both')
	plt.plot(s.trans.l[:,0],-np.asarray(s.depths)[s.trans.nn],'k',linewidth=2)
	ch.set_label('along channel velocity [m/s]')
	plt.title(leg[i])
plt.suptitle(time)

plt.subplot(2,2,4)
plt.pcolor(s.trans.l,z,ualong2-ualong)
plt.set_cmap('RdBu_r')
plt.clim((-.05,.05))
plt.title('difference')
ch=plt.colorbar(extend='both')
plt.plot(s.trans.l[:,0],-np.asarray(s.depths)[s.trans.nn],'k',linewidth=2)
plt.tight_layout()	
plt.savefig('GalperinDifference',dpi=300)



plt.figure()
plt.clf()
plt.pcolor(s.trans.l,z,ualong)
plt.title(time)
plt.set_cmap('RdBu_r')
plt.clim((-1,1))
ch=plt.colorbar(extend='both')
plt.plot(s.trans.l[:,0],-np.asarray(s.depths)[s.trans.nn],'k',linewidth=2)
ch.set_label('along channel velocity [m/s]')



