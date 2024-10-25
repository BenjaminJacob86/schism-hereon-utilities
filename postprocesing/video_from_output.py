"""
# subplot surface slab of SCHISM and Amm15 and their difference 
# iterating over timesteps outputing images for video
# including time series at given location coords
"""

#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#export OPENBLAS_NUM_THREADS=1
#
#export OMP_NUM_THREADS=4
import os, sys
import numpy as np
import netCDF4
import csv
interactive=True
import matplotlib
if interactive:
	from matplotlib import pyplot as plt
	plt.ion()
else:	
	matplotlib.use('Agg')
	#matplotlib.use('Cairo')
	
import datetime as dt
import xarray as xr
import glob
from scipy.spatial import cKDTree
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
from matplotlib import path

# @strand 

# @levante
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hereon-utilities/')

#add schism python class to path
from schism import *
import warnings
warnings.filterwarnings("ignore")	


########## settings ##################################################
# directories (have to end with '/')
# SCHIMS grid files

setupdir=['/work/gg0028/g260247/RUNV3_2D/Ghana_V1/Ghana_V2/']  #where the rundir is 
ncdir=[setupdir[0]+'outputs/']   # where the outputs is
names=['Ghana V3']

# where to store image output
outdir='/work/gg0028/g260114/RUNS/jayson/images' #setupdir[0]+'images/' #'/work/gg0028/g260114/postproc/modelcomp/comp4/'
if not os.path.exists(outdir): os.mkdir(outdir) 



#varnames=['temp','ssh','salt']	#varname ['ssh',] if only one has to have ,
varnames=['ssh',]
sdictvnames = {'temp':'temperature','ssh':'zCoordinates','salt':'salinity'}
#min_max={'ssh':(-5,5),'salt':(0,25),'temp':(5,25)}	# axis range for variables
min_max={'ssh':(-5,5),}
# considrede time periods and steps for vriables
year0=2022
year1=2023
vartimes={'ssh':{'startdate':dt.datetime(year0,11,2,1,0),'enddate':dt.datetime(year1,11,2,1,0),'step[hours]':1},\
'salt':{'startdate':dt.datetime(year0,11,2,12,0),'enddate':dt.datetime(year1,2,1,12,0),'step[hours]':24},\
'temp':{'startdate':dt.datetime(year0,11,2,12,0),'enddate':dt.datetime(year1,2,1,12,0),'step[hours]':24},\
}
# limit later by s0.dates




# coords for time series
# helgoland
if False:
    lat,lon = 43, 29
    limx=((27.3,31.8))	#((-1.14,9.84))
    limy=((41,47))	#((49.7,56.21))

    dthour=1
    ndays=25 # Run durations


# apperence
cm=plt.cm.turbo # colormap
cmDiff=plt.cm.bwr
ax_limitEQ_commonDataRange=True    # minmize ax limits to common Extend which is not Land mask in any of the models
cm.set_over(color='m', alpha=None)
cm.set_under(color='k', alpha=None)
cm.set_bad(color='gray', alpha=None)


strdate=str(dt.datetime.now())[:19]
#f=open(''.join(('model_comp',strdate.replace(' ','_'),'.log')),'w')
######

########## set default Font sizes for plots ###################
Fincrease=0.6
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
####################################

#################################################################################################






############# Functions / Classes ###############################

def datetime64_to_datetime(t):
  if len(t)==1:
    t=[t]
  return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')) for ti in t])

def schism_time_to_datetime(t, baseTime):
  if len(t)==1:
    t=[t]
  return np.asarray([ dt.timedelta(seconds=ti) + baseTime for ti in t])		

		
########################### Program start #################################################################	

os.chdir(setupdir[0])

from schism import *

p=param()
reftime=dt.datetime(np.int(p.get_parameter('start_year')),\
np.int(p.get_parameter('start_month')),\
np.int(p.get_parameter('start_day')),\
np.int(p.get_parameter('start_hour')),0,0)

	

	

i_setup=0
s=schism_setup()
s.ds=schism_outputs_by_variable(ncdir[i_setup],max_stack=14).ds
s.dates=schism_time_to_datetime(s.ds['out2d'].time.data,reftime)
synonyms={'ssh':'elevation'} #actuall variable names in the model
x=np.asarray(s.lon)
y=np.asarray(s.lat)

nt=len(s.dates)


s.wetdry = s.ds['out2d']['dryFlagElement']

# limit date range later by s0.dates
for key in list(vartimes.keys())[::1]: #vartimes.keys():
	print(key)
	
	dsi=s.ds['out2d']['elevation'] # hardcoded now
	varname=key
	vmin,vmax=min_max[varname]
	print('ploting '+ varname)
	plt.close('all')

#	
#	######### INITIAL PLOT FOR VARIABLE #######################################
#	# data1 data2 run1 run2  run2-run1

	
	

	
	outdir2=outdir+varname+'/'
	if not os.path.exists(outdir2): os.mkdir(outdir2)


	ti=0
	slab=dsi[ti].values
	elemvalues=slab[s.nvplt].mean(axis=1)    
	wetdry=np.asarray(s.wetdry[ti],bool)
	elemvalues[wetdry]=np.nan
    
    
    
	fig, axes = plt.subplots(1,1)
	fig.set_dpi(300)
	plt.suptitle( str(s.dates[0]))
	ph,ch=s.plotAtelems(elemvalues,cmap=cm,mask=None)
	plt.tight_layout()
	ph.set_clim(vmin,vmax)
	
	
	# time loop with blitting
	
	for ti in range(nt):
		print(ti)
		slab=dsi[ti].values
		elemvalues=slab[s.nvplt].mean(axis=1)    
		wetdry=np.asarray(s.wetdry[ti],bool)        
		elemvalues[wetdry]=np.nan
		plt.suptitle( str(s.dates[ti]))
		ph.set_array(np.ma.masked_array(elemvalues,mask=wetdry)) #ignore mask
		fig.canvas.draw()
		fig.canvas.flush_events()
		#print('took '+str((dt.datetime.now()-t00).total_seconds())+' s')
		plt.savefig(outdir2+'{:04d}_intercomp'.format(ti)+'_'+varname,dpi=300)
		#print('inc write took '+str((dt.datetime.now()-t00).total_seconds())+' s')
		