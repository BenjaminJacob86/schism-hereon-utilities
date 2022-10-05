"""
Perform harmonic analysis over subselection of node time series
and save to netcdf
"""


# make new own enviroment with utide


# ok this creates conflicts with utide
##conda create --name schism --clone base # make new environment with all packages #instead create new with foloowing packages:
#conda create --name harmana ipython numpy netcdf4 xarray matplotlib dask
#conda install dask	# then also install dask for xarray
#conda install -c conda-forge utide	# then also install utide

##conda activate harmana
#conda install ipython
##conda install -c conda-forge utide  - install utide in environment
##conda activate harmana


#conda create -n harmana anaconda

# conda env remove -n schism


consts=['M2','S2','M4','M6']	


import os
import sys
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/') # add path to schism.py
from schism import *

import xarray as xr
import utide # load harmonic analysis module (needs install on levante)
from matplotlib.dates import date2num

rundir='/work/gg0028/g260181/copy/'
ncdir='/work/gg0028/g260181/copy/Combi/'

os.chdir(rundir)
s=schism_setup()
ds=schism_output2(ncdir).nc # get dartaset access


# surface layer tidal decomp
# number of horizontal grid nodes to be calculated for by this instance
nde0=0    #
nde1=100  #


varname='elev'    # variable to do harmonic analysis for
#varname='hvel'    # variable to do harmonic analysis for
			# vertical level to select - always surface
			

# select nodes to iterate above			
ds_sub=ds.sel( nSCHISM_hgrid_node=slice(nde0,nde1))			


# load  surface dataData			
latitudes=np.asarray(s.lat)[nde0:nde1]

timax=1440
time=ds_sub['time'][:timax].values
if varname=='elev': #2d			
	u=ds_sub['elev'].values
	v=None
elif varname=='hvel': # velosity vecotr
	hvel=ds_sub['hvel'][:timax,:,-1,:].values
	u=hvel[:,:,0]
	v=hvel[:,:,1]
else: # other 3d scalars 
	u=ds_sub[varname][:timax,:,-1].values
	v=None

time=date2num(time) #change time format for utide
	
# harmonic analyis
for inode in range(len(latitudes)):
	utout=utide.solve(time, u=u[:,inode], v=v[:,inode], lat=latitudes[inode])			
	harmx=utide.reconstruct(time,utout) # tidal reconstruction
	#uresidual=u[:,inode]-(harmx['u']-utout['umean'])
	#vresidual=u[:,inode]-(harmx['v']+utout['vmean'])
	uresidual=u[:,inode]-(harmx['u'])#-utout['umean'])
	vresidual=u[:,inode]-(harmx['v'])#+utout['vmean'])
	
	
utout['name'[	
g phase
'tehta' direction
'Lsmaj'  #major axis
'Lsmin'  # minor axis
'theta' direction


ind=np.where(utout.name==const)[0]
if len(ind)>0:
	Amps['TG'][const][i]=utout.A[ind]
	Phas['TG'][const][i]=utout.g[ind]


Amps=dict.fromkeys(sources)
Phas=dict.fromkeys(sources)

name=names[0]
