"""
 append tracers to existing hotstart
 here only sediments with zero initialisation
"""

import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *a
from netCDF4 import Dataset

s=schism_setup()

# append tracer sediment all zero
# tracer concentrations on sides and elements will be interpolated
# tr_nd.shape has to be (nelements,nvrt,ntracers)



#!rm GB_hot_AMM1520170102nm+sediment.nc

ncfilein='GB_hot_AMM1520170102nm.nc'
ncfileout=ncfilein[:ncfilein.rindex('.')] + '+sediment.nc'

# read file
nc=Dataset(ncfilein)
old=nc['tr_nd'][:]

## add sediments
nclasses=8 # sediment add tracers with 0 as initial value
SED=np.zeros((old.shape[0],old.shape[1],nclasses))
tr_nd=np.concatenate((old,SED),axis=2)    # append traces interpolation to elements will be done by schism funxtion

s.znum=old.shape[1]
s.write_hotstart(tr_nd,filename=ncfileout)


# addd aditional dimensions and fields neeeded for sediments
nc0=Dataset(ncfileout,'a')

# additonal dimensions
dims={'two':2,'four':4,'five':5,'six':6,'seven':7,'eight':8,'nine':9,'SED_MBEDP':3,'SED_Nbed':1,
'SED_ntr':nclasses,'one_new':1}
for key in dims.keys():
	nc0.createDimension(key,dims[key]) # cannot alter sizesnc

	
# roughness not zero 
	
# additonal variables	
varnames={'SED3D_dp':('node'),'SED3D_rough':('node')} # varname;('dimensions')
values={'SED3D_dp':np.asarray(s.depths),'SED3D_rough':0} # varname;('dimensions')

for key in varnames.keys():
	newvar= nc0.createVariable(key,'float64',varnames[key])
	nc0[key][:]=values[key]

newvar= nc0.createVariable('SED3D_bed','float64',('elem','SED_Nbed','SED_MBEDP'))	# bed oncentration?
nc0['SED3D_bed'][:]=0

# working only with one bed layer
newvar= nc0.createVariable('SED3D_bedfrac','float64',('elem','SED_Nbed','SED_ntr')) # read from bedfrac ic
for i in range(dims['SED_ntr']):
	nodevalues=np.loadtxt('bed_frac_{:d}.ic'.format(i+1),skiprows=2,max_rows=s.nnodes)[:,-1]
	elemvalues=np.asarray([nodevalues[np.asarray(elem)-1].mean() for elem in s.nv])
	nc0['SED3D_bedfrac'][:,0,i]=elemvalues   #and average to element						

nc0.sync()
nc0.close()
	
# ntracers(10), one(1),  two(2), four(4), five(5), six(6), seven(7), eight(8), nine(9), SED_MBEDP(3), SED_Nbed(1), SED_ntr(8), one_new(1)
#variables(dimensions): int32 idry_e(elem), int32 idry_s(side), int32 idry(node), float64 eta2(node), float64 we(elem,nVert), float64 tr_el(elem,nVert,ntracers), float64 su2(side,nVert), float64 sv2(side,nVert), float64 tr_nd(node,nVert,ntracers), float64 tr_nd0(node,nVert,ntracers), float64 q2(node,nVert), float64 xl(node,nVert), float64 dfv(node,nVert), float64 dfh(node,nVert), float64 dfq1(node,nVert), float64 dfq2(node,nVert), f, float64 time(one_new), int32 iths(one_new), int32 ifile(one_new)

