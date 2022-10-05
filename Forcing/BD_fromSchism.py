#!/usr/bin/env python
"""
Extract bd forcing for SCHISM from SCHISM
different SCHISM run using nn interpolation
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


# setup to take forcing from
rundir='/gpfs/work/kochw/schism-routine-BS/'
ncdir='/gpfs/work/kochw/schism-routine-BS/outputs/'


# Setup to be forced
#rundir1='//gpfs/work/jacobb/data/SETUPS/NWBlackSea/ShelfNoRiver/setup/'
#ncdir2=rundir+'combined/'
rundir1='/gpfs/work/jacobb/data/SETUPS/NWBlackSea/bscut/'

#ls /gpfs/work/kochw/schism-routine-BS/outputs/schout_20??????.nc
#sed 's/#nr#/1/g' interpolate_variables.in-template > test.file


#I only need four points for the periods: 2021.06.01~2021.09.30, with time interval 1200 s
#Time    (5968200, 379270)    (5968200, 384510)    (5954200, 384510)    (5954200, 379270)

# coordinates # here utm as schism grid
# x,y pairs
#coords=[(5968200, 379270),    (5968200, 384510),    (5954200, 384510),    (5954200, 379270)] #corner coordinates from xbeach

#pairs of x,y coordinates
#coords=[(5968200, 379270)[::-1],    (5968200, 384510)[::-1],    (5954200, 384510)[::-1],    (5954200, 379270)[::-1]] #corner coordinates from

########################################################################################

# load schism grid and netcdf acess
os.chdir(rundir)
s=schism_setup()
#s.ds=schism_output2(ncdir).nc # load handle
s.init_node_tree(latlon=True) # initialie ckd tree 

#schismfiles=[] 
schismfiles=np.sort(glob(ncdir+'schout_'+'????????.nc')) # here date named files
s.ds= xr.open_mfdataset(schismfiles)


###
os.chdir(rundir1)
s1=schism_setup()


bdnodes=np.asarray(s1.bdy_segments[0])-1
plt.plot(np.asarray(s1.lon)[bdnodes],np.asarray(s1.lat)[bdnodes],'r+')

lon=np.asarray(s1.lon)
lat=np.asarray(s1.lat)
coords=list(zip(lon[bdnodes],lat[bdnodes]))
Pnn=s.node_tree_latlon.query(coords)[1] 
s.dsnn=s.ds.sel(nSCHISM_hgrid_node=Pnn) ## sele nearest neighbours within data set





plt.figure()
s.plot_domain_boundaries(latlon=True)
s.plot_mesh()
s1.plot_domain_boundaries(latlon=True,append=True)
plt.plot(lon[bdnodes],lat[bdnodes],'ko',markersize=15)
plt.plot(np.asarray(s.lon)[Pnn],np.asarray(s.lat)[Pnn],'ro',markersize=10)







#s.write_bdy_netcdf('elev2D'+timetag+'.th.nc',timesq,bdprofiles2['ssh'],frcbdnodes=frcbdnodes)



with open('bctides.in') as f:
	openbd_segs=[]
	bd=0
	for line in f.readlines()[4:]:
		print(line)
		splitted=line.split()
		bd+=np.sum([val.isdigit() for val in splitted[:5]])==5 #is openboundary 
		if (splitted[1:5]==['4', '4', '4', '4']) | (splitted[1:5]==['5', '5', '5', '5']):
			openbd_segs.append(bd-1)
print('create forcing for open boundaries '+str(openbd_segs))	

if len(openbd_segs)>0:
	frcbdnodes=[]
	for seg in openbd_segs:
		frcbdnodes+=s1.bdy_segments[seg]
		bdyvgrid = np.asarray([s1.vgrid[ii].filled(-1.) for ii in frcbdnodes ])
else:
	frcbdnodes=s1.bdy_nodes




#/gpfs/home/jacobb/git/schism-master/schismpure_evap/bin/interpolate_variables7
#/gpfs/home/jacobb/git/schism-master/src/Utility/OneWayNestScripts


bddepthsIn=np.asarray([s.vgrid[inode].filled(-1)*s.depthsdict[inode] for inode in Pnn+1])
bddepths=np.asarray([s1.vgrid[inode].filled(-1)*s1.depthsdict[inode] for inode in frcbdnodes])

#plt.figure()
#plt.plot(bddepthsIn[:,0])
#plt.plot(bddepths[:,0])
#
#plt.clf()
#plt.plot(bddepthsIn[:,:],'b')
#plt.plot(bddepths[:,:],'r')
#
#
#
#plt.plot(bddepthsIn[:,:],'.-')
#plt.plot(bddepths[:,:],'x-')
#

ibd=np.hstack([np.asarray(s1.bdy_segments[bdseg])-1 for bdseg in openbd_segs])		
d=-np.asarray(s1.depths)[ibd] # ant positive (here negative later multiplired with sigma)) values

#ibdin=np.hstack([np.asarray(s.bdy_segments[bdseg])-1 for bdseg in [np.asarray(Pnn,int)+1,]])		
ibdin=Pnn
din=-np.asarray(s.depths)[ibdin] # ant positive (here negative later multiplired with sigma)) values

xin=np.asarray(s.x)[ibdin]
yin=np.asarray(s.y)[ibdin]

xout=np.asarray(s1.x)[ibd]
yout=np.asarray(s1.y)[ibd]

nnodes,nz=bddepths.shape


t0=s.ds['time'][0].values
t1=s.ds['time'][-1].values
timetag=str(t0)[:10]+'-'+str(t1)[:10]

tout=(s.ds['time'].values-s.ds['time'][0].values)/np.timedelta64(1,'s')
tout+=tout[1]
tout=np.hstack((0,tout))

# output time step
#dt_output=tout[:2].diff()
dt_output=np.diff(tout[:2])[0]
deltaT=s.ds['time'][1]-s.ds['time'][0]

# load input
elev=s.dsnn['elev'].values
elev=np.vstack((elev[0,:],elev))
elev=elev.reshape(elev.shape+(1,1))

salt=s.dsnn['salt'].values
salt=np.vstack((salt[0,:].reshape((1,)+salt.shape[1:]),salt))
salt=salt.reshape(salt.shape+(1,))

temp=s.dsnn['temp'].values
temp=np.vstack((temp[0,:].reshape((1,)+temp.shape[1:]),temp))
temp=temp.reshape(temp.shape+(1,))

hvel=s.dsnn['hvel'].values
hvel=np.vstack((hvel[0,:].reshape((1,)+hvel.shape[1:]),hvel))



tempout=np.ones(temp.shape[:2]+(nz,1))*np.nan
saltout=np.ones(temp.shape[:2]+(nz,1))*np.nan
hvelout=np.ones(temp.shape[:2]+(nz,2))*np.nan
nt=temp.shape[0]
timesq=[ti*dt_output for ti in range(nt)]



# sort out nans in input
for inode in range(nnodes):
	inan=np.isnan(temp[0,inode,:,0])
	valid_start=inan.sum()
	temp[:,inode,inan,0]=np.tile(temp[:,inode,valid_start,0],(valid_start,1)).swapaxes(0,1)
	salt[:,inode,inan,0]=np.tile(salt[:,inode,valid_start,0],(valid_start,1)).swapaxes(0,1)
	hvel[:,inode,inan,:]=np.tile(hvel[:,inode,valid_start,:],(valid_start,1,1)).swapaxes(0,1)

	



#imaxdepin=len(depthin)-1	
for inode in range(nnodes):


	#dep=s.vgrid[ibd[inode]-1].filled(-1)*d[inode]
	dep=s1.vgrid[ibd[inode]+1].filled(-1)*d[inode]
	depthin=s.vgrid[ibdin[inode]].filled(-1)*din[inode]
	print('depth interp profiles node {:d}/{:d}'.format(inode,nnodes))
	for vind,di in enumerate(dep):	
		#iabove =  (di <= depthin).sum()-1   # schism starts from bottom
		ibelow =  (di <= depthin).sum()-1   # schism starts from bottom
		iabove=ibelow+1
		iabove-=(ibelow+1==len(depthin))   # last index len -1
		
		#ibelow=np.minimum(imaxdepin,ibelow) # extrapolate if not deep enough
		
		dztop=np.abs(di-depthin[iabove]) 
		dzbotm=depthin[ibelow]-di
		
		if dzbotm==dztop==0:
			wtop=wbotm=0.5
		elif	dztop==0:
			wtop=1
			wbotm=0
		else:	
			wtop=1-dztop/(dztop+dzbotm)  
			wbotm=1-wtop
		
		
		
		saltout[:,inode,vind,0]=salt[:,inode,ibelow,0]*wbotm+salt[:,inode,iabove,0]*wtop
		
		if np.nanmax(saltout) > 40:
			break
	if np.nanmax(saltout) > 40:
		break	
		
		tempout[:,inode,vind,0]=temp[:,inode,ibelow,0]*wbotm+temp[:,inode,iabove,0]*wtop
		
		hvelout[:,inode,vind,0]=hvel[:,inode,ibelow,0]*wbotm+hvel[:,inode,iabove,0]*wtop
		hvelout[:,inode,vind,1]=hvel[:,inode,ibelow,1]*wbotm+hvel[:,inode,iabove,1]*wtop



s1.write_bdy_netcdf('elev2D'+timetag+'.th.nc',tout,elev,frcbdnodes=frcbdnodes)
s1.write_bdy_netcdf('SAL_3D'+timetag+'.th.nc',tout,saltout,frcbdnodes=frcbdnodes)
s1.write_bdy_netcdf('TEM_3D'+timetag+'.th.nc',tout,tempout,frcbdnodes=frcbdnodes)
#hvel=hvel.reshape(hvel.shape+(1,))
s1.write_bdy_netcdf('uv3D'+timetag+'.th.nc',tout,hvelout,frcbdnodes=frcbdnodes)
# interpolate vertically to grid

try:
	os.symlink('elev2D'+timetag+'.th.nc','elev2D.th.nc')
	os.symlink('SAL_3D'+timetag+'.th.nc','SAL_3D.th.nc')
	os.symlink('TEM_3D'+timetag+'.th.nc','TEM_3D.th.nc')
	os.symlink('uv3D'+timetag+'.th.nc','uv3D.th.nc')
except:
		print('links already set')


bdprofiles2={'elev':elev,'salt':saltout,'temp':tempout,'uv':hvelout}
plot_bd=True
if plot_bd:
	plt.figure()
	#names={'ssh':name_ssh,'salt':name_salt,'temp':name_temp}

	# control plots for
	istart=0
	iend=0
	#plt.figure()
	plt.clf()
	for i in openbd_segs: #range(len(s.bdy_segments)):
		iend+=len(s1.bdy_segments[i])
		inds=np.asarray(range(istart,iend))
		istart=iend
		bdnodes=np.asarray(frcbdnodes)[inds]#-1	

		bddepths=np.asarray([s1.vgrid[inode].filled(-1)*s1.depthsdict[inode] for inode in bdnodes])
		xcoord=np.tile(inds,[len(bddepths[0]),1]).T

		label=['start','end']
		for ti in range(0,-2,-1):
			plt.clf()
			plt.subplot(2,2,1)
			plt.pcolor( xcoord,bddepths,bdprofiles2['uv'][ti,inds,:,0])
			ch=plt.colorbar()
			plt.xlabel('bdy node')
			plt.ylabel('bddepths')
			ch.set_label('u [m/s]')
		
			plt.subplot(2,2,2)
			plt.pcolor( xcoord,bddepths,bdprofiles2['uv'][ti,inds,:,1])
			ch=plt.colorbar()
			plt.xlabel('bdy node')
			plt.ylabel('bddepths')
			ch.set_label('v [m/s]')

			plt.subplot(2,2,3)
			plt.pcolor( xcoord,bddepths,bdprofiles2['temp'][ti,inds,:,0])
			ch=plt.colorbar()
			plt.xlabel('bdy node')
			plt.ylabel('bddepths')
			ch.set_label('t [deg C]')
			
			plt.subplot(2,2,4)
			plt.pcolor( xcoord,bddepths,bdprofiles2['salt'][ti,inds,:,0])
			ch=plt.colorbar()
			plt.xlabel('bdy node')
			plt.ylabel('bddepths')
			ch.set_label('s [g/kg]')
			plt.suptitle(str(t0+ti*deltaT))	
			plt.tight_layout()
			
			plt.savefig('bd_%i_%s'%(i,label[ti]),dpi=400)
			#plt.show()
			plt.close()
