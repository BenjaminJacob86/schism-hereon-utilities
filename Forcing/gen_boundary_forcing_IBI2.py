""" Generate SCHISM boundary forcing from gcoast or amm15 data"""

# eventually failure. return order for weight smight be needing to be  w1 * B[:,0,0] +w2*B[:,1,0] + w3*B[:,0,1] + w4 * B[:,1,1]    
# but this was testet  for one node dim. then it is correct.

# add IBI
# for IBI u and V on diffrent grid

# /work/gg0028/g260098/VERIBLUE/IBI_NAT_PHY/nrt
# R20230417 #start


### !! Dont use gcoast, since Sebastian saves the Depth in a seprate file the depth reading
#### From a different file and the depth in the individual files is mean depth on not related to the depth at the nodes
import numpy as np
from scipy.interpolate import griddata
import os
import sys
import pandas as pd


#Timestamp('2017-11-08 17:50:00')

sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities')
sys.path.insert(0,'/work/gg0028/g260114/RUNS/GermanBight/GB_HR_Ballje/')
sys.path.insert(0,'/gpfs/work/jacobb/data/RUNS/GB_template/schism-hzg-utilities/')

import matplotlib
#matplotlib.use('AGG')
from schism import *
from scipy.spatial import cKDTree
import xarray as xr
from glob import glob
import datetime as dt
from scipy.interpolate import interp1d
import time
plt.ion()

tstart=time.time()
# before enerting python and execution call export OMP_NUM_THREADS=1

# incorporate
#Do your spatial and temporal indexing (e.g. .sel() or .isel()) early in the pipeline, especially before calling resample() or groupby(). Grouping and rasampling triggers some computation on all the blocks, which in theory should commute with indexing, but this optimization hasn’t been implemented in dask yet. (See dask issue #746).
#Save intermediate results to disk as a netCDF files (using to_netcdf()) and then load them again with open_dataset() for further computations. For example, if subtracting temporal mean from a dataset, save the temporal mean to disk before subtracting. Again, in theory, dask should be able to do the computation in a streaming fashion, but in practice this is a fail case for the dask scheduler, because it tries to keep every chunk of an array that it computes in memory. (See dask issue #874)
#Specify smaller chunks across space when using open_mfdataset() (e.g., chunks={'latitude': 10, 'longitude': 10}). This makes spatial subsetting easier, because there’s no risk you will load chunks of data referring to different chunks (probably not necessary if you follow suggestion 1).
#

# Requires that depth in forcing model is from top to surface (this is the case for amm15 and gcoast)

# TODOS:
# !!!! distances and weights for bilinear interpolation are calculated in lon/lat
#      To make it more accurate transition into cartesian framework would be good


####### S E T T I N G S ##############################################
print('starting program')
# before program start call  'export OMP_NUM_THREADS=4' to stabalize xarrray # check if 'export OMP_NUM_THREADS=1' faster

# time
#t0=dt.datetime(2018,1,2)
#t1=dt.datetime(2018,12,31)
#dt_output=3600    # timestep for output [seconds]


year='2023'						# to limit IBI files
t0=dt.datetime(2023,1,1)		#  time start WK
t1=dt.datetime(2023,2,1)		#  time end WK
#dt_output=86400    # timestep for output [seconds] take as input
dt_output=3600         # if zero take same temporal resolution as input data
					# tine is counted from t0 on

# forcing source  amm15 gcoast.
#schismdir='/gpfs/work/jacobb/data/SETUPS/SNS_Wei/' # setup directory
schismdir=os.getcwd()  # local
#schismdir=schismdir[:schismdir.rindex('/')]+'/' #assume 
#schismdir=rundir='/gpfs/work/jacobb/data/SETUPS/GB_template/'

frocingtype='IBI'#'gcoast' # IBI # CMEMS
frcdir='/gpfs/work/jacobb/data/routines/blacksea_routine/Download/download_cmems_BS/medsea/'
frcdir='/work/gg0028/g260098/VERIBLUE/IBI_NAT_PHY/nrt/'
plot_bd=True # make control plots of boundaries

#openbd_segs=[0,8]  # open boundary segments Determine from bctides in boundary configuration corresponding to
					# 4 4 4 4 or 5 5 5 5	no need for manual specification anymore


# varnames - edit this sections and the suffixes 
# in the access handles part to configure for additionl models
name_lon='lon'
name_lat='lat'
name_time='time'
name_depth='depth'
if frocingtype=='gcoast':
	name_ssh='sossheig'
	name_salt='vosaline'
	name_temp='votemper'
	name_u='vozocrtx'
	name_v='vomecrty'
	ssh_pattern=''
	ssh_files=''
	temp_files=''
	uv_pattern=''

elif (frocingtype=='amm15') or (frocingtype=='cmems')   or (frocingtype=='IBI'):
	name_ssh='zos'
	name_salt='so'
	name_temp='thetao'
	name_u='uo'
	name_v='vo'
	ssh_pattern=''
	ssh_files=''
	temp_files=''
	uv_pattern=''

if frocingtype=='cmems':
	ssh_pattern='ASLV'
	temp_pattern='TEP'
	salt_pattern='PSAL'
	uv_pattern='RFVL'	
########################################################################

	

######################## Functions #############
def datetime64_to_datetime(t):
	if len(t)==1:
		t=[t]
	return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00Z'))/ np.timedelta64(1, 's')) for ti in t])
	
def get_4neighbours(x2d,y2d,bdx,bdy):
	""" get 4 neighbours for bilinear interlation for each node on schism boundary """
	bdcoords=np.asarray(list(zip(bdx,bdy)))
	tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))	
	dd,nn=tree2d.query(bdcoords,k=1)
	ii0,jj0=np.unravel_index(nn,lat2d.shape)

	indsx=[]
	indsy=[]

		
	for i in range(len(ii0)):
		xgt0=(x2d[ii0[i],jj0[i]]-bdx[i])>0
		ygt0=(y2d[ii0[i],jj0[i]]-bdy[i])>0
		
		if xgt0 & ygt0:  # UR
			#indx=[[1, 2],[3, 4]]
			indx=[[ii0[i]-1, ii0[i]],[ii0[i]-1, ii0[i]]]
			indy=[[jj0[i]-1, jj0[i]-1],[jj0[i], jj0[i]]]
		elif (xgt0==True) & (ygt0==False):  # LR
			indx=[[ii0[i], ii0[i]+1],[ii0[i], ii0[i]+1]]
			indy=[[jj0[i]-1, jj0[i]-1],[jj0[i], jj0[i]]]
		
		elif (xgt0==False) & (ygt0==True):  # UL
			indx=[[ii0[i]-1, ii0[i]],[ii0[i]-1, ii0[i]]]
			indy=[[jj0[i], jj0[i]],[jj0[i]+1, jj0[i]+1]]
		
		elif (xgt0==False) & (ygt0==False):  # LL
			indx=[[ii0[i], ii0[i]+1],[ii0[i], ii0[i]+1]]
			indy=[[jj0[i], jj0[i]],[jj0[i]+1, jj0[i]+1]]

		indsx.append(np.asarray(indx).T)	
		indsy.append(np.asarray(indy).T)	
		
		
		#check
		if (np.sign(x2d[np.asarray(indx).T,np.asarray(indy).T]-bdlon[i])==np.asarray([[-1,1],[-1,1]]))[0].min() & (np.sign(y2d[np.asarray(indx).T,np.asarray(indy).T]-bdlat[i])==np.asarray([[-1,-1],[1,1]]))[0].min():
			pass
		else:
			print(i)
			break
			
	#increase performance
	indsx=np.asarray(indsx)
	indsy=np.asarray(indsy)
		
		
	return indsx,indsy

	
def get_4neighbours(x2d,y2d,bdx,bdy):
	""" get 4 neighbours for bilinear interlation for each node on schism boundary """
	bdcoords=np.asarray(list(zip(bdx,bdy)))
	tree2d=cKDTree(list(zip(x2d.flatten(),y2d.flatten())))	
	dd,nn=tree2d.query(bdcoords,k=1)
	ii0,jj0=np.unravel_index(nn,y2d.shape)

	indsx=[]
	indsy=[]

	print(len(ii0))	
	for i in range(len(ii0)):
		xgt0=(x2d[ii0[i],jj0[i]]-bdx[i])>0
		ygt0=(y2d[ii0[i],jj0[i]]-bdy[i])>0
		
		if xgt0 & ygt0:  # UR
			#indx=[[1, 2],[3, 4]]
			indx=[[ii0[i]-1, ii0[i]],[ii0[i]-1, ii0[i]]]
			indy=[[jj0[i]-1, jj0[i]-1],[jj0[i], jj0[i]]]
		elif (xgt0==True) & (ygt0==False):  # LR
			indx=[[ii0[i], ii0[i]+1],[ii0[i], ii0[i]+1]]
			indy=[[jj0[i]-1, jj0[i]-1],[jj0[i], jj0[i]]]
		
		elif (xgt0==False) & (ygt0==True):  # UL
			indx=[[ii0[i]-1, ii0[i]],[ii0[i]-1, ii0[i]]]
			indy=[[jj0[i], jj0[i]],[jj0[i]+1, jj0[i]+1]]
		
		elif (xgt0==False) & (ygt0==False):  # LL
			indx=[[ii0[i], ii0[i]+1],[ii0[i], ii0[i]+1]]
			indy=[[jj0[i], jj0[i]],[jj0[i]+1, jj0[i]+1]]

		indsx.append(np.asarray(indx).T)	
		indsy.append(np.asarray(indy).T)	
		
		
		#check
		if (np.sign(x2d[np.asarray(indx).T,np.asarray(indy).T]-bdlon[i])==np.asarray([[-1,1],[-1,1]]))[0].min() & (np.sign(y2d[np.asarray(indx).T,np.asarray(indy).T]-bdlat[i])==np.asarray([[-1,-1],[1,1]]))[0].min():
			pass
		else:
			print(i)
			#from IPython import embed; embed()
			pass	
			#break
			
	#increase performance
	indsx=np.asarray(indsx)
	indsy=np.asarray(indsy)
		
		
	return indsx,indsy	
	
# das gibt richtiges	
def bilin3(x,y,v,indsx,indsy,bdx,bdy,calcw=True,w1=None,w2=None,w3=None,w4=None):
	""" bilinear interpolation """
	#interped=[]
	npbd=len(bdx)
	interped=np.zeros(npbd)#[]
	
	if (calcw==True):
		dtype=np.float64
		div_dxdy=np.zeros(npbd,dtype)#[]
		dxi=np.zeros((npbd,2),dtype)#[]
		dyi=np.zeros((npbd,2),dtype)#[]
		for i in range(len(indsx)):
			indx=indsx[i]
			indy=indsy[i]

			div_dxdy[i]=(1/((x[indx[-1,-1],indy[-1,1]]-x[indx[0,0],indy[0,0]])*(y[indx[-1,-1],indy[-1,1]]-y[indx[0,0],indy[0,0]])))

			dxi[i]=([x[indx[-1,-1],indy[-1,1]]-bdx[i],bdx[i]-x[indx[0,0],indy[0,0]]])
			
			dyi[i]=(np.asarray(([y[indx[-1,-1],indy[-1,-1]]-bdy[i],bdy[i]-y[indx[0,0],indy[0,0]]])))
			
			interped[i]=(div_dxdy[i]*np.matmul(dyi[i],np.matmul(v[indx,indy],dxi[i])))	
			
			
		# rewire such that no transpose needed
		
		# wrong order

		# above weights are shifte from what I thougt	
		#control
		#indx,indy=indsx[0],indsy[0]
		#try making float 64 yes incrases accuracy to -15, however matmul above manages exact agreement, model error shuld anyways be larger than float32 accuracy
		#x=np.asarray(x,np.float64)
		#y=np.asarray(y,np.float64)
		#v=np.asarray(v,np.float64)
		
		# check bilin interp 
		#dx=np.abs((x[indx,indy]-bdx[0]))
		#dxs=dx.sum(axis=1)
		#wx=1-dx/dxs
		#(wx*x[indx,indy]).sum(axis=1)[0]-bdx[0]
		#dy=np.abs((y[indx,indy]-bdy[0]))
		#dys=dy.sum(axis=0)
		#wy=1-dy/dys
		#W=wx*wy
		#(W*(x[indx,indy])).sum() -bdx[0]
		#(W*(x[indx,indy])).sum() -a
		
		
		#W1 W2   +_+
		#W3 W4   + +
		
		# shift calc results o above
		w1=div_dxdy* dxi[:,0]*dyi[:,0] 
		w2=div_dxdy*dxi[:,1]*dyi[:,0]
		w3=div_dxdy*dxi[:,0]*dyi[:,1]
		w4=div_dxdy*dxi[:,1]*dyi[:,1]

		#q=v[indx,indy]
		#a=w1[i] * q[0,0] + w2[i] * q[0,1] +w3[i] *q[1,0] +w4[i]*q[1,1]
#nterped[i]=(div_dxdy[i]*np.matmul(dyi[i],np.matmul(v[indx,indy],dxi[i])))	
		return interped,w1,w2,w3,w4
	else: #use precalculated
		#return w1*v[indsx,indsy][:,0,1]+w2*v[indsx,indsy][:,1,0] + w3 * v[indsx,indsy][:,0,0] + w4 * v[indsx,indsy][:,1,1]
		return w1*v[indsx,indsy][:,0,0]+w2*v[indsx,indsy][:,0,1] + w3 * v[indsx,indsy][:,1,0] + w4 * v[indsx,indsy][:,1,1]
		
		
#def bilin3UseWeights(x,y,v,indsx,indsy,bdx,bdy,w1=None,w2=None,w3=None,w4=None):
#	""" bilinear interpolation """
#	return w1*v[indsx,indsy][:,0,0]+w2*v[indsx,indsy][:,0,1] + w3 * v[indsx,indsy][:,1,0] + w4 * v[indsx,indsy][:,1,1]
#
# faster to do the subindexing once and save
def bilin3UseWeights(x,y,v,indsx,indsy,bdx,bdy,w1=None,w2=None,w3=None,w4=None):
	""" bilinear interpolation """
	q=v[indsx,indsy]
	return w1*q[:,0,0]+w2*q[:,0,1] + w3 * q[:,1,0] + w4 * q[:,1,1]

# fastest version so far:
def bilin3UseWeights_stack(x,y,v,indsx,indsy,Wstack):
	""" bilinear interpolation """
	q=v[indsx,indsy]
	#return np.multiply(q,Wstack).sum(axis=(1,2))
	return np.einsum('ijk,ijk->i', q, Wstack) # encode above in einstein summation	
	

# depth extraploate weights from  bilin3	
## try speeding up
## do over all depths
#def bilin3depth(x,y,v,indsx,indsy,bdx,bdy,calcw=True,w1=None,w2=None,w3=None,w4=None):
#	""" bilinear interpolation """
#	interped=[]
#	if (calcw==True):
#		div_dxdy=[]
#		dxi=[]
#		dyi=[]
#		for i in range(len(indsx)):
#			indx=indsx[i]
#			indy=indsy[i]
#			div_dxdy.append(1/((x[indx[-1,-1],indy[-1,1]]-x[indx[0,0],indy[0,0]])*(y[indx[-1,-1],indy[-1,1]]-y[indx[0,0],indy[0,0]])))
#			dxi.append([x[indx[-1,-1],indy[-1,1]]-bdx[i],bdx[i]-x[indx[0,0],indy[0,0]]])
#			dyi.append(np.asarray(([y[indx[-1,-1],indy[-1,-1]]-bdy[i],bdy[i]-y[indx[0,0],indy[0,0]]])))
#			interped.append(div_dxdy[i]*np.matmul(dxi[i],np.matmul(v[indx,indy].T,dyi[i])))	
#			# rewire such that no transpose needed
#		dxi=np.asarray(dxi)
#		dyi=np.asarray(dyi)
#
#		# old
#		div_dxdy=np.asarray(div_dxdy)
#		w1=div_dxdy*dxi[:,0]*dyi[:,0] 
#		w2=div_dxdy*dxi[:,0]*dyi[:,1]
#		w3=div_dxdy*dxi[:,1]*dyi[:,0]
#		w4=div_dxdy*dxi[:,1]*dyi[:,1]
#		
#		return interped,w1,w2,w3,w4
#	else: #use precalculated
#		return w1*v[:,indsx,indsy][:,:,0,1]+w2*v[:,indsx,indsy][:,:,1,0] + w3 * v[:,indsx,indsy][:,:,0,0] + w4 * v[:,indsx,indsy][:,:,1,1]	


# try speeding up
# do over all depths
def bilin3depth(x,y,v,indsx,indsy,bdx,bdy,calcw=True,w1=None,w2=None,w3=None,w4=None):
	""" bilinear interpolation """
	interped=[]
	if (calcw==True):
		div_dxdy=[]
		dxi=[]
		dyi=[]
		for i in range(len(indsx)):
			indx=indsx[i]
			indy=indsy[i]
			div_dxdy.append(1/((x[indx[-1,-1],indy[-1,1]]-x[indx[0,0],indy[0,0]])*(y[indx[-1,-1],indy[-1,1]]-y[indx[0,0],indy[0,0]])))
			dxi.append([x[indx[-1,-1],indy[-1,1]]-bdx[i],bdx[i]-x[indx[0,0],indy[0,0]]])
			dyi.append(np.asarray(([y[indx[-1,-1],indy[-1,-1]]-bdy[i],bdy[i]-y[indx[0,0],indy[0,0]]])))
			interped.append(div_dxdy[i]*np.matmul(dxi[i],np.matmul(v[indx,indy].T,dyi[i])))	
			# rewire such that no transpose needed
		dxi=np.asarray(dxi)
		dyi=np.asarray(dyi)

		# old
		div_dxdy=np.asarray(div_dxdy)
		w1=div_dxdy*dxi[:,0]*dyi[:,0] 
		w2=div_dxdy*dxi[:,0]*dyi[:,1]
		w3=div_dxdy*dxi[:,1]*dyi[:,0]
		w4=div_dxdy*dxi[:,1]*dyi[:,1]
		
		return interped,w1,w2,w3,w4
	else: #use precalculated
		return w1*v[:,indsx,indsy][:,:,0,1]+w2*v[:,indsx,indsy][:,:,1,0] + w3 * v[:,indsx,indsy][:,:,0,0] + w4 * v[:,indsx,indsy][:,:,1,1]	
	
## pure calculation function avoiding if check		
#def bilin3depthUseWeights(x,y,v,indsx,indsy,bdx,bdy,w1=None,w2=None,w3=None,w4=None):
#	""" bilinear interpolation """
#	return w1*v[:,indsx,indsy][:,:,0,1]+w2*v[:,indsx,indsy][:,:,1,0] + w3 * v[:,indsx,indsy][:,:,0,0] + w4 * v[:,indsx,indsy][:,:,1,1]			
#		
# pure calculation function avoiding if check		
def bilin3depthUseWeights(x,y,v,indsx,indsy,bdx,bdy,w1=None,w2=None,w3=None,w4=None):
	""" bilinear interpolation """
	q=v[:,indsx,indsy]
	return w1*q[:,:,0,0]+w2*q[:,:,0,1] + w3 * q[:,:,1,0] + w4 * q[:,:,1,1]					

# in 2D case einstein summation seems faster, however not for depth, maybe to many implicit loops with indices	
# pure calculation function avoiding if check		
def bilin3depthUseWeights_stack(x,y,v,indsx,indsy,WstackD):
	""" bilinear interpolation """
	q=v[:,indsx,indsy]
	return np.einsum('ijkl,ijkl->ij', q, WstackD) # encode above in einstein summation	
#############################################################

def get_file_start_dates(files):
	start_dates=[]
	for i,file in enumerate(files):
		if i==0:
			b4date='h-m_'
			i0=file.index(b4date,1)+4 #where datums tarts in string
			i1=i0+8
			startdate=file[i0:i1]
		start_dates.append(pd.to_datetime(file[i0:i1], format='%Y%m%d')) #%H%M
		#np.datetime64(file[i0:i1]))
	return(start_dates)
############## actual Program start ############################################


############## load schism info ###############
os.chdir(schismdir)
s=schism_setup()
lon=np.asarray(s.lon)
lat=np.asarray(s.lat)


# detect at which boundaries forcing has to be applied
with open('bctides.in') as f:
	openbd_segs=[]
	bd=0
	for line in f.readlines()[4:]:
		print(line)
		splitted=line.split()
		bd+=np.sum([val.isdigit() for val in splitted[:5]])==5 #is openboundary 
		if (splitted[1:5]==['4', '4', '4', '4']) | (splitted[1:5]==['5', '5', '5', '5']) | (splitted[1:5]==['5', '5', '4', '4']):
			openbd_segs.append(bd-1)
print('create forcing for open boundaries '+str(openbd_segs))	
ibd=np.hstack([np.asarray(s.bdy_segments[bdseg])-1 for bdseg in openbd_segs])		

d=-np.asarray(s.depths)[ibd] # ant positive (here negative later multiplired with sigma)) values to allow mor effective interpolation
bdlon=lon[ibd]
bdlat=lat[ibd]
depths=np.asarray(s.depths)
bdcoords=np.asarray(list(zip(bdlon,bdlat)))


timetag=(str(t0)+'_'+str(t1)).replace(' ','_').replace('-','')
day_prev_year=t0-dt.timedelta(days=1)
day_follow_year=t1+dt.timedelta(days=1)

############# End Load Schism info ##############################


############## access handles fo gcoast data ############################
## derive from date directories, get filename patterns of variables files
# checking for specified variable names ('sossheig' etc. )


if frocingtype=='gcoast':
	folder1=frcdir+day_prev_year.strftime('%Y%m%d')   # last day previous year
	folder2=frcdir+day_prev_year.strftime('%Y')+'????'   # year
	folder3=frcdir+day_follow_year.strftime('%Y%m%d') # 1s day next year
	folders= np.unique([folder1] + glob(folder2) +[folder3])
	date_prefix='_'#
	suffix='_'
	files=glob(folders[0]+'/*.nc')
elif frocingtype=='gcoast': # all files in one folder	 # Format NWS from Sebastian differen for other cmems
	folders=[frcdir,]
	date_prefix='_hi'# string to search date name fpr
	suffix='_b'
	files=glob(folders[0]+'/*.nc')
elif frocingtype=='IBI':
	folders=[frcdir,]
	files=glob(folders[0]+'*/*.nc')	
	date_prefix='_'#
	suffix='_'	
	wild_elev_files='IBI36NWS_ng_1h-m_*_2DT-oce_*.nc' # * varying infprmation	
	wild_veloU_files='IBI36NWS_ng_25h-m_*_3DU-uo_*.nc'	
	wild_veloV_files='IBI36NWS_ng_25h-m_*_3DV-vo_*.nc'	
	wild_temp_files='IBI36NWS_ng_25h-m_*_3DT-thetao_*.nc'	
	wild_salt_files='IBI36NWS_ng_25h-m_*_3DT-so_*.nc'	
	
else: # all files in one folder	 # Format NWS from Sebastian differen for other cmems
	folders=[frcdir,] # default from cmems ftp download for black sea
	date_prefix='/'# string to search date name fpr starting with date
	suffix='_'
	files=glob(folders[0]+'/*.nc')
	

	
# determine file name pattern for variable files
if frocingtype=='IBI':

	

	ssh_files=np.sort(glob(frcdir+'*/'+wild_elev_files))
	i0_i1=np.where([ year in file for file in ssh_files ])[0][[0,-1]]
	ssh_files=ssh_files[i0_i1[0]-1:i0_i1[1]+2]
	
	salt_files=np.sort(glob(frcdir+'*/'+wild_salt_files))
	i0_i1=np.where([ year in file for file in salt_files ])[0][[0,-1]]
	salt_files=salt_files[i0_i1[0]-1:i0_i1[1]+2]
	
	temp_files=np.sort(glob(frcdir+'*/'+wild_temp_files))	
	i0_i1=np.where([ year in file for file in temp_files ])[0][[0,-1]]
	temp_files=temp_files[i0_i1[0]-1:i0_i1[1]+2]

	u_files=np.sort(glob(frcdir+'*/'+wild_veloU_files))	
	i0_i1=np.where([ year in file for file in u_files ])[0][[0,-1]]
	u_files=u_files[i0_i1[0]-1:i0_i1[1]+2]

	v_files=np.sort(glob(frcdir+'*/'+wild_veloV_files))	
	i0_i1=np.where([ year in file for file in v_files ])[0][[0,-1]]
	v_files=v_files[i0_i1[0]-1:i0_i1[1]+2]

	
	
	
	
	
else:

	try:
		if frocingtype == 'cmems':
			dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex(date_prefix)+len(date_prefix)+file[file.rindex(date_prefix):].index(suffix)-1]) for file in files])
			def date_from_name(files):
				return np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex(date_prefix)+len(date_prefix)+file[file.rindex(date_prefix):].index(suffix)-1]) for file in files])
		else: #amm 15	
			dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in files])
			def date_from_name(files):
				returnnp.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in files])
		
	except:	
		date_prefix='_' # if cmems download from motu and not NWS
		suffix='_'
		dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in files])
		def date_from_name(files):
			return np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in files])
	isort=np.argsort(dates)





	if False:
		## presroting
		# working for gcoast
		ssh_dates,salt_dates,temp_dates,uv_dates=[],[],[],[]
		checkdates=isort[dates[isort]==dates[isort[0]]] # added for autocheck different file names
		for i in isort[dates[isort]==dates[isort[0]]]:
			dsi=xr.open_dataset(files[i])
			alt=checkdates[checkdates!=i][0]
			#f=np.asarray(files[i][files[i].rindex('/')+1:],str)
			#f2=np.asarray(files[alt][files[alt].rindex('/')+1:],str)
			f=files[i][files[i].rindex('/')+1:]
			f2=files[alt][files[alt].rindex('/')+1:]
			pattern=''
			for ii in np.where([ (f[j] != f2[j]) or (f[j+1] != f2[j+1]) for j in range(len(f)-1)])[0]:
				pattern+=f[ii]
			pattern=pattern[1:]
			for iname, name in enumerate([name_ssh,name_salt,name_temp,name_u]):
				iname,name
				if name in list(dsi.variables.keys()):
					if iname==0: # ssh
						np.sort(files)[isort][:4]
						ssh_pattern=pattern#files[i][files[i].rindex('/')+1:files[i].rindex(suffix)]
						ssh_files=np.sort([glob(folder+'/*'+ssh_pattern+'*') for folder in folders])[0]
						ssh_dates.append(dates[isort[0]])
					elif iname==1: # salt/temp
						salt_pattern=pattern#files[i][files[i].rindex('/')+1:files[i].rindex(suffix)]
						salt_files=np.sort([glob(folder+'/*'+salt_pattern+'*') for folder in folders])[0]
						salt_dates.append(dates[isort[0]])
					elif iname==2: # salt/temp
						temp_pattern=pattern#files[i][files[i].rindex('/')+1:files[i].rindex(suffix)]
						temp_files=np.sort([glob(folder+'/*'+temp_pattern+'*') for folder in folders])[0]
						temp_dates.append(dates[isort[0]])
					elif iname==3: # /uv
						 # added for autocheck different file names
						 # added for autocheck different file names
						
						uv_pattern=pattern#files[i][files[i].rindex('/')+1:files[i].rindex(suffix)]
						uv_files=np.sort([glob(folder+'/*'+uv_pattern+'*') for folder in folders])[0]
						uv_dates.append(dates[isort[0]])
		##########

	######### reduce file list for loading intp xarray ####################
	datestr0=day_prev_year.strftime('%Y%m%d')
	datestr1=day_follow_year.strftime('%Y%m%d')
	#ssh_dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in ssh_files])
	#salt_dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in salt_files])
	#temp_dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in temp_files])
	#uv_dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in uv_files])

	ssh_dates=date_from_name(ssh_files)
	salt_dates=date_from_name(salt_files)
	temp_dates=date_from_name(temp_files)
	uv_dates=date_from_name(uv_files)



	# use strings to limit file list in folders to acutally avialable range
	# works for gcoast not AMM15
	# limit file consideration for xarray
	def limit_file_list(datestr0,datestr1,files,dates):
		ind0=np.max((np.int(datestr0)>=dates).sum()-1,0)
		ind1=(np.int(datestr1)>dates).sum()+2
		return files[ind0:ind1]

	limit_files=False
	if limit_files: # buggy with AMM15 non gcoast:
		ssh_files=limit_file_list(datestr0,datestr1,ssh_files,ssh_dates)	
		salt_files=limit_file_list(datestr0,datestr1,salt_files,salt_dates)	
		temp_files=limit_file_list(datestr0,datestr1,temp_files,temp_dates)	
		uv_files=limit_file_list(datestr0,datestr1,uv_files,uv_dates)	

		
if frocingtype!='IBI':		
		
	if np.diff([(len(liste)) for liste in [ssh_files,salt_files,temp_files,uv_files]]).max() >0:
		print('warning different number of files available for variables')

#	
if frocingtype=='IBI':

	def reduce_file_list(files):
		""" return file list and start date for IBI files """
		npt0,npt1=np.asarray([t0,t1],np.datetime64)		
		file_dates=np.asarray(get_file_start_dates(files))
		#file0=int(np.maximum(np.where(file_dates>=npt0)[0][0]-1,0)) # has day before included
		file0=int(np.maximum(np.where(file_dates>=npt0)[0][0],0)) # has day before included, this not
		if len(np.where(file_dates>=npt1)[0]) > 0:
			file1=int(np.minimum(np.where(file_dates>=npt1)[0][0],len(files)))
		else:
			file1=len(files)		
		files=files[file0:file1]
		file_dates=file_dates[file0:file1]
		return files, file_dates
		
	# crashes   - loop over files
	ssh_files,ssh_file_dates=reduce_file_list(ssh_files)
	temp_files,temp_file_dates=reduce_file_list(temp_files)
	salt_files,salt_file_dates=reduce_file_list(salt_files)
	u_files,u_file_dates=reduce_file_list(u_files)
	v_files,v_file_dates=reduce_file_list(v_files)
#
#	
#	for count,type_file in enumerate(zip(ssh_files,temp_files,salt_files,u_files,v_files)):
#		ssh_file,temp_file,salt_file,u_file,v_file=type_file
#		#only 2d is hourly
#		
#		if count==0: #determine grids
#
#			dsi=xr.open_dataset(ssh_file)
#			ssh_dates=ssh_file_dates[count]+np.arange(len(dsi['ssh']))*np.timedelta64(1,'h')
#			LON=dsi.nav_lon.values		
#			LAT=dsi.nav_lat.values	
#		break	
#		
		

# Full hour half hour?		
		

	## vis test	
	#for k in range(LON.shape[0]):
	#	plt.plot(bdlon[k],bdlat[k],'ro')
	#	plt.plot(LON[indsx[k,:],indsy[k,:]],LAT[indsx[k,:],indsy[k,:]],'k+')
    #
	## Need to adapt for curvi
	#indsx,indsy=get_4neighbours(LON,LAT,bdlon,bdlat)	
	#
	#data=LON
	#bdx,bdy=bdlon,bdlat
	#vint,bw1,bw2,bw3,bw4=bilin3(LON,LAT,data,indsx,indsy,bdx,bdy,True) #get bilinear weights	
    #
	##vis test	
	#check=(LON >= bdlon.min()-1) & (LON <= bdlon.max()+1)  & (LAT >= bdlat.min() -1) & (LAT <= bdlat.max()+1)
    #
	#a=(dsi.ssh[0,:]>-90000)*check
	#plt.figure()
	#plt.pcolormesh(LON,LAT,a)

	

    #
	## dont trust this -> nearest neighbour
	#plt.figure()
	#plt.plot(bdx,'.-')
	#plt.plot(vint,'.-')
    #
	#data=LAT
	#vint,bw1,bw2,bw3,bw4=bilin3(LON,LAT,data,indsx,indsy,bdx,bdy,True) #get bilinear weights
    #
	#plt.figure()
	#plt.plot(bdy,'.-')
	#plt.plot(vint,'.-')

name_lon='nav_lon'
name_lat='nav_lat'
name_depth='deptht'
name_time='time_counter'	
name_ssh='ssh'
#	
for count,type_file in enumerate(zip(ssh_files[:1],temp_files[:1],salt_files[:1],u_files[:1],v_files[:1])):
	ssh_file,temp_file,salt_file,u_file,v_file=type_file
	#only 2d is hourly
	
	# make grids
	varnames={'ssh':name_ssh,'salt':name_salt,'temp':name_temp,'u':name_u,'v':name_v}
	var_files={'ssh':ssh_file,'salt':salt_file,'temp':temp_file,'u':u_file,'v':v_file}
	
	grid=dict.fromkeys(varnames.keys())
	for key in grid.keys():
		dsi=xr.open_dataset(var_files[key])
		grid[key]={}
		for var in [name_lon,name_lat,name_depth,name_time]:
			try:
				grid[key][var]=dsi[var].values#dsi[varnames[key]].values
			except:
				pass
		grid[key]['dimensions']=dsi[varnames[key]].dims # oder of dimentions



		#t=ds['temp']['time'][:].values
#	dates={}	
#	for key in varnames.keys():	
#		t=ds[key][name_time].values
#		try:
#			dates[key]=datetime64_to_datetime(t)
#		except:	
#			dates[key]=np.asarray([np.datetime64(ti) for ti in t])
#	
#	
	## temp
	#dsi=xr.open_dataset(ssh_file)
	#ssh_dates=ssh_file_dates[count]+np.arange(len(dsi['ssh']))*np.timedelta64(1,'h')
    #
	## daily
	#dsi=xr.open_dataset(salt_file)
	#salt_dates=salt_file_dates[count]+np.arange(len(dsi['so']))*np.timedelta64(1,'h')

		
		
# unsure starting 0 0.5 or 1 h		

##########################################################

## check time to not produce nan in interpolation
#if True: # workaround for cmems monthly mean
#	t0data=ds['ssh'].time[0].values
#	t0data=datetime64_to_datetime([t0data])[0]
#else:
#	t0data=np.asarray(ds['ssh'].time[0].values,dt.datetime) # not working for lengths 1
#
#if t0-t0data< dt.timedelta(0):
#	raise ValueError(t0,'selected starttime {:s} for forcing output predates first availlable  time {:s} in input data'.format(str(t0),str(t0data)))
 
 
############# input model get nearest neighbours of boundary coordinates # dimension order lat lon
#if np.prod(np.asarray(grid['ssh']['dimensions'])[np.asarray(grid['ssh']['dimensions']) != name_time]==np.asarray([name_lat,name_lon])):
#	lon2d,lat2d=np.meshgrid(ds['ssh'][name_lon][:].values,ds['ssh'][name_lat][:].values)
#	tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))
#else:
#	lat2d,lon2d=np.meshgrid(ds['ssh'][name_lat][:].values,ds['ssh'][name_on][:].values)
# nearest node tree	


lon2d=grid['ssh']['nav_lon']
lat2d=grid['ssh']['nav_lat']
tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))	

#tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))	
#dd,nn=tree2d.query(bdcoords,k=1)
#ii0,jj0=np.unravel_index(nn,lat2d.shape)
x2d=lon2d
y2d=lat2d
bdx=bdlon
bdy=bdlat
x=lon2d
y=lat2d
		

################# overwrite nan indices with closest non-nan ones
# replaces entire 4 neighbour set

# ssh # T and S same grid
# u and v each different grid


# loop over all variables

ssh_file,temp_file,salt_file,u_file,v_file
ds_temp=xr.open_dataset(temp_file)
data=ds_temp[name_temp][0,0,:].values  #0 surfae -1: bottom
# get inverse distance weights for bilinear interpolation
# replace nan from forcing grid with it nn non-nan


def get_neighbours_and_weights(lon2d,lat2d,data,ndepth=50):

	indsx,indsy=get_4neighbours(lon2d,lat2d,bdlon,bdlat)		
	vint,bw1,bw2,bw3,bw4=bilin3(x,y,data,indsx,indsy,bdx,bdy,True) #get bilinear weights
	ivalid=np.where(np.isnan(vint)==False)
	bdtree=cKDTree(bdcoords[ivalid])
	inan=np.isnan(vint)
	getfrom=bdtree.query(bdcoords[inan])[1]  #excahnge with b1 - bw4

	# position switch exchange with nearest non nan parents
	#plt.figure()
	#ds['temp'][name_temp][0,0,:].plot()
	#s.plot_domain_boundaries(latlon=True,append=True)
	#plt.plot(bdx[inan],bdy[inan],'k+')
	#plt.plot(x2d,y2d,'k.')
	for var in indsx,indsy,bw1,bw2,bw3,bw4:
		for ifrom,ito in zip(ivalid[0][getfrom],np.where(inan)[0]):
			var[ito]=var[ifrom]
	vint=bilin3(x,y,data,indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)
	if np.isnan(vint).sum()>0:
		print('vint contains still nans !!!')
	
	ndepth=len(ds_temp.deptht)

	# depth extended weights
	bw1D=np.tile(bw1,(ndepth,1))
	bw2D=np.tile(bw2,(ndepth,1))
	bw3D=np.tile(bw3,(ndepth,1))
	bw4D=np.tile(bw4,(ndepth,1))

	# stack weights for faster interpolation function
	Wstack=np.zeros((len(bdx),2,2))
	Wstack[:,0,0]=bw1
	Wstack[:,0,1]=bw2
	Wstack[:,1,0]=bw3
	Wstack[:,1,1]=bw4

	return indsx,indsy,bw1,bw2,bw3,bw4,bw1D,bw2D,bw3D,bw4D,Wstack


# assign interpolation info
lon2d=grid['ssh']['nav_lon']
lat2d=grid['ssh']['nav_lat']

ds_temp=xr.open_dataset(temp_file)
data=ds_temp[name_temp][0,0,:].values  #0 surfae -1: bottom
	
indsx,indsy,bw1,bw2,bw3,bw4,bw1D,bw2D,bw3D,bw4D,Wstack=get_neighbours_and_weights(lon2d,lat2d,data,ndepth=50)	
	
class interp():
	def __init__(self,indsx,indsy,w1,w2,w3,w4):
		self.indsx=indsx.copy()
		self.indsy=indsy.copy()
		self.w1=w1.copy()
		self.w2=w2.copy()
		self.w3=w3.copy()
		self.w4=w4.copy()

ssh_interp=interp(indsx,indsy,bw1,bw2,bw3,bw4)		
salt_interp=interp(indsx,indsy,bw1D,bw2D,bw3D,bw4D)		
temp_interp=interp(indsx,indsy,bw1D,bw2D,bw3D,bw4D)				


# assign interpolation info

## u ######
lon2d=grid['u']['nav_lon']
lat2d=grid['u']['nav_lat']
dsi=xr.open_dataset(u_file)
data=dsi['uo'][0,0,:].values  #0 surfae -1: bottom

indsx,indsy,bw1,bw2,bw3,bw4,bw1D,bw2D,bw3D,bw4D,Wstack=get_neighbours_and_weights(lon2d,lat2d,data,ndepth=50)	

u_interp=interp(indsx,indsy,bw1D,bw2D,bw3D,bw4D)		

## v ######
lon2d=grid['v']['nav_lon']
lat2d=grid['v']['nav_lat']
dsi=xr.open_dataset(v_file)
data=dsi['vo'][0,0,:].values  #0 surfae -1: bottom

indsx,indsy,bw1,bw2,bw3,bw4,bw1D,bw2D,bw3D,bw4D,Wstack=get_neighbours_and_weights(lon2d,lat2d,data,ndepth=50)	
v_interp=interp(indsx,indsy,bw1D,bw2D,bw3D,bw4D)		



#		
##ph=plt.pcolormesh(lon2d,lat2d,data)
##vmin,vmax=ph.get_clim()
##plt.scatter(bdlon,bdlat,14,vint,vmin=vmin,vmax=vmax)
##data=ds['temp'][name_temp][0,1,:].values
##vint=bilin3(x,y,data,indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)
#

#
## depth extended weights
#bw1D=np.tile(bw1,(ndepth,1))
#bw2D=np.tile(bw2,(ndepth,1))
#bw3D=np.tile(bw3,(ndepth,1))
#bw4D=np.tile(bw4,(ndepth,1))
#
## stack weights for faster interpolation function
#Wstack=np.zeros((len(bdx),2,2))
#Wstack[:,0,0]=bw1
#Wstack[:,0,1]=bw2
#Wstack[:,1,0]=bw3
#Wstack[:,1,1]=bw4
#
#salt_interp=interp(indsx,indsy,bw1D,bw2D,bw3D,bw4D)		
#temp_interp=interp(indsx,indsy,bw1D,bw2D,bw3D,bw4D)		
#
#
#


ndepth=len(ds_temp.deptht)		
### same grid
if dt_output == 0:
	#dt_output=(np.diff(ds['ssh']['time'][:2])/np.timedelta64(1,'s'))[0]
	dt_output=np.timedelta64((np.diff(ds['ssh']['time'][:2])[0]))/np.timedelta64(1,'s')
deltaT=dt.timedelta(seconds=dt_output)	
tout=np.arange(t0,t1,deltaT)
nt=np.int((t1-t0)/deltaT)+1
#dschism=np.asarray(s.depths)
nnodes=len(d)


## mainpulate indices to get data for interpolation so thatthe next non nan producing 
## indices along boundary are taken - geht bestimmt irgendwie weniger kompliziert
## u and v at different sides, need different nns parially
#data0s={}
#data1s={}
#di=0
#name=name_u
#i0,i1=0,1
#key='uv'
#data0s[name]=ds[key][name][i0,0,:].values
##uint=bilin(x,y,data0s[name][di,:],indsx,indsy,bdx,bdy,div_dxdy,dxi,dyi)
#uint=bilin3(x,y,data0s[name],indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)  # get weight for bilinear
#name=name_v
#data0s[name]=ds[key][name][i0,0,:].values
#vint=bilin3(x,y,data0s[name],indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)  # get weight for bilinear
##vint=bilin(x,y,data0s[name][di,:],indsx,indsy,bdx,bdy,div_dxdy,dxi,dyi)
#
#uvindsx={name_u:indsx.copy(),name_v:indsx.copy()}
#uvindsy={name_u:indsy.copy(),name_v:indsy.copy()}
#ws1={name_u:bw1.copy(),name_v:bw1.copy()}
#ws2={name_u:bw2.copy(),name_v:bw2.copy()}
#ws3={name_u:bw3.copy(),name_v:bw3.copy()}
#ws4={name_u:bw4.copy(),name_v:bw4.copy()}
#
## exchange interpwars with bw1 - bw4
#intpvars={name_u:uint,name_v:vint}
#for key in name_u,name_v:
#	ivalid=np.where(np.isnan(intpvars[key])==False)
#	bdtree=cKDTree(bdcoords[ivalid])
#	inan=np.isnan(intpvars[key])
#	getfrom=bdtree.query(bdcoords[inan])[1]
#	for var in ws1[key],ws2[key],ws3[key],ws4[key],uvindsx[key],uvindsy[key]:
#		for ifrom,ito in zip(ivalid[0][getfrom],np.where(inan)[0]):
#			var[ito]=var[ifrom]
#
## depth extended weights			
#ws1D={name_u: np.tile(ws1[name_u],(ndepth,1)),name_v: np.tile(ws1[name_v],(ndepth,1))}
#ws2D={name_u: np.tile(ws2[name_u],(ndepth,1)),name_v: np.tile(ws2[name_v],(ndepth,1))}
#ws3D={name_u: np.tile(ws3[name_u],(ndepth,1)),name_v: np.tile(ws3[name_v],(ndepth,1))}
#ws4D={name_u: np.tile(ws4[name_u],(ndepth,1)),name_v: np.tile(ws4[name_v],(ndepth,1))}
#



# dates

## Extract Profiles from forcing and time interpolate


# maybe I need to shift one hour
nt_file=24
ntfiles=nt_file*len(ssh_files)
dates_out=ssh_file_dates[0]+np.arange(ntfiles)*np.timedelta64(1,'h')+30*np.timedelta64(1,'m')
nt=len(dates_out)

# initialize output fields
bdprofile={
'ssh':np.zeros((nt,nnodes,1,1)),
'salt':np.zeros((nt,ndepth,nnodes)),
'temp':np.zeros((nt,ndepth,nnodes)),
name_u:np.zeros((nt,ndepth,nnodes)),
name_v:np.zeros((nt,ndepth,nnodes))
}

# schism depth
bdprofiles2={
'ssh':np.zeros((nt,nnodes,1,1)),
'salt':np.zeros((nt,nnodes,len(s.vgrid[1]),1)),
'temp':np.zeros((nt,nnodes,len(s.vgrid[1]),1)),
'uv':np.zeros((nt,nnodes,len(s.vgrid[1]),2))
}

control_plot=False

# fill bd day coordinates
npt0,npt1=np.asarray([t0,t1],np.datetime64)
count=0
for count,type_file in enumerate(zip(ssh_files,temp_files,salt_files,u_files,v_files)):
	ssh_file,temp_file,salt_file,u_file,v_file=type_file
	#only 2d is hourly
	
	dsi=xr.open_dataset(ssh_file)
	ssh_dates=ssh_file_dates[count]+np.arange(len(dsi['ssh']))*np.timedelta64(1,'h')

	#datain=ds['ssh'].interp(time=tt)[name_ssh].load().values
	datain=dsi['ssh'].values

	t00=count*nt_file
	
	print('load day file ' + str(count))
	
	for ti in range(nt_file):
		#bdprofiles2['ssh'][t00+ti,:,0,0]=bilin3UseWeights_stack(x,y,datain[ti,:],ssh_interp.indsx,ssh_interp.indsy,Wstack)
		bdprofiles2['ssh'][t00+ti,:,0,0]=bilin3UseWeights(x,y,datain[ti,:],ssh_interp.indsx,ssh_interp.indsy,bdx,bdy,w1=ssh_interp.w1,w2=ssh_interp.w2,w3=ssh_interp.w3,w4=ssh_interp.w4)
		
	if control_plot:
		plt.clf()
		vmin,vmax=-2,2
		plt.pcolormesh(x2d,y2d,datain[0,:],vmin=vmin,vmax=vmax)
		#s.plot_domain_boundaries(append=True)
		plt.xlim((bdlon.min()-0.1,bdlon.max()+0.1))
		plt.ylim((bdlat.min()-.1,bdlat.max()+0.1))
		ph1=plt.scatter(bdx,bdy,s=14,c=bdprofiles2['ssh'][t00,:],vmin=vmin,vmax=vmax,edgecolors='k',linewidth=0.1)
		plt.colorbar()
		plt.savefig('control_plot_ssh_{:02d}.png'.format(count))
		
	dsi.close()
	
	# daily mean - temp
	#temp
	dsi=xr.open_dataset(temp_file)
	datain=dsi['thetao'].values[0,:]
	interp_daymean_temp=bilin3depthUseWeights(x,y,datain,temp_interp.indsx,temp_interp.indsy,bdx,bdy,w1=temp_interp.w1,w2=temp_interp.w2,w3=temp_interp.w3,w4=temp_interp.w4)
	dsi.close()
	#salt
	dsi=xr.open_dataset(salt_file)
	datain=dsi['so'].values[0,:]
	interp_daymean_salt=bilin3depthUseWeights(x,y,datain,salt_interp.indsx,salt_interp.indsy,bdx,bdy,w1=salt_interp.w1,w2=salt_interp.w2,w3=salt_interp.w3,w4=salt_interp.w4)
	dsi.close()
	#u
	dsi=xr.open_dataset(u_file)
	datain=dsi['uo'].values[0,:]
	interp_daymean_u=bilin3depthUseWeights(x,y,datain,u_interp.indsx,u_interp.indsy,bdx,bdy,w1=u_interp.w1,w2=u_interp.w2,w3=u_interp.w3,w4=u_interp.w4)
	dsi.close()
	#v
	dsi=xr.open_dataset(v_file)
	datain=dsi['vo'].values[0,:]
	interp_daymean_v=bilin3depthUseWeights(x,y,datain,v_interp.indsx,v_interp.indsy,bdx,bdy,w1=v_interp.w1,w2=v_interp.w2,w3=v_interp.w3,w4=v_interp.w4)
	dsi.close()
	
	# dailiy mean each time step
	for ti in range(nt_file):
		bdprofile['salt'][t00+ti,:]=interp_daymean_salt
		bdprofile['temp'][t00+ti,:]=interp_daymean_temp
		bdprofile['uo'][t00+ti,:]=interp_daymean_u
		bdprofile['vo'][t00+ti,:]=interp_daymean_v
# no tide information in current

		
###############################
### get out nans in vertical levels # add BJ 14.02.2022
## take values from aobve layer if nan in current one
# depth extrapolate values

depthin=ds_temp.deptht.values
#depthin=ds['salt']['depth'].values	
for name in 'salt','temp',name_u, name_v:
	for idep in range(1,len(depthin)): # 0 is surfance (Nemo and GCOAST)
		inan=np.isnan(bdprofile[name][:,:,idep])
		bdprofile[name][:,:,idep][inan]=bdprofile[name][:,:,idep-1][inan]
################# vertical interpolation to wanted levels for schism #####################################################

# transpose again		
for key in list(bdprofile.keys())[1:]:
	key
	bdprofile[key]=bdprofile[key].swapaxes(1,2)


if len(openbd_segs)>0:
	frcbdnodes=[]
	for seg in openbd_segs:
		frcbdnodes+=s.bdy_segments[seg]
		bdyvgrid = np.asarray([s.vgrid[ii].filled(-1.) for ii in frcbdnodes ])
else:
	frcbdnodes=s.bdy_nodes

# cmems  top to bottom
if np.diff(np.abs(depthin))[:1]<0:
	print('depth sorting of input data is not top to bottom which will make problem')
imaxdepin=len(depthin)-1	
for inode in range(nnodes):
	#dep=s.vgrid[ibd[inode]-1].filled(-1)*d[inode]
	dep=s.vgrid[ibd[inode]+1].filled(-1)*d[inode]
	print('depth interp profiles node {:d}/{:d}'.format(inode,nnodes))
	for vind,di in enumerate(dep):	
		iabove =  (di >= depthin).sum()-1
		ibelow=iabove+1
		iabove+=(ibelow==0)
		
		ibelow=np.minimum(imaxdepin,ibelow) # extrapolate if not deep enough
		
		dztop=np.abs(di-depthin[iabove]) 
		dzbotm=depthin[ibelow]-di
		wtop=1-dztop/(dztop+dzbotm)  
		wbotm=1-wtop
		
		bdprofiles2['salt'][:,inode,vind,0]=bdprofile['salt'][:,inode,ibelow]*wbotm+bdprofile['salt'][:,inode,iabove]*wtop
		bdprofiles2['temp'][:,inode,vind,0]=bdprofile['temp'][:,inode,ibelow]*wbotm+bdprofile['temp'][:,inode,iabove]*wtop
		
		bdprofiles2['uv'][:,inode,vind,0]=bdprofile[name_u][:,inode,ibelow]*wbotm+bdprofile[name_u][:,inode,iabove]*wtop
		bdprofiles2['uv'][:,inode,vind,1]=bdprofile[name_v][:,inode,ibelow]*wbotm+bdprofile[name_v][:,inode,iabove]*wtop

# fill nan top to bottom
for name in 'salt','temp','uv':
	for inode in range(nnodes):
		for idep in reversed(range(len(dep)-1)):
			print(idep)
			if np.isnan(bdprofiles2[name][0,inode,idep,0]):
				bdprofiles2[name][:,inode,idep,:]=bdprofiles2[name][:,inode,idep+1,:]

name='uv' #still had nan in v
for inode in range(nnodes):
	for idep in reversed(range(len(dep)-1)):
		print(idep)
		if np.isnan(bdprofiles2[name][0,inode,idep,1]):
			bdprofiles2[name][:,inode,idep,1]=bdprofiles2[name][:,inode,idep+1,1]
				
				
for key in bdprofiles2.keys():
    print(key)
    print(bdprofiles2[key].max())				

	
# houlry mean centered at half hour?
	
# export	
#timesq=[ti*dt_output for ti in range(nt)]
timesq=(dates_out-dates_out[0])/np.timedelta64(1,'s')
s.write_bdy_netcdf('elev2D'+timetag+'.th.nc',timesq,bdprofiles2['ssh'],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('TEM_3D'+timetag+'.th.nc',timesq,bdprofiles2['temp'],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('SAL_3D'+timetag+'.th.nc',timesq,bdprofiles2['salt'],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('uv3D'+timetag+'.th.nc',timesq,bdprofiles2['uv'],frcbdnodes=frcbdnodes)

# at half houlr assume hour
timesq=dates_out
s.write_bdy_netcdf('elev2D_ath'+timetag+'.th.nc',timesq[:-1],(bdprofiles2['ssh'][:-1,:]+bdprofiles2['ssh'][1:,:])/2,frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('TEM_3D_ath'+timetag+'.th.nc',timesq[:-1],(bdprofiles2['temp'][:-1,:]+bdprofiles2['temp'][1:,:])/2,frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('SAL_3D_ath'+timetag+'.th.nc',timesq[:-1],(bdprofiles2['salt'][:-1,:]+bdprofiles2['salt'][1:,:])/2,frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('uv3D_ath'+timetag+'.th.nc',timesq[:-1],(bdprofiles2['uv'][:-1,:]+bdprofiles2['uv'][1:,:])/2,frcbdnodes=frcbdnodes)

# link under final names
for fpattern in 'elev2D','TEM_3D','SAL_3D','uv3D':
	os.system('ln -s {:s}.th.nc {:s}.th.nc'.format(fpattern+timetag,fpattern))


# 1h offset
# export	
#timesq=[ti*dt_output for ti in range(nt)]
# shift 1h seems better for IBI
timesq=(dates_out-dates_out[0])/np.timedelta64(1,'s')
s.write_bdy_netcdf('elev2D'+timetag+'_m1hshift.th.nc',timesq[:-1],bdprofiles2['ssh'][1:],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('TEM_3D'+timetag+'_m1hshift.th.nc',timesq[:-1],bdprofiles2['temp'][1:],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('SAL_3D'+timetag+'_m1hshift.th.nc',timesq[:-1],bdprofiles2['salt'][1:],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('uv3D'+timetag+'_m1hshift.th.nc',timesq[:-1],bdprofiles2['uv'][1:],frcbdnodes=frcbdnodes)

# link under final names
for fpattern in 'elev2D','TEM_3D','SAL_3D','uv3D':
	os.system('ln -sf {:s}.th.nc {:s}.th.nc'.format(fpattern+timetag+'_m1hshift',fpattern))



if plot_bd:
	plt.figure()
	names={'ssh':name_ssh,'salt':name_salt,'temp':name_temp}
	

	for varname in 'ssh','salt','temp':
		name=names[varname]
		plt.clf()
		ti=np.int(nt/2)
		t=t0+ti*deltaT
		tt=t.strftime("%Y-%m-%dT%H:%M:%S")
		
		fidx=np.where(ssh_file_dates > np.asarray(t,np.datetime64))[0][0]-1
		
		
		ilarger=np.where((np.asarray(t,np.datetime64)-dates_out)/np.timedelta64(1,'h')<0)[0][0]
		dates_out[ilarger]
		ti_in_file=ilarger%24
		
		
		ds={'ssh':xr.open_dataset(ssh_files[fidx]),'salt':xr.open_dataset(salt_files[fidx]),'temp':xr.open_dataset(temp_files[fidx])}
		
		if varname=='ssh':
			datain=ds['ssh']['ssh'][ti_in_file,:].values	
		else:
			datain=ds[varname][varnames[varname]][0,:].values	
		
		
		
		profil=bdprofiles2[varname][ti,:,:]
		if len(datain.shape)==3:
			datain=datain[0,:]
			profil=profil[:,-1,:]
		vmin=np.nanmin(datain)
		vmax=np.nanmax(datain)
		#vmin=bdprofiles2['ssh'][ti,:,:].min()
		#vmax=bdprofiles2['ssh'][ti,:,:].max()
		plt.pcolormesh(x2d,y2d,datain,vmin=vmin,vmax=vmax)
		#s.plot_domain_boundaries(append=True)
		plt.xlim((bdlon.min()-0.1,bdlon.max()+0.1))
		plt.ylim((bdlat.min()-.1,bdlat.max()+0.1))
		ph1=plt.scatter(bdx,bdy,s=12,c=profil,vmin=vmin,vmax=vmax,edgecolors='k',linewidth=0.1)
		ch=plt.colorbar()
		ch.set_label(varname)
		plt.suptitle(str(t))
		plt.title('surface '+varname)
		plt.legend([ph1],['interpolated bd points'])
		plt.tight_layout()
		plt.savefig('halftime_surface_'+varname+'.png',dpi=300)

	# control plots for
	istart=0
	iend=0
	#plt.figure()
	plt.clf()
	for i in openbd_segs: #range(len(s.bdy_segments)):
		iend+=len(s.bdy_segments[i])
		inds=np.asarray(range(istart,iend))
		istart=iend
		bdnodes=np.asarray(frcbdnodes)[inds]#-1	

		bddepths=np.asarray([s.vgrid[inode].filled(-1)*s.depthsdict[inode] for inode in bdnodes])
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
			#plt.suptitle(str(t0+timesq[ti]/timesq[1]*deltaT))	
			plt.suptitle(str(t))
			plt.tight_layout()
			
			plt.savefig('bd_%i_%s'%(i,label[ti]),dpi=400)
			#plt.show()
			plt.close()
#plt.plot(bdprofiles2['ssh'][:,0,0,0])

trun=time.time()-tstart
print('program finisched after {:f} minutes'.format(trun/60))


# amm15 is to fill hour