""" Generate SCHISM boundary forcing from gcoast or amm15 data"""

# eventually failure. return order for weight smight be needing to be  w1 * B[:,0,0] +w2*B[:,1,0] + w3*B[:,0,1] + w4 * B[:,1,1]    
# but this was testet  for one node dim. then it is correct.

### !! Dont use gcoast, since Sebastian saves the Depth in a seprate file the depth reading
#### From a different file and the depth in the individual files is mean depth on not related to the depth at the nodes
import numpy as np
from scipy.interpolate import griddata
import os
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities')
sys.path.insert(0,'/work/gg0028/g260114/RUNS/GermanBight/GB_HR_Ballje/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hereon-utilities/')

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

t0=dt.datetime(2017,11,14)		#  time start WK
t1=dt.datetime(2017,11,20)		#  time end WK
#dt_output=86400    # timestep for output [seconds] take as input
#dt_output=3600 #43200         # if zero take same temporal resolution as input data
dt_output=0 #43200  					# tine is counted from t0 on

# forcing source  amm15 gcoast.
#schismdir='/gpfs/work/jacobb/data/SETUPS/SNS_Wei/' # setup directory
schismdir=os.getcwd()  # local
#schismdir=schismdir[:schismdir.rindex('/')]+'/' #assume 
#schismdir=rundir='/gpfs/work/jacobb/data/SETUPS/GB_template/'

frocingtype='amm15'#'gcoast'
frcdir='forcing_by_file/' #'/work/gg0028/g260114/mengyao/forcing/'
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

elif (frocingtype=='amm15') or (frocingtype=='cmems') :    # change here the file patterns to adapt @ Mengyao
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
elif frocingtype=='gcoast': # all files in one folder	 # Format NWS from Sebastian differen for other cmems
	folders=[frcdir,]
	date_prefix='_hi'# string to search date name fpr
	suffix='_b'
else: # all files in one folder	 # Format NWS from Sebastian differen for other cmems
	folders=[frcdir,] # default from cmems ftp download for black sea
	date_prefix='/'# string to search date name fpr starting with date
	suffix='_'
	

# determine file name pattern for variable files
files=glob(folders[0]+'/*.nc')
try:
	if frocingtype == 'cmems':
		dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex(date_prefix)+len(date_prefix)+file[file.rindex(date_prefix):].index(suffix)-1]) for file in files])
		def date_from_name(files):
			return np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex(date_prefix)+len(date_prefix)+file[file.rindex(date_prefix):].index(suffix)-1]) for file in files])
	else: #amm 15	
		dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')].replace('-','')) for file in files])
		def date_from_name(files):
			returnnp.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in files])
	
except:	
	date_prefix='_' # if cmems download from motu and not NWS
	suffix='_'
	dates=np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')].replace('-','')) for file in files])
	def date_from_name(files):
		return np.asarray([ np.int(file[file.rindex(date_prefix)+len(date_prefix):file.rindex('.')]) for file in files])
isort=np.argsort(dates)


# working for gcoast
ssh_dates,salt_dates,temp_dates,uv_dates=[],[],[],[]
checkdates=isort[dates[isort]==dates[isort[0]]] # added for autocheck different file names
for i in isort[dates[isort]==dates[isort[0]]]:
	print(i)
	dsi=xr.open_dataset(files[i])
	alt=checkdates[checkdates!=i][0]
	#f=np.asarray(files[i][files[i].rindex('/')+1:],str)
	#f2=np.asarray(files[alt][files[alt].rindex('/')+1:],str)
    #
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

if np.diff([(len(liste)) for liste in [ssh_files,salt_files,temp_files,uv_files]]).max() >0:
	print('warning different number of files available for variables')

ds={'ssh':xr.concat([xr.open_dataset(file).chunk() for file in ssh_files], dim='time'),
'salt':xr.concat([xr.open_dataset(file).chunk() for file in salt_files], dim='time'),
'uv':xr.concat([xr.open_dataset(file).chunk() for file in uv_files], dim='time')
}
if frocingtype=='gcoast':
	ds['temp']=ds['salt'] # temperature and salinity stored in same file
else:
	ds['temp']=xr.concat([xr.open_dataset(file).chunk() for file in temp_files], dim='time')
###############################################################################


############ get time and space grids by variable ########################################
# atm only considering tempoaral inhomoginities between variables not spatial
varnames={'ssh':name_ssh,'salt':name_salt,'temp':name_temp,'uv':name_u}
grid=dict.fromkeys(['ssh','salt','temp','uv'])
for key in grid.keys():
	grid[key]={}
	for var in [name_lon,name_lat,name_depth,name_time]:
		try:
			grid[key][var]=ds[key][var].values
		except:
			pass
	grid[key]['dimensions']=ds[key][varnames[key]].dims # oder of dimentions

t=ds['temp']['time'][:].values
dates={}	
for key in varnames.keys():	
	t=ds[key][name_time].values
	try:
		dates[key]=datetime64_to_datetime(t)
	except:	
		dates[key]=np.asarray([np.datetime64(ti) for ti in t])

############### Defined reduced data set in necessary 


###############- only load necsseray domain rectangle and depth range from forcing model
ds['ssh']=ds['ssh'].sel(lon=slice(lon.min()-1,lon.max()+1),lat=slice(lat.min()-1,lat.max()+1))
try:
	ndepth=((ds['salt']['depth']>np.max(depths[ibd]))==False).sum()+1
	dmaxto=ds['salt']['depth'][ndepth].values
except:	
	ndepth=np.minimum(((ds['salt']['depth']>np.max(depths[ibd]))==False).sum(),len(ds['salt']['depth'])-1)
	dmaxto=ds['salt']['depth'][ndepth].values
if frocingtype=='gcoast':
	ds['temp']=ds['salt']=ds['salt'].sel(lon=slice(lon.min()-1,lon.max()+1),lat=slice(lat.min()-1,lat.max()+1),depth=slice(-1,dmaxto))
else:	
	ds['salt']=ds['salt'].sel(lon=slice(lon.min()-1,lon.max()+1),lat=slice(lat.min()-1,lat.max()+1),depth=slice(-1,dmaxto))
	ds['temp']=ds['temp'].sel(lon=slice(lon.min()-1,lon.max()+1),lat=slice(lat.min()-1,lat.max()+1),depth=slice(-1,dmaxto))
ds['uv']=ds['uv'].sel(lon=slice(lon.min()-1,lon.max()+1),lat=slice(lat.min()-1,lat.max()+1),depth=slice(-1,dmaxto))
ndepth=len(ds['salt']['depth'])
##########################################################

# check time to not produce nan in interpolation
if True: # workaround for cmems monthly mean
	t0data=ds['ssh'].time[0].values
	t0data=datetime64_to_datetime([t0data])[0]
else:
	t0data=np.asarray(ds['ssh'].time[0].values,dt.datetime) # not working for lengths 1

if t0-t0data< dt.timedelta(0):
	raise ValueError(t0,'selected starttime {:s} for forcing output predates first availlable  time {:s} in input data'.format(str(t0),str(t0data)))
 
############# input model get nearest neighbours of boundary coordinates # dimension order lat lon
if np.prod(np.asarray(grid['ssh']['dimensions'])[np.asarray(grid['ssh']['dimensions']) != name_time]==np.asarray([name_lat,name_lon])):
	lon2d,lat2d=np.meshgrid(ds['ssh'][name_lon][:].values,ds['ssh'][name_lat][:].values)
	tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))
else:
	lat2d,lon2d=np.meshgrid(ds['ssh'][name_lat][:].values,ds['ssh'][name_on][:].values)
# nearest node tree	
tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))	
dd,nn=tree2d.query(bdcoords,k=1)
#ii0,jj0=np.unravel_index(nn,lat2d.shape)
x2d=lon2d
y2d=lat2d
bdx=bdlon
bdy=bdlat
x=lon2d
y=lat2d
		
# check
#data=ds['temp']['votemper'][0,0,:].values

#data=np.asarray(x2d,np.float64)
#x=np.asarray(x,np.float64)
#y=np.asarray(y,np.float64)
#indsx,indsy=get_4neighbours(lon2d,lat2d,bdlon,bdlat)		
#tt0=time.time()
#for i in range(30):
#	vint,bw1,bw2,bw3,bw4=bilin3(x,y,data,indsx,indsy,bdx,bdy,True) #get bilinear weights
#dt0=time.time()-tt0
#
#tt0=time.time()
#for i in range(30):
#	vint2=bilin3(x,y,data,indsx,indsy,bdx,bdy,False,w1=bw1,w2=bw2,w3=bw3,w4=bw4) #get bilinear weights
#dt1=time.time()-tt0
#
#tt0=time.time()
#for i in range(30):
#	vint3=bilin3UseWeights(x,y,data,indsx,indsy,bdx,bdy,w1=bw1,w2=bw2,w3=bw3,w4=bw4)
#dt2=time.time()-tt0
#
#tt0=time.time()
#for i in range(30):
#	vint4=bilin3UseWeightsb(x,y,data,indsx,indsy,bdx,bdy,w1=bw1,w2=bw2,w3=bw3,w4=bw4)
#dt3=time.time()-tt0


#tt0=time.time()
#for i in range(350):
#	vint3=bilin3UseWeights(x,y,data,indsx,indsy,bdx,bdy,w1=bw1,w2=bw2,w3=bw3,w4=bw4)
#dt2=time.time()-tt0
#
#tt0=time.time()
#for i in range(350):
#	vint4=bilin3UseWeights_stack(x,y,data,indsx,indsy,Wstack)
#dt3=time.time()-tt0

#
#
#vint3[:5]-vint4[:5]
#vint5[:5]-vint4[:5]
#
#
#np.dot(q,Wstack)
#np.cross(Wstack,q)
#
#ivalid=np.where(np.isnan(vint)==False)
#bdtree=cKDTree(bdcoords[ivalid])
#inan=np.isnan(vint)
#getfrom=bdtree.query(bdcoords[inan])[1]
#vint[inan]=vint[ivalid][getfrom]		
#
##ax=plt.axis()
#plt.pcolormesh(lon2d,lat2d,data)
#plt.axis(ax)
#plt.colorbar()
#plt.clim((0,6))
#plt.scatter(xint,yint,1,np.asarray(vint))
#plt.clim((0,6))
##################################################
#

################# overwrite nan indices with closest non-nan ones
# replaces entire 4 neighbour set
data=ds['temp'][name_temp][0,0,:].values
# get inverse distance weights for bilinear interpolation
# replace nan from forcing grid with it nn non-nan
indsx,indsy=get_4neighbours(lon2d,lat2d,bdlon,bdlat)		
vint,bw1,bw2,bw3,bw4=bilin3(x,y,data,indsx,indsy,bdx,bdy,True) #get bilinear weights
ivalid=np.where(np.isnan(vint)==False)
bdtree=cKDTree(bdcoords[ivalid])
inan=np.isnan(vint)
getfrom=bdtree.query(bdcoords[inan])[1]  #excahnge with b1 - bw4

# position switch
#plt.figure()
#ds['temp'][name_temp][0,0,:].plot()
#s.plot_domain_boundaries(latlon=True,append=True)
#plt.plot(bdx[inan],bdy[inan],'k+')
#plt.plot(x2d,y2d,'k.')

for var in indsx,indsy,bw1,bw2,bw3,bw4:
	for ifrom,ito in zip(ivalid[0][getfrom],np.where(inan)[0]):
		var[ito]=var[ifrom]
vint=bilin3(x,y,data,indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)

#data=ds['temp'][name_temp][0,1,:].values
vint=bilin3(x,y,data,indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)

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




# stack weights for alernative interpolation function -- seems to get slower opposed to the 2d case
#WstackD=np.tile(Wstack,(ndepth,1,1,1))
#datain=np.tile(x2d,(ndepth,1,1))
#indsx,indsy=get_4neighbours(lon2d,lat2d,bdlon,bdlat)		
#tt0=time.time()
#for i in range(30):
#	vint4=bilin3depthUseWeights(x,y,datain,indsx,indsy,bdx,bdy,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D)
#dt0=time.time()-tt0
#
#tt0=time.time()
#for i in range(30):
#	vint5=bilin3depthUseWeights_stack(x,y,datain,indsx,indsy,WstackD)
#dt1=time.time()-tt0



		
### same grid
if dt_output == 0:
	#dt_output=(np.diff(ds['ssh']['time'][:2])/np.timedelta64(1,'s'))[0]
	dt_output=np.timedelta64((np.diff(ds['ssh']['time'][:2])[0]))/np.timedelta64(1,'s')
deltaT=dt.timedelta(seconds=dt_output)	
tout=np.arange(t0,t1,deltaT)
nt=np.int((t1-t0)/deltaT)+1
#dschism=np.asarray(s.depths)
nnodes=len(d)


# initialize output fields
# input depths
# alt
#bdprofile={
#'ssh':np.zeros((nt,nnodes,1,1)),
#'salt':np.zeros((nt,nnodes,ndepth)),
#'temp':np.zeros((nt,nnodes,ndepth)),
#name_u:np.zeros((nt,nnodes,ndepth)),
#name_v:np.zeros((nt,nnodes,ndepth))
#}
#neu -bypass transponate in loop# later transponate once to get above layout for interpolation to schism profiles
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



# mainpulate indices to get data for interpolation so thatthe next non nan producing 
# indices along boundary are taken - geht bestimmt irgendwie weniger kompliziert
# u and v at different sides, need different nns parially
data0s={}
data1s={}
di=0
name=name_u
i0,i1=0,1
key='uv'
data0s[name]=ds[key][name][i0,0,:].values
#uint=bilin(x,y,data0s[name][di,:],indsx,indsy,bdx,bdy,div_dxdy,dxi,dyi)
uint=bilin3(x,y,data0s[name],indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)  # get weight for bilinear
name=name_v
data0s[name]=ds[key][name][i0,0,:].values
vint=bilin3(x,y,data0s[name],indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)  # get weight for bilinear
#vint=bilin(x,y,data0s[name][di,:],indsx,indsy,bdx,bdy,div_dxdy,dxi,dyi)

uvindsx={name_u:indsx.copy(),name_v:indsx.copy()}
uvindsy={name_u:indsy.copy(),name_v:indsy.copy()}
ws1={name_u:bw1.copy(),name_v:bw1.copy()}
ws2={name_u:bw2.copy(),name_v:bw2.copy()}
ws3={name_u:bw3.copy(),name_v:bw3.copy()}
ws4={name_u:bw4.copy(),name_v:bw4.copy()}

# exchange interpwars with bw1 - bw4
intpvars={name_u:uint,name_v:vint}
for key in name_u,name_v:
	ivalid=np.where(np.isnan(intpvars[key])==False)
	bdtree=cKDTree(bdcoords[ivalid])
	inan=np.isnan(intpvars[key])
	getfrom=bdtree.query(bdcoords[inan])[1]
	for var in ws1[key],ws2[key],ws3[key],ws4[key],uvindsx[key],uvindsy[key]:
		for ifrom,ito in zip(ivalid[0][getfrom],np.where(inan)[0]):
			var[ito]=var[ifrom]

# depth extended weights			
ws1D={name_u: np.tile(ws1[name_u],(ndepth,1)),name_v: np.tile(ws1[name_v],(ndepth,1))}
ws2D={name_u: np.tile(ws2[name_u],(ndepth,1)),name_v: np.tile(ws2[name_v],(ndepth,1))}
ws3D={name_u: np.tile(ws3[name_u],(ndepth,1)),name_v: np.tile(ws3[name_v],(ndepth,1))}
ws4D={name_u: np.tile(ws4[name_u],(ndepth,1)),name_v: np.tile(ws4[name_v],(ndepth,1))}


## Extract Profiles from forcing and time interpolate
#for ti in range(nt):
#	t=t0+ti*deltaT
#	#print(ti)
#	if ti%24==0:
#		print('doing day ' + str(ti/24))
#
#	#datets working with gcoast but not amm15 therefore convert to string to fit both :-(	
#	# probbly makes things slower
#	# ssh
#	tt=t.strftime("%Y-%m-%dT%H:%M:%S")
#	
#	datain=ds['ssh'].interp(time=tt)[name_ssh].load().values
#	#bdprofiles2['ssh'][ti,:,0,0]=bilin3(x,y,datain,indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)
#	#bdprofiles2['ssh'][ti,:,0,0]=bilin3UseWeights(x,y,datain,indsx,indsy,bdx,bdy,bw1,bw2,bw3,bw4)
#	bdprofiles2['ssh'][ti,:,0,0]=bilin3UseWeights_stack(x,y,datain,indsx,indsy,Wstack)
#
#	# temp	
#	datain=ds['temp'].interp(time=tt)[name_temp].load().values
#	#bdprofile['temp'][ti,:]=bilin3depth(x,y,datain,indsx,indsy,bdx,bdy,calcw=False,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D).T
#	bdprofile['temp'][ti,:]=bilin3depthUseWeights(x,y,datain,indsx,indsy,bdx,bdy,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D).T
#	ds['temp'].interp(time=tt)[name_temp].load().values
#	
#	#salt
#	datain=ds['salt'].interp(time=tt)[name_salt].load().values
#	#bdprofile['salt'][ti,:]=bilin3depth(x,y,datain,indsx,indsy,bdx,bdy,calcw=False,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D).T
#	bdprofile['salt'][ti,:]=bilin3depthUseWeights(x,y,datain,indsx,indsy,bdx,bdy,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D).T
#	
#	# uv
#	for name in name_u,name_v:   # faster doing tim einterp once in xarray
#		#data0s[name]=(w1*ds[key][name][i0,:]+w2*ds[key][name][i1,:]).values # time interp in xarray
#		datain=ds['uv'].interp(time=tt)[name].load().values
#		#bdprofile[name][ti,:]=bilin3depth(x,y,datain,uvindsx[name],uvindsy[name],bdx,bdy,False,ws1D[name],ws2D[name],ws3D[name],ws4D[name]).T
#		bdprofile[name][ti,:]=bilin3depthUseWeights(x,y,datain,uvindsx[name],uvindsy[name],bdx,bdy,ws1D[name],ws2D[name],ws3D[name],ws4D[name]).T

steps_perday=86400/dt_output

# Extract Profiles from forcing and time interpolate
for ti in range(nt):
	t=t0+ti*deltaT
	#print(ti)
	if ti%steps_perday==0:
		print('doing day ' + str(ti/steps_perday))

	#datets working with gcoast but not amm15 therefore convert to string to fit both :-(	
	# probbly makes things slower
	# ssh
	tt=t.strftime("%Y-%m-%dT%H:%M:%S")
	
	datain=ds['ssh'].interp(time=tt)[name_ssh].load().values
	#bdprofiles2['ssh'][ti,:,0,0]=bilin3(x,y,datain,indsx,indsy,bdx,bdy,False,bw1,bw2,bw3,bw4)
	#bdprofiles2['ssh'][ti,:,0,0]=bilin3UseWeights(x,y,datain,indsx,indsy,bdx,bdy,bw1,bw2,bw3,bw4)
	bdprofiles2['ssh'][ti,:,0,0]=bilin3UseWeights_stack(x,y,datain,indsx,indsy,Wstack)

	# temp	
	datain=ds['temp'].interp(time=tt)[name_temp].load().values
	#bdprofile['temp'][ti,:]=bilin3depth(x,y,datain,indsx,indsy,bdx,bdy,calcw=False,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D).T
	bdprofile['temp'][ti,:]=bilin3depthUseWeights(x,y,datain,indsx,indsy,bdx,bdy,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D)
	ds['temp'].interp(time=tt)[name_temp].load().values
	
	#salt
	datain=ds['salt'].interp(time=tt)[name_salt].load().values
	#bdprofile['salt'][ti,:]=bilin3depth(x,y,datain,indsx,indsy,bdx,bdy,calcw=False,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D).T
	bdprofile['salt'][ti,:]=bilin3depthUseWeights(x,y,datain,indsx,indsy,bdx,bdy,w1=bw1D,w2=bw2D,w3=bw3D,w4=bw4D)
	
	# uv
	for name in name_u,name_v:   # faster doing tim einterp once in xarray
		#data0s[name]=(w1*ds[key][name][i0,:]+w2*ds[key][name][i1,:]).values # time interp in xarray
		datain=ds['uv'].interp(time=tt)[name].load().values
		#bdprofile[name][ti,:]=bilin3depth(x,y,datain,uvindsx[name],uvindsy[name],bdx,bdy,False,ws1D[name],ws2D[name],ws3D[name],ws4D[name]).T
		bdprofile[name][ti,:]=bilin3depthUseWeights(x,y,datain,uvindsx[name],uvindsy[name],bdx,bdy,ws1D[name],ws2D[name],ws3D[name],ws4D[name])
		

###################################
### get out nans in vertical levels # add BJ 14.02.2022
## take values from aobve layer if nan in current one
# depth extrapolate values
depthin=ds['salt']['depth'].values	
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

#positive('')
#depthin=ds['salt']['depth'].values	
# inversorder gcoats top -> to bottom
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
			if np.isnan(bdprofiles2[name][0,inode,idep,0]):
				bdprofiles2[name][:,inode,idep,:]=bdprofiles2[name][:,inode,idep+1,:]

for key in bdprofiles2.keys():
    print(bdprofiles2[key].max())				
		
# export	
timesq=[ti*dt_output for ti in range(nt)]
s.write_bdy_netcdf('elev2D'+timetag+'.th.nc',timesq,bdprofiles2['ssh'],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('TEM_3D'+timetag+'.th.nc',timesq,bdprofiles2['temp'],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('SAL_3D'+timetag+'.th.nc',timesq,bdprofiles2['salt'],frcbdnodes=frcbdnodes)
s.write_bdy_netcdf('uv3D'+timetag+'.th.nc',timesq,bdprofiles2['uv'],frcbdnodes=frcbdnodes)


# link under final names
for fpattern in 'elev2D','TEM_3D','SAL_3D','uv3D':
	os.system('ln -s {:s}.th.nc {:s}.th.nc'.format(fpattern+timetag,fpattern))


if plot_bd:
	plt.figure()
	names={'ssh':name_ssh,'salt':name_salt,'temp':name_temp}
	

	for varname in 'ssh','salt','temp':
		name=names[varname]
		plt.clf()
		ti=np.int(nt/2)
		t=t0+ti*deltaT
		tt=t.strftime("%Y-%m-%dT%H:%M:%S")
		datain=ds[varname].interp(time=tt)[name].load().values
		profil=bdprofiles2[varname][ti,:,:]
		if len(datain.shape)==3:
			datain=datain[0,:]
			profil=profil[:,-1,:]
		vmin=np.nanmin(datain)
		vmax=np.nanmax(datain)
		#vmin=bdprofiles2['ssh'][ti,:,:].min()
		#vmax=bdprofiles2['ssh'][ti,:,:].max()
		plt.pcolormesh(x2d,y2d,datain,vmin=vmin,vmax=vmax)
		#plt.scatter(bdx,bdy,s=15,c=profil,vmin=vmin,vmax=vmax)#,edgecolors=None,linewidth=0.01)
		ph1=plt.scatter(bdx,bdy,s=8,c=profil,vmin=vmin,vmax=vmax,edgecolors='k',linewidth=0.1)
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
			plt.suptitle(str(t0+timesq[ti]/timesq[1]*deltaT))	
			plt.tight_layout()
			
			plt.savefig('bd_%i_%s'%(i,label[ti]),dpi=400)
			#plt.show()
			plt.close()
#plt.plot(bdprofiles2['ssh'][:,0,0,0])

trun=time.time()-tstart
print('program finisched after {:f} minutes'.format(trun/60))