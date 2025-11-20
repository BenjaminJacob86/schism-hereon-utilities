# prepreare SCHISM boundary forcing based on cmems

import os
import sys
import glob
sys.path.append(os.environ['HOME']+'/code/python/')
sys.path.append('/pf/g/g260114/git/hzg/schism-hzg-utilities/')
from schism import schism_setup
import matplotlib
#matplotlib.use('Agg') # nin interactive backend
from matplotlib import pyplot as plt
from pylab import *
import numpy as np
import netCDF4
from scipy.spatial import cKDTree
from cftime import utime
from scipy.interpolate import interp1d
import scipy.io as io
import datetime as dt



####### Settings ###########################

# needs version of grid with only the open boundaries using forcing i.e. no rivers
openbd_segs=[0] 										  # boundary segments to create forcing fore (start counting from 0)
#setupdir='/gpfs/work/jacobb/data/SETUPS/Europe/HighRes_repair2/' # schism run direcorty containing hgrids and vgrid
#sourcedir='/gpfs/work/jacobb/data/DATA/Forcing/mount_cmems2/' # a folder containing all netcdf files to be use for 

setupdir='/work/gg0028/g260114/RUNS/Europe/c1/flow_tweak/' # schism run direcorty containing hgrids and vgrid
#sourcedir='/gpfs/work/jacobb/data/DATA/Forcing/mount_cmems2/' # a folder containing all netcdf files to be use for 
sourcedir='/gpfs/work/jacobb/data/RUNS/BlackSea/RUN24d/Forcing/data/my.cmems-du.eu/Core/MEDSEA_MULTIYEAR_PHY_006_004/all/'

#forcing
zlvls=[]
startdate=dt.datetime(2008,9,1)								  # simulation start date	
#enddate=dt.datetime(2012,7,31)								  # simulation endate date	
enddate=dt.datetime(2009,12,30)								  # simulation endate date	
output_dt=86400 										  # output timestep to interpolate to
dointerp=1													  # interpolate forcing ()  0:keep input data time step 		
#input files start date
different_files_for_vars=True								  #  tem/salt etc in different files?
 
############### end settings ##################################




######## Begin code ####################


year=startdate.year
mon0=startdate.month

# code snippets from having forcing nc in subfolders seperated by year
# check last/first day of year before/after for interpolation
#yb4=os.path.isdir(sourcedir0+str(year-1))
#yaftr=os.path.isdir(sourcedir0+str(year+1))
#files=(glob.glob(sourcedir+'*'+str(year)+'*.nc'))
#if yb4:
#	files=glob.glob(sourcedir0+str(year-1)+'/*'+str(year-1)+'*1231.nc')+files
#else:	
#	files=[sort(files)[0]]+files

# different files for varables

prevdate=startdate-dt.timedelta(days=1) # one day before model start to linptd data for interpolation
nextdate=enddate+dt.timedelta(days=1)
nextdate_string='{:d}{:02d}{:02d}'.format(nextdate.year,nextdate.month,nextdate.day)
#inputdata
startdate_string='{:d}{:02d}{:02d}'.format(prevdate.year,prevdate.month,prevdate.day)
enddate_string='{:d}{:02d}{:02d}'.format(enddate.year,enddate.month,enddate.day)
if different_files_for_vars:
	type=['ssh','sal','tem','cur'] # filenames for different variables
	files=dict.fromkeys(type)
	ds=dict.fromkeys(type)
	varnames=dict.fromkeys(type)
	for key in files.keys(): 
		filelist=np.sort(glob.glob(sourcedir+'*'+key+'*.nc'))
		i0 = [i for i, s in enumerate(filelist) if startdate_string in s][0]
		i1 = [i for i, s in enumerate(filelist) if nextdate_string in s]
		if len(i1)==0:
			i1=len(filelist)
		else:
			i1=i1[0]
	
		files[key]=	filelist[i0:i1+1]
else:
	files=np.sort(glob.glob(sourcedir+'*.nc'))
	i0 = [i for i, s in enumerate(files) if startdate_string in s][0]
	i1 = [i for i, s in enumerate(files) if nextdate_string in s]
	if len(i1)==0:
		i1=len(filels)
	else:
		i1=i1[0]
	files=files[i0:i1+1]

	
	
#if yaftr:
#	files+=glob.glob(sourcedir0+str(year+1)+'/*'+str(year+1)+'*0101.nc')
#else:	
#	files+=[sort(files[-1])]

###################################################################	


files=list(np.sort(files))
cwd=os.getcwd()

############### timie #########################################
#offstart=datetime.date(year,mon0,1)-datetime.date(year,1,1)
#offset=int(offstart.days)
#files=files[offset:]
nrdays=(enddate-startdate).days


# SCHISM Time reference
ut = utime('seconds since {0}-{1:02d}-01 00:00:00'.format(year,mon0)) 
# reference time - input nc files
nc={}
nc[0]=netCDF4.Dataset(files[0])
nc[1]=netCDF4.Dataset(files[1])
ut0=utime(nc[0]['time'].units)
date0=ut0.num2date(nc[0]['time'][0])
t0=(date0-ut.num2date(0)).total_seconds()
if len(nc[0]['time'])==1:
	date1=ut0.num2date(nc[1]['time'][0])
else:	
	date1=ut0.num2date(nc[0]['time'][1])
Dt=(date1-date0).total_seconds()
nc[0].close()
nc[1].close()


# data properties
t0=(date0-startdate).total_seconds() # start time of first data timestep relative to precribedmodel start
if t0>0:
	print('forcing data contains no time step before desired simulation start')
timesq=np.arange(0,(nrdays)*Dt+Dt,Dt) # output timestep for model in sconds since model start
nfiles=nrdays+1
################################################################


############# check input files complete
missing=[]
print('check for missing files')
fileprefix=('cmems_TS_0083deg_')
if len(files)< nrdays+2:
	print('there are input files missing')
	# check file list
	for day in range(0,int(nrdays)):
		searchstr=fileprefix+str(ut.num2date(day*86400))[0:10].replace('-','')+'.nc'
		ind=np.where(files == sourcedir+searchstr)
		if len(ind[0])==0:
			missing.append(searchstr)
	print(str(len(missing)) + ' fielse are missing (saved into missingForcing.txt)')
	f=open('missingForcing.txt','w')
	for i in missing:
		f.write(i+'\n')
	f.close()
	print(missing)
	print('exiting program')
	exit()
else:
	print('all needed files exist')
###############################################################
class cmems_one_file():

	def __init__(self,file):
		snc = netCDF4.Dataset(file)
		print(file)	
		sv = snc.variables
		self.lon = sv['longitude'][:]#[lonslice]
		self.lat = sv['latitude'][:]#[latslice]
		self.d = -sv['depth'][:]
		self.time = sv['time'][:]
		self.timeunits = sv['time'].units
		self.tidx = 0
		self.s = sv['so'][:,:,:,:]#[:,:,latslice,lonslice]
		self.t = sv['thetao'][:,:,:,:]#[:,:,latslice,lonslice]
		self.lon2,self.lat2 = meshgrid(self.lon,self.lat)
		self.u = sv['uo'][:,:,:,:]#[:,:,latslice,lonslice]
		self.v = sv['vo'][:,:,:,:]#[:,:,latslice,lonslice]
		self.ssh = sv['zos'][:,:,:]#[:,:,latslice,lonslice]
		snc.close()
		self.use_depth_slices=False
	
		print('initialize look up trees')	
		vlon = self.lon2[where(self.ssh.mask[self.tidx]==False)]
		vlat = self.lat2[where(self.ssh.mask[self.tidx]==False)]
		self.ssh_tree = cKDTree(list(zip(vlon,vlat)))
		self.ssh_var = self.ssh[self.tidx][where(self.ssh.mask[self.tidx]==False)].flatten()
	
		if not(self.use_depth_slices):
			self.d3,self.lat3,self.lon3 = meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
			vlon3 = self.lon3[where(self.s.mask[self.tidx]==False)]
			vlat3 = self.lat3[where(self.s.mask[self.tidx]==False)]
			vd3 = self.d3[where(self.s.mask[self.tidx]==False)]
			self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
		
			vlon3 = self.lon3[where(self.t.mask[self.tidx]==False)]
			vlat3 = self.lat3[where(self.t.mask[self.tidx]==False)]
			vd3 = self.d3[where(self.t.mask[self.tidx]==False)]
			self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
		
			vlon3 = self.lon3[where(self.u.mask[self.tidx]==False)]
			vlat3 = self.lat3[where(self.u.mask[self.tidx]==False)]
			vd3 = self.d3[where(self.u.mask[self.tidx]==False)]
			self.u_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
		
		
			vlon3 = self.lon3[where(self.v.mask[self.tidx]==False)]
			vlat3 = self.lat3[where(self.v.mask[self.tidx]==False)]
			vd3 = self.d3[where(self.v.mask[self.tidx]==False)]
			self.v_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
		
		
			self.s_var = self.s[self.tidx][where(self.s.mask[self.tidx]==False)].flatten()
			self.t_var = self.t[self.tidx][where(self.t.mask[self.tidx]==False)].flatten()
			self.u_var = self.u[self.tidx][where(self.u.mask[self.tidx]==False)].flatten()
			self.v_var = self.v[self.tidx][where(self.v.mask[self.tidx]==False)].flatten()
		
	
		else:
			#build trees
			self.s_tree={}
			self.t_tree={}
			self.u_tree={}
			self.v_tree={}
			self.svar={}
			self.tvar={}
			self.uvar={}
			self.vvar={}
	
	
	
			for ik,d in enumerate(self.d):
				print('	 build trees for depth %0.2f'%d)
				vlon = self.lon2[where(self.s.mask[self.tidx,ik]==False)]
				vlat = self.lat2[where(self.s.mask[self.tidx,ik]==False)]
				self.s_tree[ik] = cKDTree(zip(vlon,vlat))
				vlon = self.lon2[where(self.t.mask[self.tidx,ik]==False)]
				vlat = self.lat2[where(self.t.mask[self.tidx,ik]==False)]
				self.t_tree[ik] = cKDTree(zip(vlon,vlat))
				vlon = self.lon2[where(self.u.mask[self.tidx,ik]==False)]
				vlat = self.lat2[where(self.u.mask[self.tidx,ik]==False)]
				self.u_tree[ik] = cKDTree(zip(vlon,vlat))
				vlon = self.lon2[where(self.v.mask[self.tidx,ik]==False)]
				vlat = self.lat2[where(self.v.mask[self.tidx,ik]==False)]
				self.v_tree[ik] = cKDTree(zip(vlon,vlat))
				self.svar[ik] = self.s[self.tidx,ik][where(self.s.mask[self.tidx,ik]==False)].flatten()
				self.tvar[ik] = self.t[self.tidx,ik][where(self.t.mask[self.tidx,ik]==False)].flatten()
				self.uvar[ik] = self.u[self.tidx,ik][where(self.u.mask[self.tidx,ik]==False)].flatten()
				self.vvar[ik] = self.v[self.tidx,ik][where(self.v.mask[self.tidx,ik]==False)].flatten()


	def update(self,file):
		snc = netCDF4.Dataset(file)
		print(file)	
		sv = snc.variables
		self.lon = sv['longitude'][:]#[lonslice]
		self.lat = sv['latitude'][:]#[latslice]
		self.d = -sv['depth'][:]
		self.time = sv['time'][:]
		self.timeunits = sv['time'].units
		self.tidx = 0
		self.s = sv['so'][:,:,:,:]#[:,:,latslice,lonslice]
		self.t = sv['thetao'][:,:,:,:]#[:,:,latslice,lonslice]
		self.lon2,self.lat2 = meshgrid(self.lon,self.lat)
		self.u = sv['uo'][:,:,:,:]#[:,:,latslice,lonslice]
		self.v = sv['vo'][:,:,:,:]#[:,:,latslice,lonslice]
		self.ssh = sv['zos'][:,:,:]#[:,:,latslice,lonslice]
		snc.close()
		self.ssh_var = self.ssh[self.tidx][where(self.ssh.mask[self.tidx]==False)].flatten()
		if not(oa.use_depth_slices):
				oa.s_var = oa.s[oa.tidx][where(oa.s.mask[oa.tidx]==False)].flatten()
				oa.t_var = oa.t[oa.tidx][where(oa.t.mask[oa.tidx]==False)].flatten()
				oa.u_var = oa.u[oa.tidx][where(oa.u.mask[oa.tidx]==False)].flatten()
				oa.v_var = oa.v[oa.tidx][where(oa.v.mask[oa.tidx]==False)].flatten()
		else:
				for ik,d in enumerate(self.d):
					oa.svar[ik] = oa.s[oa.tidx,ik][where(oa.s.mask[oa.tidx,ik]==False)].flatten()
					oa.tvar[ik] = oa.t[oa.tidx,ik][where(oa.t.mask[oa.tidx,ik]==False)].flatten()
					oa.uvar[ik] = oa.u[oa.tidx,ik][where(oa.u.mask[oa.tidx,ik]==False)].flatten()
					oa.vvar[ik] = oa.v[oa.tidx,ik][where(oa.v.mask[oa.tidx,ik]==False)].flatten()
	
		
	def interpolate(self,depths,nodelon,nodelat,bidx=1):
		# start
		t = zeros((len(depths),))
		s = zeros((len(depths),))
		u = zeros((len(depths),))
		v = zeros((len(depths),))
		ssh=zeros(1)	
		
		dist,inds = self.ssh_tree.query((nodelon,nodelat),k=4)
		w = 1 / dist
		ssh[0] = np.sum(w*self.ssh_var[inds],axis=0) / np.sum(w,axis=0)
		
		
		for ik,ndepth in enumerate(depths[bidx-1:]):
			# find vertical layer in climatology
			if self.use_depth_slices:
				didx = np.abs(self.d - ndepth).argmin()
			
				#didx = int(where(ddiff==ddiff.min())[0][0])
				#vlon = self.lon2[where(self.s.mask[self.tidx,didx]==False)]
				#vlat = self.lat2[where(self.s.mask[self.tidx,didx]==False)]
				#tree = cKDTree(zip(vlon,vlat))
				
				#svar = self.s[self.tidx,didx][where(self.s.mask[self.tidx,didx]==False)].flatten()
				#tvar = self.t[self.tidx,didx][where(self.t.mask[self.tidx,didx]==False)].flatten()
			
				dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
				w = 1 / dist
				s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)
			
				dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
				w = 1 / dist
				t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
		
			# ASSUME SAME GRID # >> awakawa C
				u[bidx-1+ik] = np.sum(w*self.uvar[didx][inds],axis=0) / np.sum(w,axis=0)
				v[bidx-1+ik] = np.sum(w*self.vvar[didx][inds],axis=0) / np.sum(w,axis=0)
		
			else:
				dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
				w = 1 / dist
				s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)
		
				#dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
				#w = 1 / dist
				t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)
		
				dist,inds = self.u_tree.query((nodelon*100.,nodelat*100.,ndepth))
				u[bidx-1+ik] = np.sum(w*self.u_var[inds],axis=0) / np.sum(w,axis=0)
				dist,inds = self.v_tree.query((nodelon*100.,nodelat*100.,ndepth))
				v[bidx-1+ik] = np.sum(w*self.v_var[inds],axis=0) / np.sum(w,axis=0)
		
		
		return (t,s,u,v,ssh)

	def add_v_tree(self):
		print('rebuilding v_tree')
		vlon3 = self.lon3[where(self.v.mask[self.tidx]==False)]
		vlat3 = oa.lat3[where(self.v.mask[self.tidx]==False)]
		vd3 = oa.d3[where(oa.v.mask[self.tidx]==False)]
		self.v_tree2=cKDTree(zip(vlon3,vlat3,vd3))
		print("done building v_tree")



	def interpolate_adapt(self,depths,nodelon,nodelat,bidx=1):
		# within cmems the mask has change vor v velocity only in
		# 20090614	# here re build the look uptry temporrary
		t = zeros((len(depths),))
		s = zeros((len(depths),))
		u = zeros((len(depths),))
		v = zeros((len(depths),))
		ssh=zeros(1)	
	
		dist,inds = self.ssh_tree.query((nodelon,nodelat),k=4)
		w = 1 / dist
		ssh[0] = np.sum(w*self.ssh_var[inds],axis=0) / np.sum(w,axis=0)


		for ik,ndepth in enumerate(depths[bidx-1:]):
			# find vertical layer in climatology
			if self.use_depth_slices:
				didx = np.abs(self.d - ndepth).argmin()
			
				#didx = int(where(ddiff==ddiff.min())[0][0])
				#vlon = self.lon2[where(self.s.mask[self.tidx,didx]==False)]
				#vlat = self.lat2[where(self.s.mask[self.tidx,didx]==False)]
				#tree = cKDTree(zip(vlon,vlat))
				
				#svar = self.s[self.tidx,didx][where(self.s.mask[self.tidx,didx]==False)].flatten()
				#tvar = self.t[self.tidx,didx][where(self.t.mask[self.tidx,didx]==False)].flatten()
		
				dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
				w = 1 / dist
				s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)
		
				dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
				w = 1 / dist
				t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
		
				# ASSUME SAME GRID # >> awakawa C
				u[bidx-1+ik] = np.sum(w*self.uvar[didx][inds],axis=0) / np.sum(w,axis=0)
				v[bidx-1+ik] = np.sum(w*self.vvar[didx][inds],axis=0) / np.sum(w,axis=0)
	
			else:
	
				dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
				w = 1 / dist
				s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)
		
				#dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
				#w = 1 / dist
				t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)
	
		dist,inds = self.u_tree.query((nodelon*100.,nodelat*100.,ndepth))
		u[bidx-1+ik] = np.sum(w*self.u_var[inds],axis=0) / np.sum(w,axis=0)
		dist,inds = self.v_tree2.query((nodelon*100.,nodelat*100.,ndepth))
		v[bidx-1+ik] = np.sum(w*self.v_var[inds],axis=0) / np.sum(w,axis=0)
	
	
		return (t,s,u,v,ssh)


class amm():
  # sigma lvels surface, midwater column, near bed	
  def __init__(self,SSHfile=None,Sfile=None,Tfile=None,UVfile=None,domain_rect=[-9999, 9999, -9999, 9999]):
    self.tnc = netCDF4.Dataset(Tfile)
    self.tv = self.tnc.variables
    self.snc = netCDF4.Dataset(Sfile)
    self.sv = self.snc.variables
    self.sshnc = netCDF4.Dataset(SSHfile)
    self.sshv = self.sshnc.variables
    self.uvnc = netCDF4.Dataset(UVfile)
    self.uv = self.uvnc.variables
	
    ########### get variable file names #############################
    self.internal_names=['salt','temp','ssh','u','v']
    self.files={'salt':Sfile,'temp':Tfile,'ssh':SSHfile,'u':UVfile,'v':UVfile}
    self.ncs={}					
    for key,item in zip(self.files.keys(),self.files.items()):
        if item[1]==None:
            self.internal_names.remove(key)
        else:
            self.ncs[key]=Dataset(item[1])
	
    file=self.files[self.internal_names[0]]
    self.dir=file[:file.rfind('/')+1]
			
    # get varnames from standard names
    long_names={'salt':['salinity','salt'],'temp':'temperature','ssh':'sea surface height','u':['eastward velocity','zonal current'],'v':['northward velocity','meridional current']}
    self.varnames={} #{'salt':'so','temp':'thetao','ssh':'zos'}
    self.fname_parts={}
	
    for intern in self.internal_names: 
        long=long_names[intern]
        for key in self.ncs[intern].variables.keys():
            try: 
                if self.ncs[intern][key].long_name.lower() in long:
                    #if ('ssh' in intern ) or ('depth' in self.ncs[intern][key].dimensions):
                    self.varnames[intern]=key
            except:
                print('no match for'+ key)
                pass
        file=self.files[intern]
        file=file[file.rfind('/')+1:]
	########################### End get filenames
	
    try:
	    ilon0=find(self.sv['lon'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(self.sv['lon'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(self.sv['lon'][:])
    try:
	    ilat0=find(self.sv['lat'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(self.sv['lat'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(self.sv['lat'][:])


    self.latslice=slice(ilat0,ilat1)
    self.lonslice=slice(ilon0,ilon1)


    self.lon = self.sv['lon'][self.lonslice]
    self.lat = self.sv['lat'][self.latslice]
    self.d = -self.sv['depth'][:] # change sigma depth to minus to be in line with schism
    self.time = self.sv['time'][:]
    self.tidx = 0
    self.ssh = self.sshv[self.varnames['ssh']][:,self.latslice,self.lonslice]
    self.s = self.sv[self.varnames['salt']][:,:,self.latslice,self.lonslice]
    self.t = self.tv[self.varnames['temp']][:,:,self.latslice,self.lonslice]-273.15*np.float ('k' in  self.ncs['temp'][self.varnames['temp']].units.lower()) # into celcius
    self.u = self.uv[self.varnames['u']][:,:,self.latslice,self.lonslice]
    self.v = self.uv[self.varnames['v']][:,:,self.latslice,self.lonslice]
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    
	
    self.mask_tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.mask=self.t.mask.flatten()

    vlon = self.lon2[np.where(self.ssh.mask[self.tidx]==False)]
    vlat = self.lat2[np.where(self.ssh.mask[self.tidx]==False)]
    self.ssh_tree = cKDTree(list(zip(vlon,vlat)))
    self.ssh_var = self.ssh[self.tidx][np.where(self.ssh.mask[self.tidx]==False)].flatten()
    
    self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
    vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
    vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
    vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
    self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
    
    vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
    vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
    vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
    self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
    
    vlon3 = self.lon3[np.where(self.u.mask[self.tidx]==False)]
    vlat3 = self.lat3[np.where(self.u.mask[self.tidx]==False)]
    vd3 = self.d3[np.where(self.u.mask[self.tidx]==False)]
    self.u_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
    
    
    vlon3 = self.lon3[np.where(self.v.mask[self.tidx]==False)]
    vlat3 = self.lat3[np.where(self.v.mask[self.tidx]==False)]
    vd3 = self.d3[np.where(self.u.mask[self.tidx]==False)]
    self.v_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))
    
    self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
    self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()
    self.u_var = self.u[self.tidx][np.where(self.u.mask[self.tidx]==False)].flatten()
    self.v_var = self.v[self.tidx][np.where(self.v.mask[self.tidx]==False)].flatten()	
	

		
  def ismasked(self,nodelon,nodelat):
    "check wether nextneighbour index is masked"	
    return self.mask[self.mask_tree.query((nodelon,nodelat))[1]]
    
  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))
    u = np.zeros((len(depths),))
    v = np.zeros((len(depths),))
    ssh=np.zeros(1)	
	
    dist,inds = self.ssh_tree.query((nodelon,nodelat),k=4)
    w = 1 / dist
    ssh[0] = np.sum(w*self.ssh_var[inds],axis=0) / np.sum(w,axis=0)
	

    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)
		
        dist,inds = self.u_tree.query((nodelon*100.,nodelat*100.,ndepth))		
        w = 1 / dist
        u[bidx-1+ik] = np.sum(w*self.u_var[inds],axis=0) / np.sum(w,axis=0)
        dist,inds = self.v_tree.query((nodelon*100.,nodelat*100.,ndepth))		
        w = 1 / dist
        v[bidx-1+ik] = np.sum(w*self.v_var[inds],axis=0) / np.sum(w,axis=0)		

    return (t,s,u,v,ssh)		

  def update(self,SSHfile=None,Sfile=None,Tfile=None,UVfile=None):
    """update netcdf ids for data acces"""
    self.tnc = netCDF4.Dataset(Tfile)
    self.tv = self.tnc.variables
    self.snc = netCDF4.Dataset(Sfile)
    self.sv = self.snc.variables
    self.sshnc = netCDF4.Dataset(SSHfile)
    self.sshv = self.sshnc.variables
    self.uvnc = netCDF4.Dataset(UVfile)
    self.uv = self.uvnc.variables
	
    self.time = self.sv['time'][:]
    self.ssh = self.sshv[self.varnames['ssh']][:,self.latslice,self.lonslice]
    self.s = self.sv[self.varnames['salt']][:,:,self.latslice,self.lonslice]
    self.t = self.tv[self.varnames['temp']][:,:,self.latslice,self.lonslice]-273.15*np.float ('k' in  self.ncs['temp'][self.varnames['temp']].units.lower()) # into celcius
    self.u = self.uv[self.varnames['u']][:,:,self.latslice,self.lonslice]
    self.v = self.uv[self.varnames['v']][:,:,self.latslice,self.lonslice]
    self.ssh_var = self.ssh[self.tidx][np.where(self.ssh.mask[self.tidx]==False)].flatten()
    self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
    self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()
    self.u_var = self.u[self.tidx][np.where(self.u.mask[self.tidx]==False)].flatten()
    self.v_var = self.v[self.tidx][np.where(self.v.mask[self.tidx]==False)].flatten()	

# program start ---------------------------------------------------------------

# init
os.chdir(setupdir)
nws = schism_setup()
file=files[0]
oa = cmems_one_file(file)#cmems(file)
#oa.update(files[1])
climut = utime(oa.timeunits) # Sourcedata Time reference
times=ut.date2num(climut.num2date(oa.time)) # times in netcdf file

# specified forcing boundaries
sigma=False
#vgridtype="SZ"
vgridtype="LSC2"
if sigma==True:
		# zcor will be used to get values set dry values to depth
		#file='/work/gg0028/g260114/RUNS/GermanBight/GB2k12_realisticRiver/GB2k12neu/combined/schout_1.nc'
		#nc=netCDF4.Dataset(file)
		#ncv=nc.variables
		# improvised depth
		# guess depth from equally spaced sigma layers
	with open('vgrid.in') as f:
		ivcor=int(f.readline().split()[0])
		line=f.readline().split()[:3]
		nws.znum=int(line[0])
		nws.hs=float(line[-1])
		slvls=[]
		#zlvls=[]
		if ivcor==2: # create vgrid
			passed=False
			for i,line in enumerate(f.readlines()):
				if 'S' in line and (passed == False):
					passed=True
				elif passed and ('theta' not in line):
					slvls.append([ int(line.split()[0]), float(line.split()[1]) ])
				elif ('S' not in line) and ('Z' not in line) and ('theta' not in line):
					zlvls.append([ float(line.split()[-1]) ])					

			slvls=np.ma.masked_array(np.asarray(slvls))
			# build vgrid
			nws.vgrid={node:slvls[:,1] for node in range(1,nws.nnodes+1)}
			nws.z=np.asarray(zlvls[:-1])
			nws.bidx={node:0 for node in range(1,nws.nnodes+1)} #np.zeros(nws.nnodes,'int')



if len(openbd_segs)>0:
	frcbdnodes=[]
	for seg in openbd_segs:
		frcbdnodes+=nws.bdy_segments[seg]
		bdyvgrid = asarray([nws.vgrid[ii].filled(-1.) for ii in frcbdnodes ])
else:
	frcbdnodes=nws.bdy_nodes
#bdyvgrid = asarray([nws.vgrid[ii].filled(-1.) for ii in nws.bdy_nodes ])

timesout=[]
# loop files
for ifile,file in enumerate(files[0:nfiles+1]):
	print ('loading file ' + file)

	tbdy=[]
	sbdy=[]
	ubdy=[]
	vbdy=[]
	sshbdy=[]
	dbdy=[]
	ibdy=[]

	# update variable values
	oa.update(file)

	#==============================================	update

	# interpolate to boundary
	times=ut.date2num(climut.num2date(oa.time))
	
	nws.z=np.asarray(zlvls[:-1])
	try:
		for itime,time in enumerate(times.astype('float32')):
	
			for i,inode in enumerate(frcbdnodes):
				if (i%100) == 0:
					print('	 interpolate i = %d'%i)
				bdylon = nws.londict[inode]
				bdylat = nws.latdict[inode]
				#depths = nws.vgrid[inode].filled(-1)*nws.depthsdict[inode]
				
				if vgridtype != "SZ" or (len(nws.z)==0):
					depths = nws.vgrid[inode].filled(-1)*nws.depthsdict[inode]
				else:
					if np.abs(nws.depthsdict[inode]) < nws.hs: # shallow case sigm from local depth to 0
						depths = np.concatenate((nws.z[:,0],nws.vgrid[inode].filled(-1)*nws.depthsdict[inode]))
					else: # combine zlevels 
						depths = np.concatenate((nws.z[:,0],nws.vgrid[inode].filled(-1)*(-nws.hs)))

				
				t,s,u,v,ssh = oa.interpolate(depths,bdylon,bdylat,bidx=1)
				tbdy.append(t)
				sbdy.append(s)
				ubdy.append(u)
				vbdy.append(v)
				sshbdy.append(ssh)
				dbdy.append(depths)
				ibdy.append(i*ones(depths.shape))


			# append 3D array
			timesout.append(time)	
			if (ifile==0 and itime==0):
				tbdy_time=tbdy
				sbdy_time=sbdy
				ubdy_time=ubdy
				vbdy_time=vbdy
				sshbdy_time=sshbdy
			else:
				tbdy_time=np.dstack((tbdy_time,tbdy))
				sbdy_time=np.dstack((sbdy_time,sbdy))
				ubdy_time=np.dstack((ubdy_time,ubdy))
				vbdy_time=np.dstack((vbdy_time,vbdy))
				sshbdy_time=np.dstack((sshbdy_time,sshbdy))

	except:
		pass
print('done extracting')


# output - interpolated cmemsdata
nt=len(timesout)
m,n,p=shape(tbdy_time)

nt=len(timesq)
uvintp=np.zeros((nt,m,n,2))
tbdy_time_back=tbdy_time.copy()
sbdy_time_back=sbdy_time.copy()
sshbdy_time_back=sshbdy_time.copy()
ubdy_time_back=ubdy_time.copy()
vbdy_time_back=vbdy_time.copy()

# interpolate
if dointerp:
	tagstart=str(ut.num2date(timesq[0])).replace(' ','_')
	tagend=str(ut.num2date(timesq[-1])).replace(' ','_')
	timetag='from'+tagstart+'_to_'+tagend

	timesout=np.asarray(timesout)
	fintp=interp1d(timesout,sshbdy_time,axis=2)
	sshbdy_time=fintp(timesq)
	fintp=interp1d(timesout,tbdy_time,axis=2)
	tbdy_time=fintp(timesq)
	fintp=interp1d(timesout,sbdy_time,axis=2)
	sbdy_time=fintp(timesq)
	fintp=interp1d(timesout,ubdy_time,axis=2)
	ubdy_time=fintp(timesq)
	fintp=interp1d(timesout,vbdy_time,axis=2)
	vbdy_time=fintp(timesq)
else:
	timetag='test'
# output data	
timetag+='neu'
m,n,p=shape(tbdy_time)
nws.write_bdy_netcdf('elev2D'+timetag+'.th.nc',timesq,sshbdy_time.swapaxes(0,2).swapaxes(1,2).reshape(nt,m,1,1),frcbdnodes=frcbdnodes)
nws.write_bdy_netcdf('TEM_3D'+timetag+'.th.nc',timesq,tbdy_time.swapaxes(1,2).swapaxes(0,1).reshape(nt,m,n,1),frcbdnodes=frcbdnodes)
nws.write_bdy_netcdf('SAL_3D'+timetag+'.th.nc',timesq,sbdy_time.swapaxes(1,2).swapaxes(0,1).reshape(nt,m,n,1),frcbdnodes=frcbdnodes)

uvintp[:,:,:,0]=ubdy_time.swapaxes(1,2).swapaxes(0,1).reshape(nt,m,n)
uvintp[:,:,:,1]=vbdy_time.swapaxes(1,2).swapaxes(0,1).reshape(nt,m,n)
nws.write_bdy_netcdf('uv3D'+timetag+'.th.nc',timesq,uvintp,frcbdnodes=frcbdnodes)

print ('done writing extracted forcing data')
os.chdir(cwd)

#plt.figure()
plt.clf()
plt.pcolormesh(oa.lon,oa.lat,oa.t[0,0,:])
plt.colorbar()
vmin,vmax=oa.t[0,0,:].min(),oa.t[0,0,:].max()
plt.scatter(np.asarray(nws.lon)[np.asarray(nws.bdy_nodes)-1],np.asarray(nws.lat)[np.asarray(nws.bdy_nodes)-1],s=10,c=tbdy_time.swapaxes(1,2).swapaxes(0,1).reshape(nt,m,n,1)[-1,:,-1,0],vmin=vmin,vmax=vmax,edgecolors='w',linewidth=1)
plt.xlim((np.min(nws.lon),np.max(nws.lon)))
plt.ylim((np.min(nws.lat),np.max(nws.lat)))
plt.title(times[-1])
plt.savefig('last_temp.png')
