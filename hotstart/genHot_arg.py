# -*- coding: utf-8 -*-
"""
Created on Tue May 28 15:05:27 2019

@author: JacobB
 Generates hotstart.nc for schism from CMEMS AMM15, which can be combined in order of specified preference (related value of dictionary "prio"

call: python genHot.py
 
"""
import os
import sys
sys.path.append(os.environ['HOME']+'/code/python/')
import numpy as np
import netCDF4
from scipy.spatial import cKDTree
import os.path
import pickle
import matplotlib
#matplotlib.use('AGG')
#plt.ion()
from matplotlib import path
from matplotlib import pyplot as plt
from netCDF4 import Dataset

# own libraries
from schism import schism_setup


#rundir,'/Downloads/Download/download_cmems_GB/'
print('Argument List:', str(sys.argv))

if len(sys.argv) != 4:
	print('mistmatching number of commandline arguments: python genHot_arg <srcdir> <rundir> <date>')
	exit()

srcdir=sys.argv[1]
rundir=sys.argv[2]
date=sys.argv[3]
#################### Settings ########################################

# output setup - create a hot start for this setup
destination_dir=rundir
hotoutname='GB_hot_{:s}.nc'.format(date)    # name of hotstart.nc file which will be created as output 

os.chdir(rundir)
sigma=False   # True: Sigma  Flase: LSC2
nsigma=21

plothot=True

# use existent schism hotstart to interpolate for  the hotstart bein created here 
source_dir='/gpfs/work/jacobb/data/SETUPS/BlackSea/'
schism_hotstart_file_input='/gpfs/work/jacobb/data/SETUPS/BlackSea/blacksea_hot_jun2016.nc'#hotstart.nc0'
prio={'schism':3}  #  prio,value := priority of dataset to use in interpolaton. 0:= dont use this type. n=1,2,...N: if node of schism setup to interpolate to falls into domain
                     # of data set whith prio n use values from this data set; else continue look in dataset with prio n+1 .... i.e. use High priorities (= small values) for regional, higher resolution data,
sigma_in=False   # True: Sigma  Flase: LSC2
nsigma_in=21
                                                                                                                                                                 # lower priority for larger coarse data sets

# amm7 
#3 layer salinity and temoperature file
prio['amm7_3']=0 # use amm 7 3 layer with prio
fileS='/mnt/lustre01/work/gg0028/g260114/Data/Forcing/MyOcean2012neu/metoffice_foam1_amm7_NWS_PSAL_b20120101_hi20111231.nc'
fileT='/mnt/lustre01/work/gg0028/g260114/Data/Forcing/MyOcean2012neu/metoffice_foam1_amm7_NWS_TEMP_b20120101_hi20111231.nc'
# 24 z layer
prio['amm7_24']=0
fileS='/work/gg0028/g260099/AMM7/RAW_20170101/metoffice_foam1_amm7_NWS-HR_PSAL_b20160108_hi20160107.nc'
fileT='/work/gg0028/g260099/AMM7/RAW_20170101/metoffice_foam1_amm7_NWS-HR_TEMP_b20160108_hi20160107.nc'

# use cmems
cmemsfile='/work/gg0028/g260114/Data/Forcing/cmems_dayly/2007/cmems_global_0083deg_20070601.nc'

cmemsfile='/gpfs/work/jacobb/data/DATA/Forcing/mount_cmems2/cmems_global_0083deg_20120601.nc'
prio['cmems']=0

# use jansen
jannsen_file='/gpfs/work/jacobb/storage/SetupGeneration/1_hotstart/janssen_climatology.nc'
#/work/gg0028/g260114/SetupGeneration/1_hotstart/janssen_climatology.nc'
jannsen_month=6 # month to use from jannsen climatology
prio['jannsen']=0


prio['amm15']=1
fileS='{:s}metoffice_foam1_amm15_NWS_TEM_{:s}.nc'.format(srcdir,date)
fileT='{:s}metoffice_foam1_amm15_NWS_SAL_{:s}.nc'.format(srcdir,date)


prio['amm']=0  # replacement for amm15, amm7 #hopefully handles all
#fileS='/gpfs/work/jacobb/data/DATA/Download/motuclient-python-v1.8.4/blacksea/bs_mean_sal_20160601.nc'
#fileT='/gpfs/work/jacobb/data/DATA/Download/motuclient-python-v1.8.4/blacksea/bs_mean_tem_20160601.nc'

###############################################################################



# Classes for different hotstart data sources
class jannsen():

  def __init__(self,file,domain_rect,imon):
    tnc = netCDF4.Dataset(file)
    tv = tnc.variables
    snc = netCDF4.Dataset(file)
    sv = snc.variables

    imon=imon-1

    try:
	    ilon0=find(sv['lon'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(sv['lon'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(sv['lon'][:])
    try:
	    ilat0=find(sv['lat'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(sv['lat'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(sv['lat'][:])

    ilon0=max(ilon0,0)
    ilat0=max(ilat0,0)


    latslice=slice(ilat0,ilat1)
    lonslice=slice(ilon0,ilon1)

    self.lon = sv['lon'][lonslice]
    self.lat = sv['lat'][latslice]
    self.d = -sv['zax'][:]
    self.time = sv['time'][:]
    self.tidx = 0

    self.s = np.reshape(sv['salt'][imon,:,latslice,lonslice], (1,)+ np.shape(sv['salt'][imon,:,latslice,lonslice]))
    self.t = np.reshape(tv['temp'][imon,:,latslice,lonslice], (1,)+ np.shape(sv['salt'][imon,:,latslice,lonslice]))
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)

    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)



class cmems():

  def __init__(self,file,domain_rect):
    tnc = netCDF4.Dataset(file)
    tv = tnc.variables
    snc = netCDF4.Dataset(file)
    sv = snc.variables

    try:
	    ilon0=find(sv['longitude'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(sv['longitude'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(sv['longitude'][:])
    try:
	    ilat0=find(sv['latitude'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(sv['latitude'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(sv['latitude'][:])


    latslice=slice(ilat0,ilat1)
    lonslice=slice(ilon0,ilon1)


    self.lon = sv['longitude'][lonslice]
    self.lat = sv['latitude'][latslice]
    self.d = -sv['depth'][:]
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv['so'][:,:,latslice,lonslice]
    self.t = tv['thetao'][:,:,latslice,lonslice]
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)

    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)


class woa():

  def __init__(self):
    tnc = netCDF4.Dataset(file)
    tv = tnc.variables
    snc = netCDF4.Dataset(file)
    sv = snc.variables
    latslice=slice(528,624)
    lonslice=slice(648,860)
    self.lon = sv['lon'][lonslice]
    self.lat = sv['lat'][latslice]
    self.d = -sv['depth'][:]
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv['s_mn'][:,:,latslice,lonslice]
    self.t = tv['t_mn'][:,:,latslice,lonslice]
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)

    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)


class schismhot():

  def __init__(self,file,hotsetup):
    tnc = netCDF4.Dataset(file)
    tv = tnc.variables
    snc = netCDF4.Dataset(file)
    sv = snc.variables

    self.lon = np.asarray(hotsetup.lon)
    self.lat = np.asarray(hotsetup.lat)
    self.vgrid = hotsetup.vgrid
    self.d = hotsetup.depths
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv['tr_nd0'][:,:,1] # node vert tracer
    self.t = tv['tr_nd0'][:,:,0] # node vert tracer



  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))

    inn=np.argmin( (nodelon-self.lon)**2 + (nodelat-self.lat)**2) # horizontal next neighbour

    d=self.vgrid[inn+1]*self.d[inn]
#   d=bsc.vgrid[inn+1]*bsc.d[inn]
    #print(nodelon,nodelat,d)
    #vertical linear 

    for ik,dep in enumerate(depths[bidx-1:]):
      if np.sum(d<=dep):
        ibelow=np.where(d<=dep)[-1][-1]
        iabove=ibelow+1 #- (ibelow+1==len(d))
        dists=d[ibelow:iabove+1]-dep
        if dists[0]==0:
          w=np.double(np.logical_not(dists))
        else:
          w=1/np.abs(dists)/np.sum(1/np.abs(dists))
        s[ik]=np.sum(w*self.s[inn,ibelow:iabove+1])
        t[ik]=np.sum(w*self.t[inn,ibelow:iabove+1])
      else:
        s[ik]=np.nan
        t[ik]=np.nan


    return (t,s)
	
class amm():
  # sigma lvels surface, midwater column, near bed	
  def __init__(self,Tfile,Sfile,domain_rect,SSHfile=None,UVfile=None):
    tnc = netCDF4.Dataset(Tfile)
    tv = tnc.variables
    snc = netCDF4.Dataset(Sfile)
    sv = snc.variables

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
    #long_names=['Sea Water Salinity','Sea Water Potential Temperature','Sea surface height above geoid','Eastward Current Velocity in the Water Column','Northward Current Velocity in the Water Column']
    long_names={'salt':'salinity','temp':'temperature','ssh':'Sea surface height','u':'eastward velocity','v':'northward velocity'}
    self.varnames={} #{'salt':'so','temp':'thetao','ssh':'zos'}
    self.fname_parts={}
	
    for intern in self.internal_names: 
        long=long_names[intern]
        print(intern)
        for key in self.ncs[intern].variables.keys():
            try: 
                if long in self.ncs[intern][key].long_name.lower():
                    if ('ssh' in intern ) or ('depth' in self.ncs[intern][key].dimensions):
                        self.varnames[intern]=key
            except:
                pass
        file=self.files[intern]
        file=file[file.rfind('/')+1:]
	########################### End get filenames
	
    try:
	    ilon0=find(sv['lon'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(sv['lon'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(sv['lon'][:])
    try:
	    ilat0=find(sv['lat'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(sv['lat'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(sv['lat'][:])


    latslice=slice(ilat0,ilat1)
    lonslice=slice(ilon0,ilon1)


    self.lon = sv['lon'][lonslice]
    self.lat = sv['lat'][latslice]
    self.d = -sv['depth'][:] # change sigma depth to minus to be in line with schism
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv[self.varnames['salt']][:,:,latslice,lonslice]
    self.t = tv[self.varnames['temp']][:,:,latslice,lonslice]-273.15*np.float ('k' in  self.ncs['temp'][self.varnames['temp']].units.lower()) # into celcius
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    
	
    self.mask_tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.mask=self.t.mask.flatten()
	
    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3))) # 3 nextneighbour indey

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

		
  def ismasked(self,nodelon,nodelat):
    "check wether nextneighbour index is masked"	
    return self.mask[self.mask_tree.query((nodelon,nodelat))[1]]
    
  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)		
	
class amm7_3l():
  # sigma lvels surface, midwater column, near bed	
  def __init__(self,fileT,fileS,domain_rect):
    tnc = netCDF4.Dataset(fileT)
    tv = tnc.variables
    snc = netCDF4.Dataset(fileS)
    sv = snc.variables

    try:
	    ilon0=find(sv['lon'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(sv['lon'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(sv['lon'][:])
    try:
	    ilat0=find(sv['lat'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(sv['lat'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(sv['lat'][:])


    latslice=slice(ilat0,ilat1)
    lonslice=slice(ilon0,ilon1)


    self.lon = sv['lon'][lonslice]
    self.lat = sv['lat'][latslice]
    self.d = -sv['depth'][:] # change sigma depth to minus to be in line with schism
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv['vosaline'][:,:,latslice,lonslice]
    self.t = tv['votemper'][:,:,latslice,lonslice]-273.15 # into celcius
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    
    self.mask_tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.mask=self.t.mask.flatten()
	
    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3))) # 3 nextneighbour indey

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

		
  def ismasked(self,nodelon,nodelat):
    "check wether nextneighbour index is masked"	
    return self.mask[self.mask_tree.query((nodelon,nodelat))[1]]
    
  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)		
	
class amm7:
  # sigma lvels surface, midwater column, near bed	
  def __init__(self,fileT,fileS,domain_rect):
    tnc = netCDF4.Dataset(fileT)
    tv = tnc.variables
    snc = netCDF4.Dataset(fileS)
    sv = snc.variables

    try:
	    ilon0=find(sv['lon'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(sv['lon'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(sv['lon'][:])
    try:
	    ilat0=find(sv['lat'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(sv['lat'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(sv['lat'][:])


    latslice=slice(ilat0,ilat1)
    lonslice=slice(ilon0,ilon1)


    self.lon = sv['lon'][lonslice]
    self.lat = sv['lat'][latslice]
    self.d = -sv['depth'][:] # change sigma depth to minus to be in line with schism
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv['vosaline'][:,:,latslice,lonslice]
    self.t = tv['votemper'][:,:,latslice,lonslice]-273.15 # into celcius
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    
    self.mask_tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.mask=self.t.mask.flatten()
	
    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3))) # 3 nextneighbour indey

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

		
  def ismasked(self,nodelon,nodelat):
    "check wether nextneighbour index is masked"	
    return self.mask[self.mask_tree.query((nodelon,nodelat))[1]]
    
  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)			
	
	
	
	
	
class amm15():
  # sigma lvels surface, midwater column, near bed	
  def __init__(self,fileT,fileS,domain_rect):
    tnc = netCDF4.Dataset(fileT)
    tv = tnc.variables
    snc = netCDF4.Dataset(fileS)
    sv = snc.variables

    try:
	    ilon0=find(sv['lon'][:] >=   domain_rect[0])[0]-1
    except:
	    ilon0=0
    try:
	    ilon1=find(sv['lon'][:] >=   domain_rect[1])[0]+1
    except:
	    ilon1=len(sv['lon'][:])
    try:
	    ilat0=find(sv['lat'][:] >=   domain_rect[2])[0]-1
    except:
	    ilat0=1
    try:
	    ilat1=find(sv['lat'][:] >=   domain_rect[3])[0]+1
    except:
	    ilat1=len(sv['lat'][:])


    latslice=slice(ilat0,ilat1)
    lonslice=slice(ilon0,ilon1)


    self.lon = sv['lon'][lonslice]
    self.lat = sv['lat'][latslice]
    self.d = -sv['depth'][:] # change sigma depth to minus to be in line with schism
    self.time = sv['time'][:]
    self.tidx = 0
    self.s = sv['so'][:,:,latslice,lonslice]
    self.t = tv['thetao'][:,:,latslice,lonslice]#-273.15 # into celcius
    self.lon2,self.lat2 = np.meshgrid(self.lon,self.lat)
    
	
    self.mask_tree = cKDTree(list(zip(self.lon2.flatten(),self.lat2.flatten())))
    self.mask=self.t.mask.flatten()
	
    self.use_depth_slices=False

    if not(self.use_depth_slices):
      self.d3,self.lat3,self.lon3 = np.meshgrid(self.d,self.lat*100.,self.lon*100.,indexing='ij')
      vlon3 = self.lon3[np.where(self.s.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.s.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.s.mask[self.tidx]==False)]
      self.s_tree=cKDTree(list(zip(vlon3,vlat3,vd3))) # 3 nextneighbour indey

      vlon3 = self.lon3[np.where(self.t.mask[self.tidx]==False)]
      vlat3 = self.lat3[np.where(self.t.mask[self.tidx]==False)]
      vd3 = self.d3[np.where(self.t.mask[self.tidx]==False)]
      self.t_tree=cKDTree(list(zip(vlon3,vlat3,vd3)))

      self.s_var = self.s[self.tidx][np.where(self.s.mask[self.tidx]==False)].flatten()
      self.t_var = self.t[self.tidx][np.where(self.t.mask[self.tidx]==False)].flatten()

    else:
      #build trees
      self.s_tree={}
      self.t_tree={}
      self.svar={}
      self.tvar={}
      for ik,d in enumerate(self.d):
        print('  build trees for depth %0.2f'%d)
        vlon = self.lon2[np.where(self.s.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.s.mask[self.tidx,ik]==False)]
        self.s_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        vlon = self.lon2[np.where(self.t.mask[self.tidx,ik]==False)]
        vlat = self.lat2[np.where(self.t.mask[self.tidx,ik]==False)]
        self.t_tree[ik] = cKDTree(list(zip(vlon,vlat)))
        self.svar[ik] = self.s[self.tidx,ik][np.where(self.s.mask[self.tidx,ik]==False)].flatten()
        self.tvar[ik] = self.t[self.tidx,ik][np.where(self.t.mask[self.tidx,ik]==False)].flatten()

		
  def ismasked(self,nodelon,nodelat):
    "check wether nextneighbour index is masked"	
    return self.mask[self.mask_tree.query((nodelon,nodelat))[1]]
    
  def interpolate(self,depths,nodelon,nodelat,bidx=1):
    # start
    t = np.zeros((len(depths),))
    s = np.zeros((len(depths),))


    for ik,ndepth in enumerate(depths[bidx-1:]):
      # find vertical layer in climatology
      if self.use_depth_slices:
        didx = np.abs(self.d - ndepth).argmin()

      #didx = int(np.where(ddiff==ddiff.min())[0][0])
      #vlon = self.lon2[np.where(self.s.mask[self.tidx,didx]==False)]
      #vlat = self.lat2[np.where(self.s.mask[self.tidx,didx]==False)]
      #tree = cKDTree(zip(vlon,vlat))
    
      #svar = self.s[self.tidx,didx][np.where(self.s.mask[self.tidx,didx]==False)].flatten()
      #tvar = self.t[self.tidx,didx][np.where(self.t.mask[self.tidx,didx]==False)].flatten()

        dist,inds = self.s_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.svar[didx][inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree[didx].query((nodelon,nodelat),k=4)
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.tvar[didx][inds],axis=0) / np.sum(w,axis=0)
      else:
        dist,inds = self.s_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        s[bidx-1+ik] = np.sum(w*self.s_var[inds],axis=0) / np.sum(w,axis=0)

        dist,inds = self.t_tree.query((nodelon*100.,nodelat*100.,ndepth))
        w = 1 / dist
        t[bidx-1+ik] = np.sum(w*self.t_var[inds],axis=0) / np.sum(w,axis=0)

    return (t,s)		
	
#-----------------------------------------------------------------------------------


#-------------------------- program start -----------------------------------

# setup to create hotstart for
os.chdir(destination_dir)
stp = schism_setup()
domain_rect=[min(stp.lon), max(stp.lon), min(stp.lat), max(stp.lat)]  # rectangular domain boundingbox to subselect structured grid models for hotstart

# triangle only grid -> split quads
faces2=[]
for nr,elem in enumerate(stp.nv):
	if len(elem)==3:
		faces2.append(elem[:3])
	else: # split quad into tris
		faces2.append([elem[0],elem[1],elem[2]])
		faces2.append([elem[0],elem[2],elem[3]])
stp.faces=np.array(faces2)-1                  

## !!!! Depths has to be downward negative
if sigma==True:
	# zcor will be used to get values set dry values to depth
	#file='/work/gg0028/g260114/RUNS/GermanBight/GB2k12_realisticRiver/GB2k12neu/combined/schout_1.nc'
	#nc=netCDF4.Dataset(file)
	#ncv=nc.variables
	# improvised depth
	# guess depth from equally spaced sigma layers
    with open('vgrid.in') as f:
        ivcor=int(f.readline().split()[0])
        stp.znum=int(f.readline().split()[0])
        lvls=[]
        if ivcor==2: # create vgrid
                passed=0
                for i,line in enumerate(f.readlines()):
                        if line.split()[0]=='S':
                                passed=i+1
                        if (i>passed>0):
                                lvls.append([ int(line.split()[0]), float(line.split()[1]) ])
                lvls=np.ma.masked_array(np.asarray(lvls))
                # build vgrid
                stp.vgrid={node:lvls[:,1] for node in range(1,stp.nnodes+1)}
                stp.bidx={node:0 for node in range(1,stp.nnodes+1)} #np.zeros(stp.nnodes,'int')    
    
#d=-np.asarray(stp.depths)
#d[d>0]=0
#zcor=np.zeros((stp.nnodes,nsigma))
#for i in range(stp.nnodes):
#	zcor[i,:]=np.linspace(d[i],0,nsigma)	
#stp.vgrid=zcor
#stp.vgrid=zcor
#stp.znum=nsigma
#stp.bidx=np.zeros(stp.nnodes,'int')    
# load hotstarts
# schism setup to create hotstart from

srcdata=[]
uselist=[]
if prio['schism']>0:
    os.chdir(source_dir)
    stpin = schism_setup()

    ## !!!! Depths has to be downward negative
    if sigma_in==True:
        	# zcor will be used to get values set dry values to depth
        	#file='/work/gg0028/g260114/RUNS/GermanBight/GB2k12_realisticRiver/GB2k12neu/combined/schout_1.nc'
        	#nc=netCDF4.Dataset(file)
        	#ncv=nc.variables:w
        	# improvised depth
        	# guess depth from qually spaced sigma layers

        	with open('vgrid.in') as f:
        		ivcor=int(f.readline().split()[0])
        		stpin.znum=int(f.readline().split()[0])
        		lvls=[]
        		if ivcor==2: # create vgrid
        			passed=0
        			for i,line in enumerate(f.readlines()):
        				print(line)                      
        				if line.split()[0]=='S':
        					passed=i+1
        				if (i>passed>0):
        					lvls.append([ int(line.split()[0]), float(line.split()[1]) ])
        			lvls=np.ma.masked_array(np.asarray(lvls))
        			# build vgrid
        			stpin.vgrid={node:lvls[:,1] for node in range(1,stpin.nnodes+1)}
        			stpin.bidx={node:0 for node in range(1,stpin.nnodes+1)} #np.
            
        	#d=-np.asarray(stp.depths)
        	#d[d>0]=0
        	#zcor=np.zeros((stp.nnodes,nsigma))
        	#for i in range(stp.nnodes):
        	#	zcor[i,:]=np.linspace(d[i],0,nsigma)	
        	#stpin.vgrid=zcor
        	#stpin.vgrid=zcor
        	#stpin.znum=nsigma_in
        	#stpin.bidx=np.zeros(stpin.nnodes,'int')    
    schismhot_in=schismhot(schism_hotstart_file_input,stpin) #schism
    srcdata+=[schismhot_in]
    uselist+=['schism']
if prio['cmems']>0:
    cmem = cmems(cmemsfile,domain_rect)# cmems 
    srcdata+=[cmem]
    uselist+=['cmems']            
    
if prio['amm7_3']>0:    
    amm7_3=amm7_3l(fileT,fileS,domain_rect)# amm7 3 layer
    srcdata+=[amm7_3]
    uselist+=['amm7_3']
if prio['amm7_24']>0:    
    amm7=amm7(fileT,fileS,domain_rect)# amm7 3 layer
    srcdata+=[amm7]
    uselist+=['amm7_24']
    
if prio['jannsen']>0:        
    jan = jannsen(jannsen_file,domain_rect,imon=6)# jansen
    srcdata+=[jan]
    uselist+=['jannsen']
	
if prio['amm15']>0:    
    amm15=amm15(fileT,fileS,domain_rect)# amm7 3 layer
    srcdata+=[amm15]
    uselist+=['amm15']

if prio['amm']>0:    
    amm=amm(fileT,fileS,domain_rect)
    srcdata+=[amm]
    uselist+=['amm']
	
	
#poly=np.loadtxt('/work/gg0028/g260114/RUNS/Europe3.0/hotstarts/nbs_bound.xy')[:,1:]  # limit for jannsen usage area
#p=path.Path([(poly[i,0],poly[i,1]) for i in range(poly.shape[0])])


# sort acces to input data by decreasing priority
prios=[]
for key in uselist:
    prios.append(prio[key])
srcdata=list(np.array(srcdata)[np.argsort(prios)])
uselist=list(np.array(uselist)[np.argsort(prios)])
    
# area box
for src in srcdata:
    src.lonR=[min(src.lon), max(src.lon)]
    src.latR=[min(src.lat), max(src.lat)]


# case 1 new hotstart	
os.chdir(destination_dir)
#######  interpolate source data for hotstart   #################################
# write t,s on nodes
s = {}
t = {}
useind=np.zeros(stp.nnodes,int) 


#i=stp.inodes[ind]
#nodelon=stp.lon
#nodelat=stp.lat
#d=stp.depths[ind]

# create t,s fields:
for i,nodelon,nodelat,d in zip(stp.inodes,stp.lon,stp.lat,stp.depths):
    if (i%10000) == 0:
        print('  interpolate i = %d'%i)

    depths = stp.vgrid[i].filled(-1)*d
    bidx = stp.bidx[i]
    # interpolate from domains with decreasing priority
    #if useind[i-1]==2:
    #    src=srcdata[-1]
    #    t[i],s[i] = src.interpolate(depths,nodelon,nodelat,bidx=1)	
	#	schismhot_in.interpolate(depths,nodelon,nodelat,bidx=1)	
    #else:
    #    break		
    for isrc,src in enumerate(srcdata):
        if (src.lonR[0] <= nodelon and nodelon <= src.lonR[1] and src.latR[0] <= nodelat and nodelat <= src.latR[1]):
            t[i],s[i] = src.interpolate(depths,nodelon,nodelat,bidx=1)
            useind[i-1]=isrc+1
            break     

			
if 0:			
	hot0=srcdata[1]
	t=hot0.t*1.0
	s=hot0.s*1.0
	# use orignal hotstart and ad other data
	# create t,s fields:
	src=srcdata[0]
	isrc=0
	for i,nodelon,nodelat,d in zip(stp.inodes,stp.lon,stp.lat,stp.depths):
		if (i%10000) == 0:
			print('  interpolate i = %d'%i)

		depths = stp.vgrid[i].filled(-1)*d
		bidx = stp.bidx[i]
	   
		# interpolate from domains with decreasing priority
		#if useind[i-1]==2:
		#    src=srcdata[-1]
		#    t[i],s[i] = src.interpolate(depths,nodelon,nodelat,bidx=1)	
		#	schismhot_in.interpolate(depths,nodelon,nodelat,bidx=1)	
		#else:
		#    break		
		#for isrc,src in enumerate(srcdata):
		
		# check if larger area is masked - to limit interpolation
		inds=src.mask_tree.query((nodelon,nodelat),25)[1] 
		if (src.mask[inds]==False).sum()>0:
			if (src.lonR[0] <= nodelon and nodelon <= src.lonR[1] and src.latR[0] <= nodelat and nodelat <= src.latR[1]):
				t[i-1,:],s[i-1,:] = src.interpolate(depths,nodelon,nodelat,bidx=1)
				useind[i-1]=isrc+1
				#break    

				
				
if 0:				
	# overinterpolate river values with schism hotstart			
	# add mask for ocean land mask to avoid obscurely large interpolations distances
	# and schism hotsart exist as alternative (assuming data is close enouh)
	# this mostly concernc interpolation in rivers when points are in boundinbox of coarse data
	# but also within its land mask

	if prio['schism']!=0:
		i_overinterpolate=np.zeros(stp.nnodes,bool)		
		for isrc,src in enumerate(srcdata):
			if str(type(src))!="<class '__main__.schismhot'>":			
				src.masktree=cKDTree(list(zip(src.lon2.flatten(),src.lat2.flatten())))
				src.mask2d=src.s[0,0].mask
				src.mask1d=src.mask2d.flatten()
				src.mask1d[inds].sum() 
				inds=src.masktree.query(list(zip(stp.lon,stp.lat)),4)[1]
				i_overinterpolate+=src.mask1d[inds].sum(axis=1)==4		
			else:
				srcnr=isrc
				schismsrc=src	
		print('over interpolate river nodes')		
		inds=np.asarray(np.where(i_overinterpolate))		
		src=schismsrc		
		#for i,nodelon,nodelat,d in zip(list(np.asarray(stp.inodes)[inds]),list(np.asarray(stp.lon)[inds]),list(np.asarray(stp.lat)[inds]),list(np.asarray(stp.depths)[inds])):
		for j in inds[0]:
			i,nodelon,nodelat,d =stp.inodes[j],stp.lon[j],stp.lat[j],stp.depths[j]
			if (j%10000) == 0:
				print('  interpolate i = %d'%i)

			depths = stp.vgrid[i].filled(-1)*d
			bidx = stp.bidx[i]
		   
			# interpolate from domains with decreasing priority
			t[i],s[i] = src.interpolate(depths,nodelon,nodelat,bidx=1)
			useind[i-1]=srcnr+1

				
				
	##########################################################################################
t=np.asarray(list(t.values()))
s=np.asarray(list(s.values()))


# replace nan with value above
for i in range(len(s)):
	naninds=np.isnan(t[i])
	ibtm=np.sum(naninds)
	t[i][naninds]=t[i][ibtm]
	s[i][naninds]=s[i][ibtm]
######################################################
	
############## write hotstart
tr_nd=np.stack((t,s),axis=2) 
stp.write_hotstart(tr_nd,filename=hotoutname,time=0.0,iths=0,ifile=0,elev=0.0)
##############################################################


############ make plots of hotstart ##########################################
def plotschism(s,nodedata,cmap=plt.cm.jet):
	#plt.figure()
	ph=plt.tripcolor(s.lon,s.lat,s.faces,np.asarray(nodedata),shading='flat',cmap=cmap)
	ch=plt.colorbar()
	return ch
	ph,#ch.set_label('ssh / m')
	#show()

if plothot:
	ch=plotschism(stp,useind)
	ch.set_ticks(1+np.arange(len(uselist)))
	ch.set_ticklabels(uselist)
	plt.savefig('datasource',dpi=300)
	plotschism(stp,stp.depths)
	plt.plot(np.asarray(src.lonR)[[0,1,1,0,0]],np.asarray(src.latR)[[0,0,1,1,0]],'r')

	suffix=['surface', 'bottom' ]
	inds=[-1, 0]

	# plot surf bottmom layer
	for i, lab in zip(inds,suffix): 
		ph,ch=stp.plotAtnodes(s[:,i])
		ch.set_label('salinity g/kg')
		plt.title(lab)
		plt.savefig('salinity_'+lab,dpi=600)
		plt.close()

		#ch=schism_plotAtnodes(stp,t[:,i])
		ph,ch=stp.plotAtnodes(t[:,i])
		ch.set_label('t / degc')
		plt.title(lab)
		plt.savefig('temp_'+lab,dpi=600)
		#plt.show()
		plt.close()


###########

#plt.pcolormesh(srcdata[0].lon,srcdata[0].lat,srcdata[0].t[0,0,:,:],cmap=plt.cm.jet)
plt.colorbar()

