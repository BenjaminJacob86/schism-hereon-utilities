#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Different clases for Model and data output for loading and interpolating
"""

__author__ = "Benjamin Jacob, HZG"
__version__ = "1.0.0"


##
# 
import numpy as np
import netCDF4
from netCDF4 import Dataset,MFDataset
import datetime as dt
from scipy.spatial import cKDTree
from matplotlib import path
import datetime as dt
import glob
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from cftime import utime
import xarray as xr
######################## Observation Data CLasses ####################################################################


###### BSH ########
class bsh_spec:
	""" Read BSH Wave spectra *.spec files with format
	Line 1:  Buoy-Type Station Datum+Time
	Line 2:  up to 12 Parameters
    01 Peakfrequency (Hz),
	02 Energy Maximum (m*m*s),
	03 Dir at Peak (degN), 90 degrees is from east
	04 Spread at Peak (deg),
	05 Hs (m),
	06 Tp (s),
	07 Tm1 (s),
	08 Tm2 (s),
	09 Tm-1 (s),
	10 Nspc ,
	11 Hmax (m),
	12 T(Hmax) (s)
	% Line 3 - Nspc+3: f (Hz), e (m*m*s), dir (degN), spr (deg)
	
	Input:
	filename of spectral files as str
	or chronologically ordered list of files of bouy data:
	(bsh_spec(file=list(np.sort(glob.glob('<path>/*<station>*.spec')))   )
	
	Output class with fields:
	.name := station name as specified in file e.g. FN3
	.B.type :=  Buoy-Type - wave rider etc.
	.dates := list of dates of measurements
	.integral_header := header of integral parameters sotred in respective according numpy array (.integral_values)
	.integral_values := time times header item numpy array of integra parameters e.g. 
						bouy.integral_values[:,bouy.integral_header.index('Hs')] containing HS time series
						according to bouy.dates	
	.spectal_header := ['f', 'e', 'dir', 'spr'] header of spectral parmeters within .spectra
	.spectra := list with spectra occoring to .dates each item being an numpy array
				with spectral freuqencies in rows anbd colums according to spectal_header
	"""
	# not always all frequencies in spectra spectra
	def __init__(self,file):
	
		if type(file) == list:
			files=file
			file=files[0]
			nfiles=len(files)
		else:
			nfiles=1
			
		with open(file) as f:
			names=['Peakfrequency','Energy_Maximum','Dir_at_Peak','Spread_at_Peak',
			'Hs','Tp','Tm1','Tm2','Tm-1','Nspc','Hmax','T_Hmax']
			self.integral_header=names
			self.spectal_header=[ 'f', 'e' , 'dir', 'spr']
			lines=f.readlines()
			self.type,self.name,date=lines[0].rstrip('\n').split(' ')
			self.dates=[dt.datetime.strptime(date,'%Y%m%d%H%M')]
			self.i_nspec=names.index('Nspc')
			vals=[ float(val) for val in lines[1].rstrip('\n').split()]
			self.integral_values=[vals]
			count=2
			nspec=int(vals[self.i_nspec]) 
			
			spectrum=[[ float(item) for item in lines[count].rstrip('\n').split()]]
			for i in range(count+1,count+nspec):
				spectrum.append([ float(item) for item in lines[i].rstrip('\n').split()])
			self.spectra=[np.asarray(spectrum)]
			
			count+=nspec
			
			self.append_data(count,lines)

		if nfiles > 1:
			for file in files[1:]:
				#print(file)
				count=0
				try:
					with open(file) as f:
						self.append_data(count,f.readlines())
				except:
					print('error loading file '+ file)
					pass
		
		# convert lists to numpy arrays	
		self.integral_values=np.asarray(self.integral_values)
		
		
	def append_data(self,count,lines):
		while count < len(lines):
			date=lines[count].rstrip('\n').split(' ')[-1]
			self.dates.append(dt.datetime.strptime(date,'%Y%m%d%H%M'))
			
			count+=1 # read integral parameters
			self.integral_values.append([ float(val) for val in lines[count].rstrip('\n').split()])
			nspec=int(self.integral_values[-1][self.i_nspec]) 
			
			count+=1 # read spectra
			spectrum=[[ float(item) for item in lines[count].rstrip('\n').split()]]
			for i in range(count+1,count+nspec):
				spectrum.append([ float(item) for item in lines[i].rstrip('\n').split()])
			self.spectra.append(np.asarray([spectrum]))
			count+=nspec
		
	def get_parameter(self,parameter='Hs'):
		""" get buoy data """
		return self.integral_values[:,self.integral_header.index(parameter)]
		


class sentinel_dataset:
	""" Sentinel Dataset. Scan top directory and list dates for sat A and B
	get list of files an derives dates from folder names ala
	/S1A_IW_OCN__2SDV_20171230T172545_20171230T172619_019935_021F10_C616.SAFE"""

	def __init__(self,sentinel_dir,prefix="S1A_IW_OCN__2SDV_",suffix=".SAFE",refdate=dt.datetime(2017,1,1)):

		# start and end date extracted from folder name
		self.directories=[]
		self.start_dates=[]
		self.end_dates=[]
		# time since ref date
		self.start_days=[]  
		self.end_days=[]
		self.center_days=[]
		self.refdate=refdate
	
		# extract dates
		for fname in glob.glob(sentinel_dir+"*"+suffix):

			# check sat 
			if fname.find(prefix)>0:
				print("checking "+ prefix + fname)
				self.directories.append(fname)
				start_date,end_date=fname[fname.find(prefix)+len(prefix):].split('_')[:2]
				print(start_date + ' ' + end_date)
				start_date=dt.datetime.strptime(start_date,'%Y%m%dt%H%M%S')
				end_date=dt.datetime.strptime(end_date,'%Y%m%dt%H%M%S')
				self.start_dates.append(start_date)
				self.end_dates.append(end_date)
				self.start_days.append((start_date-refdate).total_seconds()/86400)
				self.end_days.append((end_date-refdate).total_seconds()/86400)
				self.center_days.append((self.start_days[-1]+self.end_days[-1])/2)

		# sort chronologically
		sorted_inds=np.argsort(self.start_days)
		self.directories=np.asarray(self.directories)[sorted_inds]
		self.start_dates=np.asarray(self.start_dates)[sorted_inds]
		self.end_dates=np.asarray(self.end_dates)[sorted_inds]
		self.start_days=np.asarray(self.start_days)[sorted_inds]  
		self.end_days=np.asarray(self.end_days)[sorted_inds]
		self.center_days=np.asarray(self.center_days)[sorted_inds]

		# put into pandas 
		self.dataset=pd.Series([self.directories,self.start_dates,self.end_dates,self.start_days,
		self.end_days,self.center_days],index=['dir','date0','date1','day0','day1','dayc'])






class ostia:
	""" Read Ostia netCDF  sst  Dataset into class object.
	    Get horizontal time slab in self.slab querriny
	    via datetime object. 	 
	 """
	import numpy as np
	def __init__(self, file=None):
		nc=Dataset(file)	
		self.t0=dt.datetime.strptime(nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(nc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((nc['time'][:2]))[0]))*np.arange(len(nc['time']))
		self.lon=nc['lon'][:]
		self.lat=nc['lat'][:]
		self.nc=nc
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		coords=[]
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	
		self.bbox=[(self.lon[0],self.lat[0]),(self.lon[-1],self.lat[0]),(self.lon[-1],self.lat[-1]),(self.lon[0],self.lat[-1]),(self.lon[0],self.lat[0])]
		self.bndPath=path.Path(self.bbox)
		for key in ['analysed_sst', 'sst']:
			if key in self.nc.variables.keys():
				self.varname=key
				break
	def get_slab(self,time):
		nn=np.argmin(np.abs(self.ts-time))
		self.deltaT=self.ts[nn]-time
		self.slab=self.nc[self.varname][nn,:]
		self.t=self.ts[nn]
		if 'k' in self.nc[self.varname].units:
			self.slab-=273.15
	def close(self):
		self.nc.close()



######################################     model output ###############################################################################

class myocean:
	def __init__(self, file=None):
		self.dir=file[:file.rfind('/')+1]
		self.nc=Dataset(file)
		self.t0=dt.datetime.strptime(self.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.nc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.nc['time'][:2]))[0]))*np.arange(len(self.nc['time']))
		self.lon,self.lat=self.nc['lon'][:],self.nc['lat'][:]
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		coords=[]
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	

	def get_slab(self,time):		
		nn=np.argmin(np.abs(self.ts-time))
		self.deltaT=self.ts[nn]-time
		self.T=self.nc['thetao'][nn,:]
		self.t=self.ts[nn]
	# my ocean
	def gname(self,date):	
		date2=date+dt.timedelta(days=1)	
		return 'metoffice_foam1_amm7_NWS_TEMP_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day)
	def close(self):
		self.nc.close()
	def update(self,date):
		file=self.dir+self.gname(date)
		self.nc.close()
		self.nc=Dataset(file)
		self.t0=dt.datetime.strptime(self.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.nc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.nc['time'][:2]))[0]))*np.arange(len(self.nc['time']))


class cmems:  # load any Amm product, hopefully
	"""Loading and Data acces class for copernics Amm model output.
	Call with specified Temperature, Salinity and SSH file as arguement
	testet working with amm15 and amm7. Derives filenaming convention from inpute argumenens (*file)
	and loads dataset for time (cmems.update(dt.datetime(yyyy,mm,dd,hh,MM,ss)).
	cmems.get_slab(self,time,level=0,varname='salt') loads temporal slice for layer.
	 get_hov(self,coords,varname) extracts hovmoeller data for variable.
	"""
	def __init__(self,Sfile=None, Tfile=None,SSHfile=None,UVfile=None):
		
		self.internal_names=['salt','temp','ssh','u','v']
		self.files={'salt':Sfile,'temp':Tfile,'ssh':SSHfile,'u':UVfile,'v':UVfile}
		self.ncs={}					
		for key,item in zip(self.files.keys(),self.files.items()):
			if item[1]==None:
				self.internal_names.remove(key)
				#del self.files[key]
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
			self.fname_parts[intern]=self.scan_file_name(file)
		
		key=self.internal_names[0]
		try:
			self.depths=self.ncs[key]['depth'][:]
			self.izsurf=np.argmin(np.abs(self.ncs[key]['depth'][:]))
		except:
			pass
		self.nc=self.ncs[key]
		nrs=['0','1','2','3','4','5','6','7','8','9']
		isnr=np.zeros(len(self.nc['time'].units),bool)
		for i,letter in enumerate(self.nc['time'].units):
			if letter in nrs:
				isnr[i]=True
		nrloc=np.where(isnr)[0]	
		self.ts=netCDF4.num2date(self.nc['time'][:],self.nc['time'].units)	
		self.t0=self.ts[0]
		#self.t0=dt.datetime.strptime(self.nc['time'].units[nrloc[0]:nrloc[-1]+1],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.nc['time'][0]))
		#self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.nc['time'][:2]))[0]))*np.arange(len(self.nc['time']))
		#self.lon,self.lat=self.nc['lon'][:],self.nc['lat'][:]
		try:
			self.lon,self.lat=self.nc['lon'][:],self.nc['lat'][:]
		except:	
			self.lon,self.lat=self.nc['longitude'][:],self.nc['latitude'][:]
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		
		coords=[]
		#self.slab={'salt':0,'temp':0,'ssh':0}
		self.profile={'salt':0,'temp':0}
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	
		
	def scan_file_name(self,file):
		""" function derives file name convetion from initilizing filenames to enable name generation by date"""
		nrs=['0','1','2','3','4','5','6','7','8','9']
		isnr=np.zeros(len(file),bool)
		for i,letter in enumerate(file):
			if letter in nrs:
				isnr[i]=True
		# year is coded by 8 nrs		
		# detect year blocks
		nrpos=np.where(isnr)[0]

		yearblock=[]
		for pos in nrpos:
			if (pos in nrpos) & (pos+7 in nrpos) & (pos+8 not in nrpos):
				yearblock.append((pos,pos+8))
		prefix=file[:yearblock[0][0]]
		if len(yearblock)==2: #two date name pattern (e.g. simulated date and date of simulation performed)
			midfix=file[yearblock[0][1]:yearblock[1][0]]		
			suffix=file[yearblock[1][1]:]		
		else:
			midfix=''
			suffix=file[yearblock[0][1]:]	
		return prefix,midfix,suffix
 
	def get_slab(self,time,level=0,varname='salt'):		
		nn=np.argsort(np.abs(self.ts-time))[:2] # interpolated linearly
		self.deltaT=np.abs(self.ts[nn]-time)
		dts=np.asarray([dti.total_seconds() for dti in self.deltaT])
		if len(self.ts)>1:
			w=1/dts/(1/dts).sum()
			w[np.isnan(w)]=1.0
			self.t=time 
			if self.ncs[varname][self.varnames[varname]].shape==self.ncs['salt'][self.varnames['salt']].shape:
				self.slab=self.ncs[varname][self.varnames[varname]][nn[0],level,:]*w[0]+self.ncs[varname][self.varnames[varname]][nn[1],level,:]*w[1]
			else:
				self.slab=self.ncs[varname][self.varnames[varname]][nn[0],:]*w[0]+self.ncs[varname][self.varnames[varname]][nn[1],:]*w[1]
		else:
			self.t=self.ts
			self.slab=self.ncs[varname][self.varnames[varname]][nn[0],:]
		self.slab=self.slab[level,:]		
	def get_hov(self,coords,varname):		
		""" Generate hovmoeller data for nextneighbour coordinates of coords """
		s.nn_self=self.tree.query(coords)[1]
		ii,jj=np.unravel_index(nn_self,self.LON.shape)
		if len(self.ncs[varname][self.varnames[varname]].shape)==4:
			self.profile[varname]=self.ncs[varname][self.varnames[varname]][:][::ii,jj]
		
		
	# my ocean
	def gname(self,date):	#TEM  #SAL SSH
		date2=date+dt.timedelta(days=0)	 # ...<file generation date>_<file_content_date>
		self.fnames={}
		for key in self.internal_names:
			prefix,midfix,suffix=self.fname_parts[key]
			if midfix=='':
				self.fnames[key]=self.dir+'{:s}{:04}{:02}{:02}{:s}{:s}'.format(prefix,date2.year,date2.month,date2.day,midfix,suffix)
			else:
				self.fnames[key]=self.dir+'{:s}{:04}{:02}{:02}{:s}{:04}{:02}{:02}{:s}'.format(prefix,date2.year,date2.month,date2.day,midfix,date.year,date.month,date.day,suffix)
		return self.fnames

        #def gnameExtern(date):	#TEM  #SAL SSH
	#	""" Name generation function for call outside of created class instance
	#	"""
	#	date2=date+dt.timedelta(days=1)	 # ...<file generation date>_<file_content_date>
	#	return ['metoffice_foam1_amm15_NWS_{:s}_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(var,date2.year,date2.month,date2.day,date.year,date.month,date.day) for var in ['TEM','SAL','SSH']]
		
	def close(self):
		self.snc.close()
	def update(self,date):
		self.fnames=self.gname(date)
		for key in self.internal_names:
			self.files[key]=self.fnames[key] 
		#{'salt':Sfile,'temp':Tfile,'ssh':SSHfile}
		#Tfile=tname #self.dir+
		#Sfile=sname
		#SSHfile=shname
		#self.files={'salt':Sfile,'temp':Tfile,'ssh':SSHfile}
		#Sfile,Tfile,SSHfile
		try:
			for key in self.internal_names:
				self.ncs[key].close()
				self.ncs[key]=Dataset(self.files[key])
		except:
			pass
		#self.ncs={'salt':Dataset(Sfile),'temp':Dataset(Tfile),'ssh':Dataset(SSHfile)}
		self.nc=self.ncs[self.internal_names[0]]# ['salt']
		i0=self.nc['time'].units.rindex('since')+len('since')+1
		i1=i0+19
		self.t0=dt.datetime.strptime(self.nc['time'].units[i0:i1],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.nc['time'][0]))
		timediff=np.diff((self.nc['time'][:2]))
		ut=utime(self.nc['time'].units)
		self.ts=ut.num2date(self.nc['time'][:] )
			
		
class Amm15:
	def __init__(self, Tfile=None, Sfile=None,SSHfile=None):
		self.dir=Tfile[:Tfile.rfind('/')+1]
		self.ncs={'salt':Dataset(Sfile),'temp':Dataset(Tfile),'ssh':Dataset(SSHfile)}
		self.varnames={'salt':'so','temp':'thetao','ssh':'zos'}
		self.nc=self.ncs['salt']
		self.t0=dt.datetime.strptime(self.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.nc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.nc['time'][:2]))[0]))*np.arange(len(self.nc['time']))
		self.lon,self.lat=self.nc['lon'][:],self.nc['lat'][:]
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		self.izsurf=np.argmin(np.abs(self.nc['depth'][:]))
		coords=[]
		#self.slab={'salt':0,'temp':0,'ssh':0}
		self.profile={'salt':0,'temp':0}
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	
 
	def get_slab(self,time,level=0,varname='salt'):		
		nn=np.argsort(np.abs(self.ts-time))[:2] # interpolated linearly
		self.deltaT=np.abs(self.ts[nn]-time)
		dts=np.asarray([dti.total_seconds() for dti in self.deltaT])
		w=1/dts/(1/dts).sum()
		w[np.isnan(w)]=1.0
		self.t=time 
		if self.ncs[varname][self.varnames[varname]].shape==self.ncs['salt']['so'].shape:
			self.slab=self.ncs[varname][self.varnames[varname]][nn[0],level,:]*w[0]+self.ncs[varname][self.varnames[varname]][nn[1],level,:]*w[1]
		else:
			self.slab=self.ncs[varname][self.varnames[varname]][nn[0],:]*w[0]+self.ncs[varname][self.varnames[varname]][nn[1],:]*w[1]
		
	def get_hov(self,coords,varname):		
		""" Generate hovmoeller data for nextneighbour coordinates of coords """
		s.nn_moc=self.tree.query(coords)[1]
		ii,jj=np.unravel_index(nn_moc,self.LON.shape)
		if len(self.ncs[varname][self.varnames[varname]].shape)==4:
			self.profile[varname]=self.ncs[varname][self.varnames[varname]][:][::ii,jj]
		
		
	# my ocean
	def gname(self,date):	#TEM  #SAL SSH
		date2=date+dt.timedelta(days=1)	 # ...<file generation date>_<file_content_date>
		return ['metoffice_foam1_amm15_NWS_{:s}_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(var,date2.year,date2.month,date2.day,date.year,date.month,date.day) for var in ['TEM','SAL','SSH']]

	def gnameExtern(date):	#TEM  #SAL SSH
		""" Name generation function for call outside of created class instance
		"""
		date2=date+dt.timedelta(days=1)	 # ...<file generation date>_<file_content_date>
		return ['metoffice_foam1_amm15_NWS_{:s}_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(var,date2.year,date2.month,date2.day,date.year,date.month,date.day) for var in ['TEM','SAL','SSH']]
		
	def close(self):
		self.snc.close()
	def update(self,date):
		tname,sname,sshname=self.gname(date)
		Tfile=self.dir+tname
		Sfile=self.dir+sname
		SSHfile=self.dir+sshname
		try:
			for key in self.ncs.keys():
				self.ncs[key].close()
		except:
			pass
		self.ncs={'salt':Dataset(Sfile),'temp':Dataset(Tfile),'ssh':Dataset(SSHfile)}
		self.nc=self.ncs['salt']
		self.t0=dt.datetime.strptime(self.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.nc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.nc['time'][:2]))[0]))*np.arange(len(self.nc['time']))





class schism_output():
    import netCDF4
    nc = None

    def __init__(self,filename):
      """
      read output filename and initialize grid
      """
      import netCDF4
      from netcdftime import utime
      self.nc = netCDF4.Dataset(filename)
      self.ncv = self.nc.variables
      self.lon = self.ncv['SCHISM_hgrid_node_x'][:]
      self.lat = self.ncv['SCHISM_hgrid_node_y'][:]
      self.nodeids = np.arange(len(self.lon))
      self.nv = self.ncv['SCHISM_hgrid_face_nodes'][::3]-1
      self.time = self.ncv['time'][:] # s
      self.ut = utime(self.ncv['time'].units)
      self.dates = self.ut.num2date(self.time)
      self.node_tree_latlon = None

    def init_node_tree(self,latlon=True):
      """
      build a node tree using cKDTree
      for a quick search for node coordinates
      """
      from scipy.spatial import cKDTree
      if latlon:
        self.node_tree_latlon = cKDTree(zip(self.lon,self.lat))
      else:
        self.node_tree_xy = cKDTree(zip(self.x,self.y))

    def find_nearest_node(self,x,y,latlon=True):
      """
      find nearest node for given coordinate,
      returns the node id
      """
      ridx=-1
      if latlon:
        if self.node_tree_latlon==None:
          self.init_node_tree(latlon=True)
        d,idx = self.node_tree_latlon.query((x,y),k=1)
        ridx = self.nodeids[idx]
      else:
        if self.node_tree_latlon==None:
           self.init_node_tree(latlon=False)
        d,idx = self.node_tree_latlon.query((x,y),k=1)
        ridx = self.inodes[idx]
      return ridx


if __name__ == '__main__':

    from pylab import *
    setup = schism_setup('hgrid.gr3',ll_file='hgrid.gr3')
    # plot domain
    setup.plot_domain_boundaries()

    #triplot(setup.x,setup.y,asarray(setup.nv)-1)

    show()

    if False:
      # read elevation boundaries
      t,e = setup.bdy_array('elev2D.th')
      figure()
      plot(t[:],e[:50])
      show()		

	  
#######  Wave Watch 3   ######
class WW3_mesh():
	import numpy as np
	from matplotlib import pyplot as plt
	from matplotlib.collections import PolyCollection
	import xarray as xr
	def __init__(self,file):
		with open(file) as f:
			
			lines=f.readlines()

			nodestart=lines[lines.index('$Nodes\n')+1]
			nodestart=lines.index(nodestart)
			nnodes=np.int(lines[lines.index('$Nodes\n')+1].split(' ')[0])

			elemstart=lines[lines.index('$Elements\n')+1]
			elemstart=lines.index(elemstart)
			nelem=np.int(lines[lines.index('$Elements\n')+1].split(' ')[0])

			self.nodes=np.loadtxt(file,skiprows=nodestart+1,max_rows=nnodes)[:,1:]
			self.elems=np.loadtxt(file,skiprows=elemstart+1,max_rows=nelem)[:,6:]-1
			self.elems=np.asarray(self.elems,np.int)
			self.x,self.y,self.d=self.nodes[:,0],self.nodes[:,1],self.nodes[:,2]
			
			from scipy.spatial import cKDTree
			self.node_tree = cKDTree(list(zip(self.x,self.y)))


	def find_nearest_node(self,x,y):
		"""
		find nearest node for given coordinate,
		returns the node index (zero based)
		"""
		ridx=-1
		d,idx = self.node_tree.query((x,y),k=1)
		#ridx = self.nodeids[idx]
		return idx	
			
		  # plot functions 
	def plotAtnodes(self,nodevalues,cmap=plt.cm.jet,mask=None):
		"""
		visualisation routine plotting triangles at nodes (quads are splitted)
		"""
		ph=plt.tripcolor(self.x,self.y,self.elems,nodevalues,shading='flat',mask=mask,cmap=cmap)	  
		ch=plt.colorbar()

		return ph,ch		
		
	def plot_mesh(self,tri_color='k',quad_color='m',linewidth=0.2):	  
		xy=np.c_[self.x,self.y]
		tripc = PolyCollection(xy[self.elems,:3],facecolors='none',edgecolors=tri_color,linewidth=linewidth) #, **kwargs)
		plt.gca().add_collection(tripc)
		
	def load_points(self,fname,sep=' '):	  
		with open(fname) as f:
			x,y,names=[],[],[]
			for line in f.readlines():
				line=line.split(' ')
				xi,yi,name=[li for li in line if li != '']
				x.append(np.float(xi))
				y.append(np.float(yi))
				names.append(name.split('\n')[0].replace("'",''))
			self.px=np.asarray(x)
			self.py=np.asarray(y)
			self.pnames=names
			
	def query_station_lon_lat(self,sation_name):	  
		""" returns station lon lat and index from pointlist """
		ind=self.pnames.index(sation_name)			
		return self.px[ind],self.py[ind],ind
		
	def plot_points(self,showname=False):	  
		plt.plot(self.px,self.py,'ko')
		if showname:
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,' '+ name,color='w')
		else:	
			count=0
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,str(count))		
				count+=1
				
	def open_point_output(self,file):	  
		""" initiate xarray netecd acces in ww3.ds """
		self.specds=xr.open_mfdataset(file)			
		print('calcualting Hs from point output spectra')		
		self.f=self.specds['frequency'].values
		self.directions=self.specds['direction'].values
		ddir=np.abs(np.diff(self.directions)[0])
		ddir_rad = ddir*np.pi/180
		df=np.diff(self.f)
		self.spec2d=self.specds['efth'].values
		self.spec1d=self.spec2d.sum(axis=-1)*ddir_rad # Integral over directions
		#  <n2>integral over freqencies trapez
		self.n2=np.sum(df*0.5*(self.spec1d[:,:,:-1]+self.spec1d[:,:,1:]),axis=2)
		self.Hs=4*np.sqrt(self.n2)
		self.tspec=self.specds['time'].values
		

	def open_grid_output(self,file):	  
		""" initiate xarray netecd acces in ww3.ds """
		self.gridds=xr.open_mfdataset(file)
	
	def plot_grid_timeseries(self,parameter='hs',igrid=[0],iscoordinates=False):
		if iscoordinates==True:
			igrid=self.node_tree.query(igrid)[1]
		plt.plot(self.gridds['time'],self.gridds[parameter][:,igrid])
		plt.ylabel(parameter +' '+self.gridds[parameter].units)
		plt.gcf().autofmt_xdate()

	def get_grid_timeseries(self,parameter='hs',igrid=[0],iscoordinates=False):
		if iscoordinates==True:
			igrid=self.node_tree.query(igrid)[1]
		return self.gridds['time'].values,self.gridds[parameter][:,igrid].values
		
		
		
		
		
#######  Wave Watch 3   ######
class WW3_mesh():
	import numpy as np
	from matplotlib import pyplot as plt
	from matplotlib.collections import PolyCollection
	import xarray as xr
	def __init__(self,file):
		with open(file) as f:
			
			lines=f.readlines()

			nodestart=lines[lines.index('$Nodes\n')+1]
			nodestart=lines.index(nodestart)
			nnodes=np.int(lines[lines.index('$Nodes\n')+1].split(' ')[0])

			elemstart=lines[lines.index('$Elements\n')+1]
			elemstart=lines.index(elemstart)
			nelem=np.int(lines[lines.index('$Elements\n')+1].split(' ')[0])

			self.nodes=np.loadtxt(file,skiprows=nodestart+1,max_rows=nnodes)[:,1:]
			self.elems=np.loadtxt(file,skiprows=elemstart+1,max_rows=nelem)[:,6:]-1
			self.elems=np.asarray(self.elems,np.int)
			self.x,self.y,self.d=self.nodes[:,0],self.nodes[:,1],self.nodes[:,2]
			
			from scipy.spatial import cKDTree
			self.node_tree = cKDTree(list(zip(self.x,self.y)))


	def find_nearest_node(self,x,y):
		"""
		find nearest node for given coordinate,
		returns the node index (zero based)
		"""
		ridx=-1
		d,idx = self.node_tree.query((x,y),k=1)
		#ridx = self.nodeids[idx]
		return idx	
			
		  # plot functions 
	def plotAtnodes(self,nodevalues,cmap=plt.cm.jet,mask=None):
		"""
		visualisation routine plotting triangles at nodes (quads are splitted)
		"""
		ph=plt.tripcolor(self.x,self.y,self.elems,nodevalues,shading='flat',mask=mask,cmap=cmap)	  
		ch=plt.colorbar()

		return ph,ch		
		
	def plot_mesh(self,tri_color='k',quad_color='m',linewidth=0.2):	  
		xy=np.c_[self.x,self.y]
		tripc = PolyCollection(xy[self.elems,:3],facecolors='none',edgecolors=tri_color,linewidth=linewidth) #, **kwargs)
		plt.gca().add_collection(tripc)
		
	def load_points(self,fname,sep=' '):	  
		with open(fname) as f:
			x,y,names=[],[],[]
			for line in f.readlines():
				line=line.split(' ')
				xi,yi,name=[li for li in line if li != '']
				x.append(np.float(xi))
				y.append(np.float(yi))
				names.append(name.split('\n')[0].replace("'",''))
			self.px=np.asarray(x)
			self.py=np.asarray(y)
			self.pnames=names
			
	def query_station_lon_lat(self,sation_name):	  
		""" returns station lon lat and index from pointlist """
		ind=self.pnames.index(sation_name)			
		return self.px[ind],self.py[ind],ind
		
	def plot_points(self,showname=False):	  
		plt.plot(self.px,self.py,'ko')
		if showname:
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,' '+ name,color='w')
		else:	
			count=0
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,str(count))		
				count+=1
				
	def open_point_output(self,file):	  
		""" initiate xarray netecd acces in ww3.ds """
		self.specds=xr.open_mfdataset(file)			
		print('calcualting Hs from point output spectra')		
		self.f=self.specds['frequency'].values
		self.directions=self.specds['direction'].values
		ddir=np.abs(np.diff(self.directions)[0])
		ddir_rad = ddir*np.pi/180
		df=np.diff(self.f)
		self.spec2d=self.specds['efth'].values
		self.spec1d=self.spec2d.sum(axis=-1)*ddir_rad # Integral over directions
		#  <n2>integral over freqencies trapez
		self.n2=np.sum(df*0.5*(self.spec1d[:,:,:-1]+self.spec1d[:,:,1:]),axis=2)
		self.Hs=4*np.sqrt(self.n2)
		self.tspec=self.specds['time'].values
		

	def open_grid_output(self,file):	  
		""" initiate xarray netecd acces in ww3.ds """
		self.gridds=xr.open_mfdataset(file)
	
	def plot_grid_timeseries(self,parameter='hs',igrid=[0],iscoordinates=False):
		if iscoordinates==True:
			igrid=self.node_tree.query(igrid)[1]
		plt.plot(self.gridds['time'],self.gridds[parameter][:,igrid])
		plt.ylabel(parameter +' '+self.gridds[parameter].units)
		plt.gcf().autofmt_xdate()

	def get_grid_timeseries(self,parameter='hs',igrid=[0],iscoordinates=False):
		if iscoordinates==True:
			igrid=self.node_tree.query(igrid)[1]
		return self.gridds['time'].values,self.gridds[parameter][:,igrid].values
		

file='Blankenese_WGMN!Truebung_(Tagesmaxima).txt'		
#######  Wave Watch 3   ######
class tide_portal_ascii():
	import numpy as np
	import pandas as pd
	def __init__(self,file):
		with open(file) as f:
		
		#read header
		self.meta={}
		with open(file) as f
			#lines=f.readlines()
			header=True
			for linenr,line in enumerate(lines):
			
				if (line[0] != '=') & header:
					outline=line.split('\n')[0].split(': ')
					tag,value=outline[0].strip(' '),outline[1].strip()
					self.meta[tag]=value
						
				elif (line[0] == '=') :
					if header:
						skiprows=linenr
						header=False
						#df=p.
					else:
						max_rows=linenr
						break
				else:
					pass
					#_dfs = [
					#	pandas.DataFrame([line.split(' ')], columns=columns, dtype=float),
					#pandas.read_table(td, sep='\t', header=None, names=columns)
				#]
					#df = pd.concat([df, _df])
				#df = pandas.concat(_dfs, ignore_index=True)		
				  #stop tag
		#lines			
		self.ds=pd.read_table(file,skiprows=skiprows+1,nrows=max_rows-skiprows,header=0,parse_dates=True)

	