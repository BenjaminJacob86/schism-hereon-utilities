#export OMP_NUM_THREADS=1 # call before python
"""
Script for validation of schism output against cmems tide gauge data and optional intercomparison with amm15
validates SCHISM vs tide gauge stations (for 1 year) from EMOD net 
stored by SG in /work/gg0028/g260099/OSTIA/VALIDATION/VAL_TIDE/EMOD_TIDEG_DATA/
as monthly files.
1 - search data in subfolders by name, which is in specified <year>
2 - second sub select tide gauge stations within <dtol> distance from
closest schism grid nodes. 
3 - plot 
4 - derive skill parameter (bias,rmse,correlation) after interpolation of model data
to ovebservational time steps - the temporal resolution is reduced to the model time step
Comparison is done with zcor/elev data from SCHISMs closest grid point
Uses schism class stored in /pf/g/g260114/Programs/python/scripts/schism.py

call after 
module load python3
if mistralpp busy, performance better when called from run_sript

Output:
 pdf with validation images
'errorstats.csv' - csv file with error statistics
"""

__author__ = "Benjamin Jacob"
__copyright__   = "Copyright 2020 - 03\2021 Helmholtz-Zentrum Geesthacht GmbH"
__copyright__   = "Copyright 04\2021 - Helmholtz-Zentrum Hereon GmbH"

__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"

import os
import netCDF4
import sys
import csv
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
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
# own and 3d party libraries
#sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities')
#sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')

sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')

#levante
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/Lib/')
#sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
#sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
#from techit import * # import latex script
from schism import * # import schism functions
#from techit2 import * # import latex script
from TaylorDiagram import * # import taylordiagram
#from data_and_model_classes import cmems
import pandas as pd
#import utide
from matplotlib.dates import date2num
import xarray as xr
#from numba import jit
import time
plt.ion()

########## settings #################################
#tgdir='/work/bg1186/g260099/OSTIA/VALIDATION/VAL_TIDE/EMOD_TIDEG_DATA/EMOD_DATA_201903/' #tgdir # directories (have to end with '/')

cluster='levante' # levante.


# levante
if cluster=='levante':
	#tgdir='/work/gg0028/ksddata/insitu/CMEMS/NorthWestShelf/TG/'
	tgdir='/work/gg0028/ksddata/insitu/CMEMS/NorthWestShelf/nrt.cmems-du.eu/Core/INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036/cmems_obs-ins_nws_phybgcwav_mynrt_na_irr/history/TG/'
	#tgdir='/work/gg0028/ksddata/insitu/TideGauge/MyOcean/'
elif cluster=='strand':	
	tgdir='/gpfs/work/ksddata/observation/insitu/CMEMS/NorthWestShelf/TG/'


oceandir='/work/gg0028/g260114/RUNS/GermanBight/GB_IBI_2023_AMM15_control/AMM/' #amm15 directory
ibi_dir='/work/gg0028/g260098/VERIBLUE/IBI_NAT_PHY/nrt/'			   # IBI top dir
	
# schism setup(s) for cross validation TG stations are selected based on coverage of first setup
#setupdir: schism run directory, containts hgrid.* etc.

setupdir=['/work/gg0028/g260114/RUNS/GermanBight/GB_IBI_2023_AMM15_control/','/work/gg0028/g260114/RUNS/GermanBight/GB_IBI_2023/']
#setupdir=['/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/',
#'/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/',
#'/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_max/',
#'/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_HE/',
#'/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_LE/'][:3]

setup_names=['FRC-AMM','FRC-IBI']


#ncdir=[setupdir[i] + 'outputs_all/' for i in range(len(setup_names))] 		  #   directory of schism nc output or 
ncdir=[setupdir[i] + 'outputs/' for i in range(len(setup_names))] 		  #   directory of schism nc output or 

max_stack=30 #350   #  -1:= use all stacks. maximum stack number to load for newio outputs


# True: use station output False: extract from netcdf (slow)
#use_station_in=[True,True] #,False,False]					  
use_station_in=[False]*len(setup_names) 
######################
 

outdir=setupdir[0]+'TGvalid_AMM_IBI4/' 	    # output directory where images will be stored
year=2023							# year to be analyesed if year=None use first year of model run 	 
dtol=0.4    #0.05 					# distance tolerance in degree lon/lat towards tg stations 


t00=np.datetime64('2023-01-03')#('2013-10-27') #=False		  # only use time period after t00 if t00=False use all avaialbale overlapping data	
#t11=np.datetime64('2013-11-04') #=False					  # only use time period until t00 if t11=False use all not in place yet

remove_mean=True  					  # remove temporal mean from Data and Model to compare 2
use_amm=True						  # compare against amm15
use_ibi=True
shift_dt=np.timedelta64(1,'h')		  # timestep correction if necessary


#--- what to plot True/False:			
overview_map=True												
satistic_maps=False
full_timeseries=False 				      # for flull year data hard to see anything
first_two_Weeks=True 				      # zomed period 			
running_avg_window=np.timedelta64(25,'h') # running mean filter in hours
monthly_subplots=False
taylor_diag=True
consts=['M2','S2','M4','M6']			# tidal constituents to compute - needs utide package
tidal_bar_plots=False
tidal_difference_maps=False


# FIlter Criteria Exclude bad statons or add known offsets
exclude_stations=['BallumTG','BronsTG','DelfzijlTG','HusumTG','MandoTG','NesTG','NieuweStatenzijlTG','RibeTG','VidaaTG','TerschellingNoordzeeTG','HarlingenTG'] #['Buesum','Delfzijl','Husum','Havneby','Wittduen']#,'Huibertgat'            #'HusumTG', ,'Havneby'
offset_stations={}			   #{'WittduenTG':-5,'BuesumTG':-5} # add value to stations 
dC_dt_max=1				       # maximum allowed change in sea level per hour to be considered outlier
							   # check consecutive +- deviation	
dmin=2.					       # minimum depths for schism nns 	   

###################

# images ############################
pic_format='.png'
dpivalue=300
cmap='jet'
Rcircle=150 							# radius of colorcoded cricles used in scatter maps
limit_to_data=True   					# limit map plots to bounding box of data.


station_style='k.'  #'k.-'  #'k-' # linestyle color for stations.

# Font sizes
Fincrease=1.0
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



### Latex output
put_pics_to_texdoc=True    										# images will be put in tex document
latexname='GB_IBI_valid.tex'									
latextitle='GB HRballje'
latextext='Tide Gauge validation of GB ' + str(year) +'n against Emodnet data.'
#############################################################

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

#  functions
def techit(fname='test.tex',title='manuscript',text='',figures=None,captions=None):

        from os import listdir
        from subprocess import call


        file=open(fname,'w')

        file.write("\\documentclass[a4paper]{article}\n")
        file.write("\\usepackage[a4paper]{geometry}\n")
        file.write("\\parindent0pt \n \\textwidth16cm \n \\textheight22cm \n")

        #file.write("\\usepackage{graphicx}\n")
        file.write("\\usepackage[pdftex]{graphicx}\n")
        file.write("\\usepackage{epsfig}\n")
        file.write("\\usepackage{epstopdf}\n")


        file.write("\\begin{document}\n")

        file.write("\\title{"+title+"}\n")
        file.write("\\maketitle\n")
        
        file.write("\n"+text+"\n")

        count=1
        files=sorted(listdir('./'))
        for figure,caption in zip(figures,captions):
                count=count+1
                print(figure)
                file.write("\n")
                file.write("\\begin{figure} \n \\centering \n")
                file.write(" \\includegraphics[width=1.00\\textwidth]{%s} \n" %(figure))
                file.write("\\caption{%s}\n" %(caption))
                file.write("\\end{figure}\n")
                if count%2==0:
                        file.write("\\clearpage\n")

        file.write("\\end{document}")


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
		
######### load SCHISM setup   ##################################
cwd=os.getcwd()
setups={}
output={}
access={}

newio=[]
for i,folder in enumerate(setupdir):
	os.chdir(folder)
	s=schism_setup()
	lon,lat,D=np.asarray(s.lon),np.asarray(s.lat),np.asarray(s.depths)
	lon[D>dmin]=9999
	s.nntree = cKDTree(list(zip(lon,lat))) 
	setups[i]=s
	if use_station_in[i]:
		if os.path.exists('station.in'): #& os.path.exists(setupdir[i]+'staout_1'):
			staout=schism_station_output()
			output[i]=staout
		else:
			print('station output does not exist, exiting program')
			exit()
	else:
		schismfiles=[] 
		if len(glob.glob(ncdir[i]+'schout_*.nc'))==0:
			#new io
			#s.nc=schism_outputs_by_variable(ncdir[i])
			testfile=np.sort(glob.glob(ncdir[i]+'*out2d*.nc'))[0]
			if len(testfile)>0:
				dsi=xr.open_dataset(testfile)
				use_elev= 'elevation' in list(dsi.keys())
			else:	
				use_elev=False
			if use_elev:
				access[i]=schism_outputs_by_variable(ncdir[i],varlist=['out2d'],max_stack=max_stack)
			else:	
				access[i]=schism_outputs_by_variable(ncdir[i],varlist=['zCoordinates'],max_stack=max_stack)
			newio.append(True)
		else:	#old io schout_
			for iorder in range(6): # check for schout_nc files until 99999
				schismfiles+=glob.glob(ncdir[i]+'schout_'+'?'*iorder+'.nc')
			nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
			schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
			nrs=list(np.asarray(nrs)[np.argsort(nrs)])
			s.nc=xr.open_mfdataset(schismfiles)
			access[i]=s.nc
			newio.append(False)
	
# use cmems TG files  # pre select in domain
pattern='SSH'							 							# pattern to check in amm15 files 		

if use_station_in[0]:
	T=staout.time
	T=np.asarray([np.datetime64(ti) for ti in T])
else:
	if newio[0]:
		os.chdir(setupdir[0])
		p=param()
		reftime=np.asarray(dt.datetime(np.int(p.get_parameter('start_year')),\
	np.int(p.get_parameter('start_month')),\
	np.int(p.get_parameter('start_day')),\
	np.int(p.get_parameter('start_hour')),0,0),np.datetime64)
		T=reftime+access[0].ds['out2d']['time'].values *np.timedelta64(1,'s')+shift_dt
	else:
		#T=s.nc['time'].values
		T=access[0]['time'].values

# contains year folders or station name folders
#tgdir='/gpfs/work/ksddata/observation/insitu/CMEMS/NorthWestShelf/TG/'  # year folders withs everal stations
#tgdir='/gpfs/work/ksddata/observation/insitu/TideGauge/MyOcean/'  # station folders

# files should be loaded for differetn situations
#a) monthly folders withs several stations (aquire from CMEMS)
#b) station folders with files for different dates emodnet via Sebastian

stations={'coord':[],'names':[],'TG':{}}
sources=setup_names.copy()



if use_ibi:
	sources.append('IBI')
	stations['IBI']={'nc':[],'time':[],'nn':[],'coord':[],'zeta':[],'names':[]}
	#setup_names.append('IBI-nn')

if use_amm:
	sources.append('amm')
	#setup_names.append('IBI-nn')
	
for key in sources:	
	stations[key]={'nc':[],'time':[],'nn':[],'coord':[],'zeta':[],'names':[]}

	
	
if  'history' in tgdir:	
	files=glob.glob(tgdir+'/*.nc')
	nfolders=1
	#foldertype='time'
	foldertype='history'
	folders=tgdir
else:
	
	if len(glob.glob((tgdir+'20??*'))) > 0: #ear folders
		foldertype='time'
		folders=np.sort(glob.glob(tgdir+'{:d}*'.format(year)))
		files=glob.glob(folders[0]+'/*.nc')
		nfolders=1
	else:
		foldertype='station'
		folders=glob.glob(tgdir+'/*/')
		nfolders=len(folders)
		files=glob.glob(folders[0]+'/*.nc')

	
# select stations withing domain (tolerance distance) # and constructe file acess and nearest neighbours
names=[]	
for i,setup_name in enumerate(setup_names):
	s=setups[i]
	print('find neighbours for ' + setup_name )
	if use_station_in[i]:
		staout=output[i]
		
	for folder in folders[:nfolders]:
		nfiles=1 if foldertype == 'station' else len(files) # files are eother timer or station
		if foldertype=='station':  # read new stations otherwise stations are in date folder
			#files=glob.glob(folder+'/*NO*.nc')
			files=glob.glob(folder+'/*.nc')
		for file in files[:nfiles]: 
			#if 'Helgo' in file:
			a=xr.open_dataset(file,decode_times=False)
			coord=np.float(a.geospatial_lon_min),np.float(a.geospatial_lat_min)
			if use_station_in[i]:
				D=np.sqrt((staout.coords[:,0]-coord[0])**2+(staout.coords[:,1]-coord[1])**2)
				nn=np.argmin(D)
				lldists=D[nn]
			else:
				lldists,nn=s.nntree.query(coord)		
			# in domain
			if lldists < dtol:
				#print('match')
				#break
				if foldertype=='station':
					name=folder[folder[:-1].rindex('/')+1:-1]
				else:
					name=file[file.rindex('/')+1:].split('_')[-1].split('.')[0]
					#name=file[file.rindex('/')+1:].split('_')[-2]
				#statioprint('using ' + name)
				
				if True: # (i==0):
					if nfiles > 1:
						if  foldertype=='history':
							ncfiles=[file]						
						else:						
							ncfiles=[file.replace('{:d}01'.format(year),'{:d}{:02d}'.format(year,mon+1)) for mon in range(12) if os.path.exists(file.replace('{:d}01'.format(year),'{:d}{:02d}'.format(year,mon+1)))]
					else:	
						ncfiles=np.sort(files)
						
					try:
						ncs=MFDataset(ncfiles,aggdim='TIME')
					except:
						continue
					t=(ncs['TIME'][:])
					try:
						if np.ma.isMaskedArray(ncs['SLEV_QC'][:,0]):
							igood=(ncs['SLEV_QC'][:] < 3)[:,0] & ncs['SLEV_QC'][:,0].mask==False
						else:
							igood=(ncs['SLEV_QC'][:] < 3)[:,0] # select good data (ncs['TIME_QC'][:]<3) and
					except:
						try:
							if np.ma.isMaskedArray(ncs['SLEV'][:,0]):
								igood=ncs['SLEV'][:,0].mask==False # np.ones(len(t),bool)
							else:	
								igood=np.ones(len(t),bool)
						except:
							#pass
							continue
					
					t0=dt.datetime.strptime(ncs['TIME'].units[11:30],'%Y-%m-%dT%H:%M:%S')
					timeunit=ncs['TIME'].units[:ncs['TIME'].units.index(' ')]
					t=t[igood]
					dates=np.asarray([t0+eval('dt.timedelta({:s}=float(ti))'.format(timeunit)) for ti in t],np.datetime64)
					
					# temporal overlap
					#shaere interval
					keep=(dates>T[0]-np.timedelta64(2,'h')) & (dates<T[-1]+np.timedelta64(2,'h'))
					if keep.max() and not (name in exclude_stations):
						print('adding station '+name)
						keep=(dates>T[0]-np.timedelta64(2,'h')) & (dates<T[-1]+np.timedelta64(2,'h'))
						t=t[keep]
						dates=dates[keep]
						zeta=ncs['SLEV'][:][igood][keep][:,0]
						
						dates,uinds=np.unique(dates,return_index=True)
						zeta=zeta[uinds]
						
						# deal with potential doubles
						
						#filter outliers by rate of change
						dC_dt=np.abs((zeta[1:]-zeta[:-1])/((dates[1:]-dates[:-1])/np.timedelta64(1,'h')))
						toolarge=dC_dt > dC_dt_max
						iout=np.where(toolarge[:-1] & toolarge[1:])[0]+1
						zeta=np.delete(zeta,iout)
						dates=np.delete(dates,iout)
						#zetastd=np.std(zeta)
						#zetamean=np.mean(zeta)
						#ivalid=[(zeta>zetamean-3*zetastd)	& (zeta<zetamean+3*zetastd)]
						stations['TG'][name]={'time':dates,'zeta':zeta}
						if i==0:
							names.append(name)
							stations['coord'].append(coord)
						stations[setup_name]['names'].append(name)	
						if use_station_in[i]:
							stations[setup_names[i]]['coord'].append((staout.coords[nn,0],staout.coords[nn,1]))
						else:
							stations[setup_names[i]]['coord'].append((s.lon[nn],s.lat[nn]))
						stations[setup_names[i]]['nn'].append(nn)
					else:
						print('not adding station '+name)
print('done selecting and loading tide gauges, found: '+str(len(stations['coord'])) + ' stations meeting criterea')
		


# load data from schism next neighbours to TG statsions
modeldt=[]
dates=[]
for i,setup_name in enumerate(setup_names):
	print(setup_name)
	s=setups[i]
	#name=names[i]
	#ind=stations[setup_name]['names'].index(name)
	if use_station_in[i]:
		staout=output[i]
		modeldt.append(np.timedelta64(staout.time[1]-staout.time[0]))
		dates.append(np.asarray(staout.time,np.datetime64))   # convert to np.datetime64
		stations[setup_name]['zeta']=staout.station_out['elev'][:,stations[setup_name]['nn']]
		modeldt[i]=np.timedelta64(modeldt[i])
		stations[setup_name]['time']=dates[i]#np.asarray(dates,np.datetime64)
		lons,lats=staout.coords[:,0],staout.coords[:,1]
	else:

		lons,lats=np.asarray(s.lon),np.asarray(s.lat)
		print('load SCHISM sea level time series at TG next neighbouring nodes')
		
		if newio[i]:
			s.nc=access[i].get('elevation')
			#s.nc.time[:]=T
			s.nc2=s.nc.sel(nSCHISM_hgrid_node=stations[setup_name]['nn'])		
			#,nSCHISM_vgrid_layers=-1
			stations[setup_name]['time']=T
			#stations[setup_name]['zeta']=s.nc2['elevation'].values
			stations[setup_name]['zeta']=s.nc2[:].values
		else:
			#break
			s.nc=access[i]
			use_elev='elev' in s.nc.variables.keys()
			s.nc2=s.nc.sel(nSCHISM_hgrid_node=stations[setup_name]['nn'],nSCHISM_vgrid_layers=-1)
			n=len(schismfiles)	
			print(str(n) + 'files')
			if use_elev:	
				stations[setup_name]['zeta']=s.nc2['elev'].values
			else: # use zcor
				stations[setup_name]['zeta']=s.nc2['zcor'].values
			stations[setup_name]['time']=s.nc2['time'].values
		modeldt.append(stations[setup_name]['time'][1]-stations[setup_name]['time'][0])
		dates.append(stations[setup_name]['time'])
	print('load model Tide Gauge Data')
		
		


######### load IBI setup ######################################
if use_ibi:
	wild_elev_files='IBI36NWS_ng_1h-m_*_2DT-oce_*.nc' # * varying infprmation	
	ssh_files=np.sort(glob.glob(ibi_dir+'*/'+wild_elev_files))
	#i0_i1=np.where([ year in file for file in ssh_files ])[0][[0,-1]]
	#ssh_files=ssh_files[i0_i1[0]-1:i0_i1[1]+2]
	t0=stations[setup_names[0]]['time'][0]-2*np.timedelta64(1,'h')
	t1=stations[setup_names[0]]['time'][-1]
	def reduce_file_list(files):
		""" return file list and start date for IBI files """
		npt0,npt1=np.asarray([t0,t1],np.datetime64)		
		file_dates=np.asarray(get_file_start_dates(files))
		#file0=int(np.maximum(np.where(file_dates>=npt0)[0][0]-1,0)) # has day before included
		file0=int(np.maximum(np.where(file_dates>=npt0)[0][0],0)) # has day before included, this not
		if len(np.where(file_dates>=npt1)[0]) > 0:
			file1=int(np.minimum(np.where(file_dates>=npt1)[0][0],len(files)))+1 # to be save
		else:
			file1=len(files)		
		files=files[file0:file1]
		file_dates=file_dates[file0:file1]
		return files, file_dates	
	
	ssh_files, file_dates=reduce_file_list(ssh_files)
	
	
	temp=xr.open_dataset(ssh_files[0])
	isnan=np.isnan(temp.ssh[0,:].values)
	nt_file=len(temp.time_counter)
	temp.close()
	dsIBI=xr.combine_nested([ xr.open_dataset(file)  for file in ssh_files],concat_dim='time_counter',data_vars=['ssh',])
	ntfiles=nt_file*len(ssh_files)
	dates_out=file_dates[0]+np.arange(ntfiles)*np.timedelta64(1,'h')+30*np.timedelta64(1,'m')
	lon2d,lat2d=dsIBI['nav_lon'].values,dsIBI['nav_lat'].values
	
	i_invalid=isnan#((lon2d==-1) & (lat2d==-1) ()
	lon2d[i_invalid]=-9999
	lat2d[i_invalid]=-9999
	tree2d=cKDTree(list(zip(lon2d.flatten(),lat2d.flatten())))		
	
	lldists3,nn3=tree2d.query(stations['coord'])		
	ii2,jj2=np.unravel_index(nn3,lon2d.shape)
	stations['IBI']['zeta']=np.asarray([dsIBI['ssh'][i,:].values[ii2,jj2] for i in range(len(dsIBI['time_counter']))])
	stations['IBI']['time']=dates_out
	setup_names.append('IBI-nn')
	print('done load IBI Data ')




## remove not share stations from schism files
# shared stations between setups  needs to be done flexible
#s.plot_domain_boundaries()
#stations[setup_names[i]]['coord']

#for i in range(len(a)):
#	print('dwd: ' + a[i] + 'era: ' + b[i] )

a=np.asarray(list(stations[setup_names[0]]['names']))
if len(setup_names)>1:
	b=np.asarray(list(stations[setup_names[1]]['names']))

	shared=[]
	notshared=[]
	remove=[]
	for name in a:
		if name in b:
			shared.append(name)
		else:	
			notshared.append(name)
			remove.append(np.where(a==name)[0][0])
	for key in ['nn','coord','names']:
		delete=[stations[setup_names[0]][key][ind] for ind in remove]
		for deltei in delete:
			stations[setup_names[0]][key].remove(deltei)  


	delete=[stations['coord'][ind] for ind in remove]
	for deltei in delete:
		stations['coord'].remove(deltei)  
	delete=[names[ind] for ind in remove]
	for deltei in delete:
		names.remove(deltei)
	for deltei in delete:
		del stations['TG'][deltei]	
	################################################################

	a=np.asarray(list(stations[setup_names[0]]['names']))
	b=np.asarray(list(stations[setup_names[1]]['names']))
	for i in range(len(a)):
		print('TG: ' + list(stations['TG'].keys())[i] + 'dwd: ' + a[i] + ' era: ' + b[i] )




coords=np.asarray(stations['coord'])
######### load AMM15 setup ######################################
if use_amm:
	print('load cmems Data ')
	files=np.sort(glob.glob(oceandir+'*'+pattern+'*'+str(year)+'*'))
	amm15=cmems(SSHfile=files[0])
	amm15.varnames={'ssh':'zos'}
	lldists2,nn2=amm15.tree.query(stations['coord'])		
	ii,jj=np.unravel_index(nn2,amm15.LON.shape)
	amm15.nc=xr.open_mfdataset(files)
	#stations['amm']['time']=amm15.nc.sel(time=slice("{:d}-12-31".format(year-1), "{:d}-01-01".format(year+1))).indexes['time'].to_datetimeindex()
	stations['amm']['time']=amm15.nc['time'].values
	stations['amm']['zeta']=np.asarray([amm15.nc[amm15.varnames['ssh']][i,:].values[ii,jj] for i in range(len(amm15.nc['time']))])
	print('done load cmems Data ')
	setup_names.append('AMM')
####################################




	
# limit to matching time range
nt=np.min([len(stations[key]['zeta']) for key in list(stations.keys())[3:]  if key !='TG_mean'])+1
for key in list(stations.keys())[3:]:
	if (key !='TG_mean') & (key !='IBI'):
		stations[key]['time']=stations[key]['time'][:nt]
		stations[key]['zeta']=stations[key]['zeta'][:nt,:]
	
#if remove_mean:
#		latextext+=' Time Series were mean removed and so is rmse then. In Taylor it is always mean removed'
# plt.ion() has to be off to work in background mode

# offset
for key in offset_stations.keys():	
	if key in names:
		stations[key][name]['zeta']+=offset_stations[key]
		latextext+=' station ' + name + ' was offseted by ' + str(offset_stations[key])
		
		
if 'amm' in stations.keys():		
	stations['amm']['time']=np.asarray(stations['amm']['time'],np.datetime64)
	
		
if 'IBI' in stations.keys():		
	stations['IBI']['time']=np.asarray(stations['IBI']['time'],np.datetime64)
	
### remove mean
if remove_mean:
	print('removing mean')
	for key in list(stations.keys())[2:]:
		if key=='TG':
			stations['TG_mean']={}
			for name in names:
				stations['TG_mean'][name]=stations['TG'][name]['zeta'].mean()
				stations['TG'][name]['zeta']-=stations['TG_mean'][name]
		else:	
			stations[key]['zeta_mean']=np.nanmean(stations[key]['zeta'],axis=0)
			stations[key]['zeta']-=stations[key]['zeta_mean']
	latextext+=' tide gauges and models where mean removed for every step excluding bias computation'
	
		
####### temporal interpolation of model to data  and error statistics ##########
print('interploating to common time steps of data and calculating error statistics')
sources=['TG']+sources
for i,key in enumerate(sources):
	stations[key]['interp']={name:{'time':0,'zeta':0} for name in names}
	stations[key]['stats']={'std':[],'bias':[],'rmse':[],'cor':[],'std_rel':[]}
	for subkey in stations[key]['stats'].keys():
		stations[key]['stats'][subkey]={name:0 for name in names}

n=len(names)
time_schism=stations[setup_names[0]]['time'] # dirret acces seems faster
for i,name in enumerate(names):

	
	start = time.time()
	print(str(i)+'/'+str(n)+' '+name)
	stations['TG'][name]['time']=np.asarray(stations['TG'][name]['time'],np.datetime64)	
	Todates=stations['TG'][name]['time']
	#ilast=(Todates<=time_schism[-1]).sum()
	#ifirst=(Todates<=time_schism[0]).sum()#-1?
	#ifirst=max(ifirst,0)
	#Todates=Todates[ifirst:ilast]
	iuse=(Todates >= time_schism[0]) & (Todates <= time_schism[-1])
	Todates=Todates[ iuse]
	
	if len(Todates)==0:
		continue

	# interpolate to model time step but at maximum the temporla resolution of the refrence model	
	Todates2=[Todates[0]]
	for date in Todates:
		if (date-Todates2[-1]) >= modeldt[0]:
			Todates2.append(date)
	Todates2=np.asarray(Todates2)
	# test
	if t00!=False:
		Todates2=Todates2[Todates2>=t00]
	
	stations['TG']['interp'][name]['time']=Todates2
	a,ainb,bina=np.intersect1d(Todates, Todates2, assume_unique=False, return_indices=True)
	#stations['TG']['interp'][name]['zeta']=stations['TG'][name]['zeta'][ifirst:ilast][ainb]
	stations['TG']['interp'][name]['zeta']=stations['TG'][name]['zeta'][iuse][ainb]
	
	std0=np.std(stations['TG']['interp'][name]['zeta'])
	stations['TG']['stats']['std'][name]=std0
	end = time.time()
	print("prepate TG Elapsed = %s" % (end - start))
	
	for key in sources[1:]:	

		#start = time.time()
		tin=(stations[key]['time']-time_schism[0])/(Todates2[1]-Todates2[0])
		tout=(Todates2-time_schism[0])/(Todates2[1]-Todates2[0])
		
		
		
		# datetime64
		#tin=stations[key]['time']
		#tout=Todates2
		stations[key]['interp'][name]['time']=Todates2
		fintp=interp1d(tin, stations[key]['zeta'][:,i])
		stations[key]['interp'][name]['zeta']=fintp(tout)
		#end = time.time()
		#print("interp Elapsed = %s" % (end - start))

		#start = time.time()
			
		# calculate error  stats ###########################  SLOW	
		Zeta_tg=stations['TG']['interp'][name]['zeta']#[:,0]
		Zetamod=stations[key]['interp'][name]['zeta']
		stations[key]['stats']['bias'][name]=stations[key]['zeta_mean'][i]-stations['TG_mean'][name]
		
		diff=Zeta_tg-Zetamod
		rmse=np.sqrt((diff**2).mean())
		std=np.std(stations[key]['interp'][name]['zeta']) # mean removed
		if np.isnan(np .nanmax(stations[key]['interp'][name]['zeta']) ):
			cor=np.nan
		else:	
			cor=np.corrcoef(Zetamod,Zeta_tg)[0,1]
		
		stations[key]['stats']['rmse'][name]=rmse # mean removed
		stations[key]['stats']['std'][name]=std
		stations[key]['stats']['cor'][name]=cor
		stations[key]['stats']['std_rel'][name]=std/std0
		#stations=calc_error_stats(stations,key)
		#end = time.time()
		#print("calc stats Elapsed = %s" % (end - start))


	
########## Plotting #####################################
print('plotting results')
if not os.path.exists(outdir): os.mkdir(outdir) 
##### Plot selected stations
coords=np.asarray(stations[setup_names[0]]['coord'])
xmin,ymin=np.min(coords,axis=0)
xmax,ymax=np.max(coords,axis=0)
x,y=coords[:,0],coords[:,1]
figures=[]
captions=[]
fname='0_TideGuageStationLocations.eps'
fname='0_TideGuageStationLocations2.png'
caption='Overview of Tide Gauge Stations.'

figures.append(fname)
lons,lats=np.asarray(setups[0].lon),np.asarray(setups[0].lat)
captions.append(caption)
if overview_map:
	plt.clf()
	s.plot_domain_boundaries()
	fig=plt.gcf()
	fig.set_size_inches(11,8,forward=True)
	plt.plot(lons[stations[setup_names[0]]['nn']],lats[stations[setup_names[0]]['nn']],'bo')
	for coord,name in zip(stations['coord'],names):
		lon,lat=coord
		plt.plot(lon,lat,'r+')
		plt.text(lon,lat,' '+name[:3],rotation=50, rotation_mode='anchor')
	plt.title('TG Stations')
	
	if limit_to_data:
		plt.xlim((xmin-1,xmax+1))
		plt.ylim((ymin-1,ymax+1))
	plt.tight_layout()	
	plt.savefig(outdir+fname,dpi=dpivalue)	
	plt.close()
###################################	


#name='HelgolandTG'
#ind=list(names).index(name)
#for key in setup_names:
#	plt.plot(stations[key]['time'],stations[key]['zeta'][:,ind],'.-')


####### Plot statistics #######################
names=np.asarray(names)
shortnames=np.asarray([name[:3] for name in names])
ivalids=[np.asarray(np.ones(len(coords)),bool)]*len(setup_names)

n=len(sources)+len(sources)%2
widths=np.asarray([ n*0.4/n/2 * (i+1) for i in range(np.int(n/2))])
widths=np.concatenate((-widths[::-1],widths))

colors=plt.cm.tab10(range(len(sources)))  #['b','r']
if use_amm:
	ivalid=~np.isnan(stations['amm']['zeta_mean'])
	ivalids.append(ivalid)	
plt.close()

Rcircles=[Rcircle*0.25**np.float(i) for i in range(len(sources)-1)]
phs={}
if satistic_maps:
	labels={'bias':'bias [m]','rmse':'rmse [m]','cor':'correlation [-]','std_rel':'relative std [-]'}
	for key in list(stations[setup_names[0]]['stats'].keys())[1:]:
		break
		data=[np.asarray(list(stations[setup_names[0]]['stats'][key].values()))]
		ilarger=np.abs(data[0])>3 # station offset
		data[0][ilarger]-=np.floor(data[0][ilarger])
		vmin=data[0].min()
		vmax=data[0].max()
		if len(sources[1:]) > 1:
			for i,stp in enumerate(sources[1:]):
				k=i+1
				data.append(np.asarray(list(stations[stp]['stats'][key].values())))
				data[k]=np.ma.masked_array(data[1],mask=np.isnan(data[1]))
				ilarger=np.abs(data[k])>3 # station offset
				data[k][ilarger]-=np.floor(data[1][ilarger])
		vmin=np.nanquantile(np.concatenate(data),0.1)
		vmax=np.nanquantile(np.concatenate(data),0.99)
		if key=='cor':
			vmin=0
			vmax=1
			

		plt.clf()
		s.plot_domain_boundaries()
		fig=plt.gcf()
		fig.set_size_inches(11,8,forward=True)
		for nr,model in enumerate(sources[1:]):
			phs[nr]=plt.scatter(x[ivalids[nr]],y[ivalids[nr]],s=Rcircles[nr],c=data[nr][ivalids[nr]]+nr,vmin=vmin,vmax=vmax)			
		ph=plt.colorbar()
		ph.set_label(labels[key])
		plt.legend(phs.values(),sources[1:],loc='upper center',ncol=2,frameon=False)

		for coord,name in zip(stations['coord'],names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		plt.tight_layout()
		if limit_to_data:
			plt.xlim((xmin-.2,xmax+.2))
			plt.ylim((ymin-.2,ymax+.2))

		fname='1b_'+labels[key][:labels[key].index(' ')]+'_schism_amm'+pic_format
		plt.savefig(outdir+fname,dpi=dpivalue)	
		caption='Map of '+key+ ' between ' + str(Todates2[0]) + ' and ' + str(Todates2[-1])
		figures.append(fname)
		captions.append(caption)
		plt.close()

		# bar charts
		plt.clf()
		for nr,model in enumerate(sources[1:]):
			if nr==0:
				phs[nr]=plt.bar(np.where(ivalids[nr])[0],data[nr][ivalids[nr]],widths[nr],align='edge',label=model,color=colors[nr],tick_label=shortnames[ivalids[nr]])
				plt.grid()
			else:
				phs[nr]=plt.bar(np.where(ivalids[nr])[0],data[nr][ivalids[nr]],widths[nr],align='edge',label=model,color=colors[nr])
		plt.xticks(rotation=45)
		plt.ylabel(labels[key])
		plt.legend(phs.values(),sources[1:],loc='upper center',ncol=2,frameon=False)
		fname='1c_'+labels[key][:labels[key].index(' ')]+'_bar_plot.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)	
		plt.close()
		caption='Barplot of '+key+ ' between ' + str(Todates2[0]) + ' and ' + str(Todates2[-1])
		figures.append(fname)
		captions.append(caption)
###################################	
	
#plt.figure()	
#plt.plot(stations['TG'][name]['time'],stations['TG'][name]['zeta'])
#plt.plot(stations[key]['time'],stations[key]['zeta'][:,i],'-')
#plt.title(name)
	
###### Plot first two weeks of data #####
if first_two_Weeks: 
	print('making plot of first two weeks') # of intersection
	#tmin=stations[setup_names[0]]['time'][0]
	#tmax=tmin+np.timedelta64(14,'D')
	#inds=stations[setup_names[0]]['time']<tmax
	#for key in sources[1:]:
	#	inds.append( (stations[key]['time']>=tmin) & (stations[key]['time']<=tmax))
	
	for i,name in enumerate(names):
		tmin=stations[setup_names[0]]['time'][0] if (stations[setup_names[0]]['time'][0] > stations['TG'][name]['time'][0]) else stations['TG'][name]['time'][0] 
		if t00!=False:
			tmin=np.maximum(tmin,t00)
		#tmax=tmin+np.timedelta64(14,'D')
		tmax=tmin+np.timedelta64(7,'D')
		#inds=stations[setup_names[0]]['time']<tmax
		inds=[]
		for key in sources[1:]:
			inds.append( (stations[key]['time']>=tmin) & (stations[key]['time']<=tmax))

		inds2=(stations['TG'][name]['time']>=tmin) & (stations['TG'][name]['time']<=tmax)
		#plt.clf()
		plt.close()
		plt.plot(stations['TG'][name]['time'][inds2],stations['TG'][name]['zeta'][inds2],station_style)
		for nr,key in enumerate(sources[1:]):
			plt.plot(stations[key]['time'][inds[nr]],stations[key]['zeta'][:,i][inds[nr]],'-')
		plt.legend(sources,ncol=4,loc='upper center',frameon=False)					
		plt.gcf().autofmt_xdate()
		plt.title(name)
		plt.ylabel('$\zeta $ [m]')
		plt.grid()
		fname='3_TG'+name+pic_format
		plt.tight_layout()
		plt.savefig(outdir+fname,dpi=dpivalue)
		caption='Time Series at Station ' + name
		figures.append(fname)
		captions.append(caption)
################################################	

## filter longterm TS
for key in sources[1:]:
	dates=stations[key]['time']
	zeta=stations[key]['zeta']
	zeta_filter=np.zeros(zeta.shape)
	for i,date in enumerate(dates):
		ind=np.abs(dates-date)<= running_avg_window/2	
		zeta_filter[i,:]=zeta[ind,:].mean(axis=0)
	stations[key]['zeta_filter']=zeta_filter.copy()

###slow 
#for i, name in enumerate(names):
#	start = time.time()
#	print('running average for ' + name)	
#	dates=stations['TG'][name]['time']
#	zeta=stations['TG'][name]['zeta']
#	# running mean
#	stations['TG'][name]['zeta_filter']=np.asarray([zeta[((dates[ti]-running_avg_window/2) <= dates) & (dates <= dates[ti]+running_avg_window/2)].mean() for ti in range(len(dates))])
#	end = time.time()
#	print("prepate TG Elapsed = %s" % (end - start))
#
###slow 
#wndw_half=running_avg_window/2
#for i, name in enumerate(names):
#	start = time.time()
#	print('running average for ' + name)	
#	dates=stations['TG'][name]['time']
#	zeta=stations['TG'][name]['zeta']
#	# running mean
#	zeta_filter=np.asarray([zeta[((dates[ti]-wndw_half) <= dates) & (dates <= dates[ti]+wndw_half)].mean() for ti in range(len(dates))])
#	stations['TG'][name]['zeta_filter']=zeta_filter
#	end = time.time()
#	print("prepate TG Elapsed = %s" % (end - start))

wndw_half=running_avg_window/2	
#@jit(nopython=True)
def runnin_mean(dates,zeta,wndw_half):
	zeta_filter=np.zeros(zeta.shape)
	for ti in range(len(dates)):
		zeta_filter[ti]=np.mean(zeta[((dates[ti]-wndw_half) <= dates) & (dates <= dates[ti]+wndw_half)])	
	return zeta_filter
#np.abs(dates-dates[ti])<= windw_half
	
##slow 
wndw_half=running_avg_window/2
for i, name in enumerate(names):
	start = time.time()
	print('running average for ' + name)	
	dates=stations['TG'][name]['time']
	zeta=stations['TG'][name]['zeta']
	stations['TG'][name]['zeta_filter']=runnin_mean(dates,zeta,wndw_half)
	end = time.time()
	print("prepate TG Elapsed = %s" % (end - start))	
	
	
	
###### compare data	non normalized 
if monthly_subplots:
	print('making monthly subplots of complete time series')
	
	for i, name in enumerate(names):
		plt.clf()
		
		data=[stations['TG'][name]['zeta_filter']]#[:,0]
		vmin=data[0].min()
		vmax=data[0].max()
		if len(sources[1:]) > 0:
			k=0
			for i,stp in enumerate(sources[1:]):
				k=i+1
				data.append( stations[stp]['zeta_filter'][:,i])
				data[k]=np.ma.masked_array(data[k],mask=np.isnan(data[k]))
		#vmin=np.nanquantile(np.concatenate(data),0.1)
		#vmax=np.nanquantile(np.concatenate(data),0.99)	
		vmin=np.nanmin(np.concatenate(data))
		vmax=np.nanmax(np.concatenate(data))	


		Times=stations['TG'][name]['time']
		plt.clf()
		for month in range(1,13): 
			plt.subplot(4,3,month)
			inds=np.asarray([date.astype(object).month ==month for date in Times ])
			plt.plot(Times[inds],data[0][inds],'b',linewidth=1.5)
			for nr,stp in enumerate(sources[1:]):
					time_model=stations[stp]['time']
					inds=np.asarray([np.int(str(date)[5:7])==month for date in time_model ])
					plt.plot(time_model[inds],data[nr+1][inds],'--')
			if month != 2:
				plt.title('month: '+str(month))
			else:
				plt.title(name)
			plt.grid()
			plt.ylim((vmin,vmax))
			if month%3 !=1:
				plt.tick_params(axis='y',labelleft=False)  
			plt.tick_params(axis='x',labelbottom=False)  	
			if month==12:
				plt.legend(sources,loc='upper center',ncol=3,frameon=False)
		fig=plt.gcf()
		fig.set_size_inches(8,8,forward=True)
		#plt.gcf().autofmt_xdate()
		#plt.ylabel('$\zeta $ [m]')
		plt.tight_layout()	
		fname='4_TG'+name+pic_format
		plt.savefig(outdir+fname,dpi=dpivalue)
		plt.close()	
		caption='Monthly timeseries at station '+str(name) +'. Time series where filtered with a running window of '+ str(running_avg_window)
		figures.append(fname)
		captions.append(caption)	

if taylor_diag:
	
	key=sources[1]
	samples=[[ [stations[key]['stats']['std_rel'][name],stations[key]['stats']['cor'][name],name[:3]] for name in names ]]
	for key in sources[2:]:
		samples.append([ [stations[key]['stats']['std_rel'][name],stations[key]['stats']['cor'][name],name[:3]] for name in names ])

	for fname, extend in zip(['5_taylorA.eps','5_taylorB.eps'],[True,False]):
		plt.clf()
		dia=plotTaylor(samples[0],stdref=1,extend=extend) #negative
		for point in dia.samplePoints:
			#point.set_color(colors[0])
			point.set_markeredgecolor(colors[0])
			point.set_markersize(13)

		for nr,key in enumerate(sources[2:]):
			#Add models to Taylor diagram
			for i,(stddev, corrcoef, name) in enumerate(samples[nr+1]):
				i,stddev, corrcoef, name
				dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr+1,:], mec=colors[1+nr,:],label=name)
		phis=[plt.plot(0,0,'.',markeredgecolor=colors[nr,:]) for nr in range(len(sources[2:])+1)]		
				
		plt.legend(sources[1:],loc='lower right')
		#fname='5_taylorA.eps'	
		figures.append(fname)
		caption='Taylor diagram of full signals interpolated to Tide Gauge Timesteps. Black: Schism, Red Amm15.'
		captions.append(caption)
		plt.savefig(outdir+fname,dpi=dpivalue)
		plt.close()

	
	#plt.clf()
	#dia=plotTaylor(samples[0],stdref=1,extend=False) #negative
	#for nr,key in enumerate(sources[2:3]):
	#	#Add models to Taylor diagram
	#	for i,(stddev, corrcoef, name) in enumerate(samples[nr]):
	#		i,stddev, corrcoef, name
	#		dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
	#fname='5_taylorB.eps'	
	#figures.append(fname)
	#caption='Taylor diagram of full signals interpolated to Tide Gauge Timesteps. Black: Schism, Red Amm15.'
	#captions.append(caption)
	#plt.savefig(outdir+fname,dpi=dpivalue)
	#plt.close()

	
	ivlaid=~np.isnan(list(stations[sources[-1]]['stats']['bias'].values()))
	bias=[]
	rmse=[]
	cor=[]
	
	for stp in sources[1:]:
		bias.append(np.mean(np.abs(np.asarray(list(stations[stp]['stats']['bias'].values()))[ivlaid])))
		rmse.append(np.mean((np.asarray(list(stations[stp]['stats']['rmse'].values()))[ivlaid])))
		cor.append(np.mean((np.asarray(list(stations[stp]['stats']['cor'].values()))[ivlaid])))
	#k=0
	#plt.clf()
	#xpos=np.arange(3)
	#for k,stp in enumerate(sources[1:]):
	#	plt.bar(xpos+k/n,[bias[k],rmse[k],cor[k]],1/n,align='edge',label=stp)
	#plt.xticks(xpos+2/n,('bias','mae','cor'))	
	#fname='5.2_sation_average.eps'	
	#plt.ylabel('station averaged values')
	#figures.append(fname)
	#caption='Statiscal properties averaged over commonly covered stations'
	#captions.append(caption)
	#plt.savefig(outdir+fname,dpi=dpivalue)
	#plt.close()


	
# harmonic analysis
if tidal_bar_plots:
	#Amps={'TG':dict.fromkeys(consts),'schism':dict.fromkeys(consts),'amm':dict.fromkeys(consts)}
	#Phas={'TG':dict.fromkeys(consts),'schism':dict.fromkeys(consts),'amm':dict.fromkeys(consts)}
	Amps=dict.fromkeys(sources)
	Phas=dict.fromkeys(sources)
	
	name=names[0]
	utout=utide.solve(date2num(stations['TG'][name]['time']), u=np.asarray(stations['TG'][name]['zeta']), v=None, lat=coords[i,1])
	#consts=utout.name.copy()
	
	for key in Amps.keys():
		Amps[key]=dict.fromkeys(consts)
		Phas[key]=dict.fromkeys(consts)
	
	for const in consts:
		for const in consts:
			for key in sources: #'TG','schism','amm':
				Amps[key][const]=np.zeros(len(names))
				Phas[key][const]=np.zeros(len(names))
	
	utout_mod=dict.fromkeys(sources[1:])
	for i,name in enumerate(names):			
		print(name)
		utout=utide.solve(date2num(stations['TG'][name]['time']), u=np.asarray(stations['TG'][name]['zeta']), v=None, lat=coords[i,1])
		for src in sources[1:]:
			utout_mod[src]=utide.solve(date2num(stations[src]['time']), u=stations[src]['zeta'][:,i], v=None, lat=coords[i,1])
		#if True in np.isnan(zeta_amm[:,i]):
		#	utout_amm=utout_schism
		#	ammfactor=np.nan
		#else:
			#utout_amm=utide.solve(date2num(dates2), u=np.asarray(zeta_amm[:,i]), v=None, lat=coords[i,1])
			#ammfactor=1.0
		for const in consts:
			ind=np.where(utout.name==const)[0]
			if len(ind)>0:
				Amps['TG'][const][i]=utout.A[ind]
				Phas['TG'][const][i]=utout.g[ind]
			else:	
				Amps['TG'][const][i]=np.nan
				Phas['TG'][const][i]=np.nan

			for src in sources[1:]:
				ind=np.where(utout_mod[src].name==const)[0]
				if len(ind)>0:
					Amps[src][const][i]=utout_mod[src].A[ind]
					Phas[src][const][i]=utout_mod[src].g[ind]
				else:
					Amps[src][const][i]=np.nan
					Phas[src][const][i]=np.nan
				
			
#				ind=np.where(utout_schism.name==const)
#				Amps['schism'][const][i]=utout_schism.A[ind]
#				Phas['schism'][const][i]=utout_schism.g[ind]
#				ind=np.where(utout_amm.name==const)
#				Amps['amm'][const][i]=utout_amm.A[ind]*ammfactor
#				Phas['amm'][const][i]=utout_amm.g[ind]*ammfactor
#			
	for const in consts:
		plt.clf()
		plt.subplot(2,1,1)
		xbar=np.asarray(range(len(names)))
		phs=[]
		
		phs.append(plt.bar(xbar,Amps['TG'][const],width=0.25,align='edge',label='TideGauge',color='b',tick_label=shortnames))
		for i,src in enumerate(sources[1:]):
			phs.append(plt.bar(xbar+0.25*(i+1),Amps[src][const],width=0.25,align='edge',label=src))
			#,color='r') #,tick_label=valid_names
			#ph3=plt.bar(xbar+0.5,Amps['amm'][const],width=0.25,align='edge',label='Amm15',color='k') #,tick_label=valid_names
		plt.xticks(rotation=45)
		#plt.xticks([])
		plt.gca().set_xticklabels([])
		plt.grid()
		plt.yticks(np.round(np.linspace(0,plt.ylim()[1],4)*10)/10)
		#plt.legend([ph1,ph2,ph3],('TG','SCHISM','Amm15'),ncol=3,loc='upper center')
		plt.legend(phs,sources,ncol=3,loc='upper center',frameon=False,bbox_to_anchor=(0.5, 1.5))
		plt.ylabel('$A_{'+const+'}$ [m]')
		plt.tight_layout()
		#plt.xlabel('station')	

		plt.subplot(2,1,2)
		xbar=np.asarray(range(len(names)))
		#ph1
		phs[0]=plt.bar(xbar,Phas['TG'][const],width=0.25,align='edge',label='TideGauge',color='b',tick_label=shortnames)
		for i,src in enumerate(sources[1:]):
			phs.append(plt.bar(xbar+0.25,Phas[src][const],width=0.25,align='edge',label=src)) #,tick_label=valid_names
			#phs[i+1]=plt.bar(xbar+0.25,Phas[src][const],width=0.25,align='edge',label=src) #,tick_label=valid_names
		#ph3=plt.bar(xbar+0.5,Phas['amm'][const],width=0.25,align='edge',label='Amm15',color='k') #,tick_label=valid_names
		plt.xticks(rotation=45)
		plt.grid()
		#plt.legend([ph1,ph2,ph3],('TG','SCHISM','Amm15'),ncol=3)
		plt.ylabel('$\phi_{'+const+'}$ [�]')
		plt.tight_layout()
		#plt.xlabel('station')	
		plt.yticks(np.asarray(np.round(np.linspace(0,plt.ylim()[1],4)),int))
		caption=const+' tidal amplitude (top) and phase (bottom)'
		fname=const+'amp_phase.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)
		plt.close()			
		figures.append(fname)
		captions.append(caption)

		
if tidal_difference_maps:
	for const in ['M2','M4']:		
		plt.clf()
		s.plot_domain_boundaries()
		vmin=99
		vmax=-99
		phs=[]
		for src in sources[1:]:
			damp=Amps[src][const]-Amps['TG'][const]
			#damp2=Amps['amm'][const]-Amps['TG'][const]
			#damp2=damp
			vmin=np.max((np.min((damp.min(),vmin)),-0.5))
			vmax=np.min((np.max((damp.max(),vmax)),0.5))
		for i,src in enumerate(sources[1:]):	
			phs.append(plt.scatter(x,y,s=Rcircle*0.25**i,c=damp))
			plt.clim((vmin,vmax))

		plt.legend(phs,sources[1:],loc='upper center',ncol=2,frameon=False)
		ch=plt.colorbar()
		ch.set_label(const+ 'Amplitude difference')
		plt.set_cmap(cmap)
		for coord,name in zip(coords,names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		fname=const+'_amplitude_difference.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)
		caption='Map of '+const+ ' amplitude difference (model - data)'
		figures.append(fname)
		captions.append(caption)#


		plt.clf()
		s.plot_domain_boundaries()
		vmin=9999
		vmax=-9999
		
		for src in sources[1:]:
			dphi=Phas[src][const]-Phas['TG'][const]
			dphi=-np.sign(dphi)*((np.abs(dphi)>180)*360-np.abs(dphi))
			vmin=np.min((dphi.min(),vmin))
			vmax=np.max((dphi.max(),vmax))
		phs=[]	
		for i,src in enumerate(sources[1:]):	
			phs.append(plt.scatter(x,y,s=Rcircle*0.25**i,c=dphi))
			plt.clim((vmin,vmax))

		plt.legend(phs,sources[1:],loc='upper center',ncol=2,frameon=False)
		ch=plt.colorbar()
		ch.set_label(const+ 'phase difference')
		plt.set_cmap(cmap)
		for coord,name in zip(coords,names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		fname=const+'phase_lag.eps'
		plt.savefig(outdir+fname,dpi=dpivalue)	
		caption='Map of '+const+ ' phase lag (model -data)'
		figures.append(fname)
		captions.append(caption)

	# create latex + pdf
if put_pics_to_texdoc:
		print('generating tex doc and pdf')
		os.chdir(outdir)
		#techit(latexname,latextitle,latextext)
		techit(latexname,latextitle,latextext,figures,captions)
		if os.path.isfile(latexname[:latexname.rindex('.')]+'.pdf'):
			 os.remove(latexname[:latexname.rindex('.')]+'.pdf')
		os.system('pdflatex '+latexname)
		print('done generating '+latexname[:latexname.rindex('.')]+'.pdf')


## write error stats to csv
#M=[['station name','mean data','mean SCHISM','mean Amm','stddata','std schism','std Amm','bias schism','bias Amm','rmse schism','rmse Amm','correlation schism','correlation Amm']]
#for name in names:
#	M.append([name,ZETAmean[name],zeta_schismmean[name],zeta_ammmean[name],stddata[name],stdmodel[name],stdmodel2[name],bias[name],bias2[name],rmse[name],rmse2[name],R[name],R2[name]])
#M=np.asarray(M)
#with open(outdir+'errorstats.csv', 'w') as csvFile:
#	writer = csv.writer(csvFile)
#	for i in range(M.shape[0]):
#		writer.writerow(M[i,:])
#csvFile.close()
		
#
pickle.dump(stations,open(outdir+"stations_extracted","wb"))


## compare amm15 and schism at boundaries
obd_compare=False
if obd_compare:
	## compasrions at open boundary
	ibd=np.asarray(s.bdy_segments[0][::50])
	bdcoords=[(s.lon[i-1],s.lat[i-1]) for i in ibd]
	bdnn=amm15.tree.query(bdcoords)[1]
	ii2,jj2=np.unravel_index(bdnn,amm15.LON.shape)

	date=dates2[0]
	amm15.update(date)
	tt=amm15.nc['time'][:]
	ssh=amm15.nc[amm15.varnames['ssh']][:][:,ii2,jj2]

	for date in dates2[nt2:nt2*10+1:nt2]:
		print('loading cmems ' + str(date))
		amm15.update(date)
		tt=np.concatenate((tt,amm15.nc['time'][:]))
		ssh=np.concatenate((ssh,amm15.nc[amm15.varnames['ssh']][:][:,ii2,jj2]))
		
	ncs=MFDataset(schismfiles[:10])	
	ssh2=ncs['elev'][:,ibd]
	tt2=ncs['time'][:]


	date=dates2[0]
	amm15.update(date)
	tt=amm15.t0+dt2*np.arange(len(tt))
	tt2=reftime+dt.timedelta(seconds=modeldt)*np.arange(1,len(tt2)+1) 



	nsub=int(np.ceil(np.sqrt(len(bdnn))))
	plt.close()
	for i in range(nsub):
		plt.subplot(nsub,1,i+1)
		ph1=plt.plot(tt,ssh[:,i],'b')
		ph2=plt.plot(tt2,np.asarray(ssh2[:,i]),'r--')
		plt.legend(('Amm','SCHISM'),frameon=False)
		plt.title('bd node ' + str(ibd[i]))
		plt.ylabel('ssh')
		plt.xlabel('time')
	#plt.tight_layout()
	plt.gcf().autofmt_xdate()
	plt.savefig('bnd_comparisons',dpi=dpivalue)

