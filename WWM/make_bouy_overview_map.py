

import os
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *
from data_and_model_classes import bsh_spec, WW3_mesh
from glob import glob




def set_FS(Fincrease=1.4):
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
	###########################
set_FS(Fincrease=1.4)






bouydir='//gpfs/work/jacobb/data/validation/waves/bouy/bouys/' #bsh buoy directory containing subfolders with the bouys and in 
WWMdirs=['//gpfs/work/jacobb/data/validation/waves/wwm_veg_ref/']
keysOBS=['ELB','FNO','FN3','HEL']


os.chdir(WWMdirs[0])
s=schism_setup()


class bsh_ncbouy:
	def __init__(self,files):
		self.ds=xr.open_mfdataset(files)
		self.dates=self.ds['time'].values
	def	get_parameter(self,varname=''):
		return self.ds[varname.lower()].values


# load buoys
subdirs=glob(bouydir+'*')
subdirs=np.sort([sdir[sdir.rindex('/')+1:] for sdir in subdirs])
bouys=dict.fromkeys(subdirs)
for subdir in subdirs:
	files=list(np.sort(glob(bouydir+'{:s}/*.site'.format(subdir))))
	if len(files)>0:
		bouys[subdir]=bsh_spec(files=files)
		ncbouys=False
	else: #checheck netcdf 	
		files=list(np.sort(glob(bouydir+'{:s}/*.nc'.format(subdir))))
		bouys[subdir]=bsh_ncbouy(files=files)
			
		ncbouys=True
		
		
plt.ion()		

plt.clf()
s.plot_domain_boundaries()
for key in bouys.keys():
	B=bouys[key]
	lat=float(B.ds.StartLatitude.split()[0].replace(',','.'))
	lon=float(B.ds.StartLongitude.split()[0].replace(',','.'))
	plt.plot(lon,lat,'k.')
	plt.text(lon+0.0001,lat+0.0001,key,rotation=45)
plt.savefig('bouylocations.png',dpi=300)	
	
lats,lons=[],[]
for key in bouys.keys():
	B=bouys[key]
	lats.append(float(B.ds.StartLatitude.split()[0].replace(',','.')))
	lons.append(float(B.ds.StartLongitude.split()[0].replace(',','.'))	)
	
	