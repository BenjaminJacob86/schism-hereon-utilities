import os
import netCDF4
import sys
import csv
import matplotlib
from matplotlib import pyplot as plt
background=True
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
from numba import jit
import time
plt.ion()




setupdir=['/work/bg1186/from_Mistral/bg1186/g260094/SNS/SNSE3D_01a_r212_MSLRnsob_prc50/']
setupdir+=['/work/gg0028/g260114/RUNS/CheckAmpJohannes/'] 

ncdir=[setupdir[0] + 'outputs/'] 		  #   directory of schism nc output or 
ncdir+=[setupdir[1] + 'outputs/']

#os.chdir(setupdir[0])
#so1=schism_station_output()
#os.chdir(setupdir[1])
#so2=schism_station_output()


os.chdir(ncdir[0])
stat1=np.loadtxt('staout_1')
t1=stat1[:,0]
stat1=stat1[:,1:]
os.chdir(ncdir[1])
stat2=np.loadtxt('staout_1')
t2=stat2[:,0]
stat2=stat2[:,1:]

nt=min(len(t1),len(t2))

os.chdir(setupdir[1])

istat=0
plt.plot(t2/86400,stat1[:nt,istat])
plt.plot(t2/86400,stat2[:nt,istat],'--')
plt.legend('Mistral','Repro')
plt.title('station '+str(istat))
plt.xlabel('time [days]')
plt.ylabel('Zeta [m]')
plt.grid()
plt.savefig('station_compare_'+str(istat),dpi=300)



for istat in range(stat1.shape[1]):
	plt.clf()
	plt.subplot(2,1,1)
	plt.plot(t2/86400,stat1[:nt,istat])
	plt.plot(t2/86400,stat2[:nt,istat],'--')
	plt.legend(('Mistral','Repro'))
	plt.title('station '+str(istat))
	plt.xlabel('time [days]')
	plt.ylabel('Zeta [m]')
	plt.grid()
	plt.subplot(2,1,2)
	plt.title('Repro vs Mistral')
	plt.plot(t2/86400,stat2[:nt,istat]-stat1[:nt,istat])
	plt.savefig('station_compare_'+str(istat),dpi=300)


