import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
import xarray as xr
import os
from glob import glob
plt.ion()
s=schism_setup()


dirs=['/work/gg0028/g260114/RUNS/GermanBight/WAVE/FRCww3BetaMax1.43_BetaMax1.54/','/work/gg0028/g260114/RUNS/GermanBight/WAVE/FRCww3BetaMax1.8_BetaMax1.54/','/work/gg0028/g260114/RUNS/GermanBight/WAVE/FRCww3BetaMax1.8_BetaMax2.4/']


dir1='/work/gg0028/g260114/RUNS/GermanBight/routine_GB_wave/CodeAaron/'
dir2='/work/gg0028/g260114/RUNS/GermanBight/routine_GB_wave/CodeAronTestBmax/'

nc1=xr.open_mfdataset(list(np.sort(glob(dir1+'wwm_hist_????.nc'))))
nc2=xr.open_mfdataset(list(np.sort(glob(dir2+'wwm_hist_????.nc'))))

ti=0	
for ti in range(24,49):
	plt.clf()
	plt.subplot(2,2,1)
	s.plotAtnodes(nc1['HS'][ti,:])
	plt.suptitle(str(nc1['ocean_time'][ti].values))
	plt.clim((0,0.2))
	plt.title('WWMV Bmax 1.54')
	plt.subplot(2,2,2)
	s.plotAtnodes(nc2['HS'][ti,:])
	plt.suptitle(str(nc2['ocean_time'][ti].values))
	plt.clim((0,0.2))
	plt.title('WWMV Bmax 1.74')
	plt.subplot(2,2,3)
	s.plotAtnodes(nc2['HS'][ti,:]-nc1['HS'][ti,:])
	plt.suptitle(str(nc2['ocean_time'][ti].values))
	plt.clim((-0.1,0.1))
	plt.title('WWMV: Bmax 1.74 - Bmax1.54')
	plt.savefig('{:04d}_Bmax1.54_vsBmax174.png'.format(ti))





for dir in dirs: 
	os.chdir(dir)
	nc1=xr.open_dataset(dir+'wwm_hist_0001.nc')
	ti=2
	plt.clf()
	s.plotAtnodes(nc1['HS'][ti,:])
	plt.savefig('HS_BD_check.png')

	
	
dir='./'	
nc0=xr.open_dataset(dir+'wwm_hist_0001.nc')	
for ti in range(24):
	plt.clf()
	ph,ch=s.plotAtnodes(nc0['HS'][ti,:])
	ch.set_label('hs [m]')
	plt.clim((0,0.2))
	plt.savefig('{:04d}_hs'.format(ti),dpi=200)	

	
	
	
	
	
ti=0	
for ti in range(24):
	plt.clf()
	plt.subplot(2,1,1)
	s.plotAtnodes(nc0['HS'][ti,:])
	plt.suptitle(str(nc1['ocean_time'][ti].values))
	plt.clim((0,0.2))
	plt.title('WWM3 dt 300')
	plt.subplot(2,1,2)
	s.plotAtnodes(nc1['HS'][ti,:])
	plt.suptitle(str(nc1['ocean_time'][ti].values))
	plt.clim((0,0.2))
	plt.title('WWM5 dt 300')
	plt.savefig('{:04d}'.format(ti))

dir1='/mnt/lustre01/work/gg0028/g260114/RUNS/GermanBight/WAVE/'
#dir2='/mnt/lustre01/work/gg0028/g260114/RUNS/GermanBight/WAVE/FRCww3BetaMax1.8/'
dir2='/mnt/lustre01/work/gg0028/g260114/RUNS/GermanBight/WAVE/FRCww3BetaMax1.8_BetaMax2.4/'

dir2='/work/gg0028/g260114/RUNS/GermanBight/WAVE/test_bd_correct/'

dir2='./'
nc1=xr.open_dataset(dir1+'wwm_hist_0001.nc')
plt.ion()
ti=1
plt.clf()
s.plotAtnodes(nc1['HS'][ti,:])
plt.title('1')



nc1=xr.open_dataset(dir1+'wwm_hist_0001.nc')
nc2=xr.open_dataset(dir2+'wwm_hist_0001.nc')
plt.ion()
ti=2
plt.clf()
plt.subplot(2,2,1)
s.plotAtnodes(nc1['HS'][ti,:])
plt.title('1')
plt.subplot(2,2,2)
s.plotAtnodes(nc2['HS'][ti,:])
plt.title('2')
plt.subplot(2,2,3)
s.plotAtnodes(nc2['HS'][ti,:]-nc1['HS'][ti,:])
plt.savefig('difference_forcing_change')
plt.title('2-1')

(nc2['HS'][ti,:]-nc1['HS'][ti,:]).max()

#from netCDF4 import Dataset
#nc=Dataset('ww3.2017cat0103.spec.nc','a')
#nc['efth'][:]=0
#nc.sync()
#nc.close()


# BOundary forcing not applied in WWMV

# 
