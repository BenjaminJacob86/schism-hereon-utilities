"""
Validate Model against Sentinel Satelite
"""
import os
from glob import glob
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
plt.ion()
import datetime as dt
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *
from matplotlib.path import Path
rundir='/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_900s_onpiece/'
from validation_statistics_and_plots import interp_to_data_time,QQplot, TaylorDiagram,plotTaylor
os.chdir(rundir)

####### Settings #########
#wwmdir=#
sentdir='/gpfs/work/ksddata/observation/remote/SENTINEL/Sentinel-3/SRAL/'


##########################
s=schism_setup()


wwmfiles=list(np.sort(glob('*hist*nc')))

dswwm=xr.open_mfdataset(wwmfiles)

import xarray as xr
ds=xr.open_dataset('/gpfs/work/ksddata/observation/remote/SENTINEL/Sentinel-3/SRAL/S3B_SR_2_LAN____20210723T195657_20210723T200656_20210723T211510_0599_055_071______LN3_O_NR_004.SEN3/standard_measurement.nc')


files=list(np.sort(glob(sentdir+'*_2021*.SEN3/standard*')))

files=list(np.sort(glob(sentdir+'*_2017*.SEN3/standard*')))


bd=[]
for ocn,lnd in zip(s.bdy_segments,s.land_segments):
	bd+=ocn+lnd
bd=np.asarray(bd)-1
lon,lat=np.asarray(s.lon),np.asarray(s.lat)
Poly=Path(list(zip(lon[bd],lat[bd])))


t0=np.datetime64('20170101')
t1=np.datetime64('20170401')

s.init_node_tree()
file=files[0]

sattrans=[]
wwmtrans=[]
count=0
for file in files:
	ds=xr.open_dataset(file)
	time=ds['time_20_ku'][:].values
	if ( (t0 < time ) & ( time < t1 ) ).sum(): # time match
		satlon,satlat=ds['lon_20_ku'][:],ds['lat_20_ku'][:]
		isin=Poly.contains_points(list(zip(satlon,satlat)))
		hssat=ds['swh_ocean_20_ku'][:][isin].values
		if isin.sum()>5:
			ic=np.int(np.ceil(len(time[isin])/2))
			wwm_interp=dswwm.interp(ocean_time=time[isin][ic])
			#break
			count+=1
			plt.clf()
			plt.subplot(2,1,1)
			hs=wwm_interp['HS'].values
			vmin=0
			vmax=6
			s.plotAtnodes(hs)
			plt.clim((0,5))
			#vmin,vmax=hs.min(),hs.max()
			#s.plot_domain_boundaries(append=True)
			nn=s.node_tree_latlon.query(list(zip(satlon[isin],satlat[isin])))[1]
			plt.title(str(time[isin][ic]))
			plt.scatter(satlon[isin],satlat[isin],s=4,c=hssat,vmin=vmin,vmax=vmax,cmap=plt.cm.jet)
			plt.clim((0,5))
			#plt.plot(satlon[isin],satlat[isin],'.')
			plt.subplot(2,1,2)
			#ds['swh_ocean_20_ku'][:][isin].plot()
			plt.plot(satlat[isin],hssat)
			plt.plot(satlat[isin],wwm_interp['HS'][nn])
			plt.grid()
			plt.legend(('Sentinel 1','WWM'))
			plt.ylabel('hs [m]')
			plt.tight_layout()
			plt.savefig('sentinel_wwm_'+str(time[isin][ic])[:10],dpi=300)
			sattrans.append(hssat)
			wwmtrans.append(hs[nn])

			
sattrans=np.hstack(sattrans)
wwmtrans=np.hstack(wwmtrans)
QQplot(sattrans,wwmtrans,obsname='Sentinel 1')
