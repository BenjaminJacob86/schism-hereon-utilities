import numpy as np
from matplotlib import pyplot as plt
import sys
#sys.path.insert(0,'/p/home/jusers/jacob4/juwels/shared')
sys.path.insert(0,'/gpfs/home/jacobb/code/python')
sys.path.insert(0,'//gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/pf/g/g260114/git/schism-hzg-utilities/')


from schism import *
from glob import glob
import numpy as np
import xarray as xr

#### sites########
def load_wwm_site(file):
	with open(file) as f:
		lines=f.readline()
	header=lines.split()
	m=np.loadtxt(file,skiprows=1)
	
	data=dict.fromkeys(header)
	for i,key in enumerate(header):
		data[key]=m[:,i]
	return data

	
s=schism_setup()

#stie files	
files=glob('*.site')
sites={file.split('.')[0]:None for file in files}
for site,file in zip(sites.keys(),files):
	sites[site]=load_wwm_site(file)
# all values -999 why
dss=xr.open_dataset('ww3.2017_spec.nc')


## ww3
dss=xr.open_dataset('/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd_2017/longrun/ww3.2017_01_04.nc')


################# hist #####################
ds=xr.open_mfdataset(['wwm_hist_000{:d}.nc'.format(i) for i in range(1,6)])

ti=144
for ti in range(0,83,4):
	plt.clf()
	v=ds['HS'][ti,:].values
	ch,ph=s.plotAtnodes(v)	
	ph.set_label('HS [m]')
	plt.title(str(ds['ocean_time'][ti].values)[:19])
	plt.clim((0,4))
	plt.savefig('{:04d}_HS'.format(ti))



####
dss=xr.open_dataset('/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd_2017/longrun/ww3.2017_spec.nc')
lon=dss['longitude'].values
lat=dss['latitude'].values
plt.plot(lon[0,:],lat[0,:],'ko')




ds=xr.open_dataset('wwm_hist_0001.nc')


lon=dss['longitude'].values
lat=dss['latitude'].values
lon=lon[0,:]
lat=lat[0,:]
plt.clf()
v=ds['HS'][ti,:].values
ch,ph=s.plotAtnodes(v)	
ph.set_label('HS [m]')
plt.title(str(ds['ocean_time'][ti].values)[:19])
plt.plot(lon,lat,'ko')
plt.axis(ax)

test=xr.open_dataset('ww3.2017_fake_spec.nc')
lon2=test['longitude'].values
lat2=test['latitude'].values


for i in range(len(lon)):
	plt.text(lon[i],lat[i],str(i))




