# gen Tem_1.th
# create river temperature forcing based on atmosphere

import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/') 
from schism import *

s=schism_setup()


# get temperature
# detect at which boundaries forcing has to be applied
with open('bctides.in') as f:
	river_temp_segs=[]
	bd=0
	for line in f.readlines()[4:]:
		print(line)
		splitted=line.split()
		isbd=np.sum([val.isdigit() for val in splitted[:5]])==5
		bd+=isbd #is openboundary 
		if isbd and splitted[3]=='1':
			river_temp_segs.append(bd-1)
			
refnodes=np.zeros(len(river_temp_segs),int)
for i,ibd in enumerate(river_temp_segs):
	bdnodes=np.asarray(s.bdy_segments[ibd])-1			
	refnodes[i]=bdnodes[np.int(len(bdnodes)/2)]			
	
lon,lat=np.asarray(s.lon),np.asarray(s.lat)	
coords=list(zip(lon[refnodes],lat[refnodes]))

from scipy.spatial import cKDTree
ds=xr.open_dataset('sflux/sflux_air_1.0001.nc')


lon1d=ds.lon.values.flatten()
lat1d=ds.lat.values.flatten()

sflux_nntree=cKDTree(list(zip(lon1d,lat1d)))
nn1d=sflux_nntree.query(coords)[1]

ii,jj=nn2d=np.unravel_index(nn1d,ds.lon.shape)


p=param('.//param.nml')

# ok not working for all netcdf froamts.
#ds2=xr.open_mfdataset(glob.glob('sflux/sflux_air_1.*.nc')[:2])


reftime=dt.datetime(int(p.get_parameter('start_year')),
            int(p.get_parameter('start_month')),
            int(p.get_parameter('start_day')),
            int(p.get_parameter('start_hour')),0,0)	

ntfile=len(ds.time) #file length varies
files=np.sort(glob.glob('sflux/sflux_air_1.*.nc'))
TEM=np.zeros((len(files)*ntfile,len(nn1d)+1))
ti0=0
for i,file in enumerate(files):
	dsi=xr.open_dataset(file)
	ntfile=len(dsi.time) #file length varies
	ti1=ti0+ntfile
	tair=dsi['stmp'].values[:,ii,jj]
	tair-=273.15 #convert celsius
	TEM[ti0:ti1,1:]=np.maximum(tair,0)
	
	basedate=np.asarray(dsi.time.base_date,int)
	dates=dt.datetime(basedate[0],basedate[1],basedate[2],basedate[3])+	dsi.time.values*dt.timedelta(days=1)
	TEM[ti0:ti1,0]=np.fix((dates-reftime)/dt.timedelta(seconds=1))
	ti0+=ntfile

# first zero
istart=np.where(TEM[:,0]==0)[0][0]
TEM=TEM[istart:]
	
# too long?	
iend=np.where(TEM[1:,0]==0)[0][0]	
TEM=TEM[:iend,:]
	
a=np.loadtxt('TEM_1.th')	

plt.clf()
plt.plot(a[:,0]/86400,a[:,1:])
plt.figure()
plt.plot(TEM[:,0]/86400,TEM[:,1:])





	
# interpoalte to times of flux.th
tq=np.loadtxt('flux.th')[:,0]

# find intersecting times may not work if shift
_,ainb,bina=np.intersect1d(tq,np.fix(TEM[:,0]),return_indices=True)

#otherwise requires interp.
TEM=TEM[bina,:]

np.savetxt('TEM_1_BJ.th',TEM)


xy=np.asarray(coords)
plt.ion()
s.plot_domain_boundaries()
plt.plot(xy[:,0],xy[:,1],'k+',markersize=14)


	
