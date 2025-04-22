import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *
from netCDF4 import Dataset
import xarray as xr

#export OMP_NUM_THREADS=4 
# execute within run directory expects outputs to be in subfolder combined
def find_parent_tri(gr,xq,yq,dThresh=1000):
	""" parents,ndeweights=find_parent_tri(gr,xq,yq,dThresh=1000)
		find parent for coordinates xq,yq within triangulation tris,xun,yun.
		return: parent triangle ids and barycentric weights of triangle coordinates.
		Works only for pure triangle grids	
	"""    
	#% Distance threshold for Point distance
	dThresh=dThresh**2
	
	xun=gr.x
	yun=gr.y
	tris=np.asarray(gr.faces)
	trinr=np.arange(tris.shape[0])
	
	trisX,trisY=xun[tris],yun[tris]
	#% orthogonal of side vecotrs
	SideX=np.diff(trisY[:,[0, 1, 2, 0]],axis=1)
	SideY=-np.diff(trisX[:,[0, 1, 2, 0]],axis=1)
	
	p=np.stack((xq,yq),axis=1)
	parent=-1*np.ones(len(p),int)
	for ip in range(len(p)):

			dx1=(p[ip,0]-trisX[:,0])
			dy1=(p[ip,1]-trisY[:,0])
			subind=(dx1*dx1+dy1*dy1) < dThresh # preselection
			subtris=trinr[subind]
			
			#% dot products
			parenti=(subtris[ (dx1[subind]*SideX[subind,0] + dy1[subind]*SideY[subind,0] <= 0) \
						   & ((p[ip,0]-trisX[subind,1])*SideX[subind,1] + (p[ip,1]-trisY[subind,1])*SideY[subind,1] <= 0) \
							 & ( (p[ip,0]-trisX[subind,2])*SideX[subind,2] + (p[ip,1]-trisY[subind,2])*SideY[subind,2] <= 0) ][:])
			if len(parenti):
				parent[ip]=parenti
	
	# tri nodes
	xabc=xun[tris[parent]]
	yabc=yun[tris[parent]]
	
	# barycentric weights
	divisor=(yabc[:,1]-yabc[:,2])*(xabc[:,0]-xabc[:,2])+(xabc[:,2]-xabc[:,1])*(yabc[:,0]-yabc[:,2])
	w1=((yabc[:,1]-yabc[:,2])*(xq-xabc[:,2])+(xabc[:,2]-xabc[:,1])*(yq-yabc[:,2]))/divisor
	w2=((yabc[:,2]-yabc[:,0])*(xq-xabc[:,2])+(xabc[:,0]-xabc[:,2])*(yq-yabc[:,2]))/divisor
	w3=1-w1-w2

	return parent,np.stack((w1,w2,w3)).transpose() 


class grid_from_nc:
	def __init__(self, nc):
		self.ncv=nc.variables
		self.faces=self.ncv['SCHISM_hgrid_face_nodes'][:,:3]-1
		self.x=self.ncv['SCHISM_hgrid_node_x'][:]
		self.y=self.ncv['SCHISM_hgrid_node_y'][:]

class grid_from_xr:
	def __init__(self, nc):
		self.ncv=nc.variables
		self.x=self.ncv['SCHISM_hgrid_node_x'][0,:].values
		self.y=self.ncv['SCHISM_hgrid_node_y'][0,:].values		
		# split qudas to try
		self.faces=self.ncv['SCHISM_hgrid_face_nodes'][0,:,:4].values-1
		faces2=[]
		
		if np.nansum(self.faces[:,-1]) > 0 :
			self.origins=[] # mapping from nodes to triangles for triplot
			for nr,elem in enumerate(self.faces):
				if np.isnan(elem[-1]):  # split quad into tris
					faces2.append(elem[:3])
					self.origins.append(nr)
				else:
					faces2.append(elem[[0,1,2]])
					faces2.append(elem[[0,2,3]])
					self.origins.append(nr)
					self.origins.append(nr)
			self.faces=np.asarray(faces2,int)[:,:3]                 
		else:
			self.faces=np.asarray(self.faces[:,:3],int)
			self.origins=np.arange(self.faces.shape[0])

def get_trans_data_from_nc(gr,nc,varname='temp',ti=0):
	zcor=nc['zcor'][ti,:]
	data=nc[varname][ti,:]
	zi=np.zeros((len(parent),58))
	dataTrans=np.zeros((len(parent),58))
	for i in range(58):
		zi[:,i]=(zcor[:,i][gr.faces[parent]]*w).sum(axis=1)
		# interpolation veertical linear
		# horizontal bathymetric
	for i in range(len(parent)):
		depths=zcor[gr.faces[parent[i],:],:]
		lvldata=data[gr.faces[parent[i],:],:]
		ibttm=(nc['node_bottom_index'][gr.faces[parent[i],:]]).max()-1
		dataTrans[i,ibttm:]=scipy.interp(zi[i,ibttm:],depths[0,ibttm:],lvldata[0,ibttm:])*w[i,0] + \
		scipy.interp(zi[i,ibttm:],depths[1,ibttm:],lvldata[1,ibttm:])*w[i,1] + \
		scipy.interp(zi[i,ibttm:],depths[2,ibttm:],lvldata[2,ibttm:])*w[i,2]	
		# sigma lay interp at surf and bottom
		dataTrans[i,-1]=(lvldata[:,-1]*w[i,:]).sum() 
		dataTrans[i,ibttm]=(lvldata[:,ibttm]*w[i,:]).sum() 
	dataTrans=np.ma.masked_array(dataTrans,mask=dataTrans==0.0)
	return zi,dataTrans


	
setupdir='/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/SEDIMENT/'	
ncdir=setupdir+'combined/'
bpfile='/gpfs/work/jacobb/data/bpfiles/Elbe.bp'
kmfile='/gpfs/work/jacobb/data/bpfiles/Elbe.km'
Ems.bp

bp=np.loadtxt(bpfile,skiprows=0)
km=np.loadtxt(kmfile,skiprows=0)
xq,yq=bp[:,1],bp[:,2]



os.chdir(setupdir)
s=schism_setup()
plt.ion()

# initiate file acces
schismfiles=glob(ncdir+'*.nc')	
nrs=[]
for file in schismfiles:
	nrs.append(int(file[file.rindex('_')+1:file.rindex('.')]))
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
ds=xr.open_mfdataset(schismfiles)

gr=grid_from_xr(ds)
parents,weights=find_parent_tri(gr,xq,yq,dThresh=1000)

lon,lat=np.asarray(s.lon),np.asarray(s.lat)
lon2=np.sum(lon[gr.faces[parents]] *  weights,axis=1)
lat2=np.sum(lat[gr.faces[parents]] *  weights,axis=1)


plt.ion()
s.plot_domain_boundaries()
plt.plot(lon2,lat2,'.-')

# subslect parent nodes
punq_nodes=np.unique(gr.faces[parents])
dstrans=ds.sel(nSCHISM_hgrid_node=punq_nodes)
aunq,ainb,bina=np.intersect1d(punq_nodes, gr.faces[parents], assume_unique=False, return_indices=True)
faces2=gr.faces[parents]
for nr,node in enumerate(punq_nodes):
	faces2[faces2==node]=nr
weights2=[np.tile(np.tile(weights[:,i],(dstrans['zcor'].shape[0],1)),(dstrans['zcor'].shape[2],1,1)).swapaxes(0,1).swapaxes(1,2) for i in range(3)]

	
sediment={'SED3D_{:d}'.format(i):0 for i in range(1,9)}
for key in sediment.keys():
	sediment[key]+=dstrans[key][:,faces2[:,0],:]*weights2[0]+dstrans[key][:,faces2[:,1],:]*weights2[1]+dstrans[key][:,faces2[:,2],:]*weights2[2]
zcor=dstrans['zcor'][:,faces2[:,0],:]*weights2[0]+dstrans['zcor'][:,faces2[:,1],:]*weights2[1]+dstrans['zcor'][:,faces2[:,2],:]*weights2[2]

coords={'SCHISM_hgrid_node_x':0,'SCHISM_hgrid_node_y':0}
for key in coords.keys():
	for i in range(len(weights2)):
		coords[key]+=ds[key][0,faces2[:,i]]*weights2[i][0,:,0]

l=np.sqrt(((coords['SCHISM_hgrid_node_x'].diff(0))**2 + (coords['SCHISM_hgrid_node_y'].diff(0))**2)).cumsum()
l-=l[0]
l=l.expand_dims({"nSCHISM_vgrid_layers":zcor.shape[-1]},axis=-1)
l.name='length_along_transect'		

		

plt.pcolor(zz,zcor[ti,:],sediment[key][ti,:])

zz=b['length_along_transect'][:]  

#l.encoding['units']='m'	
#l.to_netcdf(path='test.nc', mode='w')

{"my_variable": {"dtype": "int16", "scale_factor": 0.1, "zlib": True}, ...} 
 format=None, group=None, engine=None, encoding=None, unlimited_dims=None, compute=True, invalid_netcdf=False)

# export no netcdf
outfile='transect.nc'
zcor.to_netcdf(path=outfile, mode='w')
for key in sediment.keys():
	sediment[key].to_netcdf(path=outfile, mode='a')
#ds['time'].to_netcdf(path=outfile, mode='a')
l.to_netcdf(path=outfile, mode='a')
ds['time'].to_netcdf(path='test.nc', mode='a')
	

# SURFACE MEAN	
avg=dict.fromkeys(sediment.keys())
for key in avg.keys()
	avg[key]=sediment[key].mean(dim='time')
	

trans=xr.open_dataset('transect.nc')	
trans=trans.sel(nSCHISM_vgrid_layers=-1,time=slice(trans.time[0].values+np.timedelta64(30,'D'),trans.time[-1].values))
mean=trans.mean('time')
 
x0=660
x1=780 
iuse=  (660< km[:,0]/1000) & ( km[:,0]/1000 < 780)
for i,key in enumerate(sediment.keys()):
	if sd50[i] in np.asarray([0.06,0.07,0.25]):
		plt.figure()
		plt.plot(km[:,0]/1000,mean[key])
		plt.ylabel('<c_{:.2f}mm > [g/l]'.format(sd50[i]))
		plt.xlabel('km')
		plt.title(str(trans.time[0].values)[:10] + ' - ' + str(trans.time[-1].values)[:10])
		vmin=mean[key][iuse].min()
		vmax=mean[key][iuse].max()
		plt.xlim((780,660))	
		plt.ylim((vmin,vmax))
		plt.tight_layout()
 
# manifest
km=np.tile(km,(zcor.shape[-1],1)).swapaxes(0,1)

# sediment classes

import time
a=time.time()
for key in sediment.keys():
	sediment[key]=sediment[key].values
b=time.time()
print(b-a)
#[g/l]	

zcor=zcor.values

with open('/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/SEDIMENT/sediment.nml') as f:
	for line in f.readlines():
		if 'Sd50' in line:
			break
	sd50=line.replace('\n','').split('=')[1].split(',')	
	sd50=np.asarray([np.float(d.replace('d','e')) for d in sd50])


trans=Dataset('transect.nc')

plt.clf()
ti=0
plt.pcolor(km/1000,zcor[ti,:],sediment[key][ti,:])
ch=plt.colorbar(orientation='horizontal')
ch.set_label('C'+str(sd50[ised])+'mm [g/l]' )
plt.ylim((-30,2))
plt.clim((0,4))
plt.gca().set_aspect(aspect=2)
for ti in range(len(t)):
	plt.cla()
	plt.pcolor(km/1000,zcor[ti,:],sediment[key][ti,:],vmin=0,vmax=4)
	plt.ylabel('depth [m]')
	plt.ylabel('km')
	plt.ylim((-30,2))
	plt.title(str(t[ti])[:19])
	plt.gca().set_aspect(aspect=2)
	plt.savefig('{:04d}_{:s}.png'.format(ti,key),dpi=300)
