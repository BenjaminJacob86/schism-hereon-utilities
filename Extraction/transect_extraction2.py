"""
export schism outputs interpolated to transect bp file
"""
import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
import dask
dask.delayed()
import xarray as xr
plt.ion()
import time

#export OMP_NUM_THREADS=4
########### Settings ##################################
setupdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/'	
ncdir=setupdir+'outputs01/'

interp='nn' # invers distance (inverse distance is along sigma layers and not recommended where)
varlist =[] # restrict to list of varnames. All variables from variblae files will be used more than one varialbe files ar only in out2d. to include in velocities, the X and Y variable names need to be added to the list


bpfile='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/transect_south_syslt.bp'  # if bpfile=='create' open interactive map
outname='transect_south_syslt.nc' # outputname for netcdf file
########################################################

### program start ##########
os.chdir(setupdir)
s=schism_setup()

## check out put format: schout_
if len(glob.glob(ncdir+'schout_*.nc'))>0:
	ds=schism_output2(ncdir)
	newio=False
else: # new io
	newio=True
	access=schism_outputs_by_variable(ncdir,varlist=varlist,max_stack=-1)

class bp_transect_nn():
 
    def __init__(self,s,ds,x,y,latlon=False,nn=True):
      """ bp_transect(s,ds,x,y,lonlat=False,nn=False) s:=schism_setup(), ds:= xarray handle for old io and return of schism_outputs_by_variable (access=schism_outputs_by_variable) for newio code  ,x,y coordinates
	  Constraut transect object with track tangent normal and tangents, 
	  and along across projection function. If nn = True limit input coordinates
	  to those with different nearest neighbouts, only, ie to 
	  get xarray output access as nearest neighbour for transect coordinates. """

      self.x,self.y=x,y	  
      s.init_node_tree(latlon=latlon)
      coords=list(zip(x,y))
      p=np.asarray(coords)
      dd,self.nn=s.node_tree_xy.query(coords)
      self.nn,iunq=np.unique(self.nn,return_index=True)
      self.lon,self.lat=np.asarray(s.lon)[self.nn],np.asarray(s.lat)[self.nn]
      self.x,self.y=self.x[iunq],self.y[iunq]

      if nn:
          p=p[iunq,:]
      # projection for velocties
      dl=np.sqrt((np.diff(p,axis=0)**2).sum(axis=1))
      #dl=np.sqrt(np.diff(x)**2+np.diff(y)**2)
      self.l=np.hstack((0,np.cumsum(dl)))
      
      if type(ds)==schism_outputs_by_variable:
          self.access=ds
          varlist=list(self.access.vardict.keys())+self.access.varlist
          self.io='new'
          for key in self.access.ds.keys():
            try:		  
                self.access.ds[key]=self.access.ds[key].sel(nSCHISM_hgrid_node=self.nn)
            except:	 #vecotr
                pass	
          for key in self.access.ds.keys():		
            vari_vec=key
            varX=vari_vec+'X'	  
            if  varX in varlist: #key[-1]=='Y':
                varY=vari_vec+'Y'	  
                self.access.ds[key]={vari_vec: xr.concat([self.access.ds[self.access.vardict[varX]][varX], self.access.ds[self.access.vardict[varY]][varY]], dim='ivs')}
                #from IPython import embed; embed()
          #self.ds=[ds.sel(nSCHISM_hgrid_node=self.nn) for ds in access.ds.keys]
          z=self.access.get('zCoordinates')
      else:
          self.io='old'
          self.ds=ds.sel(nSCHISM_hgrid_node=self.nn)	                      
          z=ds['zcor']

	  
      self.l=np.tile(self.l,[z.shape[-1],1]).T
      tangent=np.diff(p,axis=0)/np.tile(dl,(2,1)).T
      nor=np.vstack((tangent[:,1],-tangent[:,0] )).T
      tangent=np.vstack((tangent[0,:],tangent))
      nor=np.vstack((nor[0,:],nor))
      self.tangent=np.tile(tangent,(z.shape[-1],1,1)).swapaxes(0,1)
      self.nor=np.tile(nor,(z.shape[-1],1,1)).swapaxes(0,1)	  
	  
      if self.io=='new':
          self.tangent=self.tangent.swapaxes(1,2).swapaxes(0,1)
          self.nor=self.nor.swapaxes(1,2).swapaxes(0,1)
          self.ivs_axis=0
      else:	  
          self.ivs_axis=2	  
	
      if self.io=='new':	
	  
	  # make data array
          self.tangentX=xr.DataArray(self.tangent[0,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
          self.tangentY=xr.DataArray(self.tangent[1,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
          self.norX=xr.DataArray(self.nor[0,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
          self.norY=xr.DataArray(self.nor[1,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})#
          
          self.tangentX=self.tangentX.expand_dims({"time":len(self.access.get('zCoordinates').time)},axis=0)
          self.tangentY=self.tangentY.expand_dims({"time":len(self.access.get('zCoordinates').time)},axis=0)
          self.norX=self.norX.expand_dims({"time":len(self.access.get('zCoordinates').time)},axis=0)
          self.norY=self.norY.expand_dims({"time":len(self.access.get('zCoordinates').time)},axis=0)

		  
    def proj_hvel_along_accros(self,hvel):
      """ Projet hvel of selfect to along and across direction """	
      self.ualong=(hvel*self.tangent[:,:,:]).sum(axis=self.ivs_axis) 	
      self.uacross=(hvel*self.nor[:,:,:]).sum(axis=self.ivs_axis) 	
      return self.ualong,self.uacross	        
    def dump_to_nc(self,filename='transect.nc'):
      for i,key in enumerate(self.access.vardict.keys()):
          if (key[-1]!='X') and (key[-1]!='Y'):
            mode='w' if i==0 else 'a'
            tmp=self.access.get(key)	
            try:			
                tmp.to_netcdf(filename,mode=mode)	
            except:
                from IPython import embed; embed()
				
	  # get data for projection
      h=np.diff(self.access.get('zCoordinates') )
      h=self.access.get('zCoordinates').diff(dim='nSCHISM_vgrid_layers')
      hc=0.5*h[:,0:-1,:]+h[:,1:,:]
	  
      l=xr.DataArray(self.l,dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
      dx=l[:,:-1].diff(dim='nSCHISM_hgrid_node').expand_dims({"time":len(hc.time)},axis=0)
      Ai=dx*hc  #prism cross area #[:,:-1]				
				
      #l=np.sqrt(( (self.access.get('SCHISM_hgrid_node_x').diff(0))**2 + (self.access.get('SCHISM_hgrid_node_y').diff(0))**2)).cumsum()
      l-=l[0]
      #l=l.expand_dims({"nSCHISM_vgrid_layers":self.access.get('zCoordinates').shape[-1]},axis=-1)
      l.name='length_along_transect'		
      l.to_netcdf(filename,mode=mode)	
      l.close()
      ds = xr.Dataset(data_vars=dict(lon=(["nSCHISM_hgrid_node"], self.lon),lat=(["nSCHISM_hgrid_node"], self.lat)))
      ds.to_netcdf(filename,mode=mode)	
	  # add across along alonv velocities:
      if self.io=='new':	
          # need to change variable	  
          across=(self.norX*ux + self.norY*uy)#.values
          along=(self.tangentX*ux + self.tangentY*uy)#.values
          along.name='AllongVel'
          across.name='AcrossVel'
          along.to_netcdf(filename,mode=mode)	
          across.to_netcdf(filename,mode=mode)	
          Ai=Ai.rename({'nSCHISM_hgrid_node':'nSCHISM_prism_x','nSCHISM_vgrid_layers':'nSCHISM_prism_z'})
          Ai.name='PrismCrossArea'
          Ai.to_netcdf(filename,mode=mode)
          A=Ai.sum(axis=(1,2))
          A.name='CrossArea'
          A.to_netcdf(filename,mode=mode)
          vi=(across[:,0:-1,1:] + across[:,1:,1:] + across[:,0:-1,0:-1] + across[:,1:,0:-1])/4
          vi=vi.rename({'nSCHISM_hgrid_node':'nSCHISM_prism_x','nSCHISM_vgrid_layers':'nSCHISM_prism_z'})
          vi.name='PrismAveragedVelocity'
          #vi.to_netcdf(filename,mode=mode)
          Qi=Ai*vi
          Qi.name='PrismCrossFlow'
          Qi.to_netcdf(filename,mode=mode)
          Q=Qi.sum(axis=(1,2))#.values
          Q.name='CrossSectionNetFlow'
          Qi.to_netcdf(filename,mode=mode)
          Qp=Qi.where(Qi>0).sum(axis=(1,2))#.values
          Qp.name='CrossSectionPositiveFlow'
          Qp.to_netcdf(filename,mode=mode)
          Qn=Qi.where(Qi<0).sum(axis=(1,2))#.values
          Qn.name='CrossSectionNegativeFlow'
          Qn.to_netcdf(filename,mode=mode)

# read bpfile 
#trans=bp_file(bpfile)	
#trans=bp_transect_nn(s,access,trans.x,trans.y)
#trans.dump_to_nc(outname)

#dsh=xr.open_dataset('transect_south_syslt.nc')

#access=schism_outputs_by_variable(ncdir,varlist=varlist,max_stack=-1)
access=schism_outputs_by_variable(ncdir,varlist=varlist+['zCoordinates','horizontalVelX','horizontalVelY'],min_stack=10,max_stack=10)
trans=bp_file(bpfile)	
trans=bp_transect_nn(s,access,trans.x,trans.y,nn=True)

trans.dump_to_nc(filename='test_transect.nc')

test=xr.open_dataset('test_transect.nc')

testi=test.sel(time=test.time[0])
plt.figure()
plt.clf()
plt.pcolor(testi.length_along_transect,testi.zCoordinates,testi.zCoordinates)
plt.plot(testi.length_along_transect[:,0],testi.zCoordinates[:,0],'k',linewidth=2)
plt.colorbar()

plt.figure()
plt.clf()
plt.subplot(2,1,1)
plt.pcolor(testi.length_along_transect,testi.zCoordinates,testi.AllongVel)
plt.plot(testi.length_along_transect[:,0],testi.zCoordinates[:,0],'k',linewidth=2)
plt.colorbar()
plt.subplot(2,1,2)
plt.pcolor(testi.length_along_transect,testi.zCoordinates,testi.AcrossVel)
plt.plot(testi.length_along_transect[:,0],testi.zCoordinates[:,0],'k',linewidth=2)
plt.colorbar()

plt.figure()
plt.clf()
plt.subplot(2,1,1)
plt.pcolor(testi.length_along_transect,testi.zCoordinates,testi.horizontalVelX[0,:])
plt.plot(testi.length_along_transect[:,0],testi.zCoordinates[:,0],'k',linewidth=2)
plt.colorbar()
plt.subplot(2,1,2)
plt.pcolor(testi.length_along_transect,testi.zCoordinates,testi.horizontalVelX[1,:])
plt.plot(testi.length_along_transect[:,0],testi.zCoordinates[:,0],'k',linewidth=2)
plt.colorbar()

# create transect
plt.figure()
s.plotAtnodes(s.depths,latlon=False)
n=2
plt.title("click coordinates in Fig.1. Press ESC when finished")
print("click coordinates in Fig.1. Press ESC when finished")
coords=plt.ginput(n=2,show_clicks='True')
plt.plot(coords[0][0],coords[0][1],'r+')
plt.plot(coords[-1][0],coords[-1][1],'r+')

xy=np.asarray(coords)
dxdy=np.diff(xy,axis=0)
dr=np.sqrt((dxdy**2).sum(axis=1))
r=dxdy[:,1]/dxdy[:,0]
from tkinter import filedialog
import tkinter as tk
# launch gui
root = tk.Tk()
minavgdist=tk.simpledialog.askfloat(title='interp dx', prompt='Enter distance [hgrid.gr3 units i.e. m] for between point interpolation: [0:= keep points]',initialvalue=250,minvalue=0)				
#self.minavgdist=np.float(input('Enter distance [m ord degree if ics=2] for between point interpolation:'))
coords1=[] 
if minavgdist !=0:
	for i in range(len(coords)-1):
		dx=np.linspace(0,dxdy[i,0],int(np.floor(dr[i]/minavgdist)))
		coords1+=([(xi[0][0],xi[0][1]) for xi in zip((xy[i,:]+np.stack((dx, dx*r[i]),axis=1)))])
	coords=coords1											
else:		
	coords1=coords					

coords1=np.asarray(coords1)	
m=np.zeros((coords1.shape[0],coords1.shape[1]+2))
m[:,0]=np.arange(1,m.shape[0]+1)
m[:,1:3]=coords1
m[:,3]=1

name='mytrans.bp'
f=open(name,'w')
f.write(name+'\n')
f.write('{:d}\n'.format(m.shape[0]))
np.savetxt(f,m,fmt='%d %f %f %d')
f.close()

trans=bp_file('mytrans.bp')	
trans=bp_transect_nn(s,access,trans.x,trans.y,nn=True)
trans.dump_to_nc(filename='cross_transect.nc')




trans.access.get('horizontalVel')
l=l.expand_dims({"nSCHISM_vgrid_layers":self.access.get('zCoordinates').shape[-1]},axis=-1)


ualong,uacross=trans.proj_hvel_along_accros(hvel)

# operate 

ti=0
dx=np.diff(trans.l,axis=0)
zz=trans.access.get('zCoordinates')[ti,:].values
h=np.diff(zz,axis=1)
hvel=trans.access.get('horizontalVel')[:,ti,:].values


# UL + UR +LL +LR
vi=(uacross[0:-1,1:] + uacross[1:,1:] + uacross[0:-1,0:-1] + uacross[1:,0:-1])/4
Qi=Ai*vi
Qi.sum()


hl=h[0:-1,:]
hr=h[1:,:]
hc=0.5*hl+hr
Ai=dx[:,:-1]*hc  #prism cross area


h=np.diff(trans.access.get('zCoordinates') )
h=trans.access.get('zCoordinates').diff(dim='nSCHISM_vgrid_layers')
hc=0.5*h[:,0:-1,:]+h[:,1:,:]
l=xr.DataArray(trans.l,dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
dx=l[:,:-1].diff(dim='nSCHISM_hgrid_node').expand_dims({"time":len(hc.time)},axis=0)
Ai=hc*dx  # area

# make data array
trans.tangentX=xr.DataArray(trans.tangent[0,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
trans.tangentY=xr.DataArray(trans.tangent[1,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
trans.norX=xr.DataArray(trans.nor[0,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})
trans.norY=xr.DataArray(trans.nor[1,:],dims={'nSCHISM_hgrid_node','nSCHISM_vgrid_layers'})#

trans.tangentX=trans.tangentX.expand_dims({"time":len(trans.access.get('zCoordinates').time)},axis=0)
trans.tangentY=trans.tangentY.expand_dims({"time":len(trans.access.get('zCoordinates').time)},axis=0)
trans.norX=trans.norX.expand_dims({"time":len(trans.access.get('zCoordinates').time)},axis=0)
trans.norY=trans.norY.expand_dims({"time":len(trans.access.get('zCoordinates').time)},axis=0)
# multiply
ux=trans.access.get('horizontalVelX')
uy=trans.access.get('horizontalVelY')
across=(trans.norX*ux + trans.norY*uy)#.values
along=(trans.tangentX*ux + trans.tangentY*uy)#.values

vi=(across[:,0:-1,1:] + across[:,1:,1:] + across[:,0:-1,0:-1] + across[:,1:,0:-1])/4
Qi=Ai*vi
Q=Qi.sum(axis=(1,2)).values





# load in memory or excute ...


alon=trans.tangentX*





def proj_hvel_along_accros(self,hvel):
  """ Projet hvel of transect to along and across direction """	
  self.ualong=(hvel*self.tangent[:,:,:]).sum(axis=self.ivs_axis) 	
  self.uacross=(hvel*self.nor[:,:,:]).sum(axis=self.ivs_axis) 	

vi=  #prism mean flow  

self.ualong=(hvel*transect.tangent[:,:,:]).sum(axis=2) 	
self.uacross=(hvel*transect.nor[:,:,:]).sum(axis=2) 	
self.ualong,self.uacross	        


trans.tangent.swapaxes(1,2).swapaxes(0,1).shape




hR=zz(2:end,2:end)-zz(2:end,1:end-1);
hL=zz(1:end-1,2:end)-zz(1:end-1,1:end-1);

%hR=zz(2:end,2:end)-zz(2:end,1:end-1);
%hL=zz(1:end-1,2:end)-zz(1:end-1,1:end-1);



trans.l


% trhough flow
%dx=diff(xx(:,2:end));
%hR=zz(2:end,2:end)-zz(2:end,1:end-1);
%hL=zz(1:end-1,2:end)-zz(1:end-1,1:end-1);
%hC=0.5*(hR+hL);
%Ai=hC.*dx;
%Velocity SUmme Oben unten
%Flow
%Qi=Ai.*Vi;
%A=sum(Ai(:));
%Q=sum(Qi(:));

Ai=transect.data.Ai(:,:,ti)   % transect cross sectional areas by cell
NN2=vv(:,1:end-1)+vv(:,2:end);
Vi=(NN2(1:end-1,:)+NN2(2:end,:))/4; % acerage over transect cross sectional cells
Q=sum(nansum(Ai.*Vi))
ti=1


trans=bp_file(bpfile)	
trans=bp_transect_nn(s,access,trans.x,trans.y)


trans.dump_to_nc(outname)

ph,ch,ax=s.plotAtnodes(s.depths,latlon=False)
ax.plot(trans.x,trans.y)

plt.figure()