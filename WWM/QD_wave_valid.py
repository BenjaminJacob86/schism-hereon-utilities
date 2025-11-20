import os
from glob import glob
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
plt.ion()
import datetime as dt
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')
# own libraries 
from data_and_model_classes import bsh_spec, WW3_mesh
from validation_statistics_and_plots import interp_to_data_time,QQplot, TaylorDiagram,plotTaylor
from matplotlib.collections import PolyCollection
from schism import *
from scipy.spatial import cKDTree
########### SETTINGS #########

#bouydir='/gpfs/work/jacobb/data/VALIDATION/bouys/'
#bouydir='/gpfs/work/jacobb/data/VALIDATION/ncbouys/' # ncbouys


 
ww3PointList='points.list'


bouyfile='/gpfs/work/jacobb/data/LUCIANA/WW3/NO_TS_MO_6201083.nc'
Lon,Lat,name=8.168055,53.835000,'Weser'

dsbouy=xr.open_dataset(bouyfile)

hs_bouy=dsbouy.VAVH.values
time_bouy=dsbouy.TIME.values



ww3dir='/gpfs/work/jacobb/data/LUCIANA/WW3/'
ww3_nc='/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202412.nc'
ww3mesh='NBSext_bl.msh'
ww3=WW3_mesh(ww3dir+ww3mesh)
ww3_search_tree=cKDTree(list(zip(ww3.x,ww3.y))) # nn search tree
dg_dist,ww3_nn=ww3_search_tree.query((Lon,Lat))
ds_ww3=xr.open_dataset(ww3_nc)
ds_ww3_nn=ds_ww3.sel(node=ww3_nn) # select for node
hs_ww3=ds_ww3_nn.hs.values
time_ww3=ds_ww3_nn.time.values




 

wwmdir='/gpfs/work/ksddata/ROUTINES_personal/SCHISM/jade_bay/schism-routinetest4Luciana/'
ncdir='/gpfs/work/ksddata/ROUTINES_personal/SCHISM/jade_bay/schism-routinetest4Luciana/outputs/'

os.chdir(wwmdir)
s=schism_setup()
s.init_node_tree(latlon=True)
dg_dist_wwm,wwm_nn=s.node_tree_latlon.query((Lon,Lat))
s.ds=schism_outputs_by_variable(ncdir=ncdir,varlist='out2d').ds

time_wwm=s.ds['out2d'].time
hs_wwm=s.ds['out2d'].sigWaveHeight.sel(nSCHISM_hgrid_node=wwm_nn).values

ww3.load_points(ww3dir+ww3PointList)
 
plt.plot(time_bouy, hs_bouy,'b.',label='boouy')
plt.plot(time_ww3, hs_ww3,'r.-',label='ww3')
plt.plot(time_wwm, hs_wwm,'m.-',label='wwm')
plt.legend()


ind=-1
Ti=s.ds['out2d'].time[ind].values



ww3_tind=np.argmin(np.abs(time_ww3-Ti))
time_ww3[ww3_tind]


ds_spec=xr.open_dataset('/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202412_spec.nc')

plt.figure()
plt.subplot(2,1,1)
ph,ch,ax=s.plotAtnodes(s.ds['out2d'].sigWaveHeight[ind,:],cmap=plt.cm.jet)
plt.plot(ds_spec.longitude,ds_spec.latitude,'r+')
s.plot_domain_boundaries(latlon=True,append=True)
plt.subplot(2,1,2)
xlim=ax.get_xlim()
ylim=ax.get_ylim()
ph2,ch2,ax2=ww3.plotAtnodes(ds_ww3['hs'][ww3_tind,:].values)
(7.9816755285, 8.7359873615, 53.0820031155, 53.8640312745)
plt.axis(xlim+ylim)
plt.plot(Lon,Lat,'bo')
plt.plot(ds_spec.longitude,ds_spec.latitude,'r+')
plt.suptitle(Ti)
plt.axis(xlim+ylim)
s.plot_domain_boundaries(latlon=True,append=True)
plt.clim((0,0.6))



ww3.open_point_output('/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202412_spec.nc')

# via mean
ds_spec=xr.open_dataset('/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202412_spec.nc')
#ds_spec=xr.open_dataset('/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202411_spec.nc')
# compute spectra

# Assume `ds` is your Dataset
efth = ds_spec.efth  # (time, station, frequency, direction)

# Convert direction from degrees to radians spacing:
dtheta = np.deg2rad(np.diff(ds_spec.direction).mean())  # radians
df = np.diff(ds_spec.frequency).mean()  # Hz

# Compute zeroth moment m0:
m0 = (efth.sum(dim='direction') * dtheta).sum(dim='frequency') * df

# Significant wave height:
Hs = 4.0 * np.sqrt(m0)

# Hs will have dims (time, station)
Hs


# from gh
HS[0]
#




ds_spec = xr.open_dataset('/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202412_spec.nc')
efth = ds_spec.efth  # (time, station, frequency, direction)

# Frequency bin widths
freq = ds_spec.frequency.values
df = np.diff(freq)
df = np.append(df, df[-1])  # repeat last width for last bin

# Direction bin width (convert degrees to radians if directions in degrees)
directions = ds_spec.direction.values
# Check units — assume degrees:
dtheta = np.deg2rad(np.diff(directions).mean())  # radians

# zeroth moment m0:
m0 = (efth.sum(dim='direction') * dtheta) * df
m0 = m0.sum(dim='frequency')

# significant wave height
Hs = 4.0 * np.sqrt(m0)



plt.figure()
plt.subplot(2,1,1)
plt.plot(Hs[:,:7])
plt.subplot(2,1,2)
plt.plot(Hs[:,7:])


plt.plot(ds_spec.longitude[0,6:],ds_spec.latitude[0,6:],'r+')


plt.figure()
s.plot_domain_boundaries(latlon=True)
plt.plot(xi,yi,'ko')

dist, nn=y

for nr,coord in enumerate(zip(ds_spec.longitude[0,:].values,ds_spec.latitude[0,:].values)):
   
    xi,yi=coord
    break
    plt.figure()
    dist, nn=s.node_tree_latlon.query((xi,yi))
    plt.plot(time_wwm,s.ds['out2d'].sigWaveHeight.sel(nSCHISM_hgrid_node=nn).values,'k',linewidth=2,label='WWM')
    plt.plot(efth.time,Hs[:,nr],label='WW3')
    plt.legend()
    
    
    break
    



ww3.open_point_output('/gpfs/work/jacobb/data/LUCIANA/WW3/ww3.202412_spec.nc')    


for nr,coord in enumerate(zip(ds_spec.longitude[0,:].values,ds_spec.latitude[0,:].values)):
   
    #break
   
    plt.close('all')
    xi,yi=coord
   
    dg_dist,ww3_nn=ww3_search_tree.query((xi,yi))
    ds_ww3_nn=ds_ww3.sel(node=ww3_nn) # select for node
    hs_ww3=ds_ww3_nn.hs.values
    time_ww3=ds_ww3_nn.time.values

    dist, nn=s.node_tree_latlon.query((xi,yi))

    plt.figure()
    #plt.plot(efth.time,Hs[:,nr],label='WW3_specInt')
    plt.subplot(1,2,1)
    s.plot_domain_boundaries(latlon=True,append=True)
    
    plt.plot(xi,yi,'ko',label='station')
    plt.plot(Lon,Lat,'kx',label='bouy') 
    plt.plot(s.lon[nn],s.lat[nn],'r+',label='bouy') 
    

    plt.subplot(1,2,2)

    plt.plot(time_wwm,s.ds['out2d'].sigWaveHeight.sel(nSCHISM_hgrid_node=nn).values,'k',linewidth=4,label='WWM')
    
    plt.plot(ww3.tspec,ww3.Hs[:,nr],label='WW3_specInt')
    plt.plot(time_ww3,hs_ww3,label='WW3_gridNN')
    plt.ylabel('time')
    #plt.hlines(s.depths[nn],time_ww3[0],time_ww3[-1],label='depths')
    plt.title('depth at wwm node: ' + str(s.depths[nn]))
    plt.gcf().autofmt_xdate()
    plt.legend()

    plt.pause(1)    
    
    




#coords=plt.ginput(1)
#xi2,yi2=coords[0]


#plt.plot(xi2,yi2,'m+',label='manual_sugestion')



dist, nn=s.node_tree_latlon.query((xi2,yi2))
plt.plot(s.lon[nn],s.lat[nn],'b+',label='nn to manual')

plt.legend()    
plt.close('all')
xi,yi=coord

dg_dist,ww3_nn=ww3_search_tree.query((xi,yi))
ds_ww3_nn=ds_ww3.sel(node=ww3_nn) # select for node
hs_ww3=ds_ww3_nn.hs.values
time_ww3=ds_ww3_nn.time.values

dist, nn=s.node_tree_latlon.query((xi,yi))

plt.figure()
#plt.plot(efth.time,Hs[:,nr],label='WW3_specInt')
plt.subplot(1,2,1)
s.plot_domain_boundaries(latlon=True,append=True)

plt.plot(xi,yi,'ko',label='station')
plt.plot(s.lon[nn],s.lat[nn],'b+',label='Schism-nn')
plt.legend()    

plt.subplot(1,2,2)

plt.plot(time_wwm,s.ds['out2d'].sigWaveHeight.sel(nSCHISM_hgrid_node=nn).values,'k',linewidth=4,label='WWM')

plt.plot(ww3.tspec,ww3.Hs[:,nr],label='WW3_specInt')
plt.plot(time_ww3,hs_ww3,label='WW3_gridNN')
plt.ylabel('time')
plt.gcf().autofmt_xdate()
plt.legend()

plt.pause(2)


plt.plot(time_bouy,hs_bouy,'b.')