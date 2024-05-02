
import sys
import os
sys.path.insert(0,'/gpfs/work/routine-ksd/schism-routine/outputs4cosyna/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hereon-utilities/')
from schism import *
from matplotlib import pyplot as plt
plt.ion()

########## Settings ############################################

cluster='levante'

if cluster=='strand':
    # load setup
    os.chdir('/gpfs/work/routine-ksd/schism-routine/')	
    s=schism_setup()
    os.chdir('/gpfs/work/routine-ksd/schism-routine/outputs4cosyna/lonlat')		
else:
    # load setup
    os.chdir('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/')	
    ncdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs01/'	
    s=schism_setup()
    os.chdir('/work/gg0028/g260114/RUNS/postproc/interp/')		
    rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'



#s.ds=schism_outputs_by_variable(ncdir,varlist=['out2d']).ds
#filename='schout_20230101.nc'  # file to interpolate from

# grid and interpolation info
lon,lat=np.asarray(s.lon),np.asarray(s.lat)
dx=0.005  # 0.01  # 0.025  interpolation grid spacing in degree
x=np.arange(lon.min(),lon.max(),dx) 
y=np.arange(lat.min(),lat.max(),dx)
reload_weights=True # relaod weights 


#nomenclature: name replacements and unit additions
names={'elevation':'elevation', 'salt':'sss', 'temp':'sst', 'u':'u', 'v':'v', 'velocity_magnitude':'velocity_magnitude','depth':'depth'}
units={'elevation':'m', 'salt':'psu', 'temp':'deg C', 'u':'m/s', 'v':'m/s', 'velocity_magnitude':'m/s','depth':'m'}
std_names={'elevation':'sea_water_elevation', 'salt':'surface_sea_water_salinity', 'temp':'surface_sea_water_temperature', 'u':'surface_eastward_sea_water_velocity', 'v':'surface_northward_sea_water_velocity', 'velocity_magnitude':'velocity_magnitude','depth':'depth'}

###################################################################


# grid and interpolation info
longitude=x  # the variable name is used as actually name by xarray netcdf export
latitude=y
X,Y=np.meshgrid(x,y)
if reload_weights:
	parents=np.loadtxt('parents_4cosyna_{:s}.txt'.format(str(dx).replace('.',''))).astype(int)
	ndeweights=np.loadtxt('nodeweights_4cosyna_{:s}.txt'.format(str(dx).replace('.','')))
else: #calculate weights for interpolation
	xq=X.flatten()
	yq=Y.flatten()
	parents,ndeweights=s.find_parent_tri(xq,yq,dThresh=0.02,latlon=True)
	#np.savetxt('parents_4cosyna.txt',parents)
	#np.savetxt('nodeweights_4cosyna.txt',ndeweights)
	np.savetxt('parents_4cosyna_0005.txt',parents)
	np.savetxt('nodeweights_4cosyna_0005.txt',ndeweights)

# plot test
#parents2=parents.reshape(X.shape)
#X=np.ma.masked_array(X,mask=parents2==-1)
#Y=np.ma.masked_array(Y,mask=parents2==-1)
#plt.plot(X,Y,'k+')
#s.plot_domain_boundaries(append=True)


############# Beginn actual interpolation #########################

# open refercne file
#dsin=xr.open_dataset(filename)

# go stack by stack

ncdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all2/'	

files=np.sort(glob.glob(ncdir+'*out2d*'))
file=files[0]


files1=glob.glob(ncdir+'*out2d_?.nc')
files2=glob.glob(ncdir+'*out2d_??.nc')

files=np.sort(files1+files2)

varnames=['elevation','depth']

p=param(rundir+'/param.nml')
reftime=dt.datetime(int(p.get_parameter('start_year')),int(p.get_parameter('start_month')),int(p.get_parameter('start_day')),int(p.get_parameter('start_hour')),0,0)


iuse=(parents!=-1) # intepolation parents
for file in files[1:90]:
    fname=file.split('/')[-1]
    print(file)    
    dsin=xr.open_dataset(file)
    time=dsin.time

    reftime=np.asarray(reftime,np.datetime64)
    timeval=reftime+time.values*np.timedelta64(1,'s')

    for varname in varnames:
        print(varname)
        
        if 'time' in dsin[varname].dims:

            data=np.ones((24,)+X.shape)*np.nan
            for i in range(24):
                wet=dsin.dryFlagElement[i,:].values==0
                valid_parents=parents[iuse]
                
                ivalid=wet[valid_parents]
                valid_parents=valid_parents[ivalid]
                weights=ndeweights[iuse,:][ivalid,:]
                # in domain wet elements
                temp=dsin[varname][i,:].values
                target_inds=np.where(iuse)[0][ivalid]
                ii,jj=np.unravel_index(target_inds,X.shape)
                data[i,ii,jj]=(temp[s.nvplt[valid_parents,:]]*weights).sum(axis=1)
        else:
            data=np.ones(X.shape)*np.nan
            valid_parents=parents[iuse]
            temp=dsin[varname].values
            target_inds=np.where(iuse)[0]            
            weights=ndeweights[iuse,:]
            ii,jj=np.unravel_index(target_inds,X.shape)
            data[ii,jj]=(temp[s.nvplt[valid_parents,:]]*weights).sum(axis=1)


        if varname=='elevation':
            da = xr.DataArray(name=names[varname],data=data,dims=["time","latitude", "longitude", ], coords=dict(
                    latitude=(["latitude"],latitude),
                    longitude=(["longitude"],longitude),
                    time=(["time"],timeval),
                ),
                attrs=dict(
                    description="ssh",
                    units=units[varname],
                    standard_name=std_names[varname],
                    _FillValue=np.nan,
                ))
            #da.to_netcdf('interp'+fname,mode='w')
            da.to_netcdf('interp'+fname,mode='w',encoding={'time': {'dtype': 'i4'}})
        elif varname=='depth':
            da = xr.DataArray(name=names[varname],data=data,dims=["latitude", "longitude", ], coords=dict(
                    latitude=(["latitude"],latitude),
                    longitude=(["longitude"],longitude),
                ),
                attrs=dict(
                    description="depth",
                    units=units[varname],
                    standard_name=std_names[varname],
                    _FillValue=np.nan,
                ))
            #da.to_netcdf('interp'+fname,mode='w')
            da.to_netcdf('interp'+fname,mode='a')
            
            
    #	else:
    #		da = xr.DataArray(name=names[varname],data=data,dims=["time","latitude", "longitude", ],
    #			attrs=dict(
    #				description=varname,
    #				units=units[varname],
    #				standard_name=std_names[varname],
    #				_FillValue=np.nan,
    #			))
    #		da.to_netcdf('schism_interp_20230101c.nc',mode='a')
