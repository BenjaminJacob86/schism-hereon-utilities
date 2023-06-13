
# interpolate the schism output from the German Bight routine onto a structutre grid
# for display in the cosyna NCWMS system.
# Replace variable names by cf standar conventions
# transfer to the cosyna direcory using rsync


# connection info for cosyna
#ssh behrens1@kofserver1
#/cosyna/netcdf/PO-Flow/gb_schism_original/
#/cosyna/netcdf/PO-Flow/gb_schism_interpol/
#rsync -avLz ./schism_interp_20230101.nc behrens1@kofserver1:/cosyna/netcdf/PO-Flow/gb_schism_interpol/schism_interp_20230101.nc

#~/anaconda3/envs/schism/bin/python
#conda activate schism

import sys
import os
sys.path.insert(0,'/gpfs/work/routine-ksd/schism-routine/outputs4cosyna/')
from schism import *
from matplotlib import pyplot as plt
#plt.ion()
from glob import glob

# before0.01

########## Settings ############################################

indir='/gpfs/work/routine-ksd/schism-routine/outputs/'
outdir='/gpfs/work/routine-ksd/schism-routine/outputs4cosyna/'
# load setup
os.chdir('/gpfs/work/routine-ksd/schism-routine/')	
s=schism_setup()
os.chdir('/gpfs/work/routine-ksd/schism-routine/outputs4cosyna/')
# last three available files
files=list(np.sort(glob(indir+'schout_*.nc'))[-16:]) # last three available files

overwrite=False # overwirte interpolation files if existent


# grid and interpolation info
lon,lat=np.asarray(s.lon),np.asarray(s.lat)
dx=0.005  # 0.01  # 0.025  interpolation grid spacing in degree
lonmin=5.11608851  #lon.min()
lonmax=10.40483834 #lon.max()
latmin=53.03857668 #lat.min()	
latmax=55.62872055 #lat.max()
#x=np.arange(lon.min(),lon.max(),dx) 
#y=np.arange(lat.min(),lat.max(),dx)
x=np.arange(lonmin,lonmax,dx) 
y=np.arange(latmin,latmax,dx)
reload_weights=True # relaod weights 


#nomenclature: name replacements and unit additions
varnames=['elev', 'salt', 'temp', 'u', 'v'] #, 'velocity_magnitude']
names={'elev':'ssh', 'salt':'sss', 'temp':'sst', 'u':'u', 'v':'v', 'velocity_magnitude':'velocity_magnitude'}
units={'elev':'m', 'salt':'psu', 'temp':'deg C', 'u':'m/s', 'v':'m/s', 'velocity_magnitude':'m/s'}
std_names={'elev':'sea_water_elevation', 'salt':'surface_sea_water_salinity', 'temp':'surface_sea_water_temperature', 'u':'surface_eastward_sea_water_velocity', 'v':'surface_northward_sea_water_velocity', 'velocity_magnitude':'surface_sea_water_velocity_magnitude'}

###################################################################


# grid and interpolation info
longitude=x  # the variable name is used as actually name by xarray netcdf export
latitude=y
X,Y=np.meshgrid(x,y)
if reload_weights:
	parents=np.loadtxt('parents_4cosyna_{:s}b.txt'.format(str(dx).replace('.',''))).astype(int)
	ndeweights=np.loadtxt('nodeweights_4cosyna_{:s}b.txt'.format(str(dx).replace('.','')))
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

FillValue=-9999
############# Beginn actual interpolation #########################


# open refercne file
outfiles=[]
for infile in files:

	# output file
	fname=infile.split('/')[-1].replace('schout','schism_interp')
	outfile=outdir+fname
	outfiles.append(outfile)
	print(outfile)
	if (not os.path.exists(outfile)) or overwrite:
		dsin=xr.open_dataset(infile)
		dsin=dsin.sel(nSCHISM_vgrid_layers=-1)
		dsin['u']=dsin.hvel[:,:,0]
		dsin['v']=dsin.hvel[:,:,1]
		dsin['velocity_magnitude']=np.sqrt((dsin.hvel**2).sum(axis=-1))
		# reformat test
		time=dsin.time.values
		attrs=dsin.time.attrs
		t0=time[0]-(time[1]-time[0])
		seconds=(time-t0)/np.timedelta64(1,'s')
		t0=str(t0)[:19].replace('T',' ')
		tunit= "seconds since {:s}".format(t0) 
		attrs['units']=tunit
		attrs['_FillValue']=False
		for varname in varnames:
			print(varname)
			data=np.ones((24,)+X.shape)*FillValue
			for i in range(24):
				wet=dsin.wetdry_elem[i,:].values==0
				iuse=(parents!=-1) 
				valid_parents=parents[iuse]
				
				ivalid=wet[valid_parents]
				valid_parents=valid_parents[ivalid]
				weights=ndeweights[iuse,:][ivalid,:]
				# in domain wet elements
				temp=dsin[varname][i,:].values
				target_inds=np.where(iuse)[0][ivalid]
				ii,jj=np.unravel_index(target_inds,X.shape)
				data[i,ii,jj]=(temp[s.nvplt[valid_parents,:]]*weights).sum(axis=1)
			
			if varname=='elev':
				dslat = xr.DataArray(name="latitude",data=latitude,dims=["latitude"], coords=dict(
						latitude=(["latitude"],latitude),
					),
					attrs=dict(
						description="latitude",
						units="degrees_north",
						standard_name=std_names[varname],
						_FillValue=False,
					))
					
				dslon = xr.DataArray(name="longitude",data=longitude,dims=["longitude"], coords=dict(
						longitude=(["longitude"],longitude),
					),
					attrs=dict(
						description="longitude",
						units="degrees_east",
						standard_name=std_names[varname],
						_FillValue=False,
					))	
					
				#dstime = xr.DataArray(name="time",data=time,dims=["time"], coords=dict(
				#		time=(["time"],time),
				#	),
				#	attrs=dict(
						#description="longitude",
						#units="degrees east",
						#standard_name=std_names[varname],
						#_FillValue=None,
				#	))					
				dstime = xr.DataArray(name="time",data=seconds,dims=["time"], coords=dict(
						time=(["time"],seconds),
					),					attrs=attrs,
					)			
				da = xr.DataArray(name=names[varname],data=data,dims=["time","latitude", "longitude", ], coords=dict(
						latitude=dslat, #(["latitude"],dslat),
						longitude=dslon, #(["longitude"],dslon),
						time=dstime,
					),
					attrs=dict(
						description="ssh",
						units=units[varname],
						standard_name=std_names[varname],
						_FillValue=FillValue,
					))					
			
				#da = xr.DataArray(name=names[varname],data=data,dims=["time","latitude", "longitude", ], coords=dict(
				#		latitude=(["latitude"],latitude),
				#		longitude=(["longitude"],longitude),
				#		time=time,
				#	),
				#	attrs=dict(
				#		description="ssh",
				#		units=units[varname],
				#		standard_name=std_names[varname],
				#		_FillValue=FillValue,
				#	))
				da.to_netcdf(outfile,mode='w')
			else:
				da = xr.DataArray(name=names[varname],data=data,dims=["time","latitude", "longitude", ],
					attrs=dict(
						description=varname,
						units=units[varname],
						standard_name=std_names[varname],
						_FillValue=FillValue,
					))
				da.to_netcdf(outfile,mode='a')

# clean up old files				
all_files=list(np.sort(glob(outdir+'schism_interp_*.nc')))
delete_files=all_files[:all_files.index(outfiles[0])]
if len(delete_files) >0:
	[os.remove(delete_file) for delete_file in delete_files	]			
	
# rsync	
import subprocess
#rsync -avLz ./schism_interp_20230101.nc 
#behrens1@kofserver1:/cosyna/netcdf/PO-Flow/gb_schism_interpol/schism_interp_20230101.nc	
for outfile in outfiles:
	print(outfile)
	transfered_file='behrens1@kofserver1:/cosyna/netcdf/PO-Flow/gb_schism_interpol/'+outfile.split('/')[-1]	
	result = subprocess.run(["rsync", "-avLz", "{:s}".format(outfile), "{:s}".format(transfered_file)], capture_output=True)
	print(result)




