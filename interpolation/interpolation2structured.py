#!/gpfs/home/jacobb/miniforge3/envs/schism-env/bin/python

import sys
import os

#sys.path.insert(0, '/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine/')
#sys.path.insert(0, '//work/gg0028/SCHISM/schism-hereon-utilities/')
from schism import *
from matplotlib import pyplot as plt
from glob import glob
import datetime as dt

# before0.01

########## Settings ############################################

indir = '/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all2/'
outdir = '/work/gg0028/g260114/EDITO/REF_sim_2017/'

# load setup
os.chdir('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/')

s = schism_setup()
os.chdir(outdir)


# grid and interpolation info
overwrite = False  # overwirte interpolation files if existent
lon, lat = np.asarray(s.lon), np.asarray(s.lat)
dx = 0.005  # 0.01  # 0.025  interpolation grid spacing in degree
lonmin = 5.11608851  # lon.min()
lonmax = 10.40483834  # lon.max()
latmin = 53.03857668  # lat.min()
latmax = 55.62872055  # lat.max()
x = np.arange(lonmin, lonmax, dx)# x=np.arange(lon.min(),lon.max(),dx)
y = np.arange(latmin, latmax, dx) # y=np.arange(lat.min(),lat.max(),dx)
reload_weights = True  # relaod weights
FillValue = -9999
dThresh=0.5   # allowed distance for point pre selection [for coarer resolution might need to be incresaed to avoid gaps]

# change this for the globa data settings
global_atts={}
global_atts['title'] = 'Gridded multivariable surface data from German Bight forecast with coupled SCHISM-WMM model components of GCOAST'  # Non-empty string
global_atts['creation_date'] = f'Created on {dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'  # Non-empty string
global_atts['Conventions'] = 'CF-1.8'  # This can be empty or another relevant value        
global_atts['institution'] = 'Helmholtz-Zentrum Hereon, Institute of Coastal Systems, Germany'  # This can be empty or another relevant value        
global_atts['originator']  = 'Benjamin Jacob'
global_atts['contact']  = 'Benjamin.Jacob@hereon.de'
global_atts['department'] = 'https://www.hereon.de/institutes/coastal_systems_analysis_modeling/hydrodynamics_data_assimilation/team/index.php.en'  # This can be empty or another relevant value        
global_atts['source']  = 'SCHISM_v5.11.0'
global_atts['crs']  = 'EPSG:4326'
global_atts['field_type']  = 'instantaneous'
global_atts['forecast_type']  = 'hindcast'
global_atts['forecast_range']  = '2017'
global_atts['licence']  = 'https://spdx.org/licenses/CC-BY-4.0.html'
global_atts['domain_name']  = 'German Bight'
global_atts['history']  = f"; SCHISM output processed via interpolation4HCDC.py on {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"


# make sure out2d is always first in this list
varfiles = dict.fromkeys(['out2d', 'temperature', 'salinity', 'horizontalVelX', 'horizontalVelY'])
#for key in varfiles:
#    varfiles[key] = list(np.sort(glob(indir + key + '_*.nc'))[-4:])  # last three available files

for key in varfiles:
	varfiles[key]=[]
	for ndigits in range(1,4):
		varfiles[key] += list(np.sort(glob(indir + key + '_'+'?'*ndigits+'.nc')))  # last three available files



# todaystr=dt.datetime.today().strftime('%Y%m%d')
## temp test

# combine varfiles in data set with surface from oter files

varfiles = dict.fromkeys(['out2d',])
for key in varfiles:
    varfiles[key] = list(np.sort(glob(indir + key + '_*.nc')))  # last three available file
varfiles_variables = {}


#####  Variable Definiton sections for name mapping and standard name and units (needs to be adapted for each new vairable s added)
# nomenclature: name replacements and unit additions
#varnames = ['elevation', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY',
#            'sigWaveHeight', 'meanWavePeriod', 'zeroDowncrossPeriod', 'peakPeriod', 'meanWaveDirection',
#            'dominantDirection', 'meanDirSpreading', 'peakSpreading']


varnames = ['elevation', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY',
            'sigWaveHeight','bottomStressX','bottomStressY']



names = {'depth':'depth',
         'elevation': 'ssh', 
         'salinity': 'sss', 
         'temperature': 'sst',
         'horizontalVelX': 'u', 'horizontalVelY': 'v', 'velocity_magnitude': 'velocity_magnitude',
         'sigWaveHeight': 'Hs',  # Significant wave height
         'meanWavePeriod': 'TM01',  # Mean wave period (1st moment)
         'zeroDowncrossPeriod': 'TM02',  # Zero-downcrossing wave period (2nd moment)
         'peakPeriod': 'Tp',  # Peak wave period
         'meanWaveDirection': 'Dm',  # Mean wave direction
         'dominantDirection': 'Dp',  # Dominant wave direction (often same as peak wave direction)
         'meanDirSpreading': 'sigma_dir',  # Mean directional spreading
         'peakSpreading': 'sigma_p',  # Peak directional spreading'bottomStressX': 'tau_bx',
         'bottomStressY': 'tau_bx',
         'bottomStressY': 'tau_by'
         }

units = {'depth':'m',
         'elevation': 'm', 
         'salinity': '1e-3',   # Salininity not psu but 1e-3
         'temperature': 'degC', 
         'horizontalVelX': 'm/s',
         'horizontalVelY': 'm/s',
         'velocity_magnitude': 'm/s',
         'sigWaveHeight': 'm',  # Significant wave height (meters)
         'meanWavePeriod': 's',  # Mean wave period (seconds)
         'zeroDowncrossPeriod': 's',  # Zero-downcrossing wave period (seconds)
         'peakPeriod': 's',  # Peak wave period (seconds)
         'meanWaveDirection': 'degree',  # Mean wave direction (degrees)
         'dominantDirection': 'degree',  # Dominant wave direction (degrees)
         'meanDirSpreading': 'degree',  # Mean directional spreading (degrees)
         'peakSpreading': 'degree',  # Peak directional spreading (degrees)
         'bottomStressX': 'Pa',
         'bottomStressY': 'Pa'
         }

# For some variables there might be no propper standard name (here e.g. 'meanWaveDirection': 'sea_surface_wave_mean_direction'
#so put None insead and inore writing the varialbe)
std_names = {'depth':'sea_floor_depth_below_sea_surface',
             'elevation':
             'sea_surface_elevation',
             'salinity': 'sea_water_salinity',
             'temperature': 'sea_water_temperature', 
             'horizontalVelX': 'surface_eastward_sea_water_velocity',
             'horizontalVelY': 'surface_northward_sea_water_velocity',
             'velocity_magnitude': 'surface_sea_water_velocity_magnitude',
             'sigWaveHeight': 'sea_surface_wave_significant_height',
             'meanWavePeriod': 'sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment',
             'zeroDowncrossPeriod': 'sea_surface_wave_zero_upcrossing_period', #!!!!!!!!!!!!!
             'peakPeriod': 'sea_surface_wave_period_at_variance_spectral_density_maximum',
             'meanWaveDirection': 'sea_surface_wave_to_direction', #!!!!!!!!!!!!!
             'dominantDirection': 'sea_surface_wave_from_direction_at_variance_spectral_density_maximum',
             'meanDirSpreading': 'sea_surface_wave_directional_spread', #fixedtex
             'peakSpreading': 'sea_surface_wave_directional_spread_at_variance_spectral_density_maximum',
             'bottomStressX': None,
             'bottomStressY': None
            }

long_names = {
    'depth':"sea_floor_depth_below_sea_surface",
    'elevation': 'Sea Surface Height (SSH)', 
    'salinity': 'Sea Water Salinity (SSS)', 
    'temperature': 'Sea Surface Temperature (SST)',
    'horizontalVelX': 'Eastward Surface Sea Water Velocity (U)', 
    'horizontalVelY': 'Northward Surface Sea Water Velocity (V)', 
    'velocity_magnitude': 'Surface Sea Water Velocity Magnitude',
    'sigWaveHeight': 'Significant Wave Height (Hs)',  # Significant wave height
    'meanWavePeriod': 'Mean Wave Period (TM01) (1st Moment)',  # Mean wave period (1st moment)
    'zeroDowncrossPeriod': 'Zero-Downcrossing Wave Period (TM02)',  # Zero-downcrossing wave period (2nd moment)
    'peakPeriod': 'Peak Wave Period (Tp)',  # Peak wave period
    'meanWaveDirection': 'Mean Wave Direction (Average Direction of Wave Propagation)',  # Mean wave direction
    'dominantDirection': 'Dominant Wave Direction (Peak Wave Direction)',  # Dominant wave direction (often same as peak wave direction)
    'meanDirSpreading': 'Mean Directional Spreading of Waves (sigma_dir)',  # Mean directional spreading
    'peakSpreading': 'Peak Directional Spreading of Waves (sigma_p)'  # Peak directional spreading
    'bottomStressX': 'Eastward Bottom Stress Component',
    'bottomStressY': 'Northward Bottom Stress Component'
    }

# valid ranges for all of the variables CF 1.8 requirements
# must be floats
valid_ranges = {
    'depth': [-10.0,8000.0], 
    'elevation': [-6.0,12.0], 
    'salinity': [0.0,42.0], 
    'temperature': [-2.0,28.0],
    'horizontalVelX': [-8.0,8.0], 
    'horizontalVelY': [-8.0,8.0], 
    'velocity_magnitude': [0.0,8.0],
    'sigWaveHeight': [0.0,12.0],  # Significant wave height
    'meanWavePeriod': [0.0,40.0],  # Mean wave period (1st moment)
    'zeroDowncrossPeriod':[0.0,40.0],  # Zero-downcrossing wave period (2nd moment) [0.0,26.0]
    'peakPeriod': [0.0,40.0],  # Peak wave period
    'meanWaveDirection': [0.0,360.0],  # Mean wave direction
    'dominantDirection': [0.0,360.0],  # Dominant wave direction (often same as peak wave direction)
    'meanDirSpreading': [0.0,120.0],  # Mean directional spreading
    'peakSpreading': [0.0,120.0],  # Peak directional spreading
    'bottomStressX': [-5.0, 5.0],
    'bottomStressY': [-5.0, 5.0]
    }   

###################################################################


# grid and interpolation info
longitude = x  # the variable name is used as actually name by xarray netcdf export
latitude = y
X, Y = np.meshgrid(x, y)
if reload_weights:
    parents = np.loadtxt('parents_4cosyna_{:s}.txt'.format(str(dx).replace('.', ''))).astype(int)
    ndeweights = np.loadtxt('nodeweights_4cosyna_{:s}.txt'.format(str(dx).replace('.', '')))
elif False:  # calculate weights for interpolation
    xq = X.flatten()
    yq = Y.flatten()
    parents, ndeweights = s.find_parent_tri(xq, yq, dThresh=dThresh, latlon=True)
    # np.savetxt('parents_4cosyna.txt',parents)
    # np.savetxt('nodeweights_4cosyna.txt',ndeweights)
    np.savetxt('parents_4cosyna_0005.txt', parents)
    np.savetxt('nodeweights_4cosyna_0005.txt', ndeweights)
else: #neiw hybrid mode for speed
    # hbridy finding approach
    xq = X.flatten()
    yq = Y.flatten()
    s.init_element_tree(latlon=True)
    parents, ndeweights = s.find_parent_tri_fast(xq, yq, dThresh=dThresh, latlon=True, k=10)
    # --- find missing points ---
    #missing = parents == -1
    #missing.sum()/len(missing)*100
    dist,nn=s.element_tree_latlon.query(list(zip(xq,yq)),k=1)
    missing = (parents == -1) & (dist<dThresh)
    missing.sum()/len(missing)*100

    if np.any(missing):
        print(f"Refining {np.sum(missing)} points with slow fallback...")
        parents_slow, weights_slow = s.find_parent_tri(
            xq[missing], yq[missing], dThresh=dThresh, latlon=True
        )

        # fill the missing entries
        parents[missing] = parents_slow
        ndeweights[missing] = weights_slow
    np.savetxt('parents_4cosyna_0005.txt', parents)
    np.savetxt('nodeweights_4cosyna_0005.txt', ndeweights)


iuse = (parents != -1)
valid_parents = parents[iuse]


############# Beginn actual interpolation #########################
###### write depth ####


# open refercne files go throgh by dates


outfiles = []
# for infile in files:
for idate in range(len(varfiles['out2d'])):
    infiles = [varfiles[key][idate] for key in varfiles.keys()]

    # different time in putputs
    # Open datasets
    # time miss match
    dsout2d = xr.open_dataset(infiles[0])
    ds_list = [xr.open_dataset(infiles[i]).sel(nSCHISM_vgrid_layers=-1) for i in range(1, len(infiles))]

    # Align all datasets to avoid NaNs
    dsout2d, *ds_list = xr.align(dsout2d, *ds_list, join="outer")

    # Merge datasets
    dsin = xr.merge([dsout2d] + ds_list)
    print(dsin)

    # output file
    # fname=infile.split('/')[-1].replace('schout','schism_interp')
    fname = infiles[0].split('/')[-1].replace('out2d', 'schism-wwm_interp')

    outfile = outdir + fname
    outfiles.append(outfile)
    print(outfile)
    if (not os.path.exists(outfile)) or overwrite:
        # ?? time missmatch???

        ##each day  merge out2d outputs with surface fields from other variable files
        dsout2d = xr.open_dataset(infiles[0])
        #ds2 = xr.open_dataset(infiles[1]).sel(nSCHISM_vgrid_layers=-1)  # select surface
        ds_list = [xr.open_dataset(infiles[i]).sel(nSCHISM_vgrid_layers=-1) for i in
                   range(1, len(infiles))]  # select surface

        dsin = dsmerge = xr.merge([dsout2d, ] + ds_list)  # merge date setes
        if "horizontalVelX" in dsin.keys():
            dsin['velocity_magnitude'] = np.sqrt(dsin['horizontalVelX'] ** 2 + dsin['horizontalVelY'] ** 2)

        # reformat test
        time = dsin.time.values
        attrs = dsin.time.attrs
        t0 = time[0]
        try:
            seconds = (time - t0) / np.timedelta64(1, 's')
        except:
            seconds = (time - 0)
        t0 = str(t0)[:19].replace('T', ' ')
        tunit = "seconds since {:s}".format(t0)
        attrs['units'] = tunit
        attrs['_FillValue'] = False
        for varname in varnames:
            print(varname)
            data = np.ones((24,) + X.shape) * FillValue
            for i in range(24):
                wet = dsin.dryFlagElement[i, :].values == 0
                iuse = (parents != -1)
                valid_parents = parents[iuse]

                ivalid = wet[valid_parents]
                valid_parents = valid_parents[ivalid]
                weights = ndeweights[iuse, :][ivalid, :]
                # in domain wet elements
                temperature = dsin[varname][i, :].values
                target_inds = np.where(iuse)[0][ivalid]
                ii, jj = np.unravel_index(target_inds, X.shape)
                data[i, ii, jj] = (temperature[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
                valid_data=data[data!=FillValue]

            if varname == 'elevation':
                dslat = xr.DataArray(name="latitude", data=latitude, dims=["latitude"], coords=dict(
                    latitude=(["latitude"], latitude),
                ),
                                     attrs=dict(
                                         description="latitude",
                                         units="degrees_north",
                                         standard_name="latitude",
                                         _FillValue=False,
                                     ))

                dslon = xr.DataArray(name="longitude", data=longitude, dims=["longitude"], coords=dict(
                    longitude=(["longitude"], longitude),
                ),
                                     attrs=dict(
                                         description="longitude",
                                         units="degrees_east",
                                         standard_name="longitude",
                                         _FillValue=False,
                                     ))


                dstime = xr.DataArray(name="time", data=seconds, dims=["time"], coords=dict(
                    time=(["time"], seconds),
                ), attrs=attrs,
                                      )
                da = xr.DataArray(name=names[varname], data=data, dims=["time", "latitude", "longitude", ], coords=dict(
                    latitude=dslat,  # (["latitude"],dslat),
                    longitude=dslon,  # (["longitude"],dslon),
                    time=dstime,
                ),
                                  attrs=dict(
                                      description="ssh",
                                      units=units[varname],
                                      long_name=long_names[varname],
                                      _FillValue=FillValue,
                                      actual_range = [np.min(valid_data), np.max(valid_data)],  # 
                                      valid_range = valid_ranges[varname],
                                      **({ 'standard_name': std_names[varname] } if std_names[varname] is not None else {})
                                  ))


                #da.to_netcdf(outfile, mode='w')
                ds = xr.Dataset({da.name: da})
                ds.attrs.update(global_atts)  # Set global attributes on the dataset
                ds.to_netcdf(outfile, mode='w')
                
            else:
                da = xr.DataArray(name=names[varname], data=data, dims=["time", "latitude", "longitude", ],
                                  attrs=dict(
                                      description=varname,
                                      units=units[varname],
                                      long_name=long_names[varname],
                                      _FillValue=FillValue,
                                      actual_range = [np.min(valid_data), np.max(valid_data)],  # 
                                      valid_range = valid_ranges[varname],
                                      **({ 'standard_name': std_names[varname] } if std_names[varname] is not None else {})                                      
                                  ))
                da.to_netcdf(outfile, mode='a')
                
        #add depth
        Ddata = np.ones(X.shape) * FillValue
        iuse = (parents != -1)
        valid_parents = parents[iuse]
        weights = ndeweights[iuse, :]
        target_inds = np.where(iuse)[0]
        ii, jj = np.unravel_index(target_inds, X.shape)
        Ddata[ii, jj] = (np.asarray(s.depths)[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
        valid_data=Ddata[ii, jj].flatten()
        varname='depth'                            
        da_depth = xr.DataArray(name=names[varname], data=Ddata, dims=["latitude", "longitude", ], coords=dict(
            latitude=dslat,  # (["latitude"],dslat),
            longitude=dslon,  # (["longitude"],dslon),
        ),
                          attrs=dict(
                              description="bathymetry",
                              units=units[varname],
                              long_name=long_names[varname],
                              _FillValue=FillValue,
                              actual_range = [np.min(valid_data), np.max(valid_data)],  # 
                              valid_range = valid_ranges[varname],
                              **({ 'standard_name': std_names[varname] } if std_names[varname] is not None else {})
                          ))
        
        da_depth.to_netcdf(outfile, mode='a')





#da.to_netcdf(outfile, mode='w')
ds = xr.Dataset({da.name: da})
ds.attrs.update(global_atts)  # Set global attributes on the dataset
ds.to_netcdf(outfile, mode='w')



wet = dsin.dryFlagElement[i, :].values == 0


ivalid = wet[valid_parents]
valid_parents = valid_parents[ivalid]
weights = ndeweights[iuse, :][ivalid, :]
# in domain wet elements
temperature = dsin[varname][i, :].values
target_inds = np.where(iuse)[0][ivalid]
ii, jj = np.unravel_index(target_inds, X.shape)
data[i, ii, jj] = (temperature[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
valid_data=data[data!=FillValue]
