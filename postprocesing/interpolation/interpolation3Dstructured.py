"""
Code for structured grid interpolation of SCHISM output using `schism.py`.

This script loops over a set of predefined variables, each mapped to standardized 
names, units, and CF-compliant metadata through dictionaries (`names`, `std_names`, `units`, etc.).

For 3D variables:
    - Vertical interpolation is performed first using linear interpolation 
      based on SCHISM's hybrid sigma/z-layer coordinates.
    - This is followed by horizontal interpolation using inverse distance weighting 
      over the unstructured triangular grid to map onto a structured target grid.

For 2D variables:
    - Only horizontal interpolation is applied.

To include new variables, update the dictionaries with appropriate CF-compliant 
entries (standard name, long name, units, valid range, etc.).

"""

import sys
import os

sys.path.insert(0, '/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine/')
from schism import *
from matplotlib import pyplot as plt
from glob import glob
import datetime as dt


def get_layer_weights(s, dep,ti):
    print('calculating weights for vertical interpolation')
    ibelow = np.zeros(s.nnodes, int)
    iabove = np.zeros(s.nnodes, int)
    weights = np.zeros((2, s.nnodes))

    # garbage values below ibtm different type , nan or strange values wrong values
    zcor = s.ds['zCoordinates'][ti, :,:].values  
    zcor = np.ma.masked_array(zcor, mask=s.mask3d)

    a = np.sum(zcor <= dep, 1) - 1
    # ibelow=a+np.sum(zcor.mask,1)-1
    ibelow = a + s.ibbtm - 1
    # ibelow=a+(self.nc['bottom_index_node'][:]-1)-1
    iabove = np.minimum(ibelow + 1, s.nz - 1)
    inodes = np.where(a > 0)[0]
    ibelow2 = ibelow[inodes]
    iabove2 = iabove[inodes]

    d2 = zcor[inodes, iabove2] - dep
    d1 = dep - zcor[inodes, ibelow2]
    ivalid = d1 > 0.0
    iset = d1 == 0.0
    d1 = d1[ivalid]
    d2 = d2[ivalid]
    weights[0, inodes[ivalid]] = 1 - d1 / (d1 + d2)  # 1/d1/(1/d1+1/d2)
    weights[1, inodes[ivalid]] = 1 - weights[0, inodes[ivalid]]  # 1/d2/(1/d1+1/d2)
    weights[0, inodes[iset]] = 1
    weights[:, np.sum(weights, 0) == 0.0] = np.nan

    return ibelow, iabove, weights
########## Settings ############################################

indir = '/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine/outputs/'
outdir = '/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine/outputs4HCDC/'
# load setup
os.chdir('/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine/')
s = schism_setup()
os.chdir(outdir)
# last three available files
# files=list(np.sort(glob(indir+'schout_*.nc'))[-16:]) # last three available files


# grid and interpolation info
overwrite = False  # overwirte interpolation files if existent
lon, lat = np.asarray(s.lon), np.asarray(s.lat)
dx = 0.005  # 0.01  # 0.025  interpolation grid spacing in degree
lonmin = 5.11608851  # lon.min()
lonmax = 10.40483834  # lon.max()
latmin = 53.03857668  # lat.min()
latmax = 55.62872055  # lat.max()
# x=np.arange(lon.min(),lon.max(),dx)
# y=np.arange(lat.min(),lat.max(),dx)
x = np.arange(lonmin, lonmax, dx)
y = np.arange(latmin, latmax, dx)
reload_weights = True  # relaod weights
FillValue = -9999


nt=24 #time steps per file

# vertical levels in sigma or z (depth, positive down or up depending on your convention)
nz=21  # number inerpolated vertical layers
z_levels = np.linspace(0, 50, nz)  # Example: 20 layers from surface to -30 m

# Target grid
Z, Y3D, X3D = np.meshgrid(z_levels, y, x, indexing='ij')  # Shape (nz, ny, nx)


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
global_atts['forecast_type']  = 'forecast'
global_atts['forecast_range']  = '2-day_forecast'
global_atts['licence']  = 'https://spdx.org/licenses/CC-BY-4.0.html'
global_atts['domain_name']  = 'German Bight'
global_atts['history']  = f"; SCHISM output processed via interpolation4HCDC.py on {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"


# make sure out2d is always first in this list
varfiles = dict.fromkeys(['out2d', 'temperature', 'salinity', 'horizontalVelX', 'horizontalVelY','zCoordinates'])
for key in varfiles:
    varfiles[key] = list(np.sort(glob(indir + key + '_*.nc'))[-4:])  # last three available files

# todaystr=dt.datetime.today().strftime('%Y%m%d')
## temp test

# combine varfiles in data set with surface from oter files


varfiles_variables = {}


#####  Variable Definiton sections for name mapping and standard name and units (needs to be adapted for each new vairable s added)
# nomenclature: name replacements and unit additions
varnames = ['elevation', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY',
            'sigWaveHeight', 'meanWavePeriod', 'zeroDowncrossPeriod', 'peakPeriod', 'meanWaveDirection',
            'dominantDirection', 'meanDirSpreading', 'peakSpreading']

#names = {'elevation': 'ssh', 
#         'salinity': 'sss', 
#         'temperature': 'sst',
#         'horizontalVelX': 'u', 'horizontalVelY': 'v', 'velocity_magnitude': 'velocity_magnitude',
#         'sigWaveHeight': 'Hs',  # Significant wave height
#         'meanWavePeriod': 'TM01',  # Mean wave period (1st moment)
#         'zeroDowncrossPeriod': 'TM02',  # Zero-downcrossing wave period (2nd moment)
#         'peakPeriod': 'Tp',  # Peak wave period
#         'meanWaveDirection': 'Dm',  # Mean wave direction
#         'dominantDirection': 'Dp',  # Dominant wave direction (often same as peak wave direction)
#         'meanDirSpreading': 'sigma_dir',  # Mean directional spreading
#         'peakSpreading': 'sigma_p'  # Peak directional spreading
#         }

names = {
    'elevation': 'zos',
    'salinity': 'so',
    'temperature': 'thetao',
    'horizontalVelX': 'uo',
    'horizontalVelY': 'vo',
    'velocity_magnitude': 'speed',  # Not a CMEMS var, but descriptive
    'sigWaveHeight': 'swh',
    'meanWavePeriod': 'mwp',
    'zeroDowncrossPeriod': 'tz',
    'peakPeriod': 'pp',
    'meanWaveDirection': 'mwd',
    'dominantDirection': 'dp',
    'meanDirSpreading': 'sprd_mean',
    'peakSpreading': 'sprd_peak'
}




units = {'elevation': 'm', 
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
         'peakSpreading': 'degree'  # Peak directional spreading (degrees)
         }

# For some variables there might be no propper standard name (here e.g. 'meanWaveDirection': 'sea_surface_wave_mean_direction'
#so put None insead and inore writing the varialbe)
std_names = {'elevation':
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
             'peakSpreading': 'sea_surface_wave_directional_spread_at_variance_spectral_density_maximum'}

long_names = {
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
}


# valid ranges for all of the variables CF 1.8 requirements
# must be floats
valid_ranges = {
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
    'peakSpreading': [0.0,120.0]  # Peak directional spreading
}


###################################################################


# grid and interpolation info
longitude = x  # the variable name is used as actually name by xarray netcdf export
latitude = y
X, Y = np.meshgrid(x, y)
ny, nx =  X.shape
if reload_weights:
    parents = np.loadtxt('parents_4cosyna_{:s}b.txt'.format(str(dx).replace('.', ''))).astype(int)
    ndeweights = np.loadtxt('nodeweights_4cosyna_{:s}b.txt'.format(str(dx).replace('.', '')))
else:  # calculate weights for interpolation
    xq = X.flatten()
    yq = Y.flatten()
    parents, ndeweights = s.find_parent_tri(xq, yq, dThresh=0.02, latlon=True)
    # np.savetxt('parents_4cosyna.txt',parents)
    # np.savetxt('nodeweights_4cosyna.txt',ndeweights)
    np.savetxt('parents_4cosyna_0005.txt', parents)
    np.savetxt('nodeweights_4cosyna_0005.txt', ndeweights)


############# Beginn actual interpolation #########################


# open refercne files go throgh by dates

outfiles = []
# for infile in files:
for idate in range(len(varfiles['out2d'])):
   
    infiles = [varfiles[key][idate] for key in varfiles.keys()]

    # different time in putputs
    # Open datasets
    # time miss match
    dsout2d = xr.open_dataset(infiles[0])
    #ds_list = [xr.open_dataset(infiles[i]).sel(nSCHISM_vgrid_layers=-1) for i in range(1, len(infiles))]
    ds_list = [xr.open_dataset(infiles[i]) for i in range(1, len(infiles))]

    # Align all datasets to avoid NaNs
    dsout2d, *ds_list = xr.align(dsout2d, *ds_list, join="outer")

    # Merge datasets
    dsin = xr.merge([dsout2d] + ds_list)
    print(dsin)

    # output file
    # fname=infile.split('/')[-1].replace('schout','schism_interp')
    fname = infiles[0].split('/')[-1].replace('out2d', 'schism-wwm_3Dinterp')

    outfile = outdir + fname
    outfiles.append(outfile)
    print(outfile)
    if (not os.path.exists(outfile)) or overwrite:
        # ?? time missmatch???

        ##each day  merge out2d outputs with surface fields from other variable files
        dsout2d = xr.open_dataset(infiles[0])
        ds2 = xr.open_dataset(infiles[1]) #.sel(nSCHISM_vgrid_layers=-1)  # select surface
        ds_list = [xr.open_dataset(infiles[i]) for i in
                   range(1, len(infiles))]  # select surface

        dsin = dsmerge = xr.merge([dsout2d, ] + ds_list)  # merge date setes
        dsin['velocity_magnitude'] = np.sqrt(dsin['horizontalVelX'] ** 2 + dsin['horizontalVelY'] ** 2)
        s.ds=dsin
        s.nz=dsin['zCoordinates'].shape[-1]
        
        s.ibbtm = s.ds.bottom_index_node[:].values - 1
        s.mask3d = np.zeros((s.nnodes, s.nz), bool)  # mask for 3d field at one time step
        for inode in range(s.nnodes):
            s.mask3d[inode, :s.ibbtm[inode]] = True  # controlled that corresponding z is depths
        
        # reformat test
        time = dsin.time.values
        attrs = dsin.time.attrs
        t0 = time[0]
        seconds = (time - t0) / np.timedelta64(1, 's')
        t0 = str(t0)[:19].replace('T', ' ')
        tunit = "seconds since {:s}".format(t0)
        attrs['units'] = tunit
        attrs['_FillValue'] = False
        
        for varname in varnames:
            print(varname)
            
            if 'nSCHISM_vgrid_layers' in dsin[varname].dims:
                #3D
                data = np.ones((nt, nz, ny, nx)) * FillValue
                dim3D=True
            else:            
            #2D
                data = np.ones((nt,) + X.shape) * FillValue
                dim3D=False

            iuse = (parents != -1)
            target_inds_all = np.where(iuse)[0]
            ii_all, jj_all = np.unravel_index(target_inds_all, X.shape)
            

            for t in range(nt):
            
                wet = dsin.dryFlagElement[t, :].values == 0
                valid_parents = parents[iuse]
                wet_mask = wet[valid_parents]
                
                #updated logic to avodi recalc
                # Subset weights + indices for current time step
                valid_parents_t = valid_parents[wet_mask]
                weights_t = ndeweights[iuse, :][wet_mask]
                target_inds = target_inds_all[wet_mask]
                ii = ii_all[wet_mask]
                jj = jj_all[wet_mask]
                            
                
                iuse = (parents != -1)
                valid_parents = parents[iuse]

                ivalid = wet[valid_parents]
                valid_parents = valid_parents[ivalid]
                weights = ndeweights[iuse, :][ivalid, :]
                # in domain wet elements
                target_inds = np.where(iuse)[0][ivalid]
                ii, jj = np.unravel_index(target_inds, X.shape)
                
                if dim3D:
                    for k in range(nz):  # vertical layers
                    
                        ibelow, iabove, weights_vert = get_layer_weights(s, -z_levels[k],t) # get weights for above and below depths layer
                        # interpolate slize to search interpolation z level
                        
                        s.nodeinds=np.arange(s.nnodes)
                        field3D=dsin[varname][t, :].values
                        interp_slice=weights_vert[0, :] * field3D[s.nodeinds, ibelow] + weights_vert[1, :] * \
                                          field3D[:, :][s.nodeinds, iabove]
                                          
                        data[t,k, ii, jj] = (interp_slice[s.nvplt[valid_parents, :]] * weights).sum(axis=1)                  
                                            
                else:
                    interp_slice  = dsin[varname][t, :].values
                    data[t, ii, jj] = (interp_slice[s.nvplt[valid_parents, :]] * weights).sum(axis=1)
                
            valid_data=data.flatten()[data.flatten()!=FillValue]
            
            if varname == 'elevation':
                depth=z_levels # usae varaible name in dict
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

                # Create DataArrays for depth and time
                dsdepth = xr.DataArray(
                    name="depth",
                    data=depth,
                    dims=["depth"],
                    coords={"depth": (["depth"], depth)},
                    attrs={
                        "standard_name": "depth",
                        "long_name": "Depth below sea surface",
                        "units": "m",
                        "positive": "down",
                        "axis": "Z"
                    }
                )

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
                #ds = xr.Dataset({da.name: da})
                #ds.attrs.update(global_atts)  # Set global attributes on the dataset
                #ds.to_netcdf(outfile, mode='w')
                ##Create depth variable
                #xr.Dataset({"depth_var": dsdepth}).to_netcdf(outfile, mode='a')
                
                ds = xr.Dataset({da.name: da, "depth_below_surface": dsdepth})
                ds.attrs.update(global_atts)
                ds.to_netcdf(outfile, mode="w")
                

                
            elif dim3D:
                da = xr.DataArray(name=names[varname], data=data, dims=["time","depth", "latitude", "longitude", ],
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
                
            else: #2d data
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
    
