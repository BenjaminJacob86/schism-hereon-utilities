# ~/miniforge3/envs/geo_env/bin/ipython
import xarray as xr
from glob import glob
import yaml
import datetime as dt
import numpy as np
import os


yaml_file = 'metadata.yaml'
with open(yaml_file, 'r') as f:
    meta_defs = yaml.safe_load(f)

with open('units.yaml', 'r') as f:
    var_defs = yaml.safe_load(f)


global_atts = meta_defs.get('global_attributes', {}).copy()
global_atts['creation_date'] = f'Created on {dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
global_atts['history'] = f"; SCHISM output processed via interpolation_hybrid.py on {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
global_atts['creation_date'] = f'Created on {dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

global_atts["Conventions"] = "CF-1.10, UGRID-1.0, ACDD-1.3"


crs=xr.DataArray(
    0,
    attrs={
        "grid_mapping_name": "transverse_mercator",
        "epsg_code": "EPSG:32632",
        "semi_major_axis": 6378137.0,
        "inverse_flattening": 298.257223563,
        "longitude_of_central_meridian": 9.0,
        "latitude_of_projection_origin": 0.0,
        "scale_factor_at_central_meridian": 0.9996,
        "false_easting": 500000.0,
        "false_northing": 0.0,
    },
)    





ds=xr.open_mfdataset(glob('*.nc'))

ds.attrs.update(global_atts)
ds["crs"] =crs

names = var_defs['names']
units = var_defs['units']
std_names_2d = var_defs['std_names_2d']
std_names_3d = var_defs['std_names_3d']
long_names_2d = var_defs['long_names_2d']
long_names_3d = var_defs['long_names_3d']
valid_ranges = var_defs['valid_ranges']



global_atts = var_defs.get('global_attributes', {}).copy()
# Add dynamic fields that change with each run


#### overwrite units

FillValue = np.nan

remove=['Veg_max',]

for var in list(ds.data_vars):   # use data_vars, not variables
    print(var)
    
    
    for varname in names:
        if varname in var:

            addendum = [part for part in var.split(varname) if part != ''][0]
            
            for rmv in remove:
                addendum=addendum.replace(rmv,'')
            
            new_name = names[varname] + addendum
            
            #repairs of inconsistent format
            new_name = new_name.replace('mean__','_mean').replace('mean_','_mean')
            
            long_name = addendum.replace('_', '') + ' of ' + long_names_2d[varname]
            std_name = None
            valid_data = valid_ranges[varname]

            # Rename variable at dataset level
            ds = ds.rename({var: new_name})

            # Now update attributes
            ds[new_name].attrs = dict(
                description=varname,
                units=units[varname],
                long_name=long_name,
                #_FillValue=FillValue,
                #actual_range=[np.min(valid_data), np.max(valid_data)],
                actual_range = [float(ds[new_name].min().values),   float(ds[new_name].max().values),],
                valid_range=valid_ranges[varname],
                **({'standard_name': std_name} if std_name is not None else {}),
                grid_mapping='crs',
                mesh="mesh",
                location="node"
            )
            ds[new_name].encoding["_FillValue"] = FillValue
            break





# add spatial coordinates from schism raw output reference
dsi_4_coords=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs01/out2d_1.nc')
#plt.tripcolor(x,y,ds['SCHISM_hgrid_face_nodes'].values[:,:-1]-1,y,shading='flat')


coord_attrs = {
    "SCHISM_hgrid_node_x": dict(
        standard_name="projection_x_coordinate",
        long_name="UTM easting",
        units="m",
    ),
    "SCHISM_hgrid_node_y": dict(
        standard_name="projection_y_coordinate",
        long_name="UTM northing",
        units="m",
    ),
    "SCHISM_hgrid_face_nodes": dict(
        long_name="Triangular face node connectivity",
        cf_role="face_node_connectivity",
        start_index=1,  # SCHISM is 1-based
    ),
}

for key, attrs in coord_attrs.items():
    ds[key] = dsi_4_coords[key]

    if key == "SCHISM_hgrid_face_nodes":
        # remove unused quad column if only triangles
        #ds = ds.drop_dims("nMaxSCHISM_hgrid_face_nodes")
        #ds[key] = dsi_4_coords[key].isel(nMaxSCHISM_hgrid_face_nodes=slice(0,3))

        conn = dsi_4_coords["SCHISM_hgrid_face_nodes"][:, :3]
        conn = conn.rename({"nMaxSCHISM_hgrid_face_nodes": "nSCHISM_hgrid_face_nodes"})
        ds[key] = conn
        ds["SCHISM_hgrid_face_nodes"].attrs.update(
            dict(
                cf_role="face_node_connectivity",
                start_index=1,
                long_name="Triangular face connectivity"
            )
        )
    ds[key].attrs = attrs

# add mesh variable for ugidr
ds["mesh"] = xr.DataArray(
    0,
    attrs=dict(
        cf_role="mesh_topology",
        topology_dimension=2,
        node_coordinates="SCHISM_hgrid_node_x SCHISM_hgrid_node_y",
        face_node_connectivity="SCHISM_hgrid_face_nodes",
    ),
)


ds.to_netcdf(
    "annual_statistics_veg_MAX_2090.nc",
    format="NETCDF4",
    engine="netcdf4"
)


ds.to_netcdf(
    "annual_statistics_veg_CNTRL_2090.nc",
    format="NETCDF4",
    engine="netcdf4"
)

ds.to_netcdf(
    "annual_statistics_veg_MAX_1997.nc",
    format="NETCDF4",
    engine="netcdf4"
)

ds.to_netcdf(
    "annual_statistics_veg_CNTRL_1997.nc",
    format="NETCDF4",
    engine="netcdf4"
)


