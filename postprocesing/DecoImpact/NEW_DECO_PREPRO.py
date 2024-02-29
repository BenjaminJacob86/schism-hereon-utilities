""" PreProc converter for D-ECO impact 
Load a stack of SCHISM output
Modifies it to ugrid standard using functions of Deltares
saves output in ugrid
"""

import glob
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from netCDF4 import default_fillvals
import time
import shapely
#import shapely.geometry
import xugrid as xu # pip install xugrid, adds ugrid and grid accessors to xarray datasets


import numpy as np
import netCDF4 as nc
from typing import Union, Optional
import datetime as dt
import pandas as pd
import dask
import math
import time 
plt.ion()

# try parallel execution at different cpus


# DecoImpact prproc
	
def preprocess_schism_scribedio(ds):
    """
    This preprocessing function describes the minimal changes that would have to be made 
    in the SCHISM output in order for it to be read via ugrid conventions with xugrid.
    It is probably not a complete overview.
    """
    # set varnames
    gridname = "SCHISM_hgrid"
    fnc_varn = f"{gridname}_face_nodes"
    enc_varn = f"{gridname}_edge_nodes"
    node_x = f"{gridname}_node_x"
    node_y = f"{gridname}_node_y"
    
    # set topology attributes to empty variable
    topo_attrs = {"cf_role": "mesh_topology",
                 "topology_dimension": 2, # has to be int, not str
                 "node_coordinates": f"{node_x} {node_y}",
                 "face_node_connectivity": fnc_varn,
                 "edge_node_connectivity": enc_varn,
                 }
    ds[gridname] = xr.DataArray(np.array(default_fillvals['i4'], dtype=np.int32), attrs=topo_attrs)
    
    # assign necessary attributes to connectivity variables
    fnc_attrs = {"_FillValue":-1, "start_index":1}
    ds[fnc_varn] = ds[fnc_varn].assign_attrs(fnc_attrs)
    ds[enc_varn] = ds[enc_varn].assign_attrs(fnc_attrs)
    
    # set node_x/node_y as coordinate variables instead of data_vars
    ds = ds.set_coords([node_x,node_y])
    
    # to prevent xugrid UserWarning, but this is hardcoded and it should be different for lat/lon models.
    # "UserWarning: No standard_name of ('projection_x_coordinate', 'longitude', 'projection_y_coordinate', 'latitude') in
    # ['SCHISM_hgrid_node_x', 'SCHISM_hgrid_node_y']. Using SCHISM_hgrid_node_x and SCHISM_hgrid_node_y as projected x and y coordinates."
    projected = True
    if projected:
        ds[node_x] = ds[node_x].assign_attrs({"standard_name":"projection_x_coordinate"})
        ds[node_y] = ds[node_y].assign_attrs({"standard_name":"projection_y_coordinate"})
    else:
        ds[node_x] = ds[node_x].assign_attrs({"standard_name":"longitude"})
        ds[node_y] = ds[node_y].assign_attrs({"standard_name":"latitude"})
    
    # add variable with coordinate system, optional but convenient for loading into QGIS and other tools
    if projected:
        grid_mapping_name = 'Unknown projected'
        crs_varn = 'projected_coordinate_system'
        crs_num = 0 # TODO: this is clearly incorrect
    else:
        grid_mapping_name = 'latitude_longitude'
        crs_varn = 'wgs84'
        crs_num = 4326
    crs_str = f'EPSG:{crs_num}'
    crs_attrs = {'epsg': crs_num, # epsg or EPSG_code are required for correct interpretation by QGIS
                  'EPSG_code': crs_str, # epsg or EPSG_code are required for correct interpretation by QGIS
                  'grid_mapping_name': grid_mapping_name,
                  }
    ds[crs_varn] = xr.DataArray(np.array(default_fillvals['i4'],dtype=int),dims=(),attrs=crs_attrs)
    
    # mesh attribute is required for d-ecoimpact #TODO: why?
    # valueable other attrs are "location" (node/face/edge), 
    # "standard_name", "long_name", "units", "grid_mapping"
    for varn in ds.data_vars:
        ds[varn] = ds[varn].assign_attrs({'mesh': gridname})
    
    # time requires "units" attribute to be converted by xarray and other tools
    # refdate taken from params.nml
    ds['time'] = ds['time'].assign_attrs({"units":"seconds since 2017-01-02 00:00:00"})
    ds = xr.decode_cf(ds)
    
    #TODO: set _FillValue attribute for data_vars, test dataset did not seem to have nans
    
    return ds	
	
def uda_nodes_to_faces(uda_node : xu.UgridDataArray) -> xu.UgridDataArray:
    """
    Credits : J. Veenstra (Deltares)
    https://github.com/Deltares/dfm_tools/blob/a5dd98f73e20d684f0a25a05d4e445a1b00a14b4/dfm_tools/xugrid_helpers.py#L466
    Interpolates a ugrid variable (xu.DataArray) with an node dimension to the faces by averaging the 3/4 nodes around each face.
    Since node variables are mostly defined on interfaces, it also interpolates from interfaces to layers

    Parameters
    ----------
    uda_node : xu.UgridDataArray
        DESCRIPTION.

    Raises
    ------
    KeyError
        DESCRIPTION.

    Returns
    -------
    uda_face : xu.UgridDataArray
        DESCRIPTION.

    """
    
    dimn_faces = uda_node.grid.face_dimension
    dimn_maxfn = 'nMax_face_nodes' #arbitrary dimname that is reduced anyway
    dimn_nodes = uda_node.grid.node_dimension
    fill_value = uda_node.grid.fill_value
    
    if dimn_nodes not in uda_node.sizes:
        raise KeyError(f'varname "{uda_node.name}" does not have an node dimension ({dimn_nodes})')
    
    # construct indexing array
    data_fec = xr.DataArray(uda_node.grid.face_node_connectivity,dims=(dimn_faces,dimn_maxfn))
    data_fec_validbool = data_fec!=fill_value
    data_fec = data_fec.where(data_fec_validbool,-1)
    
    print('node-to-face interpolation: ',end='')
    dtstart = dt.datetime.now()
    # for each face, select all corresponding edge values (this takes some time)
    uda_face_alledges = uda_node.isel({dimn_nodes:data_fec})
    # replace nonexistent edges with nan
    uda_face_alledges = uda_face_alledges.where(data_fec_validbool) #replace all values for fillvalue edges (-1) with nan
    # average edge values per face
    uda_face = uda_face_alledges.mean(dim=dimn_maxfn,keep_attrs=True)
    #update attrs from edge to face
    face_attrs = {'location': 'face', 'cell_methods': f'{dimn_faces}: mean'}
    uda_face = uda_face.assign_attrs(face_attrs)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return uda_face	
    
    
## load the data from OG SCHISM
def reduce_file_list(file_nc_list):
	keep_files=[]
	keep=['horizontalVelX',  'horizontalVelY',  'out2d',  'salinity',  'temperature',  'verticalVelocity', 'zCoordinates']
	keep_files=[]
	for file in file_nc_list:
		for pattern in keep:
			if pattern in file:
				keep_files.append(file)
				break
	return keep_files


#stack=1
stack0=1
stack1=5

for stack in range(stack0,stack1+1):

    print('Preprocessing Stack '+ str(stack))
    # open dataset 
    file_nc_list = glob.glob('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all/*_?.nc') + glob.glob('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all/*_1?.nc') +glob.glob('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all/*_2?.nc') +glob.glob('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all/*_30.nc')		
    a = glob.glob('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all/*_'+str(stack)+'.nc')
    a=reduce_file_list(a)

    dsa = xr.open_mfdataset(a) 
    dsa = dsa.sel(nSCHISM_vgrid_layers=-1)

    # make it a ugrid
    dsa = preprocess_schism_scribedio(dsa) # verggessen an welcher setelle
    uds=xu.UgridDataset(dsa)
    uds.ugrid.set_crs("EPSG:25832") 

    uds['SCHISM_hgrid']=dsa['SCHISM_hgrid']  # add attributes again

    #uds.ugrid.to_netcdf('test_x.nc') # writes to ugrid

    # then convert coordinates, also converts standard names from projection_x_coordinate to longitude

    ### ADd variables for Deco
    mesh2d=uds

    mesh_name="SCHISM_hgrid"   # "schism_mesh"

    mesh2d["waterdepth"] = mesh2d["elevation"] - mesh2d["depth"]
    mesh2d["waterdepth"] = mesh2d["waterdepth"].assign_attrs({
      "mesh" : mesh_name
     })

    # Creates Error
    #calculate horizontalvelocity (based on pythagoram of x and y)
    V = np.sqrt((mesh2d["horizontalVelX"])**2 + (mesh2d["horizontalVelY"])**2)
    V.attrs= {'mesh': 'SCHISM_hgrid'}
    mesh2d["horizontalVel"] = V


    # Map variables to to face
    mesh2d["elevation_face"] = uda_nodes_to_faces(mesh2d["elevation"]) # only this works
    mesh2d["waterdepth_face"] = uda_nodes_to_faces(mesh2d["waterdepth"])
    mesh2d["depth_face"] = uda_nodes_to_faces(mesh2d["depth"])
    mesh2d["horizontalVel_face"] = uda_nodes_to_faces(mesh2d["horizontalVel"]) # takes longer
    mesh2d["salinity_face"] = uda_nodes_to_faces(mesh2d["salinity"])


    file_nc_out='preproc_schout_{:02d}.nc'.format(stack)
    #mesh2d.ugrid.to_netcdf(file_nc_out)


    # maybe needs fxed order of packates



    ####### Experimental part ########

    ### test subsetting only element variables and needed coords

    #mesh2d_red=mesh2d[["elevation_face","waterdepth_face","depth_face",
    #"horizontalVel_face","salinity_face","projected_coordinate_system",
    #"SCHISM_hgrid_node_x","SCHISM_hgrid_node_y",'nSCHISM_hgrid_face',
    #'nSCHISM_hgrid_node',
    #'nSCHISM_hgrid_edge']]
    #mesh2d_red[mesh_name]=mesh2d[mesh_name]  # needs to be seperate assignment

    #mesh2d_red['nSCHISM_hgrid_edge']
    #mesh2d_red.to_netcdf('reduced_test.nc')

    # still makes probles
    temp=mesh2d
    variables=list(temp.variables)
    for var  in variables:
        if not (('SCHISM' in var) or ('grid' in var) or ('face' in var) or ('projected_coordinate_system' in var) or ('dryFlagElement' in var) or ('time' in var)  ):
            print('dropping variable: ' + var)
            temp=temp.drop(var)
    temp['SCHISM_hgrid']=mesh2d['SCHISM_hgrid'] # add again
    temp['time']= mesh2d['time']

    temp.ugrid.to_netcdf(file_nc_out)


#.to_netcdf('reduced_test.nc')
#mesh2d_red[mesh_name]=mesh2d[mesh_name]

# ncrcat files merge files 

files2=np.sort(glob.glob('preproc_schout_*.nc'))
ds1=xr.open_mfdataset(files2)
ds1.to_netcdf('merged_out.nc')