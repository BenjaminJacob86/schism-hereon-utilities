""" Overwrite temperature and salinity fields """
import xarray as xr
#a=xr.open_dataset('hotstart.nc_sed_merged')

#/gpfs/work/ksddata/ROUTINES_personal/SCHISM/gb_wave_routine_seagrass/
a=xr.open_dataset('hotstart_dummy_with_sediments.nc')

b=xr.open_dataset('hotstart_20260801.nc') #take T und S from this hotstart to overwrite in a

a.tr_nd0[:,:,0]=b.tr_nd0[:,:,0]
a.tr_nd0[:,:,1]=b.tr_nd0[:,:,1]

a.tr_nd[:,:,0]=b.tr_nd[:,:,0]
a.tr_nd[:,:,1]=b.tr_nd[:,:,1]

a.tr_el[:,:,0]=b.tr_el[:,:,0]
a.tr_el[:,:,1]=b.tr_el[:,:,1]

a.to_netcdf('hotstart.nc_sed_merged_modified')



dsx=xr.open_dataset('hotstart.nc_sed_merged_modified')