"""
Valdiation of tide portal files for Gertman Bight,
assumption is that model grid is in UTM otherwise rework for 
lon lat necessrary
"""
import os
import sys
sys.path.insert(0, '/work/gg0028/SCHISM/schism-hereon-utilities/')
from schism import *
sys.path.insert(0, '/work/gg0028/SCHISM/schism-hereon-utilities/Validation/')
from kuestendaten import *
/work/gg0028/SCHISM/schism-hereon-utilities/Validation/kuestendaten.py

from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
import datetime
import numpy as np
plt.ion()

##settings
tp_dir='/work/gg0028/g260114/PROJECTS/FOCCUS/downloads'  #tide portal files
ncdir='/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000/outputs01/'
max_dist=200 #maximum dinstance betweens tation and grid node for consideration   

os.chdir(tp_dir)

station_files=glob('*.txt')
station_names=[station.split('!')[0] for station in station_files]
stations=dict.fromkeys(station_names)

for station in station_files:
    key=station.split('!')[0]  
    print("loading " + key)
    stations[key]=KuestenDataFile(station)


# Extract coordinates (skip missing ones)
def extreact_locations(stations):
    xs = []
    ys = []
    labels = []

    for key, obj in stations.items():
        try:
            xs.append(float(obj.East))
            ys.append(float(obj.North))
            labels.append(key)
        except:
            # Some stations don't have UTM coordinates
            pass
    return xs,ys,labels



# UTM Zone 32N projection (ETRS89)
proj_utm32 = ccrs.UTM(zone=32)
proj_geo = ccrs.PlateCarree()  # for coastlines

fig = plt.figure(figsize=(8, 10))
ax = plt.axes(projection=proj_utm32)

# Add coastlines & borders
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')

# Plot points
xs,ys,labels=extreact_locations(stations)
ax.scatter(xs, ys, transform=proj_utm32, s=30)

# Add labels
for x, y, label in zip(xs, ys, labels):
    ax.text(x, y, label, fontsize=8, transform=proj_utm32)

ax.set_title("Küstenmessstationen — ETRS89 / UTM 32N")
plt.show()


#Filter  with model extends
os.chdir('/work/gg0028/g260114/PROJECTS/FOCCUS/RUNS/Ehype000')
s=schism_setup()

s.plot_domain_boundaries(append=True,plot_legend=False,latlon=False)
plt.scatter(xs, ys,s=30)
# Add labels
for x, y, label in zip(xs, ys, labels):
    plt.text(x, y, label, fontsize=8)
    
    
## FIlter station distance to remove pegel From side arms    
s.init_node_tree(latlon=False)    
station_coords=list(zip(xs,ys))
plt.close()

dists,nns=s.node_tree_xy.query(station_coords)
remove_indices=dists>max_dist

keys = list(stations.keys())

stations_filtered = {
    k: stations[k]
    for k, remove in zip(keys, remove_indices)
    if not remove
}

stations = stations_filtered

xs,ys,labels=extreact_locations(stations)

# updates nns
station_coords=list(zip(xs,ys))
dists,nns=s.node_tree_xy.query(station_coords)

# Plot of Filer stations
s.plot_domain_boundaries(append=True,plot_legend=False,latlon=False)
plt.scatter(xs, ys,s=30)
# Add labels
for x, y, label in zip(xs, ys, labels):
    plt.text(x, y, label, fontsize=8)





# load schism
acc=schism_outputs_by_variable(ncdir=ncdir,varlist=['out2d','salinity'])

# depends on schism date/time output format
inidate= acc.ds['salinity'].time[0].values.astype('datetime64[ms]').astype(datetime.datetime)
findate= acc.ds['salinity'].time[-1].values.astype('datetime64[ms]').astype(datetime.datetime)


inidate=datetime.datetime(2017,1,2,0,0,0)
findate=datetime.datetime(2017,1,2,0,0,0) + acc.ds['salinity'].time[-1].values/86400 * datetime.timedelta(days=1)
dates=inidate+datetime.timedelta(seconds=1)*acc.ds['salinity'].time.values


# filter data date range
for key, st in stations.items():
    df = st.data
    # Filter to the model time range
    mask = (df['timestamp'] >= inidate) & (df['timestamp'] <= findate)
    st.data = df.loc[mask].copy()  # assign back filtered data

# interpolate
plt.figure()
for key, st in stations.items():
    df = st.data

    # interpoalte
    mask = (df['timestamp'] >= model_times.min()) & (df['timestamp'] <= model_times.max())
    df = df.loc[mask].copy()
    df.set_index('timestamp', inplace=True)

    df_interp = df.reindex(model_times)  # introduces NaNs at model times
    df_interp['value'] = df_interp['value'].interpolate(method='time')  # linear in time
    df_interp['status'] = df_interp['status'].ffill()
    #df_interp['remark'] = df_interp['remark'].ffill()
    
    st.data=df_interp

    print(df_interp.head())
    print(df_interp.tail())
    st.data.value.plot()


# subselect at station locations by nerest neighbour
subset=acc.ds['salinity'].sel(nSCHISM_hgrid_node=nns,nSCHISM_vgrid_layers=-1)

# load model time series
values=subset.salinity.values    

plt.figure()
plt.plot(model_times,subset.salinity)


print("Mean error:", mean_error)
print("RMSE:", rmse)
print("MAE:", mae)


# Assuming stations are in the same order as the model nodes in `subset`
station_array = np.column_stack([
    stations[key].data['value'].values
    for key in stations.keys()
])  # s

model_array = subset.salinity.values  # shape: (ntime, nstations)


residuals = station_array - model_array   # same shape
abs_errors = np.abs(residuals)
rmse_per_station = np.sqrt(np.nanmean(residuals**2, axis=0))
mae_per_station = np.nanmean(abs_errors, axis=0)


overall_rmse = np.sqrt(np.nanmean(residuals**2))
overall_mae = np.nanmean(abs_errors)

s.plot_domain_boundaries(append=True,plot_legend=False,latlon=False)
# Scatter stations colored by RMSE
ph = plt.scatter(
    xs,                 # Easting (UTM X)
    ys,                 # Northing (UTM Y)
    c=rmse_per_station, # RMSE per station
    s=30,               # marker size
    vmin=0, vmax=5,     # color scale limits
    cmap='viridis',     # optional colormap
)
# Add colorbar
cb = plt.colorbar(ph)
cb.set_label("RMSE [pmill]")

plt.title("Station RMSE vs Model")
plt.xlabel("UTM Easting (m)")
plt.ylabel("UTM Northing (m)")
plt.show()

plt.figure(figsize=(12,4))
for i, key in enumerate(stations.keys()):
    plt.plot(model_times, residuals[:, i], label=key)
plt.axhline(0, color='k', linestyle='--')
plt.legend()
plt.ylabel("Residual [pmill]")
plt.xlabel("Time")
plt.show()

