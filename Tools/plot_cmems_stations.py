#/gpfs/home/jacobb/anaconda3/bin/python3 plot_station_positions.py
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from glob import glob

# cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

plt.ion()
import numpy as np


### cartopy ########
proj=ccrs.Mercator()  # define Projection
## load higher resolutions coastline assets
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',edgecolor='face',                                        facecolor=cfeature.COLORS['land'])
ocean_10m = cfeature.NaturalEarthFeature('physical', 'ocean', '10m', edgecolor='face',facecolor=[0,0,1])
area=(20,42,38,48)


tmin=np.datetime64('2021-11-01')
tmax=np.datetime64('2021-12-01')

# Data type
# NetCDF file names are as follow: GL_XX_YY_CODE_<_ZZZ>.nc 
# GL: region bigram corresponding to global area
# XX:
# TS (timeseries)
# PR (profile)
# YY:
# BO bottle data
# CT oceanographic CTD profiles
# FB ferrybox
# GL gliders
# ML mini logger for fishery observing system
# MO fixed buoys, mooring time series, fixed observations
# PF profiling floats SM Sea mammals 

types=['FB'] 
#types=['BO', 'CT', 'DB', 'DC', 'FB', 'GL', 'HF', 'ML', 'MO', 'PF', 'RF', 'SD', 'SF', 'SM', 'TG']
types=['BO','CT','DB','FB','GL','HF','ML','MO','PF','RF','TG']
year='2021'
month='06'


land_color='gray'
fill_color='none'


typedict={'BO': 'bottle data','CT': 'oceanographic CTD profiles','DB': 'drifting buoys','DC': 'drifting buoy reporting calculated sea water current','FB':'ferrybox','GL':'gliders','HF':'HF radar','ML':'mini logger for fishery observing system','MO':'fixed buoys, mooring time series, fixed observations','PF':'profiling floats','RF':'River flow','SD':'Saildrone','SF':'towed CTD data (ex: scanfish)','SM':'Sea mammals data','TG':'Tide gauge station'}
#GL  gliders
#AR  'Aritiert

mobiles=['DB','DC','FB','GL','SD']


plt.clf()
ax = plt.axes(projection=proj)
ax.set_extent(area)
#ax.redraw_in_frame(zoom_extend)
plt.plot(stati['lon'],stati['lat'],'+',color=cmap(i),transform=proj,markersize=30)
outproj=ax.projection.transform_points(ccrs.Mercator(),stati['lon'],stati['lat'])
ax.add_feature(land_10m,zorder=-2)

# Set tick labels
xticks=np.unique(np.round(np.linspace(area[0],area[1],8)))
#yticks=np.arange(np.round(lat.min()),lat.max(),0.5)
yticks=np.unique(np.round(np.linspace(area[2],area[3],8)))
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.1f',degree_symbol='',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f',degree_symbol='')
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
plt.tight_layout()


for type in types:

	#folder='nrt.cmems-du.eu/Core/INSITU_NWS_NRT_OBSERVATIONS_013_036/nws_multiparameter_nrt/monthly/'+type+'/'
	folder='./'+type+'/'
	#folders=glob(folder+'??????')
	#folders=glob(folder+year+month)
	#files=glob(folder+'*'+year+'*')
	AR=glob(folder+'*')

#for file in files:
	stations=[]
	if type in mobiles:
		for ARi in AR:
		  ds=xr.open_dataset(ARi)	
		  
		  time=ds['TIME'].values
		  if ((tmin < time) & (tmax > time)).sum():
			stations.append({'name':ARi[ARi.rindex(type)+3:ARi.rindex('_')] ,'lon':ds['LONGITUDE'][:].values,'lat':ds['LATITUDE'][:].values},'time':ds['TIME'].values)
		  ds.close()
	else:
		for ARi in AR:
		  ds=xr.open_dataset(ARi)	
		  time=ds['TIME'].values
		  if ((tmin < time) & (tmax > time)).sum():
			stations.append({'name':ARi[ARi.rindex(type)+3:ARi.rindex('_')] ,'lon':ds['LONGITUDE'][0].values,'lat':ds['LATITUDE'][0].values},'time':ds['TIME'].values))
		  ds.close()
		  
	cmap=plt.cm.get_cmap('jet',len(stations))
	if type in mobiles:
		for stati in stations:
			#m.plot(stati['lon'],stati['lat'],'r+')
			#xi,yi=m(stati['lon'],stati['lat']) 
			outproj=ax.projection.transform_points(ccrs.Geodetic(),stati['lon'],stati['lat'])
			lon,lat=outproj[:,0],outproj[:,1]
			plt.plot(lon,lat,'+',color=cmap(i))

			
			try:
				#m.plot(stati['lon'],stati['lat'],latlon=True)
				plt.plot(lon,lat,'.',color=cmap(i),transform=proj)
			except:
				#m.plot(stati['lon'][0],stati['lat'][0],latlon=True)
				plt.plot(lon,lat,'.',color=cmap(i),transform=proj)
			#plt.plot(xi[:1000:3],yi[:1000:3])
			#plt.plot(xi,yi)
			#plt.text(xi,yi,stati['name'])
	else:
		for i,stati in enumerate(stations):
			#plt.plot(stati['lon'],stati['lat'],'+',color=cmap(i),transform=proj)
			outproj=ax.projection.transform_points(ccrs.Geodetic(),stati['lon'],stati['lat'])
			lon,lat=outproj[:,0],outproj[:,1]
			plt.plot(lon,lat,'+',color=cmap(i))
			#xi,yi=m(stati['lon'],stati['lat']) 
			#plt.plot(xi,yi,'r+')
			#plt.text(xi,yi,stati['name'])
	if len(stations) < 20:		
		plt.legend([stati['name'] for stati in stations], bbox_to_anchor=(-0.4, 1), loc='upper left',ncol=1)	
	else:	
		plt.legend([stati['name'] for i,stati in enumerate(stations) if i <20 ]+['...'], bbox_to_anchor=(-0.4, 1), loc='upper left',ncol=1)	
		
	ds=xr.open_dataset(ARi)	
	vars=np.unique([ds[key].long_name for key in ds.variables.keys()])
	ds.close()
	
	txt='variables:'
	for i in range(len(vars)):
		txt+='\n'+vars[i] 
	plt.text(0.78, 0.1, txt, fontsize=7, transform=plt.gcf().transFigure)
		
	plt.savefig('sations_black'+type+year+month+'.png',dpi=300)	
	#plt.show()	  
