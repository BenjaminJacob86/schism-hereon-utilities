## import sediment data from pange via
## pangeapy


#https://doi.org/10.1594/PANGAEA.882561
#This dataset comprises results from filter weights of more than 5000 water samples taken during numerous field surveys between 1998 and 2017 in the Odra Lagoon (German Baltic Sea coast), several parts of the German Wadden Sea, the Exclusive Economic Zone of Germany in the German Bight (Southern North Sea), the Limfjorden (Denmark), the Oosterschelde (The Netherlands) and the Ria de Vigo (Spain). From the filter weights and filtered water volumes the suspended particulate matter concentrations (SPMC) and, in most cases, the fractions of organic matter were determined by combustion of the loaded filters (Loss on Ignition - LoI). Over the years, the laboratory methods and the type of filter (Whatman GF/C glass fibre filter, 47 mm diameter) were kept identical, but the sampling methods were adapted to technical demands, to the specific conditions of the sampling areas and to novel methodological insight. The samples had to undergo a number of quality checks regarding sampling time and co-ordinates and all laboratory processing steps. Depending on which test they passed they were assigned two types of quality flags, for (1) space and time information and (2) for the sample itself. They range from 1 ("good": all tests passed) to 4 ("bad and not correctable"); further, 9 is assigned for cases with missing information of time, latitude/longitude or water pressure. Further, the filter weights were corrected for filter offsets and - when Loss on Ignition was also determined - for losses of structural water.

import numpy as np
import pangaeapy.pandataset as pgd
import pandas as pd
from matplotlib import pyplot as plt
import sys
import os
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')
from schism import *
import xarray as xr
import glob

setupgdir='/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/SEDIMENT/'
ncdir=setupgdir+'combined_deflate/'
os.chdir(setupgdir)
s=schism_setup()

schismfiles=[]
for iorder in range(6): # check for schout_nc files until 99999
	schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
nrs=list(np.asarray(nrs)[np.argsort(nrs)])
s.nc=xr.open_mfdataset(schismfiles[6:])  #error 5


datasets

ds= pgd.PanDataSet(882561)

				 897292 # Riethmüller, Rolf; Flöser, Götz (2019): Suspended particulate matter concentrations and organic matter fractions from water samples - update 2017/2018. Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research, PANGAEA, https://doi.org/10.1594/PANGAEA.897292
ds= pgd.PanDataSet(897292)
date_colum='Date/Time'
date_colum='Date/time start'
#25	Suspended matter, particulate/solids	TSS	mg/l

#t0=np.datetime64('2017-01-01') 
#t1=np.datetime64('2018-01-01') 
t0=s.nc['time'][0].values
t1=s.nc['time'][-1].values
t=s.nc['time'].values
#filter=(ds.data.Area=='German Bight') & ( t0 <= ds.data['Date/Time']) & ( t1 >= ds.data['Date/Time']) & (ds.data['Depth water']< 2)

lon=np.asarray(s.lon)
lat=np.asarray(s.lat)
#882561
inDistTrehs=np.asarray([np.min((coord[0]-lon)**2+(coord[1]-lat)**2)<0.3 for coord in list(zip(ds.data['Longitude'],ds.data['Latitude']))])
filter=inDistTrehs & ( t0 <= ds.data['Date/Time']) & ( t1 >= ds.data['Date/Time']) & (ds.data['Depth water']< 2)
#filter=(ds.data.Area=='German Bight') & ( t0 <= ds.data['Date/Time']) & ( t1 >= ds.data['Date/Time']) & (ds.data['Depth water']< 2)
#ds.data=ds.data.loc(filter)
#ds.data=ds.data[filter]

#897292
ds= pgd.PanDataSet(897292)
ds.data[date_colum]=pd.to_datetime(ds.data[date_colum])
inDistTrehs=np.asarray([np.min((coord[0]-lon)**2+(coord[1]-lat)**2)<0.3 for coord in list(zip(ds.data['Longitude'],ds.data['Latitude']))])
filter=inDistTrehs & ( t0 <= ds.data[date_colum]) & ( t1 >= ds.data[date_colum]) & (ds.data['Depth water']< 2)
ds.data=ds.data[filter]



s.init_node_tree()

#for key in ds.data.head().keys():
#	print(ds.data[key][tis])

data=[]
model=[]
times=[]
depths=[]
lons=[]
lats=[]
campaigns=[]
areas=[]
deltaT=np.timedelta64(30,'m') # allowed time difference
for ti in np.unique(ds.data[date_colum]):
	print(str(ti))
	#break
	
	
	deltas=t-ti
	if (np.abs(deltas)<deltaT).sum()>0:
		tnn=np.argmin(np.abs(t-ti))
	else:
		continue
	
	
	# average
	
	tis=(ds.data[date_colum]==ti)
	
	#useforaverage=ds.data['Depth water'][tis]==np.min(q[tis])
	useforaverage=ds.data['Depth water'][tis]==ds.data['Depth water'][tis].min()
	if (tis.sum()>1) and np.diff(ds.data.Latitude[tis]).max() > 0:
		print('deviation')
		break
		
	qlon=ds.data.Longitude[tis].values[0]
	qlat=ds.data.Latitude[tis].values[0]
	vmean=ds.data.TSS[tis][useforaverage].mean() #mg/l
	#ds.data.Latitude[tis]
	#ds.data.Latitude[tis]
	#ds.data.Longitude[tis]
	#ds.data['Date/Time'][tis]
	nn=s.node_tree_latlon.query((qlon,qlat))[1]
	vschismInt=np.sum( [s.nc['SED3D_{:d}'.format(i)][tnn,nn,:].values for i in range(1,9)],axis=0)*1000  #[mg/l] #int 
	#sample_depth=-np.unique(ds.data['Depth water'][tis])
	sample_depth=-np.unique(np.round(ds.data['Depth water'][tis],1))
	
	if len(sample_depth) > 1:
		print('to much samples')
		break
	vschismInt=np.interp(sample_depth,s.nc['zcor'][tnn,nn,:],vschismInt)
	#classes
	#vschismIntb=np.sum( [s.nc['SED3D_{:d}'.format(i)][tnn,nn,-2].values for i in range(1,9)])*1000  #[mg/l] #int classes

	
	data.append(vmean)
	model.append(vschismInt)
	times.append(ti)
	depths.append(sample_depth)
	lons.append(qlon)
	lats.append(qlat)
	campaigns.append(list(ds.data['Campaign'][tis])[0])
	areas.append(list(ds.data['Area'][tis])[0])
	#s.node_tree_latlon.query(list(zip(ds.data.Longitude[tis],ds.data.Latitude[tis])))
	#ds.data['Depth water'][tis]
	#plt.clf()
	#s.plot_domain_boundaries(append=True)
	#plt.scatter(ds.data.Longitude[tis],ds.data.Latitude[tis],s=10,c=ds.data.TSS[tis],cmap=plt.cm.jet)
	#plt.title(str(t))
	#plt.colorbar()
	#ds.data.['Depth water'][tis]
 
pairs=pd.DataFrame({'time':times,'lon':lons,'lat':lats,'depth':depths,'campaign':campaigns,'SPM mes':data,'SPM sim':model}) 

 pairs['SPM sim']/pairs['SPM mes']

plt.figure()
plt.subplot(1,2,1) 
plt.scatter(pairs['SPM mes'],pairs['SPM sim'],c=range(len(pairs['SPM sim']))) #,cmap=plt.cm.jet()
plt.set_cmap('jet')
plt.xlabel('TSS ' + pairs['campaign'].values[0] + ' [mg/l]')
plt.ylabel('SPM SCHISM [mg/l]')
plt.subplot(1,2,2) 
s.plot_domain_boundaries(append=True)
#plt.plot(pairs['lon'],pairs['lat'],'ko')
plt.scatter(pairs['lon'],pairs['lat'],c=range(len(pairs['SPM sim']))) #,cmap=plt.cm.jet()
plt.set_cmap('jet')

#TSS suspended matter, particulate solits