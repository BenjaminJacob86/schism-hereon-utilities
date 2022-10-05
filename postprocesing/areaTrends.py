#area trends
import os
import netCDF4
import sys
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
import matplotlib
matplotlib.use('AGG')
from netCDF4 import Dataset
import numpy as np
from matplotlib import path
from matplotlib import pyplot as plt
import datetime as dt
import glob
from scipy.spatial import cKDTree
import xarray as xr
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
import pickle
from schism import *
import datetime as dt

# calculate spatial mean as function of time for defined polygons (kml from google earth)

########## settings #################################
# export OMP_NUM_THREADS=4 # before call
outname='BlackTrends2008.nc'
outname='EuropeTrends2012.nc'

kmlfiles=['NorthAgeanSea.kml','NorthAgeanSeaBSdomain.kml']
setupdir='/work/gg0028/g260114/RUNS/Europe/c1/flow_tweak/'
setupdir='/gpfs/work/jacobb/data/SETUPS/BlackSea/'
ncdir=setupdir+'/combined/'

ocnfile='/work/gg0028/g260114/RUNS/postproc/water_mass_balance/OceanDB.kml'  #ocean regions kml file
ocnfile='/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Tools/CalcWaterMassBalance/OceanDB.kml'  #ocean regions kml file

basins=['Mediterranean Sea','Sea of Marmara','Black Sea']
basins=['Sea of Marmara','Black Sea']
kmlfiles=['NorthAgeanSeaBSdomain.kml']

nspool = 960 
ihfskip = 9600 
nstep=np.int(ihfskip/nspool)

nstep=20 # ihfskip/nspool

cwd=os.getcwd()
def google_kml_coords(kmlfiles):
	"""
		extract coordinates from google earth polygon kml file
		as array of x,y tuples
	"""
	
	seas={}
	for kmlfile in kmlfiles:
		name=[kmlfile[:kmlfile.index('.')]][0]
		with open(kmlfile,'r') as infile:
			prevline=[]
			for line in infile:
				if '<coordinates>' in prevline:
					temp=line.replace('\t','').replace(' \n','').split(' ')
					coords=[]
					for coord in temp:
						coords.append(tuple(np.double(coord.split(',')[:2])))
					seas[name]=np.asarray(coords)
				prevline=line
	return seas

def kml_coords(kmlfile):
	"""
        extract coordinates from google earth polygon kml file
        as array of x,y tuples
	"""
	seas={}
	start = '<name>'
	end = '<\name>'
	string=''
	read_cords=False
	name=''
	
	with open(kmlfile,'r') as infile:
		prevline=[]
		for line in infile:
			
			if ('name' in line) and ('Events' not in line) and ('<Folder>' in prevline):
				name=line[line.find(start)+len(start):line.rfind(end)-len(end)-1]
			if ('<coordinates>' in prevline) and ('<Placemark>' not in prevline): #
				read_cords=True
				coords=[]
			if  ('</coordinates>' in line) and ('<Placemark>' not in line):
				seas[name]=np.asarray(coords)
				read_cords=False
			if read_cords:
				coords.append( (float(line.split(',')[0]),float(line.split(',')[1]) ) )
			prevline=line
	
	return seas
	
def extrac_basin_averaged_timeseries(hxarr,faces,A,eleminds,regname='Baltic Sea',varnames='zcor',nt=np.inf):
	""" extract spatial mean time series for define region of SCHISM grid using xarray acces """
	trinodes=faces[eleminds[regname]]
	nodes=np.unique(trinodes)
	trinodes_reduced=np.zeros(trinodes.shape,int)
	for i,node in enumerate(nodes):
		node
		trinodes_reduced[np.where(trinodes==node)]=i
	selection=hxarr.sel(nSCHISM_hgrid_node=nodes)
	
	if nt==np.inf:
		nt=len(selection['time'])
	t=selection['time'][:nt]
	
	if type(varnames) != list:
		varnames=[varnames,]
	basin_ts=dict.fromkeys(varnames)	
	for varname in varnames:
		if len(hxarr[varname].shape)==2:
			basin_ts[varname]=(np.asarray(selection[varname][:nt,:])[:,trinodes_reduced].mean(axis=2)*A[eleminds[regname]]).sum(axis=1)/A[eleminds[regname]].sum()
		else:
			basin_ts[varname]=(np.asarray(selection[varname][:nt,:,-1])[:,trinodes_reduced].mean(axis=2)*A[eleminds[regname]]).sum(axis=1)/A[eleminds[regname]].sum()
	
	return t,basin_ts

def datetime64_to_datetime(t):
	if len(t)==1:
		t=[t]
	return np.asarray([ dt.datetime.utcfromtimestamp((ti - np.datetime64('1970-01-01T00:00:00Z'))/ np.timedelta64(1, 's')) for ti in t])
	
	
s=schism_setup()	
faces=s.nvplt
x=np.asarray(s.lon)
y=np.asarray(s.lat) 
# element centers
cx=np.mean(x[faces],axis=1)
cy=np.mean(y[faces],axis=1)
elcoords=[[cx[i],cy[i]] for i in range(len(cx))] # pooint pairs of nodes
elem_nn_tree = cKDTree(elcoords) # next neighbour search tree	     

# get element indices for basins
seas=kml_coords(ocnfile)	
regions=google_kml_coords(kmlfiles)	
key=list(regions.keys())
#seas[key[0]]=regions[key[0]]
seas={**seas,**regions}


basins+=list(regions.keys())
eleminds={}
regioncolor=np.zeros(s.nnodes)	
for ifile,tag in enumerate(basins):
	areaPoly=path.Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
	eleminds[tag]=areaPoly.contains_points(elcoords)
	regioncolor[s.nvplt[eleminds[tag]]]=ifile+1

plt.close()	
ph,ch=s.plotAtnodes(regioncolor)
ch.set_label('regions')
ch.set_ticks(range(len(eleminds)+1))
#plt.set_cmap('hsv')
ch.set_ticklabels(['']+list(eleminds.keys()))
plt.tight_layout()
plt.savefig('regions.png')
plt.close()	


#quad
A=[]
for i in range(s.nvplt.shape[0]):
	nodes=s.nvplt[i,:]+1
	A.append(s.proj_area(nodes))
A=np.asarray(A)

################   net cdf access   ########################
## check for schout_nc files
schismfiles=[]
for iorder in range(6):
    iorder
    schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')

# regular schismfiles for acces
nrs=[]  
for file in schismfiles:
    nr=int(file[file.rfind('_')+1:file.index('.nc')])
    nrs.append(nr)
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
nrs=list(np.asarray(nrs)[np.argsort(nrs)])

# initialize file acces
ds={}
ds['schism']=xr.open_mfdataset(schismfiles)


s.hgridnodes=np.arange(s.nnodes)

varnames=['elev','salt','temp']
ti=0
nt=len(ds['schism']['time'])
Warea={regname:np.tile(A[eleminds[regname]],(nstep,1))/A[eleminds[regname]].sum() for regname in basins}
# extract mean for region
basin_ts=dict.fromkeys(basins)
for key in basin_ts.keys():
	basin_ts[key]={varname:[] for varname in varnames}
while(ti<nt):
	for varname in varnames:
		if len(ds['schism'][varname].shape)==2:
			data_exert=ds['schism'][varname][ti:ti+nstep,s.hgridnodes].values
		else:	
			data_exert=ds['schism'][varname][ti:ti+nstep,s.hgridnodes,-1].values
			
		for regname in basin_ts.keys():		
			trinodes=faces[eleminds[regname]]
			basin_ts[regname][varname].append((data_exert[:,trinodes].mean(axis=2)*Warea[regname]).sum(axis=1))
	if ti%(nstep*100)==0:
		print(ti/nstep*100)
	ti+=nstep
	
basin_ts[regname][varname]=np.concatenate(basin_ts[regname][varname])		
for key in basin_ts.keys():
		for varname in varnames:
			try:
				basin_ts[key][varname]=np.concatenate(basin_ts[key][varname])	
			except:
				pass
###################### other model#####################


ab=datetime64_to_datetime([ds['schism']['time'][0].values])[0]
basedate=(ab.year,ab.month,ab.day,0,0,0)

# put to nc
from netCDF4 import Dataset
from cftime import utime
nc = Dataset(outname,'w', format='NETCDF4') #'w' stands for write
nc.createDimension('time',None)
tv = nc.createVariable('time','f8',('time'))
tv.long_name = 'Time'
tv.standard_name = 'time'
tv.units = 'days since %4i-%02i-%02i %02i:%02i:%02i'%basedate
tv.base_date =  list(basedate)
tout=datetime64_to_datetime(ds['schism']['time'][:])
tv[:]=np.asarray([(ti-tout[0]).total_seconds()/86400 for ti in tout])
#ds['schism']['time'][:]
for key in basin_ts.keys():
	for varname in varnames:
		outvarname=key.replace(' ','')+'_mean_'+varname
		tv = nc.createVariable(outvarname,'f8',('time'))
		tv[:]=basin_ts[key][varname]
nc.close()


		
# calc anual mean		
year_means=ds['schism'].groupby('time.year').mean('time')
varname='elev'

for i,year in enumerate(year_means.year):
	#plt.subplot(2,2,i+1)
	plt.clf()
	ph,ch=s.plotAtnodes(year_means[varname][i,:])		
	plt.title(year)
	ch.set_label('<ssh> [m]')
	cxs=[]
	cys=[]
	for key in basin_ts.keys():
		if key=='NorthAgeanSea':
			continue
		trinodes=faces[eleminds[key]]
		cx=x[trinodes].mean()
		cy=y[trinodes].mean()
		cxs.append(cx)
		cys.append(cx)
		meanval=(year_means['elev'][i,:].values[trinodes].mean(axis=1)*Warea[key][0,:]).sum()
		plt.plot(seas[key][:,0],seas[key][:,1],'k--')
		plt.text(cx,cy,'mean : {:.2f}'.format(meanval))
	plt.xlim(np.min(cxs)-3,np.max(cxs)+3)	
	plt.xlim(np.min(cys)-3,np.max(cys)+3)	
	plt.tight_layout()
	plt.savefig('ssh_mean_year_{:d}'.format(int(year)),dpi=300)
	

	
# black Sea means
xmin,xmax=np.min(s.lon),np.max(s.lon)		
ymin,ymax=np.min(s.lat),np.max(s.lat)		
xmin=24.00605199
xmax=41.74752356	
ymin=39.00458572
ymax=47.25219778


labels=['<ssh> [m]','<sst> [degC]','<sss> [Psu]']
for ivar,varname in enumerate(varnames):
	for i,year in enumerate(year_means.year):
		plt.clf()
		if varname=='elev':
			ph,ch=s.plotAtnodes(year_means[varname][i,:])		
		else:	
			ph,ch=s.plotAtnodes(year_means[varname][i,:,-1])		
		ch.set_label(labels[ivar])
		plt.title(year)
		plt.xlim(xmin,xmax)	
		plt.ylim(ymin,ymax)	
		#plt.title('BlackSea '+ str(int(year)))
		if i==0:
			plt.clim((-0.6,-0.2))
		elif i==1:
			plt.clim((-0.5,-0.3))
		else:	
			plt.clim((-0.5,-0.1))
		plt.title('EU-BlackSea '+ str(int(year)))
		plt.tight_layout()
		#plt.savefig('sttandalone_blacksea_ssh_mean_year_{:d}'.format(int(year)),dpi=300)	
		plt.savefig('neuEurope_blacksea_'+str(varname)+'_mean_year_{:d}'.format(int(year)),dpi=300)	
	


	
import numpy as np	
import datetime as dt
from netCDF4 import Dataset
from matplotlib import pyplot as plt	

ncbs=Dataset('BlackTrends2008b.nc')
nceu=Dataset('EuropeTrends2012.nc')   

nceu=Dataset('EuropeTrends.nc')   

area='SeaofMarmara_mean_elev'	
t0bs=dt.datetime.strptime(ncbs['time'].units[11:],'%Y-%m-%d %H:%M:%S')
datesbs=t0bs+np.arange(1,len(ncbs['time'])+1)*dt.timedelta(days=np.diff(ncbs['time'][:2])[0])

t0eu=dt.datetime.strptime(nceu['time'].units[11:],'%Y-%m-%d %H:%M:%S')
dateseu=t0eu+np.arange(1,len(nceu['time'])+1)*dt.timedelta(hours=6)

bsyear=np.asarray([ti.year for ti in datesbs])
euyear=np.asarray([ti.year for ti in dateseu])

plt.tight_layout()
for year1,year2 in zip([2008,2009,2010],[2012,2013,2014]):
	plt.clf()
	plt.subplot(2,1,1)
	for varname in 'NorthAgeanSeaBSdomain_mean_elev','SeaofMarmara_mean_elev','BlackSea_mean_elev':
		plt.plot(datesbs[bsyear==year1],ncbs[varname][bsyear==year1])
	plt.grid()
	plt.ylabel('spat. avg ssh [m]')
	plt.legend(['N. Agean','Marmara','Black Sea'])	
	plt.title('Black Sea' + str(year1))	
	plt.tight_layout()
	plt.subplot(2,1,2)
	for varname in 'NorthAgeanSeaBSdomain_mean_elev','SeaofMarmara_mean_elev','BlackSea_mean_elev','MediterraneanSea_mean_elev':
		plt.plot(dateseu[euyear==year2],nceu[varname][euyear==year2])
	plt.grid()
	plt.ylabel('spat. avg ssh [m]')
	#plt.legend(['N. Agean','Marmara','Black Sea','Meditaranean'])	
	plt.legend(['N. Agean','Marmara','Black Sea','Meditaranean'],loc='upper center',ncol=4)		
	plt.title('EU - Black Sea' + str(year2))	
	plt.tight_layout()
	plt.savefig('basin_trends_'+str(year1)+'vs'+str(year2))


for year1,year2 in zip([2008,2009,2010],[2012,2013,2014]):
	plt.clf()
	plt.subplot(2,1,1)
	for varname in 'NorthAgeanSeaBSdomain_mean_salt','SeaofMarmara_mean_salt','BlackSea_mean_salt':
		plt.plot(datesbs[bsyear==year1],ncbs[varname][bsyear==year1])
	plt.grid()
	plt.ylabel('spat. avg salinity [PSU]')
	#plt.legend(['N. Agean','Marmara','Black Sea'])	
	plt.title('Black Sea' + str(year1))	
	plt.subplot(2,1,2)
	for varname in 'NorthAgeanSeaBSdomain_mean_salt','SeaofMarmara_mean_salt','BlackSea_mean_salt','MediterraneanSea_mean_salt':
		plt.plot(dateseu[euyear==year2],nceu[varname][euyear==year2])
	plt.grid()
	plt.ylabel('spat. avg salinity [PSU]')
	plt.legend(['N. Agean','Marmara','Black Sea','Meditaranean'],loc='upper center',ncol=4)		
	plt.title('EU - Black Sea' + str(year2))	
	plt.tight_layout()
	plt.savefig('basin_trends_salt_'+str(year1)+'vs'+str(year2))

for year1,year2 in zip([2008,2009,2010],[2012,2013,2014]):
	plt.clf()
	plt.subplot(2,1,1)
	for varname in 'NorthAgeanSeaBSdomain_mean_temp','SeaofMarmara_mean_temp','BlackSea_mean_temp':
		plt.plot(datesbs[bsyear==year1],ncbs[varname][bsyear==year1])
	plt.grid()
	plt.ylabel('spat. avg temperature [degC]')
	#plt.legend(['N. Agean','Marmara','Black Sea'])	
	plt.title('Black Sea' + str(year1))	
	plt.subplot(2,1,2)
	for varname in 'NorthAgeanSeaBSdomain_mean_temp','SeaofMarmara_mean_temp','BlackSea_mean_temp','MediterraneanSea_mean_temp':
		plt.plot(dateseu[euyear==year2],nceu[varname][euyear==year2])
	plt.grid()
	plt.ylabel('spat. avg temperature [degC]')
	plt.legend(['N. Agean','Marmara','Black Sea','Meditaranean'],loc='upper center',ncol=4)		
	plt.title('EU - Black Sea' + str(year2))	
	plt.tight_layout()
	plt.savefig('basin_trends_temp_'+str(year1)+'vs'+str(year2))


	
# volume changes	
varname='BlackSea_mean_elev'
hmean=np.asarray([nceu[varname][euyear==year].mean() for year in range(2013,2017)])
dh=np.diff(hmean)
dv=A[eleminds[regname]].sum()*dh   # m³
dv=A[eleminds[regname]].sum()*dh/1000**3   # km³

varname='NorthAgean_mean_elev' 
hmean2=np.asarray([nceu[varname][euyear==year].mean() for year in range(2013,2017)])

plt.plot(pyears,hmean-hmean2)
plt.grid()
plt.ylabel('Dh m')
plt.tight_layout()
plt.title(' Blacksea<ssh> - N. agea<ssh>')
plt.savefig('blacksea_develope.png',dpi=300)

#Hdecrease_in_year=np.asarray([nceu[varname][euyear==year][-28:].mean()-nceu[varname][euyear==year][:28].mean() for year in range(2013,2017)])
Hdecrease_in_year=np.asarray([nceu[varname][euyear==year][-(4*5):].mean()-nceu[varname][euyear==year][:(4*5)].mean() for year in range(2013,2017)])
dv=A[eleminds[regname]].sum()*Hdecrease_in_year/1000**3   # km³

#2013 to 2016
netsumAll=[-74.57593124,30.55830667,-129.0203637,-33.91895648]
netsumFlux=[-244.057884,-170.0869803,-292.224382215133,-213.0067796]

pyears=range(2013,2017)
plt.clf()
plt.plot(pyears,dv)
plt.plot(pyears,netsumAll)
plt.plot(pyears,netsumFlux)
plt.legend(['dv from dh','Netflux(evp,Q,strait)','strait net out'])
plt.grid()
plt.ylabel('Dv km³/a')
plt.tight_layout()
plt.savefig('blacksea_develope.png',dpi=300)


# Trend Amplifitiert in Black Sea
# Barotropische karft
# schwarzes meer sink stärker als mittelmeer
# verschiedene Gebiete Auslenekungen relative zu einander stabil
# S. Grayek, Idee



-74.57593124	-1.535489533


regname='Black Sea'






	
ncbs['NorthAgeanSeaBSdomain_mean_elev'][:]==ncbs['SeaofMarmara_mean_elev'][:]
 
basin_ts['Sea of Marmara'][varname]==basin_ts['Black Sea'][varname]
 
years=year


#datesbs=t0bs+np.arange()*dt.timedelta(seconds=np.diff(ncbs['time'][:2])[0]/86400)



dt.datet
ncbs['time'][:2]


ncbs[area]
years
ncbs[area]
ncbs[area]



BlackSea_mean_elev
NorthAgeanSeaBSdomain_mean_salt


		
reg_avg=extrac_basin_averaged_timeseries(ds['schism'],faces,A,eleminds,regname='NorthAgeanSea',varnames=['zcor','salt'],nt=np.inf)




## extract for Grayeks modell
# varnames
name_lon='lon'
name_lat='lat'
name_time='time'
name_depth='depth'

name_ssh='sossheig'
name_salt='vosaline'
name_temp='votemper'
name_u='vozocrtx'
name_v='vomecrty'


from glob import glob
folder='/work/gg0028/g260099/GCOAST_10y_OUTPUT/GCWAVE_COUPLED_02_RUN/20120631' # patterns 
folders=np.sort([glob('/work/gg0028/g260099/GCOAST_10y_OUTPUT/GCWAVE_COUPLED_02_RUN/2012*') + glob('/work/gg0028/g260099/GCOAST_10y_OUTPUT/GCWAVE_COUPLED_02_RUN/2013*') + glob( '/work/gg0028/g260099/GCOAST_10y_OUTPUT/GCWAVE_COUPLED_02_RUN/2014*')])[0]
ind=list(folders).index('/work/gg0028/g260099/GCOAST_10y_OUTPUT/GCWAVE_COUPLED_02_RUN/20120601')
folders=folders[ind:]


# gather filelist
ssh_pattern='GCWAVE_COUPLED_02_1h_OUTPUT' # shh 
ssh_files=[glob(folder+'/*'+ssh_pattern+'*') for folder in folders]
# temp / salt
temp_pattern='GCWAVE_COUPLED_02_6h_OUTPUT'
temp_files=[glob(folder+'/*'+temp_pattern+'*') for folder in folders]
salt_pattern='GCWAVE_COUPLED_02_6h_OUTPUT'
salt_files=temp_files#[glob(folder+'/*'+salt_pattern+'*') for folder in folders]
uv_pattern='GCWAVE_COUPLED_02_3h_OUTPUT' # uv
uv_files=[glob(folder+'/*'+uv_pattern+'*') for folder in folders]

# establish file access
ds={'ssh':xr.open_mfdataset(ssh_files),
'salt':xr.open_mfdataset(salt_files),
'uv':xr.open_mfdataset(uv_files)
}




varnames={'ssh':name_ssh,'salt':name_salt,'temp':name_temp,'uv':name_u}
grid=dict.fromkeys(['ssh','salt','temp','uv'])
for key in grid.keys():
	grid[key]={}
	for var in [name_lon,name_lat,name_depth,name_time]:
		try:
			grid[key][var]=ds[key][var].values
		except:
			pass
	grid[key]['dimensions']=ds[key][varnames[key]].dims # oder of dimentions


	
ds['ssh'][name_ssh]

	
	
lon2d,lat2d=np.meshgrid(ds['ssh'][name_lon][:].values,ds['ssh'][name_lat][:].values)
	
grid['ssh']['lon'].shape
grid['ssh']['lon'].shape	
	
eleminds2={}
eleminds3={}
for ifile,tag in enumerate(basins):
	areaPoly=path.Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
	eleminds2[tag]=np.where(areaPoly.contains_points(list(zip(lon2d.flatten(),lat2d.flatten()))))	
	ii,jj=np.unravel_index(eleminds2[tag],lon2d.shape)
	eleminds3[tag]={'ii':ii,'jj':jj}

for ti in range(len(grid['ssh']['time'])):
	ti
	
# extract mean for region
basin_ts2=dict.fromkeys(basins)
for key in basin_ts2.keys():
	basin_ts2[key]={varname:[] for varname in varnames}
ti=0
nstep=120 # ihfskip/nspool
while(ti<nt):
	#grid['ssh']['time']
	x1 = xr.DataArray(grid['ssh']['lon'], dims='z')
	y1 = xr.DataArray(qlat1, dims='z')
	ds['ssh'][name_ssh]

	for varname in varnames:
		if len(ds['schism'][varname].shape)==2:
			data_exert=ds['ssh'][name_ssh][ti:ti+nstep,:].values
		else:	
			data_exert=ds['schism'][varname][ti:ti+nstep,:,-1].values
			
		for regname in basin_ts.keys():		
			np.nanmean(np.squeeze(data_exert[:,eleminds3[regname]['ii'],eleminds3[regname]['jj']]),axis=1)
			basin_ts2[regname][varname].append((data_exert[:,trinodes].mean(axis=2)*Warea[regname]).sum(axis=1))
	if ti%(nstep*10)==0:
		print((ti/nstep)*100)
	ti+=nstep
	
	

	
	
# convert schism time	
t=ds['temp']['time'][:].values
dates={}	
for key in varnames.keys():	
	t=ds[key][name_time].values
	dates[key]=datetime64_to_datetime(t)
