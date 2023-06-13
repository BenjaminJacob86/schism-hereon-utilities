"""
SCHISM Wave modeling
Perform intercomparison of bsh bouys against
station output from WW3 and WWM model runs and WAM netcdf field extracted nearest neighbour timeseries.
Making plots of time series qqplot and taylordiagram
exportet as images and pdf.
"""
# imports
import os
from glob import glob
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
plt.ion()
import datetime as dt
import sys
# cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')
# own libraries 
from data_and_model_classes import bsh_spec, WW3_mesh
from validation_statistics_and_plots import interp_to_data_time,QQplot, TaylorDiagram,plotTaylor


########### SETTINGS #########
pdfname="WAM_SCHISM_HS_Validation_2017b.pdf" # name of image pdf
image_output_dir='/gpfs/work/jacobb/data/validation/waves/pics/'
skip_n_first_days=5 #skip first days of WWM where model is building up
skip_n_first_days=np.timedelta64(skip_n_first_days,'1D')

maxdate=np.datetime64('2017-10-31') # manually set maximumdate

bouydir='/gpfs/work/jacobb/data/validation/waves/bouy/bouys/' #bsh buoy directory containing subfolders with the bouys and in those a file with name <subfolder>_BSH.spec

# WWM run directories containing model station output
#'/gpfs/work/jacobb/data/RUNS/routine_GB_wave/from_mistral/Bmax1.52Code/','/gpfs/work/jacobb/data/RUNS/routine_GB_wave/from_mistral/Bmax1.3_succ_CodeCHange/',

## WWM runs
WWMdirs=['/gpfs/work/jacobb/data/validation/waves/wwm_veg_ref/',]
wwm_names=['WWM_2017',]  # names of scenario for display
wwm_descriptions=['SCHISM vegetation runs 2017/2018  WWM dt 240s']
comments='' 


## WAM runs
#WWMdirs=['/gpfs/work/jacobb/data/validation/waves/wwm_veg_ref/',]
#wam_names=['NEMO_WAM_SNS',]  # names of scenario for display
#wam_descriptions=['cuple NEMO_RC42 WAM SNS WD',]
#comments='' 

#No FIno1  wave after 2017


addWWM=True
addWam=False#True

display_names=[]
if addWWM:
	display_names+=wwm_names

if addWam:
	display_names+=wam_names
else:
	wam_names=[]
#keysWWM=['ELBE_joint', 'FINO1_joint', 'FINO3_joint', 'HELGON_joint', 'Westerland_joint']
#keysOBS=['ELB','FN1','FN3','HEL','WES']
#keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland','Westerland']


# no FIno1 after 2017
keysWWM=['ELBE_joint', 'FINO3_joint', 'HELGON_joint', 'Westerland_joint']
keysOBS=['ELB','FN3','HEL','WES']
keysWW3=['Elbe','Fino-3','Helgoland','Westerland']


#display_names=wam_names
# names of stations
#keysOBS=['ELB','FN1','FN3','HEL']
#keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland']
#keysWWM=['prep_ELBE', 'prep_FINO1', 'prep_FINO3', 'prep_HELGON']


from scipy.spatial import cKDTree
class wam_hs:
	def __init__(self,file):
		self.ds=xr.open_dataset(file)
		self.dates=self.ds['time'].values
		self.lat,self.lon=np.meshgrid(self.ds.lon.values,self.ds.lat.values)

	def	get_hs(self,lon,lat,date_slice=None):
		dsi = wam.ds.sel(lon = lonq, lat = latq, method = 'nearest')
		if date_slice != None:
			dsi=dsi.sel(time=slice(date_slice[0],date_slice[-1]))
		return dsi.time.values,dsi.hs.values
wam=wam_hs('/gpfs/work/jacobb/data/immerse/wave/WAVES_HS_2018.nc')


#keysOBS=['ELB','FN1','FN3','HEL','WES']
#keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland','Westerland']
#keysWWM=['prep_ELBE', 'prep_FINO1', 'prep_FINO3', 'prep_HELGON','prep_Westerland']


#keysOBS.remove('FN1')
#keysWWM.remove('FINO1_joint')
#keysWW3.remove('Fino-1')


## WWW3 ###########################################################################
addWW3=False

ww3dirs=['/gpfs/work/jacobb/data/RUNS/routine_GB_wave/WW4NBSbnd/','/gpfs/work/jacobb/data/RUNS/routine_GB_wave/WW4NBSbnd/GBgrid/']
ww3meshs=['NBSext_bl.msh','GB.msh']
ww3PointList='points.list'
ww3GridOutput=[ww3dirs[0]+'ww3.202108.nc',]
ww3PointOutput=[ww3dirs[0]+'ww3.2021_08_09_spec.nc',ww3dirs[1]+'ww3.202108_spec.nc']
ww3_names=['WW3Bmax1.5 DWD','WW3GBBmax1.5 DWD']#,'WW3Bmax1.6 DWD','WW3Bmax1.8 DWD'
ww3_descriptions=[': WW3 run with betamax of 1.5 forced with DWD forecast ',': WW3 run with betamax of 1.5 forced with DWD forecast on German Bight mesh.']

if not addWW3:
	ww3_names=[]


#parameter to validate
varnameWWM='HS'
varnameBSH='Hs'
unit='[m]'


##### appearance###################################################################

def set_FS(Fincrease=1.4):
	SMALL_SIZE = 8*Fincrease
	MEDIUM_SIZE = 10*Fincrease
	BIGGER_SIZE = 12*Fincrease
	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure
	###########################
set_FS(Fincrease=1.4)
#################################################



class bsh_ncbouy:
	def __init__(self,files):
		self.ds=xr.open_mfdataset(files)
		self.dates=self.ds['time'].values
	def	get_parameter(self,varname=''):
		return self.ds[varname.lower()].values


# make image output dir
if not os.path.exists(image_output_dir):
	os.mkdir(image_output_dir)		


# load buoys
subdirs=glob(bouydir+'*')
subdirs=np.sort([sdir[sdir.rindex('/')+1:] for sdir in subdirs])
bouys=dict.fromkeys(subdirs)
for subdir in subdirs:
	files=list(np.sort(glob(bouydir+'{:s}/*.spt'.format(subdir))))
	try:
		if len(files)>0:
			bouys[subdir]=bsh_spec(files=files)
			ncbouys=False
		else: #checheck netcdf 	
			files=list(np.sort(glob(bouydir+'{:s}/*.nc'.format(subdir))))
			bouys[subdir]=bsh_ncbouy(files=files)
			ncbouys=True
	except:
		bouys[subdir]=None

bouy_coords=[]#dict.from_keys()		
count=0
for key in (keysOBS):
	print('Buoys: loading '+key)
	if bouys[key]!=None:
		count+1
		lonq=float(bouys[key].ds.StartLongitude.split()[0].replace(',','.'))
		latq=float(bouys[key].ds.StartLatitude.split()[0].replace(',','.'))
		bouy_coords.append((lonq,latq))	

#######  WAM ######
if addWam:
	wam=wam_hs('/gpfs/work/jacobb/data/immerse/wave/WAVES_HS_2018.nc')
	parameters=['hs']		
	wam_stations=dict.fromkeys(wam_names)
	wam_stations_intp=dict.fromkeys(wam_names)
	
	for stp_name in wam_stations.keys():
		wam_stations[stp_name]=dict.fromkeys(keysOBS)
		wam_stations_intp[stp_name]=dict.fromkeys(keysOBS)
	
		#build dictionary
		for key in (keysOBS):
			print('WAM: loading '+key)
			if bouys[key]!=None:
				lonq=float(bouys[key].ds.StartLongitude.split()[0].replace(',','.'))
				latq=float(bouys[key].ds.StartLatitude.split()[0].replace(',','.'))
				date_slice=tuple(bouys[key].dates[[0,-1]])
				wam_stations[stp_name][key]=dict.fromkeys(['dates']+parameters)
				ti,hsi=wam.get_hs(lonq,latq,date_slice)
				wam_stations[stp_name][key]['dates']=ti
				wam_stations[stp_name][key]['hs']=hsi
			else:
				wam_stations[stp_name][key]=dict.fromkeys(['dates']+parameters)
				wam_stations[stp_name][key]['dates']=np.ma.masked_array(np.nan,mask=True)
				wam_stations[stp_name][key]['hs']=np.ma.masked_array(np.nan,mask=True)

			
#######  Wave Watch 3   ######
if addWW3:
	parameters=['hs']		
	ww3PointList='points.list'
	ww3Bouys_specInt=dict.fromkeys(ww3_names)
	ww3_stations_intp=dict.fromkeys(ww3_names)

	## grid load
	#for name,diri,Pout in zip (ww3_names,ww3dirs,ww3PointOutput,ww3GridOutput):	
	#	ww3.open_grid_output(Pout)
	#
	#	plt.figure()
	#	plt.clf()
	#	ww3.plotAtnodes(ww3.gridds['hs'][5,:])
	#	
	#	ww3.query_station_lon_lat(key)
	#for key,nr in zip(keysWW3,nrs):
	#	lon,lat=ww3.query_station_lon_lat(key)[:2]
	#	inode=ww3.find_nearest_node(lon,lat)
	# point load

	count=0
	for name,diri,Pout in zip (ww3_names,ww3dirs,ww3PointOutput):	
		ww3Bouys_specInt[name]=dict.fromkeys(keysWW3)
		ww3_stations_intp[name]=dict.fromkeys(keysWW3)
		ww3=WW3_mesh(diri+ww3meshs[count])
		ww3.load_points(diri+ww3PointList)
		ww3.open_point_output(Pout)
		nrs=[ww3.pnames.index(key) for key in keysWW3]

		count+=1
		#build dictionary
		for key,nr in zip(keysWW3,nrs):
			ww3Bouys_specInt[name][key]=dict.fromkeys(['dates']+parameters)
			ww3_stations_intp[name][key]=dict.fromkeys(['dates']+parameters)#dict.fromkeys(ww3Bouys_specInt.keys())	
			ww3Bouys_specInt[name][key]['hs']=ww3.Hs[:,nr]
			ww3Bouys_specInt[name][key]['dates']=ww3.tspec
# sometimes formatting errors 


def load_wwm_site(file):
	""" line by line reading only first to columns to circumvent issues """
	with open(file) as f:
		line=f.readline()
		header=line.split()
		#data=dict.fromkeys(header)
		data={h:[] for h in header}#dict.fromkeys(header)
		nheader=len(header)
		
		for line in f.readlines():
			line=line.replace('***************',' nan')
			ar=line.split()[:nheader]
			for key,val in zip(header,ar):
				data[key].append(val)
		for key in header[1:]:
			data[key]=np.asarray(data[key],float)
			
		import pandas as pd
		try:
			data['dates']=pd.to_datetime(data['TIME'], format='%Y%m%d.%H%M%S')		
		except:
			from IPython import embed; embed()
		#min=(data['TIME'])*100
		#min=np.asarray(np.fix(100*(min-np.fix(min))),int)
		#hour=(data['TIME'])
		#hour=np.asarray(np.round(100*(hour-np.fix(hour))),int)
		#day=(data['TIME']/100)
		#day=np.asarray(np.fix(100*(day-np.fix(day))),int)
		#mon=(data['TIME']/10000)
		#mon=np.asarray(np.fix(100*(mon-np.fix(mon))),int)
		#year=np.asarray(np.fix(data['TIME']/10000),int)
		#data['dates']=np.asarray(dt.datetime(year[0],mon[0],day[0],hour[0],0,0)+np.arange(len(year))*dt.timedelta(hours=(np.diff(min[:2])/60+np.diff(hour[:2]))[0]),np.datetime64)

		return data	

		
def load_wwm_site_multiple_files(files,write_joint_file=False):
	""" line by line reading only first to columns to circumvent issues """

	
	for i,file in enumerate(files):
		with open(file) as f:
			line=f.readline()
			
			if i==0:
				header=line.split()
				#data=dict.fromkeys(header)
				data={h:[] for h in header}#dict.fromkeys(header)
				nheader=len(header)
			
			for line in f.readlines():
				line=line.replace('***************',' nan')
				ar=line.split()[:nheader]
				for key,val in zip(header,ar):
					data[key].append(val)


    #uniqute time steps (remove double times due to restarts)						
	unqindex=np.unique(data['TIME'],return_index=True)[1]
	for key in header:
		data[key]=np.asarray(data[key],float)[unqindex]


	if write_joint_file:
		outname=files[0].split('/')[-1].replace('.','_joint.')		
		#header='\t TIME	\t HS \t TM01 \t TM02'
		#m=np.vstack((data['TIME'],data['HS'],data['tm01'],data['tm02'])).T
		m=np.vstack((data[key] for key in data.keys() )).T
		np.savetxt(outname,m,fmt='%.6f \t '*m.shape[1],header=header)

	import pandas as pd
	data['dates']=pd.to_datetime(data['TIME'], format='%Y%m%d.%H%M%S')	
	
	return data	
			

	
wwm_stations=dict.fromkeys(wwm_names)	
for dir,name in zip(WWMdirs,wwm_names):
	os.chdir(dir)	
	print(dir)
	##### work around station output formatting issues
	# error * in files.
	files=glob('*.site')
	sites={file.split('.')[0]:None for file in files}
	for file in files:
		with open(file) as file_in:
			if not os.path.exists('prep_'+file) and (not 'prep' in file):
				with open('prep_'+file,'w') as file_out:
					for nr,line in enumerate(file_in.readlines()):
						#if nr==99:
						#	break
						#line=line.replace('*','0').replace('Infinity','inf').replace('NaN            NaN','NaN')#.replace('NaN\tNaN','NaN') #.replace('NaN','-99')
						line=line.replace('*','0').replace('Infinity','').replace('NaN            NaN','NaN')#.replace('NaN\tNaN','NaN') #.replace('NaN','-99')
						file_out.write(line)
	files=np.sort(glob('prep*.site'))
	#files=np.sort(glob('*.site'))
	sites={file.split('.')[0].replace('prep_',''):None for file in files}
	for site,file in zip(sites.keys(),files):
		print('WWM: loading '+site)
		sites[site]=load_wwm_site(file)
	wwm_stations[name]=sites
	os.chdir('../')
#####################



# wam buoys
#n,m=int(np.ceil(np.sqrt(len(wam_stations)))),int(np.floor(np.sqrt(len(wam_stations))))
#wam_stations_intp=dict.fromkeys(wam_names)
#for i,name in enumerate(wam_names):
#	wam_stations_intp[name]=dict.fromkeys(wam_stations[name].keys())
#					

n,m=int(np.ceil(np.sqrt(len(wwm_stations)))),int(np.floor(np.sqrt(len(wwm_stations))))
wwm_stations_intp=dict.fromkeys(wwm_names)
for i,name in enumerate(wwm_names):
	wwm_stations_intp[name]=dict.fromkeys(wwm_stations[name].keys())


	

################## Plot ##############################################
from cycler import cycler
import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = cycler(color='rcmykbg')

varnamesBSH=['Hs']#,'tm2']
scatter={var:[] for var in varnamesBSH}
scatter_ww3={var:[] for var in varnamesBSH}
scatter_wam={var:[] for var in varnamesBSH}

ts={var:[] for var in varnamesBSH}
taylor={var:[] for var in varnamesBSH}

#hs
wwm_bias={var:[] for var in varnamesBSH}
wwm_rmse={var:[] for var in varnamesBSH}
wwm_corr={var:[] for var in varnamesBSH}


#hs
wam_bias={var:[] for var in varnamesBSH}
wam_rmse={var:[] for var in varnamesBSH}
wam_corr={var:[] for var in varnamesBSH}

## Plot map of bouys
m=np.asarray(bouy_coords)

plt.figure()
proj=projection=ccrs.PlateCarree()
ax = plt.axes(projection=proj)    
ax.set_extent((m[:,0].min()-.5,m[:,0].max()+.5,m[:,1].min()-.5,m[:,1].max()+.5))
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',edgecolor='face',facecolor=cfeature.COLORS['land'])
ax.add_feature(land_10m,zorder=-2)
for i,coord in enumerate(bouy_coords):
	ax.plot(coord[0],coord[1],'ko',transform=proj)
	ax.text(coord[0],coord[1],str(i+1),transform=proj)
fname=image_output_dir+'BouyLocs.png'
plt.savefig(fname,dpi=300)



for varnameBSH,varnameWWM,varnameWW3,unit in zip(['Hs','tm2'], ['HS','tm02'],['hs','tm2'],['[m]','[h]']):
#for varnameBSH,varnameWWM,varnameWW3,unit in zip(['Hs',], ['HS'],['hs'],['[m]']):
	# bouy time series
	for key1,key2 in zip(keysOBS,keysWWM):
		key1,key2
		
		print(key2)
		
		#if bouys[key1] == None:
		#	continue
		B=bouys[key1]
		if bouys[key1] != None:
			dataB=B.get_parameter(varnameBSH)
			date0=B.dates[0]
			date1=B.dates[-1]

			if addWam:
				for i,name in enumerate(wam_names):
					wami=wam_stations[name][key1]
					date0=np.maximum(date0,wami['dates'][0])
					date1=np.minimum(date1,wami['dates'][-1])
				date0=np.maximum(date0,wami['dates'][0]+skip_n_first_days) #ofset days
				if maxdate != []: # set manual maximum date
					date1=np.minimum(date1,maxdate)
					
			if addWWM:
				for i,name in enumerate(wwm_names):
					WWM=wwm_stations[name][key2]
					date0=np.maximum(date0,WWM['dates'][0])
					date1=np.minimum(date1,WWM['dates'][-1])
				date0=np.maximum(date0,WWM['dates'][0]+skip_n_first_days) #ofset days
				if maxdate != []: # set manual maximum date
					date1=np.minimum(date1,maxdate)

		#wam
		if addWam:
			for i,name in enumerate(wam_names):
				wami=wam_stations[name][key1]
				wam_stations_intp[name][key1]=dict.fromkeys(['dates',varnameWW3,'obs'+varnameWW3])

				if bouys[key1]!=None:
					tdata,param_intp,dataWamintp=interp_to_data_time(B.dates,dataB,wami['dates'],wami[varnameWW3],min_maxdate=[date0,date1])
				else:
					tdata=param_intp=dataWamintp=np.ma.masked_array(np.nan,mask=True)
					
				wam_stations_intp[name][key1]={'dates':tdata,varnameWW3:dataWamintp,'obs'+varnameWW3:param_intp}
				#plt.subplot(m,n,i+1)
				plt.clf()
				if len(dataWamintp)>0:
					QQplot(param_intp,dataWamintp,stationname=key1,obsname='bouy',modname=name,parameter=varnameBSH,unit=unit)

				plt.suptitle(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])	
				plt.tight_layout()
				fname=image_output_dir+'wam_'+varnameBSH+key1+'scatter_'+name+'.png'
				plt.savefig(fname,dpi=300)
				scatter_wam[varnameBSH].append(fname)

			
		if addWW3 & (varnameBSH =='Hs'):
			key3=keysWW3[keysOBS.index(key1)]	
			WW3s=[]
			for name in ww3_names:
				WW3s.append(ww3Bouys_specInt[name][key3])
				
		if addWWM:			
			for i,name in enumerate(wwm_names):
				WWM=wwm_stations[name][key2]
				wwm_stations_intp[name][key2]=dict.fromkeys(['dates',varnameWWM,'obs'+varnameWWM])

				tdata,param_intp,dataWWMintp=interp_to_data_time(B.dates,dataB,WWM['dates'],WWM[varnameWWM],min_maxdate=[date0,date1])

				wwm_stations_intp[name][key2]={'dates':tdata,varnameWWM:dataWWMintp,'obs'+varnameWWM:param_intp}
				#plt.subplot(m,n,i+1)
				plt.clf()
				if len(dataWWMintp)>0:
					QQplot(param_intp,dataWWMintp,stationname=key1,obsname='bouy',modname=name,parameter=varnameBSH,unit=unit)

				plt.suptitle(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])	
				plt.tight_layout()
				fname=image_output_dir+'wwm_'+varnameBSH+key1+'scatter_'+name+'.png'
				plt.savefig(fname,dpi=300)
				scatter[varnameBSH].append(fname)
			
		if addWW3 & (varnameBSH =='Hs'): # currently not working for HS
			for WW3,name in zip(WW3s,ww3_names):
				tdata,param_intp,dataWW3intp=interp_to_data_time(B.dates,dataB,WW3['dates'],WW3['hs'],min_maxdate=[date0,date1])
				ww3_stations_intp[name][key3]={'dates':tdata,'hs':dataWW3intp,'obs'+varnameWWM:param_intp}
				plt.clf()	
				if len(dataWW3intp)>0:
					QQplot(param_intp,dataWW3intp,stationname=key1,obsname='bouy',modname=name,parameter='Hs',unit='[m]')
				plt.suptitle(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])	
				plt.tight_layout()
				fname=image_output_dir+varnameBSH+key1+'scatter_'+name+'.png'
				scatter_ww3[varnameBSH].append(fname)
				plt.savefig(fname,dpi=300)

		# time series
		#if bouys[key1] != None:
		plt.clf()
		plt.plot(B.dates,dataB,'g.',  markersize=4,color='gray')


		if addWWM:
			for i,name in enumerate(wwm_names):
				WWM=wwm_stations[name][key2]
				plt.plot(WWM['dates'],WWM[varnameWWM],'-')
		if addWam:
			for i,name in enumerate(wam_names):
				wami=wam_stations[name][key]
				plt.plot(wami['dates'],wami[varnameWW3],'-')
			
		plt.xlim((date0,date1))
		plt.ylabel(varnameWWM+unit)
		plt.gcf().autofmt_xdate()
		plt.title(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])
		if addWW3 & (varnameBSH =='Hs'):
			for WW3,name in zip(WW3s,ww3_names):
				plt.plot(WW3['dates'],WW3['hs'],'--')
			plt.legend(['buoy']+display_names+ww3_names,loc='best',frameon=False)		
		else:	
			plt.legend(['buoy']+display_names,loc='best',frameon=False)			
		plt.grid()
		plt.tight_layout()
		fname=image_output_dir+'wwm_'+varnameBSH+key1+'timeseries_'+'_'.join(display_names)+'.png'
		ts[varnameBSH].append(fname)
		plt.savefig(fname,dpi=300)
	## Collect samplews per run for taylor diag	
	samples=[]	
	
	# interp
	
	if addWWM:
		for i,name in enumerate(wwm_names):	
			corr=np.asarray([np.corrcoef(wwm_stations_intp[name][key][varnameWWM],wwm_stations_intp[name][key]['obs'+varnameWWM])[0,1] for key in wwm_stations_intp[name].keys()])

			std_rel=np.asarray([np.std(wwm_stations_intp[name][key][varnameWWM])/np.std(wwm_stations_intp[name][key]['obs'+varnameWWM]) for key in wwm_stations_intp[name].keys()])

			rmse=np.asarray([np.sqrt(np.mean((wwm_stations_intp[name][key][varnameWWM]-wwm_stations_intp[name][key]['obs'+varnameWWM])**2)) for key in wwm_stations_intp[name].keys()])

			samples.append(list(zip(std_rel,corr,keysOBS)))

			bias=np.asarray([np.mean((wwm_stations_intp[name][key][varnameWWM]-wwm_stations_intp[name][key]['obs'+varnameWWM])) for key in wwm_stations_intp[name].keys()])

			wwm_bias[varnameBSH].append(bias)
			wwm_rmse[varnameBSH].append(rmse)
			wwm_corr[varnameBSH].append(corr)
			
			#display_names+=wwm_names

	if addWam:
		for i,name in enumerate(wam_names):	
			corr=np.asarray([np.corrcoef(wam_stations_intp[name][key][varnameWW3],wam_stations_intp[name][key]['obs'+varnameWW3])[0,1] for key in wam_stations_intp[name].keys()])

			std_rel=np.asarray([np.std(wam_stations_intp[name][key][varnameWW3])/np.std(wam_stations_intp[name][key]['obs'+varnameWW3]) for key in wam_stations_intp[name].keys()])

			rmse=np.asarray([np.sqrt(np.mean((wam_stations_intp[name][key][varnameWW3]-wam_stations_intp[name][key]['obs'+varnameWW3])**2)) for key in wam_stations_intp[name].keys()])

			samples.append(list(zip(std_rel,corr,keysOBS)))

			bias=np.asarray([np.mean((wam_stations_intp[name][key][varnameWW3]-wam_stations_intp[name][key]['obs'+varnameWW3])) for key in wam_stations_intp[name].keys()])

			wam_bias[varnameBSH].append(bias)
			wam_rmse[varnameBSH].append(rmse)
			wam_corr[varnameBSH].append(corr)			
			
	if addWW3 & (varnameBSH =='Hs'):

		ww3_bias={var:[] for var in varnamesBSH}
		ww3_rmse={var:[] for var in varnamesBSH}
		ww3_corr={var:[] for var in varnamesBSH}
		for i,name in enumerate(ww3_names):	
			corr=np.asarray([np.corrcoef(ww3_stations_intp[name][key]['hs'],ww3_stations_intp[name][key]['obs'+varnameWWM])[0,1] for key in ww3_stations_intp[name].keys()])

			std_rel=np.asarray([np.std(ww3_stations_intp[name][key]['hs'])/np.std(ww3_stations_intp[name][key]['obs'+varnameWWM]) for key in ww3_stations_intp[name].keys()])

			rmse=np.asarray([np.sqrt(np.mean((ww3_stations_intp[name][key]['hs']-ww3_stations_intp[name][key]['obs'+varnameWWM])**2)) for key in ww3_stations_intp[name].keys()])

			samples.append(list(zip(std_rel,corr,keysOBS)))

			bias=np.asarray([np.mean((ww3_stations_intp[name][key]['hs']-ww3_stations_intp[name][key]['obs'+varnameWWM])) for key in ww3_stations_intp[name].keys()])

			ww3_bias[varnameBSH].append(bias)
			ww3_rmse[varnameBSH].append(rmse)
			ww3_corr[varnameBSH].append(corr)

		display_names+=ww3_names

	from matplotlib.text import TextPath
	label = TextPath((0,0), str(1))#, linewidth=1)
	plt.close('all')
	dia=plotTaylor(samples[0],stdref=1,extend=False) #negative
	colors=plt.cm.tab10(range(len(samples)))  #['b','r']
	#Add models to Taylor diagram
	phs=[plt.plot(np.nan,np.nan,color='k',marker=label,linestyle='',markersize=14)[0]]
	for nr,sample in enumerate(samples[1:]):
		for i,(stddev, corrcoef, name) in enumerate(sample):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
		phs+=plt.plot(np.nan,np.nan,marker=label,color=colors[nr,:],linestyle='',markersize=14)
	plt.legend(phs,display_names)
	#plt.tight_layout()
	fname=image_output_dir+'wwm_'+varnameBSH+key1+'taylour_'+'_'.join(display_names)+'.png'
	taylor[varnameBSH].append(fname)
	plt.savefig(fname,dpi=300)
	plt.close()

	if addWW3 & (varnameBSH =='Hs'):
		plt.close('all')
		dia=plotTaylor(samples[0],stdref=1,extend=False) #negative
		colors=plt.cm.tab10(range(len(samples)))  #['b','r']
		#Add models to Taylor diagram
		phs=[plt.plot(np.nan,np.nan,color='k',marker=label,linestyle='',markersize=14)[0]]
		for nr,sample in enumerate(samples[1:len(wwm_names)]):
			for i,(stddev, corrcoef, name) in enumerate(sample):
				i,stddev, corrcoef, name
				dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
			phs+=plt.plot(np.nan,np.nan,marker=label,color=colors[nr,:],linestyle='',markersize=14)
		plt.legend(phs,wwm_names)
		#plt.tight_layout()
		fname=image_output_dir+varnameBSH+'wwm_only_'+key1+'taylour_'+'_'.join(display_names)+'.png'
		taylor[varnameBSH].append(fname)
		plt.savefig(fname,dpi=300)
		plt.close()

		plt.close('all')
		dia=plotTaylor(samples[len(wwm_names)],stdref=1,extend=False) #negative
		#Add models to Taylor diagram
		phs=[plt.plot(np.nan,np.nan,color='k',marker=label,linestyle='',markersize=14)[0]]
		for nr,sample in enumerate(samples[len(wwm_names)+1:]):
			for i,(stddev, corrcoef, name) in enumerate(sample):
				i,stddev, corrcoef, name
				dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
			phs+=plt.plot(np.nan,np.nan,marker=label,color=colors[nr,:],linestyle='',markersize=14)
		plt.legend(phs,ww3_names)
		#plt.tight_layout()
		fname=image_output_dir+varnameBSH+'ww3_only_'+key1+'taylour_'+'_'.join(display_names)+'.png'
		taylor[varnameBSH].append(fname)
		plt.savefig(fname,dpi=300)
		plt.close()
		

		
####### Report #############################


###### make report #######		
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
import time

chart_style = TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                          ('VALIGN', (0, 0), (-1, -1), 'CENTER')])
chart_style = TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTER')])

doc = SimpleDocTemplate(image_output_dir+pdfname,pagesize=letter,
                        rightMargin=72,leftMargin=72,
                        topMargin=72,bottomMargin=18,title='Blubb')
doc.title="WWM Validation"
					

time.ctime()
title = 'Validation of WAM run against BSH bouys'							
formatted_time = time.ctime()
full_name = "Benjamin Jacob"
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
Story=[]
Story.append(Paragraph(title, styles["title"]))
Story.append(Spacer(1, 12))
Story.append(Spacer(1, 12))
####

ptext = '%s' % formatted_time
Story.append(Paragraph(ptext, styles["Normal"]))

# Create return address
ptext = '%s' % full_name
Story.append(Paragraph(ptext, styles["Normal"]))    
Story.append(Spacer(4, 12))

ptext='Validation of WWM '
if len(ww3_names) >0:
	ptext+='and WW3 '
if len(wam_names) >0:
	ptext+='and WAM '
	
ptext+='run(s) for period ' +  str(date0)[:10] + ' - '+ str(date1)[:11] +'.'

Story.append(Paragraph(ptext, styles["Normal"]))    
if len(wwm_names) >0:
	ptext='WWM scenarios are: ' 
	Story.append(Paragraph(ptext, styles["Normal"]))    

	for name,desc in zip(wwm_names,wwm_descriptions):
		ptext=' '+name + desc +', '
		Story.append(Spacer(0.5, 5))
		Story.append(Paragraph(ptext, styles["Normal"]))    

if len(wam_names) >0:
	ptext='WAM scenarios are: ' 
	Story.append(Paragraph(ptext, styles["Normal"]))    

	for name,desc in zip(wam_names,wam_descriptions):
		ptext=' '+name + desc +', '
		Story.append(Spacer(0.5, 5))
		Story.append(Paragraph(ptext, styles["Normal"]))    

#ptext=ptext[:-1]	
if len(ww3_names) >0:
	Story.append(Spacer(1, 12))
	ptext='WW3 scenarios are ' 
	Story.append(Paragraph(ptext, styles["Normal"]))    
	for name,desc in zip(ww3_names,ww3_descriptions):
		ptext=' '+name + desc
		Story.append(Spacer(0.5, 5))
		Story.append(Paragraph(ptext, styles["Normal"]))    

ptext=comments
Story.append(Paragraph(ptext, styles["Normal"]))    
counter=1

for varnameBSH in varnamesBSH:

	Story.append(Spacer(4, 12))
	Story.append(Paragraph('Variable {:s}'.format(varnameBSH), styles["h1"]))

	Story.append(Spacer(4, 12))
	Story.append(Paragraph('Overview', styles["h1"]))
	
	ims=[]
	im = Image(image_output_dir+'BouyLocs.png',4*inch,4*inch)
	ims.append(im)
	im = Image(taylor[varnameBSH][0],4*inch,4*inch)
	ims.append(im)
	Story.append(Table([ims],colWidths=[4 * inch, 4 * inch],rowHeights=[4 * inch], style=chart_style))	
	#Story.append(Table([ims],colWidths=[inchW * inch, inchW * inch],				rowHeights=[inchH * inch], style=chart_style))	
	#Story.append(im)
	#Story.append(im)
	Story.append(Paragraph('Fig {:d}: taylor diagram (all runs) of significant wave height performance against BSH buoys for period {:s} - {:s} -All samples.'.format(counter,str(date0),str(date1)), styles["Normal"]))
	Story.append(Spacer(2, 12))

	if addWW3 & (varnameBSH =='Hs'):
		counter+=1
		Story.append(Spacer(4, 12))
		im = Image(taylor[varnameBSH][1],4.5*inch,4*inch)
		Story.append(im)
		Story.append(Paragraph('Fig {:d}: taylor diagram (WWM only) of significant wave height performance against BSH buoys for period {:s} - {:s}.'.format(counter,str(date0),str(date1)), styles["Normal"]))
		Story.append(Spacer(2, 12))

		counter+=1
		Story.append(Spacer(4, 12))
		im = Image(taylor[varnameBSH][2],4.5*inch,4*inch)
		Story.append(im)
		Story.append(Paragraph('Fig {:d}: taylor diagram (WW3 only) of significant wave height performance against BSH buoys for period {:s} - {:s}.'.format(counter,str(date0),str(date1)), styles["Normal"]))
		Story.append(Spacer(2, 12))

	ptext='The station averaged values (each station weighted equally independent of nr of observtions) correspond to simulation: bias/rmse/correlation '
	Story.append(Paragraph(ptext, styles["Normal"]))    
	Story.append(Spacer(1, 12))

	
	if (len(wam_names) >0) & addWam:
		for i,name in enumerate(wam_names):
			Story.append(Spacer(1, 12))
			ptext='\n{:s}: {:f}/{:f}/{:f}, '.format(name,wam_bias[varnameBSH][i].mean(),wam_rmse[varnameBSH][i].mean(),wam_corr[varnameBSH][i].mean())
			Story.append(Paragraph(ptext, styles["Normal"]))    
	
	if (len(wwm_names) >0) & addWWM:
		for i,name in enumerate(wwm_names):
			Story.append(Spacer(1, 12))
			ptext='\n{:s}: {:f}/{:f}/{:f}, '.format(name,wwm_bias[varnameBSH][i].mean(),wwm_rmse[varnameBSH][i].mean(),wwm_corr[varnameBSH][i].mean())
			Story.append(Paragraph(ptext, styles["Normal"]))    

		Story.append(Paragraph(ptext, styles["Normal"]))    
	if (len(ww3_names) >0) & addWW3 & (varnameBSH =='Hs'):	
		for i,name in enumerate(ww3_names):
			Story.append(Spacer(1, 12))
			ptext='\n{:s}: {:f}/{:f}/{:f}, '.format(name,ww3_bias[varnameBSH][i].mean(),ww3_rmse[varnameBSH][i].mean(),ww3_corr[varnameBSH][i].mean())
			Story.append(Paragraph(ptext, styles["Normal"]))    	

	Story.append(Spacer(4, 12))
	Story.append(Paragraph('time series', styles["h1"]))

	for i,fig in enumerate(ts[varnameBSH]):
		Story.append(Spacer(4, 12))
		#print('appending '+ts[varnameBSH][i])
		im = Image(ts[varnameBSH][i],4*inch,3*inch)
		Story.append(im)
		counter+=1
		Story.append(Paragraph('Fig {:d}: Time Series of {:s} at station {:s}.'.format(counter,varnameBSH,str(keysWW3[i])), styles["Normal"]))

		
		
	Story.append(Spacer(4, 12))
	Story.append(Paragraph('scatter plots', styles["h1"]))

	# 1 row per station
	nwam=len(wam_names)
	nwwm=len(wwm_names)*addWWM
	nww3=len(ww3_names)*addWW3

	nruns=len(display_names)
	sqrt=np.sqrt(nruns)	

	ncols=np.minimum(nruns,4)
	nrows=int(np.ceil(nruns/ncols))

	inchW=8.2/ncols*1.1
	inchH=6/ncols*1.1

	for i in range(len(keysWW3)):
		Story.append(Spacer(4, 12))
		counter+=1	
		ims=[]
		for k in range(i*(nwwm),(i+1)*nwwm):
			ims.append(Image(scatter[varnameBSH][k],inchW*inch,inchH*inch))
		
		if addWam & (varnameBSH =='Hs'):
			for k in range(i*(nwam),(i+1)*nwam):
				ims.append(Image(scatter_wam[varnameBSH][k],inchW*inch,inchH*inch))

		if addWW3 & (varnameBSH =='Hs'):
			for k in range((i)*nww3,(i+1)*nww3):
				ims.append(Image(scatter_ww3[varnameBSH][k],inchW*inch,inchH*inch))

		for irow in range(nrows):
			i0=irow*ncols
			i1=i0+ncols	
			i1=np.minimum(i1,len(ims))
			Story.append(Table([ims[i0:i1]],
				colWidths=[inchW * inch, inchW * inch],
				rowHeights=[inchH * inch], style=chart_style))				

		Story.append(Spacer(2, 0))			
		Story.append(Paragraph('Fig {:d}: Scatter plot of {:s} from top left to bottom right (counting along rows) simulated by '.format(counter,varnameBSH) +' '.join(display_names) +' against BSH buoys at station {:s}.'.format(str(keysWW3[i])), styles["Normal"]))
doc.build(Story)				