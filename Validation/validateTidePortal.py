import os
import netCDF4
import sys
import csv
import matplotlib
matplotlib.use('Agg') # backend
import datetime as dt
import glob
from scipy.spatial import cKDTree
import datetime as dt
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
# own and 3d party libraries
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
from schism import * # import schism functions
from techit import * # import latex script
from TaylorDiagram import * # import taylordiagram
from data_and_model_classes import cmems
from pyproj import Proj

# validates SCHISM vs tide gauge stations (for 1 year) from EMOD net 
# stored by SG in /work/gg0028/g260099/OSTIA/VALIDATION/VAL_TIDE/EMOD_TIDEG_DATA/
# as monthly files.
# 1 - search data in subfolders by name, which is in specified <year>
# 2 - second sub select tide gauge stations within <dtol> distance from
# closest schism grid nodes. 
# 3 - plot 
# 4 - derive skill parameter (bias,rmse,correlation) after interpolation of model data
# to ovebservational time steps - the temporal resolution is reduced to the model time step
# Comparison is done with zcor/elev data from SCHISMs closest grid point
# Uses schism class stored in /pf/g/g260114/Programs/python/scripts/schism.py
#
# call after 
# module load python3
# python validateTide.ps
# if mistralpp busy, performance better when called from run_sript
#
# Output:
#  pdf with validation images
# 'errorstats.csv' - csv file with error statistics
#
########## settings #################################
# directories (have to end with '/')
#tgdir='/work/gg0028/g260099/OSTIA/VALIDATION/VAL_TIDE/EMOD_TIDEG_DATA/EMOD_DATA_201903/' # Direcotry of 
tgdir='/work/gg0028/SCHISM/validation/tidegauge/validate_tideportal/'
tg_utc=+1  # utc reference of tide gauge data for correction
#oceandir='/work/gg0028/g260099/AMM15/2017/'
#sshfile='/work/gg0028/g260099/AMM15/2017/metoffice_foam1_amm15_NWS_SSH_b20171228_hi20171227.nc' #ssh file necessary to intialize amm data axxes
#pattern='SSH'												 # pattern reocurring in al ssh files	
															   # contains data somewhat between start 2012 and end 2018
setupdir='//work/gg0028/g260114/RUNS/GermanBight/GermanBight/'  # schism run directory, containts hgrid.* etc.
ncdir=setupdir+'combined2/' 								   # directory of schism nc output
year=2017													   # year to be analyesed 	 
dtol=0.01           										   # distance tolerance in degree lon/lat 
															   # between SCHISM closest grid node and tide gauge station
p = Proj(proj="utm",zone=32,ellps="WGS84",south=False)         # set Projection according to textfiles from
															   # tide portal	
															   
outdir='/work/gg0028/SCHISM/validation/tidegauge/GB_tideportal/'	   # output directory where images will be stored
if not os.path.exists(outdir): os.mkdir(outdir) 

remove_mean=True  # remove temporal mean from Data and Model to compare 

#--- what to plot True/False:			
overview_map=True												
satistic_maps=True
full_timeseries=False # for flull year data hard to see anything
first_two_Weeks=True # zomed period 			
monthly_subplots=True
taylor_diag=True
put_pics_to_texdoc=True    										# images will be put in tex document
latexname='GB_amm15_valid.tex'										# in alphanumerical order 
latextitle='GB 2017'
latextext='Tide Gauge validation of Europe ' + str(year) +'period of 10 year run against Emodnet data. \
\n Ocean Forcing CMEMS dayli mean + Tide \n Atmosphere forcing: Era5. Data is Blue. Model is red'
if remove_mean:
		latextext+=' Time Series were mean removed and so is rmse then. In Taylor it is always mean removed'
plt.ioff() # plt.ion() has to be off to work in background mode
exclude_by_names=['dummy1','dummy2'] # list of Station names to remove from analysis
				     # to exclude bad station in second iterarion (names have to match those in the overviewmap)



Rcircle=150 # radius of colorcoded cricles used in scatter maps

limit_to_data=True   # limit map plots to bounding box of data. If data only covers sub domain of area (Europe case)					 
#############################################################


######### load SCHISM setup   ##################################
cwd=os.getcwd()
os.chdir(setupdir)
s=schism_setup()
s.nntree = cKDTree(list(zip(s.lon,s.lat))) 


# initiate file access
schismfiles=[] 
for iorder in range(6): # check for schout_nc files until 99999
	schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
nrs=list(np.asarray(nrs)[np.argsort(nrs)])


s.nc=Dataset(schismfiles[0])
# set reference time
reftime=dt.datetime.strptime(s.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')
lons=np.asarray(s.lon)
lats=np.asarray(s.lat)
modeldt=np.diff(s.nc['time'][:2])[0]				# use model time step as resolution
nt=len(s.nc['time'])
daysPerStack=modeldt*nt/86400.0
print('done loading schism setup')

# restrict files to year
stack_startdates=reftime+dt.timedelta(seconds=modeldt) + dt.timedelta(seconds=(daysPerStack*86400.0))*np.arange(len(schismfiles))     
ilast=(stack_startdates<=dt.datetime(year+1,1,1,0,0,0)).sum()
ifirst=(stack_startdates<=dt.datetime(year,1,1,0,0,0)).sum()
schismfiles=schismfiles[ifirst:ilast]

################################################################################################


######### identify Gauge Data ##########################################
#folders=glob.glob(tgdir+'/*/')
## check if data for 2017
#print('searching coordinates for available data in year {:d}'.format(year))
#coords=[]
#folders2=[]
#for folder in folders:
#	files=glob.glob(folder+'*{:d}*'.format(year))
#	if len(files)==0:
#		#folders.remove(folder)
#		pass
#	else: # check geographical coordinates and existenc of SLEV field
#		nc=Dataset(files[0])
#		if 'SLEV' in nc.variables.keys():
#			coords.append((nc['LONGITUDE'][0],nc['LATITUDE'][0]))
#			folders2.append(folder)
#		nc.close()
######### identify Gauge Data ##########################################
files=glob.glob(tgdir+'/*.txt')
# check if data for 2017
print('searching coordinates for available data in year {:d}'.format(year))
coords=[]
folders2=[]

#def get_from_tideporal_txt(file):
#	station={}
#	station['date']=[]
#	station['zeta']=[]
#	start=False
#	with open(file) as f:
#		for i,line in enumerate(f.readlines()):
#				if (line[0]=='=') & (start == False):
#					start=True
#				if 'Datum' in line:
#					continue
#				elif (line[0]=='=') & (i > 30):
#					break
#
#				if start==False:
#					station[line.split(':')[0].replace(' ','_')]=[line.split(':')[1].split('\n')[0]]
#				elif line[0] != '=':
#					try:
#						a,b=line.split('\t')[:2]
#						station['date'].append( dt.datetime.strptime(a,'%Y-%m-%d %H:%M:%S'))
#						station['zeta'].append(np.float(b))
#					except:	
#						pass
#	station['date']=np.asarray(station['date'])
#	station['zeta']=np.asarray(station['zeta'])
#	station['lon'],station['lat']=p(np.float(station['East________________________________'][0]), np.float(station['North_______________________________'][0]), inverse=True)	
#	return station
def get_from_tideporal_txt(file,utc=0):
	station={}
	station['date']=[]
	station['zeta']=[]
	start=False
	with open(file) as f:
		for i,line in enumerate(f.readlines()):
				if (line[0]=='=') & (start == False):
					start=True
				if 'Datum' in line:
					continue
				elif (line[0]=='=') & (i > 30):
					break

				if start==False:
					station[line.split(':')[0].replace(' ','_')]=[line.split(':')[1].split('\n')[0]]
				elif line[0] != '=':
					a,b=line.split('\t')[:2]
					if len(b) == 0:
						b='nan'
					station['date'].append( dt.datetime.strptime(a,'%Y-%m-%d %H:%M:%S'))
					station['zeta'].append(np.float(b))
					
	station['date']=np.asarray(station['date'])+dt.timedelta(hours=-utc)
	station['zeta']=np.asarray(station['zeta'])
	station['lon'],station['lat']=p(np.float(station['East________________________________'][0]), np.float(station['North_______________________________'][0]), inverse=True)	
	return station
	
coords=[]
stations={}
names0=[]
for file in	files:
	station=get_from_tideporal_txt(file)
	stations[station['Stationsname________________________'][0][1:]]=station
	coords.append((station['lon'],station['lat']))
	plt.plot(station['lon'],station['lat'],'ko')
	plt.text(station['lon'],station['lat'],station['Stationsname________________________'][0][1:])
	names0.append(station['Stationsname________________________'][0][1:])
names=names0	

#----  reduce to data in domain as given by dtol -----------
		
lldists,nn=s.nntree.query(coords)		
#nn=nn[lldists<=dtol]
#stations=np.asarray(folders2)[lldists<=dtol]
#coords=np.asarray(coords)[lldists<=dtol]

#print('done selecting tide gauges, found: '+str(len(coords)) + ' stations meeting criterea')
#names=[folders2[i][:-1][folders2[i][:-1].rindex('/')+1:] for i in range(len(folders2)) ]
#names=list(stations.keys()) #this changes the order
########################################################

all_names=names0.copy()
################## remove e.g bad station from list
#for name,coord in zip(names,coords):
#	if name in exclude_by_names:
#		names.remove(name)
##################

######### initialize TG Data acces ##############################################
#ncs={name:0 for name in names}
#invalid_names=[]
#for name,folder,coord in zip(names,folders2,coords):
#	#print(folder+ ' ' + name)
#	files=np.sort(glob.glob(folder+'*{:d}*'.format(year)))
#	try:
#		ncs[name]=(MFDataset(files))
#	except:
#		#print('error loading '+ name)
#		#print('deleting ' + name+ ' from list ')
#		invalid_names.append(name)
#		#names.remove(name)    # this shifts folders with respect to name
#		#ncs.pop(name, None) # delete item from netcdf
#
## remove missing and defect data		
#keep=np.asarray([name not in invalid_names for name in all_names])		
#coords=coords[keep]		
#folders2=folders2[keep]
#for name in invalid_names:
#	names.remove(name)
#	ncs.pop(name, None)
###############################################################################





###############################

xmin,ymin=np.min(coords,axis=0)
xmax,ymax=np.max(coords,axis=0)




##### Plot selected stations
if overview_map:
	s.plot_domain_boundaries()
	fig=plt.gcf()
	fig.set_size_inches(11,8,forward=True)
	plt.plot(lons[nn],lats[nn],'bo')
	for coord,name in zip(coords,names0):
		lon,lat=coord
		name
		plt.plot(lon,lat,'r+')
		plt.text(lon,lat,' '+name,rotation=50, rotation_mode='anchor')
	plt.tight_layout()
	plt.title('TG Stations')
	
	if limit_to_data:
		plt.xlim((xmin-1,xmax+1))
		plt.ylim((ymin-1,ymax+1))
	plt.savefig(outdir+'0_TideGuageStationLocations.png',dpi=300)	
	plt.close()
###################################	







	
# load data from schism next neighbours to TG statsions
print('load SCHISM sea level time series at TG next neighbouring nodes')
if 'elev' in s.nc.variables.keys():
	use_elev=True
else:
	use_elev=False # use zcor

n=len(schismfiles)	
print(str(n) + 'files')
if use_elev:	
	file=schismfiles[0]
	nc=Dataset(file)
	t2=nc['time'][:]
	zeta_schism=nc['elev'][:,:][:,nn]
	nc.close()
	for file in schismfiles[1:n]:
		print('loading file '+ file)
		nc=Dataset(file)
		zeta_schism=np.concatenate((zeta_schism,nc['elev'][:,:][:,nn]))
		t2=np.concatenate((t2,nc['time'][:]))
		nc.close()
else: # use zcor
	file=schismfiles[0]
	nc=Dataset(file)
	t2=nc['time'][:]
	zeta_schism=nc['zcor'][:,:,-1][:,nn]
	nc.close()
	for file in schismfiles[1:n]:
		print('loading file '+ file)
		nc=Dataset(file)
		zeta_schism=np.concatenate((zeta_schism,nc['zcor'][:,:,-1][:,nn]))
		t2=np.concatenate((t2,nc['time'][:]))
		nc.close()
dates=reftime+dt.timedelta(seconds=t2[0])+dt.timedelta(seconds=modeldt)*np.arange(len(t2))
#dates=reftime+dt.timedelta(seconds=t2[1]-t2[0])*np.arange(1,len(t2)+1) # convert into datetime object
																	# has to start from 1

# load tide gauge data
print('load Tide Gauge Data')
Time={name:0 for name in names}
ZETA={name:0 for name in names}
for i, name in enumerate(names):
	print(name)
	# load tide gauge
	#t=(ncs[name]['TIME'][:])
	#try:
	#	igood= (ncs[name]['SLEV_QC'][:] < 3)[:,0] # select good data (ncs[name]['TIME_QC'][:]<3) and
	#except:
	#	igood=np.ones(len(t),bool)
	#t=t[igood]
	#t0=dt.datetime.strptime(ncs[name]['TIME'].units[11:30],'%Y-%m-%dT%H:%M:%S')
	#timeunit=ncs[name]['TIME'].units[:ncs[name]['TIME'].units.index(' ')]
	#Time[name]=[t0+eval('dt.timedelta({:s}=float(ti))'.format(timeunit)) for ti in t]
	Time[name]=stations[name]['date']
	#ZETA[name]=(ncs[name]['SLEV'][:])[igood]
	ZETA[name]=stations[name]['zeta']/100
print('done loading Tide Gauge Data')

#[len(Time[name]) == len(ZETA[name]) for name in names]

########################### CMEMS data acces #####################################
# utc + 1

### remove mean
ZETAmean={name:0 for name in names}
zeta_schismmean={name:0 for name in names}
if remove_mean:
	print('removing mean')
	for i,name in enumerate(names[:]):
		ZETAmean[name]=np.nanmean(ZETA[name])
		#zeta_schismmean[name]=zeta_schism[:,i].mean()
		ZETA[name]-=ZETAmean[name]
		#zeta_schism[:,i]-=zeta_schismmean[name]


# debug		
#pickle.dump((Time,ZETA),open(outdir+"time_zetaData","wb"))
#pickle.dump((dates,zeta_schism),open(outdir+"time_zetaSCHISM","wb"))
#pickle.dump((ZETAmean,zeta_schismmean),open(outdir+"mean_zetas","wb"))

		
####### temporal interpolation of model to data  and error statistics ##########
print('interploating to common time steps of data and calculating error statistics')
timeintp={name:0 for name in names}
zetaintp={name:0 for name in names}
zetaintp_schism={name:0 for name in names}
#zetaintp_amm={name:0 for name in names}

# error stats
bias={name:0 for name in names}
rmse={name:0 for name in names}
stddata={name:0 for name in names}
stdmodel={name:0 for name in names}
R={name:0 for name in names}      

# same for cmems against tide gauge
#bias2={name:0 for name in names}
#rmse2={name:0 for name in names}
#stddata2={name:0 for name in names}
#stdmodel2={name:0 for name in names}
#R2={name:0 for name in names}    


n=len(names)
for i,name in enumerate(names):
	print(str(i)+'/'+str(n)+' '+name)
	Todates=np.asarray(Time[name])
	ilast=(Todates<=dates[-1]).sum()
	ifirst=(Todates<=dates[0]).sum()
	Todates=Todates[ifirst:ilast]

	# interpolate to model time step but at maximum the temporla resolution of the model	
	Todates2=[Todates[0]]
	for date in Todates:
		if (date-Todates2[-1]).total_seconds() >= modeldt:
			Todates2.append(date)
	
	# make simple seconds since start time frame for interpolation (not working with date objects)
	tin=[(timei-reftime).total_seconds() for timei in dates]
	tout=[(timei-reftime).total_seconds() for timei in Todates2]
	
	timeintp[name]=Todates2
	zetaintp[name]=ZETA[name][ifirst:ilast][np.asarray([datei in Todates2 for datei in Todates])]
	fintp=interp1d(tin, zeta_schism[:,i])
	zetaintp_schism[name]=fintp(tout)
	
	#tin=[(timei-reftime).total_seconds() for timei in dates2]
	#fintp=interp1d(tin, zeta_amm[:,i])
	#zetaintp_amm[name]=fintp(tout)
	
	# calculate error  stats ###########################	
	diff=(zetaintp_schism[name]-zetaintp[name])
	bias[name]=zeta_schismmean[name]-ZETAmean[name]
	rmse[name]=np.sqrt((diff**2).mean()) # mean removed
	stddata[name]=np.std(zetaintp[name])
	stdmodel[name]=np.std(zetaintp_schism[name])
	ivalid=np.isnan(zetaintp[name])==False
	R[name]=np.corrcoef(zetaintp[name][ivalid],zetaintp_schism[name][ivalid])[1]

	#diff2=(zetaintp_amm[name]-zetaintp[name])
	#bias2[name]=zeta_ammmean[name]-ZETAmean[name]
	#rmse2[name]=np.sqrt((diff2**2).mean()) # mean removed
	#stdmodel2[name]=np.std(zetaintp_amm[name])
	#R2[name]=np.corrcoef(zetaintp[name][:,0],zetaintp_amm[name])[0,1]

	
# relative score eg model vs climatology
#def BS(data): # Brier score
#	return	((data-data.mean())**2).sum()/len(data)	
#def BSS(BSdata,BSmodel): # Brier skill score	
#	return (BSdata-BSmodel)/BSdata
#BS(ZETA[name])


#mask_amm=np.asarray(list(zeta_ammmean.values()))<-999
#ivalid=mask_amm==False
#valid_names=np.asarray(list(bias2.keys()))[ivalid]

#for name,ismasked in zip(names,mask_amm):
	#name,ismasked
	#bias2[name]=np.array(bias2[name])
	#rmse2[name]=np.array(rmse2[name])
	#stdmodel2[name]=np.array(stdmodel2[name])
	#R2[name]=np.array(R[name])

	#bias2[name]=np.ma.masked_array(bias2[name],mask=(ismasked))
	#rmse2[name]=np.ma.masked_array(rmse2[name],mask=(not ismasked))
	#stdmodel2[name]=np.ma.masked_array(stdmodel2[name],mask=(ismasked))
	#R2[name]=np.ma.masked_array(R[name],mask=(ismasked))

	
	#bias2[name]=np.ma.masked_array(bias2[name],mask=(not ismasked))
	#rmse2[name]=np.ma.masked_array(rmse2[name],mask=(not ismasked))
	#stdmodel2[name]=np.ma.masked_array(stdmodel2[name],mask=(not ismasked))
	#R2[name]=np.ma.masked_array(R[name],mask=(not ismasked))

# mask 



names=np.asarray(names)

####### Plot statistics #######################
x,y=coords[:,0],coords[:,1]
# if satistic_maps:
	# labels=['bias','rmse','correlation']
	# itdata=[bias,rmse,R]
	# for data,label in zip(itdata,labels):
		# s.plot_domain_boundaries()
		# fig=plt.gcf()
		# fig.set_size_inches(11,8,forward=True)
		# plt.scatter(x,y,s=Rcircle,c=np.asarray([data[name] for name in names]))
		# ch=plt.colorbar()
		# ch.set_label(label)
		# for coord,name in zip(coords,names):
			# xi,yi=coord
			# plt.text(xi+0.01,yi+0.01,''+name[:5],rotation=50,rotation_mode='anchor')
		# plt.tight_layout()
		# if limit_to_data:
			# plt.xlim((xmin-1,xmax+1))
			# plt.ylim((ymin-1,ymax+1))
		# plt.savefig(outdir+'1_'+label+'.png',dpi=300)	
		# plt.close()
		
	#add ration of stadnard deviatons	
	# s.plot_domain_boundaries()
	# fig=plt.gcf()
	# fig.set_size_inches(11,8,forward=True)
	# plt.scatter(x,y,s=Rcircle,c=np.asarray([stdmodel[name]/stddata[name] for name in names]))
	# ch=plt.colorbar()
	# ch.set_label('std rel. to data')
	# for coord,name in zip(coords,names):
		# xi,yi=coord
		# plt.text(xi+0.01,yi+0.01,''+name[:5],rotation=50,rotation_mode='anchor')
	# plt.tight_layout()
	# if limit_to_data:
		# plt.xlim((xmin-1,xmax+1))
		# plt.ylim((ymin-1,ymax+1))

	# plt.savefig(outdir+'1_'+'relative_std'+'.png',dpi=300)	
	# plt.close()
	
	#---------- analogous for amm 15 ##################
	# labels=['bias','rmse','correlation']
	# itdata=[bias2,rmse2,R2]
	# for data,label in zip(itdata,labels):
		# s.plot_domain_boundaries()
		# fig=plt.gcf()
		# fig.set_size_inches(11,8,forward=True)
		# plt.scatter(x[ivalid],yivalid],s=Rcircle,c=np.asarray([data[name] for name in names[ivalid]]))
		# ch=plt.colorbar()
		# ch.set_label(label)
		# for coord,name in zip(coords[ivalid],names[ivalid]):
			# xi,yi=coord
			# plt.text(xi+0.01,yi+0.01,''+name[:5],rotation=50,rotation_mode='anchor')
		# plt.tight_layout()
		# if limit_to_data:
			# plt.xlim((xmin-1,xmax+1))
			# plt.ylim((ymin-1,ymax+1))
		# plt.title(label + ' Amm')	
		# plt.savefig(outdir+'1_'+label+'amm'+'.png',dpi=300)	
		# plt.close()
		
	#add ration of stadnard deviatons	
	# s.plot_domain_boundaries()
	# fig=plt.gcf()
	# fig.set_size_inches(11,8,forward=True)
	# plt.scatter(xs,ys,s=Rcircle,c=np.asarray([stdmodel2[name]/stddata[name] for name in names[ivalid]]))
	# ch=plt.colorbar()
	# ch.set_label('std rel. to data')
	# for coord,name in zip(coords[ivalid],names[ivalid]):
		# xi,yi=coord
		# plt.text(xi+0.01,yi+0.01,''+name[:5],rotation=50,rotation_mode='anchor')
	# plt.tight_layout()
	# if limit_to_data:
		# plt.xlim((xmin-1,xmax+1))
		# plt.ylim((ymin-1,ymax+1))
	# plt.title(label + ' Amm')	
	# plt.savefig(outdir+'1_'+'relative_std'+'_amm'+'.png',dpi=300)	
	# plt.close()
###############################################	




####### Plot statistics #######################
#ivalid=mask_amm==False
#valid_names=np.asarray(list(bias2.keys()))[ivalid]
shortnames=[name[:3] for name in bias.keys()]
width=0.4
names=np.asarray(names)
x,y=coords[:,0],coords[:,1]
if satistic_maps:
	labels=['bias','rmse','correlation']
	itdata=[bias,rmse,R]
	itdata2=[bias2,rmse2,R2]
	for data,data2,label in zip(itdata,itdata2,labels):
		s.plot_domain_boundaries()
		fig=plt.gcf()
		fig.set_size_inches(11,8,forward=True)
		ph1=plt.scatter(x,y,s=Rcircle,c=np.asarray([data[name] for name in names]))
		ch=plt.colorbar()
		ch.set_label(label)
		for coord,name in zip(coords,names):
			xi,yi=coord
			plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
		plt.tight_layout()
		if limit_to_data:
			plt.xlim((xmin-1,xmax+1))
			plt.ylim((ymin-1,ymax+1))

		#ph2=plt.scatter(x[ivalid],y[ivalid],s=Rcircle*0.25,c=np.asarray([data2[name] for name in names[ivalid]]))
		plt.legend(ph1,'SCHISM',loc='upper center',ncol=2)
		plt.savefig(outdir+'1b_'+label+'_schism'+'.png',dpi=300)	
		plt.close()

		#plt.legend([ph1,ph2],('SCHISM','AMM15'),loc='upper center',ncol=2)
		#plt.savefig(outdir+'1b_'+label+'_schism_amm'+'.png',dpi=300)	
		plt.close()

		# bar charts
		val=np.ma.masked_array(list(data2.values()),mask=mask_amm)
		ph1=plt.bar(range(len(data.keys())),data.values(),-width,align='edge',label='SCHISM',color='b',tick_label=shortnames)
		#ph2=plt.bar(np.where(ivalid)[0],val[ivalid],width,align='edge',label='Amm15',color='r') #,tick_label=valid_names
		plt.xticks(rotation=45)
		plt.grid()
		#plt.legend([ph1,ph2],('SCHISM','Amm15'))
		plt.legend(ph1,ph2,'SCHISM')
		plt.ylabel(label)
		plt.savefig(outdir+'1c_'+label+'_bar_plot',dpi=300)	
		plt.close()
			
	# add ration of stadnard deviatons	
	s.plot_domain_boundaries()
	fig=plt.gcf()
	fig.set_size_inches(11,8,forward=True)
	ph1=plt.scatter(x,y,s=Rcircle,c=np.asarray([stdmodel[name]/stddata[name] for name in names]))
	ch=plt.colorbar()
	ch.set_label('std rel. to data')
	for coord,name in zip(coords,names):
		xi,yi=coord
		plt.text(xi+0.01,yi+0.01,' '+name[:3],rotation=50,rotation_mode='anchor')
	plt.tight_layout()
	if limit_to_data:
		plt.xlim((xmin-1,xmax+1))
		plt.ylim((ymin-1,ymax+1))

	#ph2=plt.scatter(x[ivalid],y[ivalid],s=Rcircle*0.25,c=np.asarray([stdmodel2[name]/stddata[name] for name in names[ivalid]]))
	#plt.legend([ph1,ph2],('SCHISM','AMM15'),loc='upper center',ncol=2)
	plt.legend([ph1,ph2],'SCHISM',loc='upper center',ncol=2)
	#plt.savefig(outdir+'1b_'+'relative_to_data_std_schism_amm'+'.png',dpi=300)	
	plt.savefig(outdir+'1b_'+'relative_to_data_std_schism'+'.png',dpi=300)	
	plt.close()

	# bar charts
	val=np.asarray([stdmodel2[name]/stddata[name] for name in names[ivalid]])
	ph1=plt.bar(range(len(data.keys())),data.values(),-width,align='edge',label='SCHISM',color='b',tick_label=shortnames)
	#ph2=plt.bar(np.where(ivalid)[0],val,width,align='edge',label='Amm15',color='r') #,tick_label=valid_names
	plt.xticks(rotation=45)
	plt.grid()
	#plt.legend([ph1,ph2],('SCHISM','Amm15'))
	plt.legend(ph1,'SCHISM')
	plt.ylabel('std rel. to data')
	plt.savefig(outdir+'1c_relative_std_bar_plot',dpi=300)	
###################################	
	
		











	
##### Plot complete available time series #####
if full_timeseries: 
	print('making plot of complete time series')
	for name in names:
		plt.clf()
		plt.plot(dates,zeta_schism[:,i],'r')
		plt.plot(np.asarray(Time[name]),np.asarray(ZETA[name]),'--')
		#plt.plot(dates2,zeta_amm[:,i],'k')
		plt.legend(('SCHISM','TG'),ncol=2,loc='upper center')					
		#plt.legend(('SCHISM','TG','Amm'),ncol=2,loc='upper center')					
		plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.title(name)
		plt.savefig(outdir+'2_TG'+name+'.png',dpi=300)
###############################################	


	
##### Plot first two weeks of data #####
if first_two_Weeks: 
	print('making plot of first two weeks')
	inds=dates<reftime+dt.timedelta(days=14)
	#inds3=dates2<reftime+dt.timedelta(days=14)
	
	for i,name in enumerate(names):
		Time[name]=np.asarray(Time[name])
		inds2=Time[name]<reftime+dt.timedelta(days=14)
		plt.clf()
		plt.plot(dates[inds],zeta_schism[:,i][inds],'r')
		plt.plot(Time[name][inds2],ZETA[name][inds2],'--') #-dt.timedelta(hours=1) correction now when reading ts
		#plt.plot(dates2[inds3],zeta_amm[:,i][inds3],'k')
		#plt.legend(('SCHISM','TG','Amm'),ncol=3,loc='upper center')					
		plt.legend(('SCHISM','TG'),ncol=2,loc='upper center')					
		plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.title(name)
		plt.grid()
		plt.savefig(outdir+'3_mutc_TG'+name.replace(' ','_').replace('(','').replace(')','').replace('-','_')+'.png',dpi=300)
###############################################	


###### compare data	non normalized 
if monthly_subplots:
	print('making monthly subplots of complete time series')
	for i, name in enumerate(names):
		plt.clf()
		#plt.suptitle(name)
		vmin=np.min((ZETA[name].min(),zeta_schism[:,i].min()))
		vmax=np.max((ZETA[name].max(),zeta_schism[:,i].max()))
		print(name)
		for month in range(1,13): 
			plt.subplot(4,3,month)
			inds=np.asarray([date.month==month for date in Time[name] ])
			plt.plot(np.asarray(Time[name])[inds],np.asarray(ZETA[name])[inds],'b.--')
			inds=np.asarray([date.month==month for date in dates ])
			plt.plot(dates[inds],zeta_schism[:,i][inds],'r')
			inds=np.asarray([date.month==month for date in dates2 ])
			plt.plot(dates2[inds],zeta_amm[:,i][inds],'k--')
									
			if month != 2:
				plt.title('month: '+str(month))
			else:
				plt.title(name)
			plt.grid()
			plt.ylim((vmin,vmax))
			if month%3 !=1:
				plt.tick_params(axis='y',labelleft=False)  
			plt.tick_params(axis='x',labelbottom=False)  					
		#plt.legend(('TG','SCHISM'))		
		fig=plt.gcf()
		fig.set_size_inches(8,8,forward=True)
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()	
		plt.savefig(outdir+'4_TG'+name+'.png',dpi=300)
	
plt.close()	

add_amm=False
if taylor_diag:
	samples=[ [stdmodel[name]/stddata[name],R[name],name[:6]] for name in names ]
	plotTaylor(samples,stdref=1,extend=True) #negative
	#plt.tight_layout()
	plt.savefig(outdir+'5_taylorA.png',dpi=300)
	plt.close()
	plotTaylor(samples,stdref=1,extend=False) #negative
	#plt.tight_layout() # add Amm7 with same nr other color
	plt.savefig(outdir+'5_taylorB.png',dpi=300)
	plt.close()

	# amm 15
	if add_amm:
		samples=[ [stdmodel2[name]/stddata[name],R2[name],name[:6]] for name in names ]
		plotTaylor(samples,stdref=1,extend=True) #negative
		#plt.tight_layout()
		plt.savefig(outdir+'6_taylorA.png',dpi=300)
		plt.close()
		plotTaylor(samples,stdref=1,extend=False) #negative
		#plt.tight_layout() # add Amm7 with same nr other color
		plt.savefig(outdir+'6_taylorB.png',dpi=300)
		plt.close()

	

	
if taylor_diag:
	samples=[ [stdmodel[name]/stddata[name],R[name],name[:3]] for name in names ]
	dia=plotTaylor(samples,stdref=1,extend=True) #negative
	#Add models to Taylor diagram
	if add_amm:
		samples2=[ [stdmodel2[name]/stddata[name],R2[name],name[:3]] for name in names ]
		for i,(stddev, corrcoef, name) in enumerate(samples2):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc='r', mec='r',label=name)

	plt.savefig(outdir+'5_taylorA.png',dpi=300)
	plt.close()

	#samples=[ [stdmodel[name]/stddata[name],R[name],name[:3]] for name in names ]
	dia=plotTaylor(samples,stdref=1,extend=False) #negative
	#Add models to Taylor diagram
	if add_amm:	
		#samples2=[ [stdmodel2[name]/stddata[name],R2[name],name[:3]] for name in names ]
		for i,(stddev, corrcoef, name) in enumerate(samples2):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc='r', mec='r',label=name)

		plt.savefig(outdir+'5_taylorB.png',dpi=300)
		plt.close()

	
#rms=np.sqrt(1-R[name]**2)*stdmodel[name]
#rmse[name]
	
# create latex + pdf
if put_pics_to_texdoc:
		print('generating tex doc and pdf')
		os.chdir(outdir)
		techit(latexname,latextitle,latextext)
		if os.path.isfile(latexname[:latexname.rindex('.')]+'.pdf'):
			 os.remove(latexname[:latexname.rindex('.')]+'.pdf')
		os.system('pdflatex '+latexname)
		print('done generating '+latexname[:latexname.rindex('.')]+'.pdf')


## write error stats to csv
M=[['station name','mean data','mean SCHISM','mean Amm','stddata','std schism','std Amm','bias schism','bias Amm','rmse schism','rmse Amm','correlation schism','correlation Amm']]
for name in names:
	M.append([name,ZETAmean[name],zeta_schismmean[name],zeta_ammmean[name],stddata[name],stdmodel[name],stdmodel2[name],bias[name],bias2[name],rmse[name],rmse2[name],R[name],R2[name]])
M=np.asarray(M)
with open(outdir+'errorstats.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    for i in range(M.shape[0]):
        writer.writerow(M[i,:])
csvFile.close()
		
# hceck shift for coralation +-1 hour	
#output all data as pickle	
analysis={'stations':{'coords':coords,'names':names},
'mes':{'data':{'time':Time,'ssh':ZETA},'data_intp':{'time':timeintp,'ssh':zetaintp}},
'schism':{'data':{'time':dates,'ssh':zeta_schism},'data_intp':{'time':timeintp,'ssh':zetaintp_schism}},
'amm':{'data':{'time':dates2,'ssh':zeta_amm},'data_intp':{'time':timeintp,'ssh':zetaintp_amm}},
'mes_stats':{'mean':ZETAmean,'std':stddata},
'schism_stats':{'mean':zeta_schismmean,'std':stdmodel},
'amm_stats':{'mean':zeta_ammmean,'std':stdmodel2},
'schism_errors_stats':{'bias':bias,'rmse':rmse,'cor':R},
'amm_errors_stats':{'bias':bias2,'rmse':rmse2,'cor':R2}}

pickle.dump(analysis,open(outdir+"analysisdata","wb"))



## compasrions at open boundary
ibd=np.asarray(s.bdy_segments[0][::50])
bdcoords=[(s.lon[i-1],s.lat[i-1]) for i in ibd]
bdnn=amm15.tree.query(bdcoords)[1]
ii2,jj2=np.unravel_index(bdnn,amm15.LON.shape)

date=dates2[0]
amm15.update(date)
tt=amm15.nc['time'][:]
ssh=amm15.nc[amm15.varnames['ssh']][:][:,ii2,jj2]

for date in dates2[nt2:nt2*10+1:nt2]:
	print('loading cmems ' + str(date))
	amm15.update(date)
	tt=np.concatenate((tt,amm15.nc['time'][:]))
	ssh=np.concatenate((ssh,amm15.nc[amm15.varnames['ssh']][:][:,ii2,jj2]))
	
ncs=MFDataset(schismfiles[:10])	
ssh2=ncs['elev'][:,ibd]
tt2=ncs['time'][:]


date=dates2[0]
amm15.update(date)
tt=amm15.t0+dt2*np.arange(len(tt))
tt2=reftime+dt.timedelta(seconds=modeldt)*np.arange(1,len(tt2)+1) 



nsub=int(np.ceil(np.sqrt(len(bdnn))))
plt.close()
for i in range(nsub):
	plt.subplot(nsub,1,i+1)
	ph1=plt.plot(tt,ssh[:,i],'b')
	ph2=plt.plot(tt2,np.asarray(ssh2[:,i]),'r--')
	plt.legend(('Amm','SCHISM'),frameon=False)
	plt.title('bd node ' + str(ibd[i]))
	plt.ylabel('ssh')
	plt.xlabel('time')
#plt.tight_layout()
plt.gcf().autofmt_xdate()
plt.savefig('bnd_comparisons',dpi=300)


