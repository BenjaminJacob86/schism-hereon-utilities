import os
import netCDF4
import sys
import matplotlib
matplotlib.use('Agg') # backend
from matplotlib import path
import datetime as dt
import glob
from scipy.spatial import cKDTree
import xarray as xr
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
import pickle
import datetime as dt
from netCDF4 import Dataset,MFDataset
import pickle
from matplotlib import gridspec  # unequal subplots
from scipy.interpolate import interp1d
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
from techit import * # import latex script
from TaylorDiagram import * # import taylordiagram
from schism import * # import schism functions

# validate agains marnet data as stored as netcdf by sebastian

# Marnet Data: netcd Variables in files
#<station_name>_year.nc
# WTxm Wasser temperature at depth x m
# WSxm Wasser Salnity  at depth x m    


#### I have added 1 day to compensate error in falsely set param.in start date


########## settings #################################
# directories
mardir='/mnt/lustre01/work/gg0028/g260099/MARNET/' 					  # myocean
setupdir='//work/gg0028/g260114/RUNS/GermanBight/GermanBight/'
ncdir=setupdir+'combined2/' #'combined_start_wrongly_1.1/'						  # schism run	
year=2017		    # currently data between 2011 and 2018				
dtol=0.01           										   # distance tolerance in degree lon/lat 

outdir='/work/gg0028/SCHISM/validation/marnet/EU2015/'	   # output directory where images will be stored
if not os.path.exists(outdir): os.mkdir(outdir) 

add_amm15=True   											# add plot comparing with Amm15 | Amm15 profiles have to be extracted with 
	     # other scripts before, because doing it internally for a full year
	     #  takes forever 		


put_pics_to_texdoc=True    										# images will be put in tex document
latexname='GBmarnetAMM.tex'										# in alphanumerical order 
latextitle='GB 2017 vs Marnet'
latextext= str(year) + ' validation against Marnet'

#plt.ion() 
######### load SCHISM setup   ##################################
cwd=os.getcwd()
os.chdir(setupdir)
s=schism_setup()
s.nntree = cKDTree(list(zip(s.lon,s.lat))) 

# initiate file acces
schismfiles=glob.glob(ncdir+'*.nc')	
nrs=[]
for file in schismfiles:
	nrs.append(int(file[file.rindex('_')+1:file.rindex('.')]))
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
s.nc=MFDataset(schismfiles)

# set reference time
reftime=dt.datetime.strptime(s.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')#+dt.timedelta(days=1)
time=reftime+dt.timedelta(days=0.5)
################################################################################################


######### initialize marnet Data in schism domain ##########################
marnet=glob.glob(mardir+'*_{:d}.nc'.format(year))	


######### initialize marnet Data located in SCHISM Domain #######
marnet=glob.glob(mardir+'*_{:d}.nc'.format(year))	

# check if data for 2017
print('searching coordinates for available data in year {:d}'.format(year))
coords=[]
stations=[]
for station in marnet:
	nc=Dataset(station)
	coords.append((nc['ppos'][1],nc['ppos'][0]))
	nc.close()
print('done selecting tide gauges, found: '+str(len(coords)) + ' stations meeting criterea')
#----  reduce to data in domain as given by dtol -----------
lldists,nn=s.nntree.query(coords)		
nn=nn[lldists<=dtol]
marnet=np.asarray(marnet)[lldists<=dtol]
coords=np.asarray(coords)[lldists<=dtol]	
	
# initiate file acces
names=[ file[file.rindex('/')+1:file.rindex('_')]  for file in marnet]
ncs={tag:Dataset(file) for tag,file in zip(names,marnet)}       # create netcdf acces
##############################################################################

	



s.plot_domain_boundaries()
coords=[]
for item in ncs.items():
	tag,nc=item
	yi,xi=nc['ppos'][:]
	coords.append((xi,yi))
	plt.plot(xi,yi,'ro')
	plt.text(xi,yi,' '+tag)
plt.title('Marnet Stations')
plt.savefig(outdir+'0_MarnetStations.png',dpi=300)




# Data put Data in TIme Depth Matrix
T={name:0 for name in names}
S={name:0 for name in names}
Ztemp={name:0 for name in names}
Zsalt={name:0 for name in names}
Time={name:0 for name in names}
for item in ncs.items():
	tag,nc=item
	ncv=nc.variables
	Ti=[]
	Zi_temp=[]
	Zi_salt=[]
	Si=[]
	nc['time'].units
	t0=dt.datetime.strptime(nc['time'].units[12:31],'%Y-%m-%d %H:%M:%S') # data start time
	timeunit=nc['time'].units[:nc['time'].units.index(' ')]
	t=nc['time'][:]
	Time[tag]=[t0+eval('dt.timedelta({:s}=float(ti))'.format(timeunit)) for ti in t]
	for key in ncv.keys():
		if 'WT' in key:
			print(key)
			Zi_temp.append(float(key[2:key.rindex('m')]))
			Ti+=[ncv[key][:]]
		if 'WS' in key:
			print(key)
			Si+=[ncv[key][:]]
			Zi_salt.append(float(key[2:key.rindex('m')]))
	Ztemp[tag]=np.asarray(Zi_temp)			
	Zsalt[tag]=np.asarray(Zi_salt)			
	T[tag]=np.asarray(Ti)			
	S[tag]=np.asarray(Si)			
	


###### SCHISM read next neighbours profiles #############################
ref=dt.datetime.now()
nn=s.nntree.query(coords)[1] # nextneighbour coordinates to marnet stations

nhours=24*360  # length of period
wndw=24#*5       # indexing windows larger faster / memory
				 # too much memory maybe 
# may be loop over files faster
reload=1				 
if reload:				 
	salt=s.nc['salt'][:wndw,:,:][:,nn,:] 
	temp=s.nc['temp'][:wndw,:,:][:,nn,:]  
	zcor=s.nc['zcor'][:wndw,:,:][:,nn,:] 

	# loop in windows
	for ti in range(wndw,nhours,wndw):
		print(ti)
		salt=np.concatenate((salt,s.nc['salt'][ti:ti+wndw,:,:][:,nn,:]),axis=0)
		temp=np.concatenate((temp,s.nc['temp'][ti:ti+wndw,:,:][:,nn,:]),axis=0)
		zcor=np.concatenate((zcor,s.nc['zcor'][ti:ti+wndw,:,:][:,nn,:]),axis=0)
	#####################################################################
	pickle.dump(salt,open("saltAtmarnet.pickle","wb"))
	pickle.dump(temp,open("tempAtmarnet.pickle","wb"))
	pickle.dump(zcor,open("zcorAtmarnet.pickle","wb"))
else:
	salt=pickle.load(open("saltAtmarnet.pickle","rb"))
	temp=pickle.load(open("tempAtmarnet.pickle","rb"))
	zcor=pickle.load(open("zcorAtmarnet.pickle","rb"))





time=reftime+dt.timedelta(hours=1)*np.arange(1,len(s.nc['time']))
TT=np.tile(time[:nhours],[21,1]).swapaxes(0,1)
nhours=8640



# Plot Mar net stations
# monthly   - salt
doplots=0
if doplots:
	plt.figure()
	for i,tag in enumerate(names):
		print(tag)
		
		# salt grid for hovmoeller N
		XX=-np.tile(Zsalt[tag].T,[nhours,1])
		YY=np.tile(Time[tag][23:23+nhours],[len(Zsalt[tag]),1]).swapaxes(0,1)
		
		
		# salt
		# range
		data=S[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([salt[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([salt[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		gs = gridspec.GridSpec(10, 1)
		plt.subplot(gs[0:2])
		plt.plot(Time[tag][23:23+nhours],data[0,23:23+nhours])
		plt.plot(time[:nhours],salt[:,i,-1])
		plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		plt.ylabel('Surf Salt [psu]')
		plt.tick_params(labelbottom=False)

		plt.subplot(gs[3:10])
		
		
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],salt[:,i,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Salinity [psu]')
		plt.ylabel('H [m]')
		plt.title(tag)
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.savefig(outdir+'1_Profile_salt_'+tag+'.png',dpi=300)
		
		
		# temp
		# salt grid for hovmoeller 
		XX=-np.tile(Ztemp[tag].T,[nhours,1])
		YY=np.tile(Time[tag][23:23+nhours],[len(Ztemp[tag]),1]).swapaxes(0,1)

		
		data=T[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([temp[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([temp[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		gs = gridspec.GridSpec(10, 1)
		plt.subplot(gs[0:2])
		plt.plot(Time[tag][23:23+nhours],data[0,23:23+nhours])
		plt.plot(time[:nhours],temp[:,i,-1])
		plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		plt.ylabel('Tsurf [deg C]')
		plt.tick_params(labelbottom=False)

		plt.subplot(gs[3:10])
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],temp[:,i,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Temperature [deg C]')
		plt.ylabel('H [m]')
		plt.title(tag)
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()

		plt.savefig(outdir+'1_Profile_temp_'+tag+'.png',dpi=300)

	
# statistics
# interp model to data
T_schism={name:0 for name in names}
S_schism={name:0 for name in names}
Ztemp_schism={name:0 for name in names}
Zsalt_schism={name:0 for name in names}
timeout={name:0 for name in names}
for i,tag in enumerate(names):
	
	# interpolate to layers
	
	# salt
	Si=[]
	Zi=[]
	for zi in Zsalt[tag]:
		diffs=zi+zcor[:,i,:]
		ibelow=(diffs<0).sum(axis=1)-1
		iabove=ibelow+1
		diffs[range(len(ibelow)),ibelow]
		d1=zcor[range(len(iabove)),i,iabove]+zi
		d2=-zcor[range(len(iabove)),i,ibelow]-zi
		w1=1/d1/(1/d1+1/d2)
		w2=1/d2/(1/d1+1/d2)
		D=w1*zcor[range(len(iabove)),i,iabove]+w2*zcor[range(len(iabove)),i,ibelow]
		Zi.append(D)
		Si.append(w1*salt[range(len(iabove)),i,iabove]+w2*salt[range(len(iabove)),i,ibelow])
	Zsalt_schism[tag]=np.asarray(Zi)
	S_schism[tag]=np.asarray(Si)
	
	# temp
	Ti=[]
	Zi=[]
	for zi in Ztemp[tag]:
		diffs=zi+zcor[:,i,:]
		ibelow=(diffs<0).sum(axis=1)-1
		iabove=ibelow+1
		diffs[range(len(ibelow)),ibelow]
		d1=zcor[range(len(iabove)),i,iabove]+zi
		d2=-zcor[range(len(iabove)),i,ibelow]-zi
		w1=1/d1/(1/d1+1/d2)
		w2=1/d2/(1/d1+1/d2)
		D=w1*zcor[range(len(iabove)),i,iabove]+w2*zcor[range(len(iabove)),i,ibelow]
		Zi.append(D)
		Ti.append(w1*temp[range(len(iabove)),i,iabove]+w2*temp[range(len(iabove)),i,ibelow])
	Ztemp_schism[tag]=np.asarray(Zi)
	T_schism[tag]=np.asarray(Ti)

	# temporal interpolation to data
	Todates=np.asarray(Time[tag])
	ilast=(Todates<=time[nhours-1]).sum()
	ifirst=(Todates<=time[0]).sum()
	Todates=Todates[ifirst:ilast]

	
	tin=[(timei-time[0]).total_seconds() for timei in time]
	tout=[(timei-time[0]).total_seconds() for timei in Todates]

	timeout[tag]=Todates	
	fintp=interp1d(tin[:nhours], S_schism[tag], axis=1)
	S_schism[tag]=fintp(tout)
	fintp=interp1d(tin[:nhours], T_schism[tag], axis=1)
	T_schism[tag]=fintp(tout)
	
	# limit data to smae range
	# Data put Data in TIme Depth Matrix
	T[tag]=T[tag][:,ifirst:ilast]
	S[tag]=S[tag][:,ifirst:ilast]
	Time[tag]=Todates
	
	#mask dumm entries
	T[tag]=np.ma.masked_array(T[tag],mask=T[tag]==-10)
	S[tag]=np.ma.masked_array(S[tag],mask=S[tag]==-10)


	

	
# Same for amm15
oceandir='/work/gg0028/g260099/AMM15/2017/' 						  # myocean
#nn=s.nntree.query(coords)[1] # nextneighbour coordinates to marnet stations

class Amm15:
	def __init__(self, Tfile=None, Sfile=None):
		self.dir=Tfile[:Tfile.rfind('/')+1]
		self.tnc=Dataset(Tfile)
		self.snc=Dataset(Sfile)
		self.t0=dt.datetime.strptime(self.tnc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.tnc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.tnc['time'][:2]))[0]))*np.arange(len(self.tnc['time']))
		self.lon,self.lat=self.tnc['lon'][:],self.tnc['lat'][:]
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		coords=[]
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	

	def get_slab(self,time):		
		nn=np.argsort(np.abs(self.ts-time))[:2] # interpolated linearly
		self.deltaT=np.abs(self.ts[nn]-time)
		dts=np.asarray([dti.total_seconds() for dti in self.deltaT])
		w=1/dts/(1/dts).sum()
		self.T=self.tnc['thetao'][nn[0],:]*w[0]+self.tnc['thetao'][nn[1],:]*w[1]
		self.t=time #self.ts[nn]
		self.S=self.snc['so'][nn[0],:]*w[0]+self.snc['so'][nn[1],:]*w[1]
		
	def get_hov(self,coords):		
		s.nn_moc=self.tree.query(coords)[1]
		ii,jj=np.unravel_index(nn_moc,moc.LON.shape)
		
		self.Sprofile=moc.snc['so'][:][:,:,ii,jj]
		self.Tprofile=moc.tnc['thetao'][:][:,:,ii,jj]
		
	# my ocean
	def gname(self,date):	
		date2=date+dt.timedelta(days=1)	
		return 'metoffice_foam1_amm15_NWS_TEM_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day), 'metoffice_foam1_amm15_NWS_SAL_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day)

	def close(self):
		self.snc.close()
	def update(self,date):
		tname,sname=self.gname(date)
		Tfile=self.dir+tname
		Sfile=self.dir+sname
		self.tnc.close()
		self.snc.close()
		self.tnc=Dataset(Tfile)
		self.snc=Dataset(Sfile)
		self.t0=dt.datetime.strptime(self.tnc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.tnc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.tnc['time'][:2]))[0]))*np.arange(len(self.tnc['time']))


# my ocean
def gname(date):	
	date2=date+dt.timedelta(days=1)	
	return 'metoffice_foam1_amm15_NWS_TEM_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day),'metoffice_foam1_amm15_NWS_SAL_b{}{:02}{:02}_hi{}{:02}{:02}.nc'.format(date2.year,date2.month,date2.day,date.year,date.month,date.day)
tname,sname=gname(time[0])
tmocean=np.unique([dt.datetime(ti.year,ti.month,ti.day) for ti in time])
moc=Amm15(oceandir+tname,oceandir+sname)
nn_moc=moc.tree.query(coords)[1]
ii,jj=np.unravel_index(nn_moc,moc.LON.shape)
# get neixhtneighbours for 

	
# plots

# temp
plt.figure()#(figsize=(25,33)) # full dina 4 widht in inch
#gcf=plt.gcf()
#gcf.set_size_inches(33,33)
#plt.savefig(outdir+'test.png',dpi=300)

compare={tag:{'temp':0,'salt':0} for tag in names}
for key in compare.keys():
	for key2 in compare[key].keys():
		compare[key][key2]={tag:0 for tag in  ['bias','rmse','cor','std1','std2']}
plt.clf()
for i,tag in enumerate(names):
	zmat=np.tile(Ztemp[tag],[nhours-1,1]).T
	plt.clf()
	plt.subplot(2,2,1)
	plt.pcolor(Time[tag],zmat,T[tag])
	plt.colorbar()
	plt.title(tag)
	vmin, vmax = plt.gci().get_clim()
	plt.tick_params(labelbottom=False)
	plt.subplot(2,2,2)
	plt.pcolor(Time[tag],zmat,T_schism[tag],vmin=vmin,vmax=vmax)
	plt.colorbar()
	plt.title('schism')
	plt.subplot(2,2,3)
	plt.pcolor(Time[tag],zmat,T[tag]-T_schism[tag])
	plt.colorbar()
	plt.title(tag + ' - schism')
	plt.gcf().autofmt_xdate()
	mng = plt.get_current_fig_manager()
	plt.tight_layout()
	plt.subplot(2,2,4)
	bias=(T[tag]-T_schism[tag]).mean(axis=1)
	rmse=np.sqrt(((T[tag]-T_schism[tag])**2).mean(axis=1))
	R=[np.corrcoef(T[tag][j,:],T_schism[tag][j,:])[0,1] for j in range(T[tag].shape[0])]
	plt.bar(Ztemp[tag]-0.5,bias, label='bias')
	plt.bar(Ztemp[tag]+0.5, rmse,color='g', label='RMST')
	plt.title('total Bias/rmse {:.2f}/{:.2f}'.format((T[tag]-T_schism[tag]).mean(),np.sqrt( ((T[tag]-T_schism[tag])**2).mean()))) 
	compare[tag]['temp']['bias']=bias		
	compare[tag]['temp']['rmse']=bias		
	compare[tag]['temp']['std1']=np.std(T[tag],axis=1)
	compare[tag]['temp']['std2']=np.std(T_schism[tag],axis=1)	
	compare[tag]['temp']['cor']=R
	for zi,biasi,rmsei in zip(Ztemp[tag],bias,rmse):
		plt.text(zi-0.6,biasi,'{:.2f}'.format(biasi),rotation=90 )
		plt.text(zi+0.6,rmsei,'{:.2f}'.format(rmsei),rotation=90 )
	plt.xlabel('depth')
	plt.legend(('bias','rmse'),frameon=False)
	plt.grid()
	plt.tight_layout()
	plt.savefig(outdir+'3_'+tag+'Temp_statistics.png',dpi=300)


# salt
plt.clf()
for i,tag in enumerate(names):
	zmat=np.tile(Zsalt[tag],[nhours-1,1]).T
	plt.clf()
	plt.subplot(2,2,1)
	plt.pcolor(Time[tag],zmat,S[tag])
	plt.colorbar()
	plt.title(tag)
	vmin, vmax = plt.gci().get_clim()
	plt.tick_params(labelbottom=False)
	plt.subplot(2,2,2)
	plt.pcolor(Time[tag],zmat,S_schism[tag],vmin=vmin,vmax=vmax)
	plt.colorbar()
	plt.title('schism')
	plt.subplot(2,2,3)
	plt.pcolor(Time[tag],zmat,S[tag]-S_schism[tag])
	plt.colorbar()
	plt.title(tag + ' - schism')
	plt.gcf().autofmt_xdate()
	mng = plt.get_current_fig_manager()
	plt.tight_layout()
	plt.subplot(2,2,4)
	bias=(S[tag]-S_schism[tag]).mean(axis=1)
	rmse=np.sqrt(((S[tag]-S_schism[tag])**2).mean(axis=1))
	plt.bar(Zsalt[tag]-0.5,bias, label='bias')
	plt.bar(Zsalt[tag]+0.5, rmse,color='g', label='RMST')
	plt.title('total Bias/rmse {:.2f}/{:.2f}'.format((S[tag]-S_schism[tag]).mean(),np.sqrt( ((S[tag]-S_schism[tag])**2).mean()))) 

	for zi,biasi,rmsei in zip(Zsalt[tag],bias,rmse):
		plt.text(zi-0.6,biasi,'{:.2f}'.format(biasi),rotation=90 )
		plt.text(zi+0.6,rmsei,'{:.2f}'.format(rmsei),rotation=90 )
	plt.xlabel('depth')
	plt.legend(('bias','rmse'),frameon=False)
	plt.grid()
	plt.tight_layout()
	plt.savefig(outdir+'3_'+tag+'Salt_statistics.png',dpi=300)

	bias=(S[tag]-S_schism[tag]).mean(axis=1)
	rmse=np.sqrt(((S[tag]-S_schism[tag])**2).mean(axis=1))
	R=[np.corrcoef(S[tag][j,:],S_schism[tag][j,:])[0,1] for j in range(S[tag].shape[0])]
	compare[tag]['salt']['bias']=bias		
	compare[tag]['salt']['rmse']=bias		
	compare[tag]['salt']['std1']=np.std(S[tag],axis=1)
	compare[tag]['salt']['std2']=np.std(S_schism[tag],axis=1)	
	compare[tag]['salt']['cor']=R
	
plt.figure()
plt.bar(Zsalt[tag]-0.1,compare[tag]['salt']['std1'])
plt.bar(Zsalt[tag]+0.1,compare[tag]['salt']['std2'])


# top botoopm
label=['top','bottom']
for nr,i in enumerate([0,-1]):
	samples=[ [compare[tag]['temp']['std2'][i]/compare[tag]['temp']['std1'][i],compare[tag]['temp']['cor'][i] ,tag+' '+str(Zsalt[tag][i])+'m'] for tag in names ]
	plotTaylor(samples=samples)
	plt.suptitle('Temperature ' +label[nr])
	plt.savefig(outdir+'5_taylor_temp'+label[nr]+'.png',dpi=300)
	plt.close()

# top botoopm
label=['top','bottom']
for nr,i in enumerate([0,-1]):
	samples=[ [compare[tag]['salt']['std2'][i]/compare[tag]['salt']['std1'][i],compare[tag]['salt']['cor'][i] ,tag+' '+str(Zsalt[tag][i])+'m'] for tag in names ]
	plotTaylor(samples=samples)
	plt.suptitle('Salinity ' +label[nr])
	plt.savefig(outdir+'5_taylor_salt'+label[nr]+'.png',dpi=300)	
	plt.close()	
	
# comparison with amm15	
if add_amm15: 
	##### Add Amm 15	
	sfiles15=np.sort(glob.glob('/work/gg0028/g260114/RUNS/GermanBight/paraExtract/*SAL*merged*'))
	tfiles15=np.sort(glob.glob('/work/gg0028/g260114/RUNS/GermanBight/paraExtract/*TEM*merged*'))
	##
	ncsS={tag:Dataset(file) for tag,file in zip(names,sfiles15)}       # create netcdf acces 
	ncsT={tag:Dataset(file) for tag,file in zip(names,tfiles15)}       # create netcdf acces 
	 
	T15={name:0 for name in names}
	S15={name:0 for name in names}
	Time15={name:0 for name in names}
	Z15=moc.snc['depth'][:]	
	for item in ncsS.items():
		tag,nc=item
		nc['time'].units
		t0=dt.datetime.strptime(nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S') # data start time
		timeunit=nc['time'].units[:nc['time'].units.index(' ')]
		t=nc['time'][:]
		Time15[tag]=[t0+eval('dt.timedelta({:s}=float(ti))'.format(timeunit)) for ti in t]
		S15[tag]=nc['so'][:,:,0,0].swapaxes(0,1)
		
	for item in ncsT.items():
		tag,nc=item
		T15[tag]=nc['thetao'][:,:,0,0].swapaxes(0,1)
		

	############### plot itner comparison
	plt.figure()
	for i,tag in enumerate(names):
		print(tag)
		
		# salt grid for hovmoeller N
		XX=-np.tile(Zsalt[tag].T,[nhours,1])
		YY=np.tile(Time[tag][23:23+nhours],[len(Zsalt[tag]),1]).swapaxes(0,1)
		
		
		# salt
		# range
		data=S[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([salt[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([salt[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		#gs = gridspec.GridSpec(10, 1)
		#plt.subplot(gs[0:2])
		#plt.plot(Time[tag][23:23+nhours],data[0,23:23+nhours])
		#plt.plot(time[:nhours],salt[:,i,-1])
		#plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		#plt.ylabel('Surf Salt [psu]')
		#plt.tick_params(labelbottom=False)

		#plt.subplot(gs[3:10])
		
		plt.subplot(2,1,1)
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],salt[:,i,:],vmin=vmin,vmax=vmax)
		#ch=plt.colorbar(orientation='horizontal')
		#ch.set_label('Salinity [psu]')
		plt.ylabel('H [m]')
		plt.title(tag + 'SCHISM')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.subplot(2,1,2)
		# plot amm15 profile
		iend=(T15[tag][:,0].mask==False).sum()
		plt.pcolormesh(Time15[tag],-Z15[:iend],S15[tag][:iend,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Salinity [psu]')
		plt.ylabel('H [m]')
		plt.title(tag + ' Amm15')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		
		plt.savefig(outdir+'Profile_salt2_'+tag+'.png',dpi=300)
		
		
		# temp
		# salt grid for hovmoeller 
		XX=-np.tile(Ztemp[tag].T,[nhours,1])
		YY=np.tile(Time[tag][23:23+nhours],[len(Ztemp[tag]),1]).swapaxes(0,1)

		
		data=T[tag]
		data=np.ma.masked_array(data,mask=data==-10) # mask  Salt
		vmin=np.min([temp[:,i,:].flatten().min(),data[:,23:23+nhours].swapaxes(0,1).flatten().min()])
		vmax=np.max([temp[:,i,:].flatten().max(),data[:,23:23+nhours].swapaxes(0,1).flatten().max()])
		
		# range
		plt.clf()
		#gs = gridspec.GridSpec(10, 1)
		#plt.subplot(gs[0:2])
		#plt.plot(Time[tag][23:23+nhours],data[0,23:23+nhours])
		#plt.plot(time[:nhours],temp[:,i,-1])
		#plt.legend((tag,'Schism'),loc='upper center',ncol=2, bbox_to_anchor=(0.5, 1.5))	
		#plt.ylabel('Tsurf [deg C]')
		#plt.tick_params(labelbottom=False)

		#plt.subplot(gs[3:10])
		plt.subplot(2,1,1)
		# plot schism profile
		plt.pcolor(TT,zcor[:,i,:],temp[:,i,:],vmin=vmin,vmax=vmax)
		#ch=plt.colorbar(orientation='horizontal')
		#ch.set_label('Temperature [deg C]')
		plt.ylabel('H [m]')
		plt.title(tag + 'SCHISM')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		
		plt.subplot(2,1,2)
		# plot amm15 profile
		plt.pcolormesh(Time15[tag],-Z15[:iend],T15[tag][:iend,:],vmin=vmin,vmax=vmax)
		ch=plt.colorbar(orientation='horizontal')
		ch.set_label('Temperature [deg C]')
		plt.ylabel('H [m]')
		plt.title(tag + ' Amm5')
		
		# add marnet data as scatter
		sampling=6
		colors=data[:,23:23+nhours].swapaxes(0,1)[::sampling,:].flatten()
		ivalid= colors.mask==False
		if sum(ivalid):
			ph=plt.scatter(YY[::sampling,:].flatten()[ivalid],XX[::sampling,:].flatten()[ivalid],30,c=colors[ivalid],vmin=vmin,vmax=vmax,linewidths=None,marker='s',edgecolors='none')
		#plt.gcf().autofmt_xdate()
		plt.tight_layout()
		plt.savefig(outdir+'Profile_temp2_'+tag+'.png',dpi=300)

	


# Statists
#salt
#b=np.vstack((compare[tag]['salt']['std1'],compare[tag]['salt']['std2']))

# create latex + pdf
if put_pics_to_texdoc:
		print('generating tex doc and pdf')
		os.chdir(outdir)
		techit(latexname,latextitle,latextext)
		if os.path.isfile(latexname[:latexname.rindex('.')]+'.pdf'):
			 os.remove(latexname[:latexname.rindex('.')]+'.pdf')
		os.system('pdflatex '+latexname)
		print('done generating '+latexname[:latexname.rindex('.')]+'.pdf')
