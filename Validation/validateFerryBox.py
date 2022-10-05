import os
import sys
import matplotlib
matplotlib.use('Agg') # backend
from matplotlib import gridspec
import datetime as dt
import glob
from scipy.spatial import cKDTree
import pickle
from netCDF4 import Dataset,MFDataset
import pickle
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')

from schism import * # import schism functions
import datetime as dt
try:
	from netcdftime import utime
except:	
	from cftime import utime

# validate SCHISM and Amm15 against netcdf Ferrybox track data as available from KSAMBA:
# \\ksamba\cosyna\netcdf\ferrybox_hzg.
# validation is done for all ferrybox nc files within folder  fbdir (below)
# Amm15 is currently linearly interpolated in time using the sptial next neighbour points
# schism is currently using nextneighbourt time step field using a barycentric interpolation in space
#
# Use:
#	module load python3
#	set adaquate setting in setting section below and run with python
#
# output:
# each carrying name of NR assoiciate with transect name
#
# png images (stored in subfolders named by ferryboxtrack nr):
# Plot Temperature and Salinity fields for SCHISM and Amm15 each overimposed with colorcoded
# values from ferrybox observations as also interpolated graphs plottet as along track timeseries  
#
# exported ascii files: 
#	1) FBvalid_Track<NR>.ascii
#     interpolated along track data (time and coordinates of ferrybox track) for fields
#     time(days since 2017-01-01 00:00:00) lon lat  obs_T[degc] obs_S[psu]  schism_T[degc] schism_S[psu]  amm_T[degc] amm_S[psu] 
#   example:
#  	time(days since 2017-01-01 00:00:00) lon lat  obs_T[degc] obs_S[psu]  schism_T[degc] schism_S[psu]  amm_T[degc] amm_S[psu] 
#	2.133182870370370381e+02 8.842829704284667969e+00 5.412107086181640625e+01 2.112454986572265625e+01 2.654088592529296875e+01 1.926928928797392615e+01 2.757626810651600735e+01 1.985830555555555321e+01 2.230683333333333351e+01
#	2.133185185185185162e+02 8.840401649475097656e+00 5.412164688110351562e+01 2.084615707397460938e+01 2.654479217529296875e+01 1.926857362027735476e+01 2.757671601187153598e+01 1.985847777777777878e+01 2.230566666666666364e+01
#	2.133187500000115620e+02 8.838033676147460938e+00 5.412221145629882812e+01 2.072570037841796875e+01 2.654428100585937500e+01 1.927040762047127842e+01 2.757693339053808046e+01 1.985865000000860903e+01 2.230449999994166177e+01
#
#	2) FBvalid_Track<NR>_schismfield.ascii
#	values of schism surface slabe at nextneighbour time step for transect median time
#   example:
#	% SCHISM field of validation for  /work/gg0028/SCHISM/validation/FerryBox/data2017/HZG-Ferrybox_Bues-HelgoNR65536.nc
#	% at 2017-08-02 08:36:20.000010 
#	%   lon lat T[degc] S[degc]
#	 7.314304199999999589e+00 5.303926374000000266e+01 1.831864166259765625e+01 1.661940752144712634e-20
#	...						  ...						...						...	
#
# 	3) FBvalid_Track{<NR>_ammfield.ascii
#	values of amm surface slab (one dimensionlaizeed (reshape(nlat,nlon) for 2D field )) at nextneighbour time step for transect median time
#   example:
#   #  Amm  field (nlon=958,nlat=1240  as row vector) of validation for  /work/gg0028/SCHISM/validation/FerryBox/data2017/HZG-Ferrybox_Bues-HelgoNR65536.nc
#    at 2017-08-02 08:36:20.000010 
#    lon lat T[degc] S[degc]
#	-1.600000000000000000e+01 4.600000000000000000e+01 1.833658333337500324e+01 3.561021111111666926e+01
#	-1.596969985961914062e+01 4.600000000000000000e+01 1.837368888893333363e+01 3.560421111111666903e+01
#	-1.593939018249511719e+01 4.600000000000000000e+01 1.841458333337500264e+01 3.560600000000000165e+01
#	-1.590909004211425781e+01 4.600000000000000000e+01 1.847347777781666878e+01 3.560800000000000409e+01


########## settings #################################
# directories
fbdir='/gpfs/work/jacobb/data/SETUPS/SNS_Wei/valid/FB/' 					  # folder containing FerryBox nc tracks 
setupdir='/gpfs/work/chenw1/SCHISM/cfgs/Control_Run/SNS_GB_Combined/'  				  # schism rundir	
ncdir='/gpfs/work/chenw1/SCHISM/cfgs/Control_Run/SNS_GB_Combined/outputs/Combined_drag_turned_1/'					  #	output dir
ncdir='/gpfs/work/jacobb/data/SETUPS/SNS_Wei/valid/combined_linked/'
#oceandir='/work/gg0028/g260099/AMM15/RAW/' 									  # amm15 directory
oceandir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/amm15/'
outdir='/gpfs/work/jacobb/data/SETUPS/SNS_Wei/valid/' # output directory where images will be stored

year=2018		  														  # year used as reference time (year,0,0,0,0,0)			
# appearance
lw=0.0																	  # linewidth around overimposed circles
ms=45																	  #  markersize	
cm=plt.cm.jet				  											  # colormap
ds=0.5						 			 # area of zoom pictures is ds large in all direction of rectangle containing ferry track
names=['FB','SCHISM','AMM15'] 												#names

# output
make_plots=True   						# True: plot each transect while extracting date. False: extract all data and plot afterwards
export_to_ascii=True    				# export transects as ascii
comment='%'								# comment symbole before e.g. % or # to load with matlab/python 
######################### / Settings ##############################################################################

################ functions #####################################
def find_parent_tri(gr,xq,yq,dThresh=1000):
	""" parents,ndeweights=find_parent_tri(gr,xq,yq,dThresh=1000)
		find parent for coordinates xq,yq within triangulation tris,xun,yun.
		return: parent triangle ids and barycentric weights of triangle coordinates.
		Works only for pure triangle grids	
	"""    
	#% Distance threshold for Point distance
	dThresh=dThresh**2
	
	xun=gr.x
	yun=gr.y
	tris=np.asarray(gr.faces)
	trinr=np.arange(tris.shape[0])
	
	trisX,trisY=xun[tris],yun[tris]
	#% orthogonal of side vecotrs
	SideX=np.diff(trisY[:,[0, 1, 2, 0]],axis=1)
	SideY=-np.diff(trisX[:,[0, 1, 2, 0]],axis=1)
	
	p=np.stack((xq,yq),axis=1)
	parent=-1*np.ones(len(p),int)
	for ip in range(len(p)):

			dx1=(p[ip,0]-trisX[:,0])
			dy1=(p[ip,1]-trisY[:,0])
			subind=(dx1*dx1+dy1*dy1) < dThresh # preselection
			subtris=trinr[subind]
			
			#% dot products
			parenti=(subtris[ (dx1[subind]*SideX[subind,0] + dy1[subind]*SideY[subind,0] <= 0) \
						   & ((p[ip,0]-trisX[subind,1])*SideX[subind,1] + (p[ip,1]-trisY[subind,1])*SideY[subind,1] <= 0) \
							 & ( (p[ip,0]-trisX[subind,2])*SideX[subind,2] + (p[ip,1]-trisY[subind,2])*SideY[subind,2] <= 0) ][:])
			if len(parenti):
				parent[ip]=parenti
	
	# tri nodes
	xabc=xun[tris[parent]]
	yabc=yun[tris[parent]]
	
	# barycentric weights
	divisor=(yabc[:,1]-yabc[:,2])*(xabc[:,0]-xabc[:,2])+(xabc[:,2]-xabc[:,1])*(yabc[:,0]-yabc[:,2])
	w1=((yabc[:,1]-yabc[:,2])*(xq-xabc[:,2])+(xabc[:,2]-xabc[:,1])*(yq-yabc[:,2]))/divisor
	w2=((yabc[:,2]-yabc[:,0])*(xq-xabc[:,2])+(xabc[:,0]-xabc[:,2])*(yq-yabc[:,2]))/divisor
	w3=1-w1-w2

	return parent,np.stack((w1,w2,w3)).transpose() 

class grid_from_nc:
	def __init__(self, nc):
		self.ncv=nc.variables
		self.faces=self.ncv['SCHISM_hgrid_face_nodes'][:,:3]-1
		self.x=self.ncv['SCHISM_hgrid_node_x'][:]
		self.y=self.ncv['SCHISM_hgrid_node_y'][:]

class Amm15:
	def __init__(self, Tfile=None, Sfile=None):
		self.dir=Tfile[:Tfile.rfind('/')+1]
		self.tnc=Dataset(Tfile)
		self.snc=Dataset(Sfile)
		self.t0=dt.datetime.strptime(self.tnc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=float(self.tnc['time'][0]))
		self.ts=self.t0+dt.timedelta(seconds=np.double(np.diff((self.tnc['time'][:2]))[0]))*np.arange(len(self.tnc['time']))
		self.lon,self.lat=self.tnc['lon'][:],self.tnc['lat'][:]
		self.LON,self.LAT=np.meshgrid(self.lon,self.lat)
		self.depths=self.tnc['depth']
		coords=[]
		for xi,yi in zip(self.LON.flatten(),self.LAT.flatten()):
			coords.append((xi,yi))
		self.tree = cKDTree(coords) 	

	def get_slab(self,time,level=0):		
		nn=np.argsort(np.abs(self.ts-time))[:2] # interpolated linearly
		self.deltaT=np.abs(self.ts[nn]-time)
		dts=np.asarray([dti.total_seconds() for dti in self.deltaT])
		w=1/dts/(1/dts).sum()
		self.T=self.tnc['thetao'][nn[0],level,:]*w[0]+self.tnc['thetao'][nn[1],level,:]*w[1]
		self.t=time #self.ts[nn]
		self.S=self.snc['so'][nn[0],level,:]*w[0]+self.snc['so'][nn[1],level,:]*w[1]
		
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
####################################/ functions #################################################



################################################ Program Start ###############################################
if not os.path.exists(outdir): os.mkdir(outdir) 

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

gr=grid_from_nc(s.nc)
gr.x=np.asarray(s.lon)
gr.y=np.asarray(s.lat)
ut2=utime(s.nc['time'].units)
schismdates=ut2.num2date(s.nc['time'][:])

nrs=[]
for file in schismfiles:
	nrs.append(int(file[file.rindex('_')+1:file.rindex('.')]))
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
s.nc=MFDataset(schismfiles)		
################################################################################################

###########  iterate over transects #################################
fields=['Nr','ctime','ship','lon','lat','time','obs','schism','amm','fname']
tslices={'schism':[],'amm':[]}
sslices={'schism':[],'amm':[]}
extract=[]

FBt=list(np.sort(glob.glob(fbdir+'*.nc')))	
for filenr,file in	enumerate(FBt): #FBt[10],: #
	print('doing ' +file)
	
	extract.append(dict.fromkeys(fields))
	
	try:
		NR=file[file.index('NR')+2:file.index('.')]
	except:	
		NR=file[file.rindex('_')+1:file.index('.')]
	outfolder=outdir+NR+'/'
	if not os.path.exists(outfolder): os.mkdir(outfolder) 

	# load Data
	fbnc=Dataset(file)
	lon=fbnc['LONGITUDE'][:]
	lat=fbnc['LATITUDE'][:]
	time=fbnc['TIME'][:]
	T=fbnc['TEMP'][:]
	S=fbnc['PSAL'][:]
	Tqc=fbnc['TEMP_QC'][:]
	Sqc=fbnc['PSAL_QC'][:]

	# get dates
	ut=utime(fbnc['TIME'].units)
	fbdates=ut.num2date(time)
	
	# calculate parent grid elements
	parents,w=find_parent_tri(gr,lon,lat,dThresh=3)

	# get track part whithin schism domain 
	iuse=np.where(parents!=-1)	# trackpoints within region
	if len(iuse[0])==0:
		continue
	cdate=fbdates[int(np.median(iuse))]  # central time to plot SCHISM
	
	# schism profiles
	Tintp=np.zeros(len(T[iuse]))
	Sintp=np.zeros(len(T[iuse]))
	tnni=np.zeros(len(T[iuse]))
	
	# amm profiles
	Tintp2=np.zeros(len(T[iuse]))
	Sintp2=np.zeros(len(T[iuse]))
	tnni2=np.zeros(len(T[iuse]))
	parents2=parents[iuse]
	
	# load amm
	tname,sname=gname(fbdates[iuse][0])
	moc=Amm15(oceandir+tname,oceandir+sname)
	dt_amm=(moc.ts[2]- moc.ts[1]).total_seconds()/3600
	lon2=lon[iuse]
	lat2=lat[iuse]
	
	# interpolate Amm
	for i,date in enumerate(fbdates[iuse]):
			tnni[i]=np.argmin(np.abs(schismdates-date))
			if (fbdates[iuse][0]+dt.timedelta(hours=dt_amm)-moc.ts[0]).days>0:
				moc.update(date)	
			moc.get_slab(date,level=0)	
			Sintp2[i]=moc.S.flatten()[moc.tree.query((lon2[i],lat2[i]))[1]]
			Tintp2[i]=moc.T.flatten()[moc.tree.query((lon2[i],lat2[i]))[1]]
	print('done interpolating AMM15')		
			
	# interpolate schism
	w2=w[iuse,:][0,:]	
	for iun in np.unique(tnni):
		inds=np.where(tnni==iun)[0]
		tslice=s.nc['temp'][iun,:,-1]
		sslice=s.nc['salt'][iun,:,-1]
		Tintp[inds]=(tslice[gr.faces[parents2[inds],:]]*w2[inds,:]).sum(axis=1)
		Sintp[inds]=(sslice[gr.faces[parents2[inds],:]]*w2[inds,:]).sum(axis=1)
	print('done interpolating SCHISM')		

	# append Data	
	extract[filenr]['fname']=file
	extract[filenr]['ctime']=cdate
	extract[filenr]['ship']=fbnc.platform_code
	extract[filenr]['lon']=lon2
	extract[filenr]['lat']=lat2
	extract[filenr]['obs']={'T':np.asarray(T[:,0][iuse]),'S':np.asarray(S[:,0][iuse])}
	extract[filenr]['schism']={'T':np.asarray(Tintp),'S':np.asarray(Sintp)}
	extract[filenr]['amm']={'T':np.asarray(Tintp2),'S':np.asarray(Sintp2)}
	extract[filenr]['time']=fbdates[iuse]
	extract[filenr]['Nr']=NR
	fbnc.close()
	
	# apply quality control
	extract[filenr]['obs']['T'][Tqc[:,0][iuse]>2]=np.nan
	extract[filenr]['obs']['S'][Sqc[:,0][iuse]>2]=np.nan	


	# horizontal slice data
	tslices['amm'].append(moc.T.flatten())
	sslices['amm'].append(moc.S.flatten())
	tslices['schism'].append(tslice)
	sslices['schism'].append(sslice)

	# make plots
	
	for var in 'salt','temp':
	
		if var=='salt':
			varname='S'
			
			pdata2=moc.S
			label='Salinity [psu]'
		else:
			varname='T'
			pdata2=moc.T
			label='T [deg C]'
	
		# stats 
		diff=extract[filenr]['schism'][varname]-extract[filenr]['obs'][varname]
		diff2=extract[filenr]['amm'][varname]-extract[filenr]['obs'][varname]
		bias=[np.nanmean(diff), np.nanmean(diff2)]
		rmse=[np.sqrt(np.nanmean((diff**2))),np.sqrt(np.nanmean((diff2**2)))]
		ivalid=((np.isnan(extract[filenr]['amm'][varname]) ==0) & (np.isnan(extract[filenr]['obs'][varname])==0) )
		R=np.corrcoef(extract[filenr]['obs'][varname][ivalid],extract[filenr]['schism'][varname][ivalid])[0,1]
		R2=np.corrcoef(extract[filenr]['obs'][varname][ivalid],extract[filenr]['amm'][varname][ivalid])[0,1]
		cor=[R,R2]	
		extract[filenr]['stats']={'label':names[1:],'bias':bias,'rmse':rmse,'cor':cor}
		
		if make_plots:
			fig = plt.figure(figsize=(14,10))
			gs = gridspec.GridSpec(2, 2)
			
			tnn=np.argmin(np.abs(schismdates-cdate))
			moc.update(cdate)	
			moc.get_slab(cdate)	

			fig = plt.figure(figsize=(14,10))
			gs = gridspec.GridSpec(2, 2)

			
			pdata=s.nc[var][tnn,:,-1]

			vmin=pdata.min()
			vmax=pdata.max()
			for key in 'obs','schism','amm':
				vmin=np.floor(np.min((vmin,np.nanmin(extract[filenr][key][varname]))))
				vmax=np.ceil(np.max((vmax,np.nanmax(extract[filenr][key][varname]))))
		
			plt.subplot(gs[0])
			ph,ch=s.plotAtnodes(pdata)
			plt.clim(vmin,vmax)
			plt.scatter(extract[filenr]['lon'],extract[filenr]['lat'],s=ms,c=extract[filenr]['obs'][varname],vmin=vmin, vmax=vmax, cmap=cm,edgecolor=None,linewidth=lw)
			plt.axis('tight')
			ch.set_label(label)
			plt.title('SCHISM @' +str(ut2.num2date(s.nc['time'][tnn])))
			plt.text(extract[filenr]['lon'][0],extract[filenr]['lat'][0], str(extract[filenr]['time'][0]))
			plt.text(extract[filenr]['lon'][-1],extract[filenr]['lat'][-1], str(extract[filenr]['time'][-1]))
			plt.tight_layout()		
			ax=plt.axis()
			
			plt.subplot(gs[1])
			plt.pcolormesh(moc.lon,moc.lat,pdata2,vmin=vmin,vmax=vmax,cmap=cm)
			ch=plt.colorbar()
			plt.clim(vmin,vmax)
			plt.scatter(extract[filenr]['lon'],extract[filenr]['lat'],s=ms,c=extract[filenr]['obs'][varname],vmin=vmin, vmax=vmax, cmap=cm,edgecolor=None,linewidth=lw)
			plt.axis('tight')
			ch.set_label(label)
			plt.title('AMM15 @' +str(moc.t))
			plt.text(extract[filenr]['lon'][0],extract[filenr]['lat'][0], str(extract[filenr]['time'][0]),rotation=0)
			plt.text(extract[filenr]['lon'][-1],extract[filenr]['lat'][-1], str(extract[filenr]['time'][-1]),rotation=0)
			plt.axis(ax)
			plt.tight_layout()		
			
			plt.subplot(gs[2:4])
			extract[filenr]['time']=np.asarray([np.datetime64(ti) for ti in extract[filenr]['time']])
			for key in 'obs','schism','amm':
				plt.plot(extract[filenr]['time'] ,extract[filenr][key][varname])
			plt.grid('on')
			plt.legend(names,ncol=3,loc='upper center',frameon=False)
			plt.ylabel(label)
			plt.gcf().autofmt_xdate()
			plt.tight_layout()
			plt.savefig(outfolder+'FB_{:d}_{:s}.png'.format(i,var),dpi=300)
			print('done plotting transect')	
			# zoom
			ax=(np.min(extract[filenr]['lon'])-ds,np.max(extract[filenr]['lon'])+ds,np.min(extract[filenr]['lat'])-ds,np.max(extract[filenr]['lat'])+ds )
			plt.subplot(gs[0])
			plt.axis(ax)
			plt.subplot(gs[1])
			plt.axis(ax)
			plt.savefig(outfolder+'FB_{:d}_{:s}_zoom.png'.format(i,var),dpi=300)
			plt.close()
	
			# statistics plots
			ph1=plt.bar([1,4],bias,color='b')
			ph2=plt.bar([2,5],rmse,color='r')
			ph3=plt.bar([3,6],cor,color='k')
			
			plt.legend((ph1,ph2,ph3),('bias','rmse','cor'),ncol=2,loc='upper center',frameon=False)
			plt.grid('on')
			plt.xticks([3, 6],('SCHISM','AMM15'))
			plt.ylabel(label)
			plt.savefig(outfolder+'FB_{:d}_{:s}_stats2.png'.format(i,var),dpi=300)
			plt.title( extract[filenr]['ship'] +' ' +str(extract[filenr]['time'][0])  +str(extract[filenr]['time'][-1]) )
			plt.close()
pickle.dump(extract,open("FerryBox.pickle","wb"))	

# export ascii data
if export_to_ascii:
	reftime=dt.datetime(year,1,1,0,0,0)
	amm_lon,amm_lat=np.meshgrid(moc.lon,moc.lat)

	for i in range(len(extract)):
		time=np.asarray([ (ti-reftime).total_seconds()/86400 for ti in extract[i]['time'] ])
	
		header=' validation for ' + file + '\n'
		header+=' ship: ' + extract[i]['ship'] + ' \n stats: \n' 
		header+='quantity {:s} {:s}'.format(names[1],names[2])
		for tag in 'bias','rmse','cor':
			header+='\n {:s} {:f} {:f} '.format(tag,extract[i]['stats'][tag][0],extract[i]['stats'][tag][1])
		header+='\n time(days since {:s}) lon lat '.format(str(reftime))
		
		m=([time,]+[ extract[i][tag] for tag in ['lon','lat'] ] )
		for tag in ['obs','schism','amm']:
			m+=[extract[i][tag]['T'],extract[i][tag]['S']]
			header+=' {:s}_T[degc] {:s}_S[psu] '.format(tag,tag)
		m=np.asarray(m).T
		np.savetxt(outdir+'FBvalid_Track{:s}.ascii'.format(extract[i]['Nr']),m,header=header)

		header=' SCHISM field of validation for  ' + file + '\n at ' + str(extract[i]['ctime']) + ' \n  '
		header+=' lon lat T[degc] S[psu]'
		m=np.vstack((s.lon,s.lat,tslices['schism'],sslices['schism'])).T
		np.savetxt(outdir+'FBvalid_Track{:s}_schismfield.ascii'.format(extract[i]['Nr']),m,header=header,comments='%')

		header=' Amm  field (nlon={:d},nlat={:d}  as row vector) of validation for  '.format(len(moc.lon),len(moc.lat)) + file + '\n at ' + str(extract[i]['ctime']) + ' \n  '
		header+=' lon lat T[degc] S[psu]'
		m=np.vstack((amm_lon.flatten(),amm_lat.flatten(),tslices['amm'],sslices['amm'])).T		
		np.savetxt(outdir+'FBvalid_Track{:s}_ammfield.ascii'.format(extract[i]['Nr']),m,header=header)
print('done validating')