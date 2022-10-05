""""
SCHISM Wave modeling

Perform intercomparison of WW3 data Model

"""

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
plt.ion()


##########S E T T I N G S###
## WWW3 ###
ww3dir='/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd_2017/longrun/'
ww3mesh='NBSext_bl.msh'
ww3PointList='GB_bd_points.list'
ww3PointOutput='ww3_out_1monthsWInd/ww3.201701_spec.nc'
ww3GridOutput='ww3_out_1monthsWInd/ww3.2017.nc'
##WWMM ###


ww3PointOutput=[ww3dir+'ww3_out_1monthsWInd/ww3.201701_spec.nc',ww3dir+'ww3_out_1monthsWInd/ww3.201702_spec.nc']


##### appearance

def set_FS(Fincrease=1.4):
	#SMALL_SIZE = 10*Fincrease
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


############################


#######  Wave Watch 3   ######
class WW3_mesh():
	import numpy as np
	from matplotlib import pyplot as plt
	from matplotlib.collections import PolyCollection
	import xarray as xr
	def __init__(self,file):
		with open(file) as f:
			
			lines=f.readlines()

			nodestart=lines[lines.index('$Nodes\n')+1]
			nodestart=lines.index(nodestart)
			nnodes=np.int(lines[lines.index('$Nodes\n')+1].split(' ')[0])

			elemstart=lines[lines.index('$Elements\n')+1]
			elemstart=lines.index(elemstart)
			nelem=np.int(lines[lines.index('$Elements\n')+1].split(' ')[0])

			self.nodes=np.loadtxt(file,skiprows=nodestart+1,max_rows=nnodes)[:,1:]
			self.elems=np.loadtxt(file,skiprows=elemstart+1,max_rows=nelem)[:,6:]-1
			self.elems=np.asarray(self.elems,np.int)
			self.x,self.y,self.d=self.nodes[:,0],self.nodes[:,1],self.nodes[:,2]
			
			from scipy.spatial import cKDTree
			self.node_tree = cKDTree(list(zip(self.x,self.y)))


	def find_nearest_node(self,x,y):
		"""
		find nearest node for given coordinate,
		returns the node index (zero based)
		"""
		ridx=-1
		d,idx = self.node_tree.query((x,y),k=1)
		#ridx = self.nodeids[idx]
		return idx	
			
		  # plot functions 
	def plotAtnodes(self,nodevalues,cmap=plt.cm.jet,mask=None):
		"""
		visualisation routine plotting triangles at nodes (quads are splitted)
		"""
		ph=plt.tripcolor(self.x,self.y,self.elems,nodevalues,shading='flat',mask=mask,cmap=cmap)	  
		ch=plt.colorbar()

		return ph,ch		
		
	def plot_mesh(self,tri_color='k',quad_color='m',linewidth=0.2):	  
		xy=np.c_[self.x,self.y]
		tripc = PolyCollection(xy[self.elems,:3],facecolors='none',edgecolors=tri_color,linewidth=linewidth) #, **kwargs)
		plt.gca().add_collection(tripc)
		
	def load_points(self,fname,sep=' '):	  
		with open(fname) as f:
			x,y,names=[],[],[]
			for line in f.readlines():
				line=line.split(' ')
				xi,yi,name=[li for li in line if li != '']
				x.append(np.float(xi))
				y.append(np.float(yi))
				names.append(name.split('\n')[0].replace("'",''))
			self.px=np.asarray(x)
			self.py=np.asarray(y)
			self.pnames=names
			
	def query_station_lon_lat(self,sation_name):	  
		""" returns station lon lat and index from pointlist """
		ind=self.pnames.index(sation_name)			
		return self.px[ind],self.py[ind],ind
		
	def plot_points(self,showname=False):	  
		plt.plot(self.px,self.py,'ko')
		if showname:
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,' '+ name,color='w')
		else:	
			count=0
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,str(count))		
				count+=1
				
	def open_point_output(self,file):	  
		""" initiate xarray netecd acces in ww3.ds """
		self.specds=xr.open_mfdataset(file)				

	def open_grid_output(self,file):	  
		""" initiate xarray netecd acces in ww3.ds """
		self.gridds=xr.open_mfdataset(file)
	
	def plot_grid_timeseries(self,parameter='hs',igrid=[0],iscoordinates=False):
		if iscoordinates==True:
			igrid=self.node_tree.query(igrid)[1]
		plt.plot(ww3.gridds['time'],ww3.gridds[parameter][:,igrid])
		plt.ylabel(parameter +' '+self.gridds[parameter].units)
		plt.gcf().autofmt_xdate()

	def get_grid_timeseries(self,parameter='hs',igrid=[0],iscoordinates=False):
		if iscoordinates==True:
			igrid=self.node_tree.query(igrid)[1]
		return ww3.gridds['time'].values,ww3.gridds[parameter][:,igrid].values
		
		
# open data		
ww3=WW3_mesh(ww3dir+ww3mesh)
ww3.load_points(ww3dir+ww3PointList)
ww3.open_point_output(ww3dir+ww3PointOutput)
ww3.open_grid_output(ww3dir+ww3GridOutput)


ww3.open_point_output(ww3PointOutput)

		

#ww3GridOutput='ww3.201702.nc'
#ww3.open_grid_output(ww3dir+ww3GridOutput)
#
#qx,qy,qind=ww3.query_station_lon_lat('Fino-1')
#igrid=ww3.find_nearest_node(qx,qy)
#plt.figure()
#ww3.plot_grid_timeseries( parameter='hs',igrid=igrid,iscoordinates=False)
#
#plt.figure()
#ti=-1
#ww3.plotAtnodes(ww3.gridds['hs'][ti,:])
#plt.title(str(ww3.gridds['time'][ti].values))
#
#ww3GridOutput='ww3.201703.nc'
#ww3.open_grid_output(ww3dir+ww3GridOutput)
#plt.figure()
#ww3.plot_grid_timeseries( parameter='hs',igrid=(-40,40),iscoordinates=True)


###### BSH ########
class bsh_spec:
	""" Read BSH Wave spectra *.spec files with format
	Line 1:  Buoy-Type Station Datum+Time
	Line 2:  up to 12 Parameters
    01 Peakfrequency (Hz),
	02 Energy Maximum (m*m*s),
	03 Dir at Peak (degN), 90 degrees is from east
	04 Spread at Peak (deg),
	05 Hs (m),
	06 Tp (s),
	07 Tm1 (s),
	08 Tm2 (s),
	09 Tm-1 (s),
	10 Nspc ,
	11 Hmax (m),
	12 T(Hmax) (s)
	% Line 3 - Nspc+3: f (Hz), e (m*m*s), dir (degN), spr (deg)
	
	Input:
	filename of spectral files as str
	or chronologically ordered list of files of bouy data:
	(bsh_spec(file=list(np.sort(glob.glob('<path>/*<station>*.spec')))   )
	
	Output class with fields:
	.name := station name as specified in file e.g. FN3
	.B.type :=  Buoy-Type - wave rider etc.
	.dates := list of dates of measurements
	.integral_header := header of integral parameters sotred in respective according numpy array (.integral_values)
	.integral_values := time times header item numpy array of integra parameters e.g. 
						bouy.integral_values[:,bouy.integral_header.index('Hs')] containing HS time series
						according to bouy.dates	
	.spectal_header := ['f', 'e', 'dir', 'spr'] header of spectral parmeters within .spectra
	.spectra := list with spectra occoring to .dates each item being an numpy array
				with spectral freuqencies in rows anbd colums according to spectal_header
	"""
	# not always all frequencies in spectra spectra
	def __init__(self,file):
	
		if type(file) == list:
			files=file
			file=files[0]
			nfiles=len(files)
		else:
			nfiles=1
			
		with open(file) as f:
			names=['Peakfrequency','Energy_Maximum','Dir_at_Peak','Spread_at_Peak',
			'Hs','Tp','Tm1','Tm2','Tm-1','Nspc','Hmax','T_Hmax']
			self.integral_header=names
			self.spectal_header=[ 'f', 'e' , 'dir', 'spr']
			lines=f.readlines()
			self.type,self.name,date=lines[0].rstrip('\n').split(' ')
			self.dates=[dt.datetime.strptime(date,'%Y%m%d%H%M')]
			self.i_nspec=names.index('Nspc')
			vals=[ float(val) for val in lines[1].rstrip('\n').split()]
			self.integral_values=[vals]
			count=2
			nspec=int(vals[self.i_nspec]) 
			
			spectrum=[[ float(item) for item in lines[count].rstrip('\n').split()]]
			for i in range(count+1,count+nspec):
				spectrum.append([ float(item) for item in lines[i].rstrip('\n').split()])
			self.spectra=[np.asarray(spectrum)]
			
			count+=nspec
			
			self.append_data(count,lines)

		if nfiles > 1:
			for file in files[1:]:
				#print(file)
				count=0
				try:
					with open(file) as f:
						self.append_data(count,f.readlines())
				except:
					print('error loading file '+ file)
					pass
		
		# convert lists to numpy arrays	
		self.integral_values=np.asarray(self.integral_values)
		
		
	def append_data(self,count,lines):
		while count < len(lines):
			date=lines[count].rstrip('\n').split(' ')[-1]
			self.dates.append(dt.datetime.strptime(date,'%Y%m%d%H%M'))
			
			count+=1 # read integral parameters
			self.integral_values.append([ float(val) for val in lines[count].rstrip('\n').split()])
			nspec=int(self.integral_values[-1][self.i_nspec]) 
			
			count+=1 # read spectra
			spectrum=[[ float(item) for item in lines[count].rstrip('\n').split()]]
			for i in range(count+1,count+nspec):
				spectrum.append([ float(item) for item in lines[i].rstrip('\n').split()])
			self.spectra.append(np.asarray([spectrum]))
			count+=nspec
		
	def get_parameter(self,parameter='Hs'):
		""" get buoy data """
		return self.integral_values[:,self.integral_header.index(parameter)]

# bash merege
"""
for folder in *
	do
	cd $folder
		 cat *_$folder_*.spec >> ${folder}_BSH.spec
	cd ..	
	done
"""

# load buoys
buoydir='/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/WAVE/bojen/'
subdirs=glob(buoydir+'*')
subdirs=np.sort([sdir[sdir.rindex('/')+1:] for sdir in subdirs])
bouys=dict.fromkeys(subdirs)
for subdir in subdirs:
	bouys[subdir]=bsh_spec(buoydir+'{:s}/{:s}_BSH.spec'.format(subdir,subdir))
B=bouys[subdirs[0]]
################
# Bojen
B.get_parameter('Hs')



### W W M ############
#### sites########

##### work around station output formatting issues
# error * in files.
files=glob('*.site')
sites={file.split('.')[0]:None for file in files}
for file in files:
	with open(file) as file_in:
		with open('prep_'+file,'w') as file_out:
			for nr,line in enumerate(file_in.readlines()):
				#if nr==99:
				#	break
				#line=line.replace('*','0').replace('Infinity','inf').replace('NaN            NaN','NaN')#.replace('NaN\tNaN','NaN') #.replace('NaN','-99')
				line=line.replace('*','0').replace('Infinity','').replace('NaN            NaN','NaN')#.replace('NaN\tNaN','NaN') #.replace('NaN','-99')
				file_out.write(line)

# eventually an extra in spaced
# occurance of double NaN ?

def load_wwm_site(file):
	with open(file) as f:
		lines=f.readline()
	header=lines.split()[:-1]
	try: 
		m=np.loadtxt(file,skiprows=1)
	except:
		m=np.loadtxt(file,skiprows=2)
	data=dict.fromkeys(header)
	for i,key in enumerate(header):
		data[key]=m[:,i]

	min=(data['TIME'])*100
	min=np.asarray(np.fix(100*(min-np.fix(min))),int)
	hour=(data['TIME'])
	hour=np.asarray(np.round(100*(hour-np.fix(hour))),int)
	day=(data['TIME']/100)
	day=np.asarray(np.fix(100*(day-np.fix(day))),int)
	mon=(data['TIME']/10000)
	mon=np.asarray(np.fix(100*(mon-np.fix(mon))),int)
	year=np.asarray(np.fix(data['TIME']/10000),int)
	data['dates']=dt.datetime(year[0],mon[0],day[0],hour[0],0,0)+np.arange(len(year))*dt.timedelta(hours=(np.diff(min[:2])/60+np.diff(hour[:2]))[0])

	return data

files=np.sort(glob('prep*.site'))
#files=glob('*.site')
sites={file.split('.')[0]:None for file in files}
for site,file in zip(sites.keys(),files):
	sites[site]=load_wwm_site(file)
# all values -999 why
#####################



################ compare Plotting ######################


 
#'Westerland',
#'Fino-3']


keysOBS=['ELB','FNO','FN3','HEL']
keysWWM=['prep_ELBE', 'prep_FINO1', 'prep_FINO3', 'prep_HELGON']
keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland']



#build dictionary
ww3Bouys=dict.fromkeys(keysWW3)
parameters=['hs',]
for key in keysWW3:
	ww3Bouys[key]=dict.fromkeys(['dates']+parameters)
	qx,qy,qind=ww3.query_station_lon_lat(key)
	for parameter in parameters:
		time,params=ww3.get_grid_timeseries(parameter=parameter,igrid=(qx,qy),iscoordinates=True)
		ww3Bouys[key][parameter]=params
	ww3Bouys[key]['dates']=time

query_station_lon_lat(self,sation_name)

pName0='hs'
pName='Hs'
pName2='HS'
unit='[m]'








def get_time_indices(date0,date1,datesOBS,datesMod):
	""" restrict to data range """
	# convert inputs to dateime 64 to avoid compatability issues
	date0,date1=np.asarray([date0,date1],np.datetime64)
	datesOBS=np.asarray(datesOBS,np.datetime64)
	datesMod=np.asarray(datesMod,np.datetime64)
	
	iuse= (date0 <= datesMod) & (date1 >= datesMod)
	nns=[]
	for i,date in enumerate(datesMod[iuse]):
		inn=np.argmin(np.abs(date-datesOBS))
		if np.abs(date-datesOBS)[inn]  < np.timedelta64(15,'m'):  #< dt.timedelta(hours=15):
			nns.append(inn)
		else:
			iuse[i]=False
	return nns,iuse

	
def interp_to_data_time(tdata,vdata,tmod,vmod,min_maxdate=None):
	""" interpolates inputs to common time of tdata
	tinterp,vdata_interp,vmod_interp=interp_to_data_time(tdata,vdata,tmod,vmod,min_maxdate=None)  """
	# convert inputs to dateime 64 to avoid compatability issues
	from scipy.interpolate import interp1d	
	tdata=np.asarray(tdata,np.datetime64)
	tmod=np.asarray(tmod,np.datetime64)

	if min_maxdate !=None:
		date0,date1=np.asarray(min_maxdate,np.datetime64)
	else:
		date0=[]
		date1=[]
	date0=np.max([tdata[0],tmod[0],date0])
	date1=np.min([tdata[-1],tmod[-1],date1])

	idata=[(tdata >= date0) & (tdata<= date1)]
	tdata_intp=tdata[idata]
	vdata_intp=vdata[idata]
	
	tin=(tmod-tmod[0])/(tmod[1]-tmod[0])
	tout=(tdata_intp-tmod[0])/(tmod[1]-tmod[0])

	fintp=interp1d(tin,vmod)
	vmod_intp=fintp(tout)

	return tdata_intp,vdata_intp,vmod_intp

	
def QQplot(dataObs,dataMod,date0,date1,stationname='',obsname='bouy',modname='SCHISM WWM',parameter='Hs',unit=' [m]'):	

	# stats
	entries=len(dataObs)
	meanR=dataObs.mean()
	meanM=dataMod.mean()
	SDR=np.std(dataObs)
	SDM=np.std(dataMod)
	RMSE=np.sqrt(np.mean((dataObs-dataMod)**2))
	E=dataMod-dataObs
	SE=np.sqrt(1/(entries-1)*np.sum((E-E.mean())**2))
	SI=SE/dataObs.mean()
	#np.sqrt(np.mean(((dataMod-meanR)-(dataObs-meanR).mean())**2))/dataObs.mean()
	#np.sqrt(np.sum(((dataMod-meanR)-(dataObs-meanR))**2)/dataObs.sum())
	bias=np.mean(dataObs-dataMod)
	cor=np.corrcoef(dataObs,dataMod)[0,1]

	txt=' Entries = {:d} \n Mean R ={:.3} {:s} \n Mean M ={:.3} {:s} \n SD R ={:.3} {:s} \n SD M ={:.3} {:s} \n RMSE ={:.3} {:s} \n SI ={:.3} \n bias ={:.3} {:s} \n CORR ={:.3}'.format(entries,meanR,unit,meanM,unit,SDR,unit,SDM,unit,RMSE,unit,SI,bias,unit,cor)


	coef = np.polyfit(dataObs,dataMod,1)
	poly1d_fn = np.poly1d(coef) ## poly1d_fn is now a function which takes in x and returns an estimate for y 
	lbl=str('{:.2f}x + {:.2f}'.format(coef[0],coef[1]))

	# scatter
	plt.clf()
	plt.scatter(dataObs,dataMod,s=0.4)
	plt.axis('square')
	plt.xlabel(obsname+ ' ' +parameter + unit)
	plt.ylabel(modname+ ' ' +parameter + unit)
	plt.title(stationname + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])
	ylim=plt.ylim()
	xlim=plt.xlim()
	plt.text(0,ylim[1]/1.8,txt)
	plt.plot([0,xlim[1]],[0,ylim[1]],'k')
	ph=plt.plot([0,xlim[1]], poly1d_fn([0,ylim[1]]), 'r')	
	plt.legend(ph,[lbl],loc='lower right',frameon=False)	

	
#####
for key1,key2,key3 in zip(keysOBS,keysWW3,keysWWM):
	key1,key2,key3

	B=bouys[key1]
	WW3=ww3Bouys[key2]	
	WWM=sites[key3]	
	
	datesWW3=WW3['dates']
	dataWW3=WW3[pName0]
	datesWWM=WWM['dates']
	dataWWM=WWM[pName2]

	B.dates=np.asarray(B.dates)
	dataB=B.get_parameter(parameter=pName)
	
	# nearest neighbours in time interval
	date0=np.maximum(B.dates[0],WWM['dates'][0])
	date1=np.minimum(B.dates[-1],WWM['dates'][-1])
	

	
	tdata,param_intp,dataWW3intp=interp_to_data_time(B.dates,dataB,WW3['dates'],WW3[pName0],min_maxdate=[date0,date1])
	QQplot(param_intp,dataWW3intp,date0,date1,obsname='bouy',modname='WW3',parameter='Hs',unit='[m]')
	plt.savefig('ww3_scatter_'+key1+'png',dpi=300)
	tdata,param_intp,dataWWMintp=interp_to_data_time(B.dates,dataB,WWM['dates'],WWM[pName2],min_maxdate=[date0,date1])
	QQplot(param_intp,dataWWMintp,date0,date1,obsname='bouy',modname='SCHISM WWM',parameter='Hs',unit='[m]')
	plt.savefig('wwm_900s_scatter_'+key1+'png',dpi=300)


	
	# time series
	plt.clf()
	plt.plot(B.dates,dataB)
	plt.plot(datesWWM,dataWWM,'-')
	plt.plot(datesWW3,dataWW3,'--')
	plt.xlim((date0,date1))
	plt.ylabel(pName+unit)
	plt.gcf().autofmt_xdate()
	plt.title(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])
	plt.legend(['buoy','SCHISM WWM dt=900','ww3'],loc='best')
	plt.grid()
	plt.savefig('wwm_900s_timeseries_'+key1+'png',dpi=300)