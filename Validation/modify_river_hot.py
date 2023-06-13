import numpy as np	
import datetime as dt	
import glob
# mark on schism map	
import sys 	
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities')
from schism import* # import schism class to read grid structure	
from matplotlib.path import Path	
from scipy.spatial import cKDTree
from netCDF4 import Dataset



ncfilein='GB_hot_AMM720111101.nc'
ncfileout='GB_hot_AMM720111101_correct_river.nc'


s=schism_setup()
if not (os.path.exists('rivers.reg')):
	exit()

riverPoly=np.loadtxt('river.reg',skiprows=3)
riverPoly=path.Path(list(zip(riverPoly[:,0],riverPoly[:,1])))
setsalt0=riverPoly.contains_points(list(zip(s.x,s.y)))


# read file
nc=Dataset(ncfilein)
tr_nd=nc['tr_nd'][:]
tr_nd[setsalt0,:,1]=0
nc.close()
os.rename(ncfilein,ncfilein+'bck_up')
s.write_hotstart(tr_nd,filename=ncfilein)

#check sohle und tahl

def get_from_tideporal_txt(file,headeronly=False):
	""" load tide portal data, and return as utc array"""
    #import pyproj	

	with open(file)	as f:
		#read header
		header={}
		line=''
		while '====' not in line:
			line=f.readline()
			if '====' not in line:
				split=line.split(':')
				header[split[0].strip()]=''.join(split[1:]).strip()
		#read data	
		line=f.readline() #skip header
		varname=header['Messwert'].split()[0]
		station={}
		station['date']=[]
		station[varname]=['varname']
		station['name']=header['Stationsname'].replace(' ','')
		if '(' in station['name']:
			station['name']=station['name'][:station['name'].index('(')]
		
		try:			
			header['East']=float(header['East'])		
			header['North']=float(header['North'])
		except:
			print('no easting weasting in header')

		station['header']=header
		station['East']=float(header['East'])
		station['North']=float(header['North'])
			
		if headeronly:
		
			return station
		
		while '====' not in line:
			try:
				line=f.readline() #ignore defect line
				if '====' not in line:
					a,b=line[:-1].split('\t')[:2]
					if len(b) == 0:
						b='nan'
					station['date'].append( dt.datetime.strptime(a,'%Y-%m-%d %H:%M:%S'))
					station[varname].append(float(b))
			except:
				continue
	
		
	if 'Zeitzone' in list(header.keys()):
		time=header['Zeitzone']
		if 'UTC' in time:
			utc=float(time[time.index('UTC')+3:-1])
		else:
			utc=0 # wrong if MEZ 
	else:		
		utc=0
	station['date']=np.asarray(station['date'])+dt.timedelta(hours=-utc)
	station[varname]=np.asarray(station[varname])	
	return station



# sort in folders by river	
files=glob.glob('*.txt')
import os
for folder in ['ems','weser','elbe']:
	if not os.path.isdir(folder):
		os.mkdir(folder)
		
for file in files:
	station=get_from_tideporal_txt(file,headeronly=True)
	zipfile=file[:-4]+'.zip'
	print(file)
	for river in ['ems','weser','elbe']:
		if river in station['header']['Gewässer'].lower():
		
			varname=station['header']['Messwert'].split()[0]
			if not os.path.isdir(river+'/'+varname):
				os.mkdir(river+'/'+varname)
			os.rename(file,river+'/'+varname+'/'+file)
			os.rename(zipfile,river+'/'+varname+'/'+zipfile)
			break

temp_stations={}	
salt_stations={}			
topdir1='/gpfs/work/jacobb/data/RUNS/GermanBight/GB_2013_xaver/download/elbe/'	
topdir2='/gpfs/work/jacobb/data/RUNS/GermanBight/GB_2013_xaver/download/weser/'	
for topdir in [topdir1,topdir2]:
	vardirs=glob.glob(topdir+'*')	

	# load stations for river	

	for vardir in vardirs:	
		files=glob.glob(vardir+'/*.txt')	
		if 'temp' in vardir:
			for file in files:
				print(file)
				stat=get_from_tideporal_txt(file)
				temp_stations[stat['name']]=stat	
		else: #salt
			for file in files:
				print(file)
				stat=get_from_tideporal_txt(file)
				salt_stations[stat['name']]=stat	

		

plt.ion()
s=schism_setup()


def load_reg(file):
	""" load region file from xmgredi"""
	return np.loadtxt(file,skiprows=3)

#elbe=load_reg('/gpfs/work/jacobb/data/RUNS/GermanBight/GB_2013_xaver/download/elbe/elbe.reg')
regions=['/gpfs/work/jacobb/data/RUNS/GermanBight/GB_2013_xaver/download/elbe/elbe.reg','/gpfs/work/jacobb/data/RUNS/GermanBight/GB_2013_xaver/download/weser/weser.reg'] # file regions for river


coords=list(zip(s.x,s.y))
s.x,s.y=np.asarray(s.x),np.asarray(s.y)	

# load values for date
date=dt.datetime(2013,11,1) # referecne date
# get values

tempinterp=np.zeros(s.nnodes)# indices for nn region
 # how many neighbouzr to ocnsider

tempinterp,tempnames=set_river_values(temp_stations,do_plot=True)
saltinterp,tempnames=set_river_values(salt_stations,do_plot=True)

def set_river_values(stations,npoints=4,do_plot=False):
	""" perform inverse distance value interpoaltion from npoint closest measurement stations"""

	k=0
	nstat=len(stations)
	x,y=np.zeros(nstat),np.zeros(nstat)
	#salt,temp=np.zeros(nstat),np.zeros(nstat)
	temp=np.zeros(nstat)
	#temp_names=[]
	#salt_names=[]
	names=[]


	for i,station in enumerate(stations.values()):
		if  len(station['date'])>0:
			tnn=np.argmin(np.abs(station['date']-date)) #next time step
			
			#plt.plot(station['East'],station['North'],'kx')
			consonants=['A','E','I','O','U','Ö','Ü']
			name=station['header']['Stationsname']
			for con in consonants:
				name=name.replace(con,'')	
			temp_names.append(name)	
			x[k],y[k]=station['East'],station['North']	
			temp[k]=station[station['header']['Messwert'].split()[0]][tnn]	
			k+=1
	x=x[:k]
	y=y[:k]
	temp=temp[:k]

	schismpoints={}
	datapoints={}
	for ifile,file in enumerate(regions):
		tag=file.split('/')[-1].split('.')[0]
		region=load_reg(file)
		areaPoly=Path(list(zip(region[:,0],region[:,1])))
		schismpoints[tag]=areaPoly.contains_points(coords)
		datapoints[tag]=areaPoly.contains_points(list(zip(x,y))) &  ~ np.isnan(temp)
	####################### #######################

		xriver=x[datapoints[tag]]
		yriver=y[datapoints[tag]]
		tempriver=temp[datapoints[tag]]
		nntree=cKDTree(list(zip(xriver,yriver))) # next neighbour search tree	
		
		for key in datapoints.keys():
			for nde in np.where(schismpoints[tag])[0]:
				nde
				riverinds=nntree.query(coords[nde],npoints)[1]
				dist=np.sqrt((s.x[riverinds]-xriver[riverinds])**2+(s.y[riverinds]-yriver[riverinds])**2)
				tempinterp[nde]= ((1/dist/(1/dist).sum())*tempriver[riverinds]).sum()
	
	if do_plot:
		plt.figure()
		ph0,ch,ax=s.plotAtnodes(tempinterp,latlon=False)
		ph=ax.scatter(x[datapoints[tag]],y[datapoints[tag]],s=15,c=temp[datapoints[tag]],edgecolor='k',clim=ph0.get_clim())
	
	return tempinterp,names

# file in data derived values	
tempinterp
saltinterp
tempind=tempinterp!=0
saltind=saltinterp!=0
for zi in range(tr_nd.shape[1]):
	tr_nd[tempind,zi,0]=tempinterp[tempind]
	tr_nd[saltind,zi,1]=saltinterp[saltind]

s.write_hotstart(tr_nd,filename=ncfileout)


# plot
suffix=['surface', 'bottom' ]
inds=[-1, 0]
# plot surf bottmom layer
for i, lab in zip(inds,suffix): 
	ph,ch,ax=s.plotAtnodes(tr_nd[:,i,1])
	ch.set_label('salinity g/kg')
	plt.title(lab)
	plt.savefig('salinity_correct'+lab,dpi=600)
	plt.close()

for i, lab in zip(inds,suffix): 
	ph,ch,ax=s.plotAtnodes(tr_nd[:,i,0])
	ch.set_label('temp deg C')
	plt.title(lab)
	plt.savefig('temp_correct'+lab,dpi=600)
	plt.close()
	