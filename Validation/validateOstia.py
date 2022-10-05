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
from matplotlib import path
# own and 3d party libraries
sys.path.insert(0,'/pf/g/g260114/Programs/python/scripts/')
sys.path.insert(0,'/work/gg0028/SCHISM/validation/lib/')
from matplotlib.path import Path
from schism import * # import schism functions
from techit import * # import latex script
from TaylorDiagram import * # import taylordiagram
import numpy as np
from data_and_model_classes import ostia

		
########## settings #################################
# directories (have to end with '/')
oceandir='/work/gg0028/g260099/AMM15/2017/'  					  # myocean
setupdir='/work/gg0028/g260114/RUNS/GermanBight/GermanBight/'     # schism run
ncdir=setupdir+'combined2/' #'combined_start_wrongly_1.1/'		  # schism run	
outdir='/work/gg0028/g260114/RUNS/GermanBight/GermanBight/ostia/'	      # image output folder
if not os.path.exists(outdir): os.mkdir(outdir) 
#plt.ion()
year=2017 # year of ostia data
ostiafile='/work/gg0028/g260114/Data/sat/Ostia/ostia_2017.nc'   # path of ostia file containong data for year


vmin=2
vmax=24
difflim=4
struct=[1,0,1]  # structured rid model

# coords for time series
lon=7.2
lat=54.5

cm=plt.cm.jet # colormap

dthours=12  # plot each 12 hours
ndays=30    # show time series for days
##############################################################################		
		

######### load setups   ##################################

# shcism
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
s.nc=Dataset(schismfiles[0])

# set reference time
reftime=dt.datetime.strptime(s.nc['time'].units[14:33],'%Y-%m-%d %H:%M:%S')#+dt.timedelta(days=1)
time=reftime+dt.timedelta(days=0.5)

# enhance schism for conveinet load
s.ncdir=ncdir		
s.reftime=reftime
s.file=s.ncdir+'schout_1.nc'
s.nc=Dataset(s.file)
s.dt=np.diff(s.nc['time'][:2])[0]
s.nt=len(s.nc['time'])
s.date=s.reftime+dt.timedelta(seconds=s.nc['time'][0]-s.dt)
def update(self,date):
        try:
            self.nc.close()
        except:
            pass
        (date-self.reftime)
        self.file=self.ncdir+'schout_'+str(int(np.ceil((date-self.reftime).total_seconds()/(self.dt*self.nt))))+'.nc'     
        self.nc=Dataset(self.file)
        self.date=self.reftime+dt.timedelta(seconds=self.nc['time'][0]-self.dt)
        self.time=self.date+dt.timedelta(seconds=self.dt)*np.arange(self.nt)          
def get_slab(self,time,varname,layer=-1):
		nn=np.argmin(np.abs(self.time-time))
		self.deltaT=self.time[nn]-time
		self.slab=self.nc[varname][nn,:,layer]
		self.ti=self.time[nn]
s.update=update
s.get_slab=get_slab


ostia=ostia(ostiafile)

############# load data
date=reftime+dt.timedelta(seconds=3600)
s.update(s,date)
#ostia.update(date)
s.get_slab(s,date,'temp')
ostia.get_slab(time=date)
# interp to structured grid

# horizontal next neighbours towards schism
distsi,nnsi=ostia.tree.query([ (xi,yi) for xi,yi in zip(np.asarray(s.lon),np.asarray(s.lat)) ]) # 
ostia.nns={'schism':nnsi}


######## make plots ###########################################
x=np.asarray(s.lon)
y=np.asarray(s.lat)


models=[ostia,s]
names=['ostia','SCHISM']
nmodels=len(models)


# derive boundary of schism model
bdnodes=[]
for ocean, land in zip(s.bdy_segments,s.land_segments):
	bdnodes+=ocean
	bdnodes+=land[1:]
bdnodes=np.asarray(bdnodes)-1
#plt.figure(),#plt.plot(x[bdnodes],y[bdnodes])
p=path.Path([(x[bdnodes][i],y[bdnodes][i]) for i in range(len(bdnodes))])
isin=p.contains_points(list(zip(ostia.LON.flatten(),ostia.LAT.flatten()))) # check Amm15 points within
isin2d=np.reshape(isin,ostia.LON.shape)  								   #schism domain
ostia.Mdiff=np.ma.masked_array(np.zeros(ostia.LON.shape),mask=(isin2d==False) |  ostia.slab.mask)
ostia.bias=np.ma.masked_array(np.zeros(ostia.LON.shape),mask=(isin2d==False) |  ostia.slab.mask)

###### nn to ostia #######################
coords=list(zip(ostia.LON.flatten()[isin],ostia.LAT.flatten()[isin]))
nn=s.nntree.query(coords)[1]
cx=x[s.nvplt].mean(axis=1)
cy=y[s.nvplt].mean(axis=1)
elemtree=cKDTree(list(zip(cx,cy)))
nn_element=elemtree.query(coords)[1]
###########################

nx=nmodels
ny=nx
# plot timeseries ######
nn1=ostia.tree.query((lon,lat))[1]
nn2=s.nntree.query((lon,lat))[1]


#time0=reftime+dt.timedelta(hours=1)
time1=reftime+dt.timedelta(hours=1)
times=[time1+i*dt.timedelta(hours=dthours) for i in range(365*int(24/dthours))]




varname='temp'
plt.close('all')
print('ploting '+ varname)

# data1 data2 run1 run2  run2-run1
fig, axes = plt.subplots(nrows=nmodels,ncols=nmodels)
plt.tight_layout()

outdir2=outdir+varname+'/'
if not os.path.exists(outdir2): os.mkdir(outdir2)

# plot values
phs=[]
for i,model in enumerate(models):
	i,model
	plt.subplot(ny,ny,i+1)
	plt.set_cmap('jet')
	if not 'schism' in str(type(model)):
		exec('ph{:d}=plt.pcolormesh(model.lon,model.lat,model.slab)'.format(i))
		plt.colorbar()
		model.get_slab(time=date)

	else:
		model.get_slab(model,time=date,varname=varname)
		exec('ph{:d},ch=s.plotAtelems(s.slab,cmap=cm)'.format(i))
	eval('phs.append(ph{:d})'.format(i))
	plt.title(names[i])
	plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
	plt.clim((vmin,vmax))

# plot differences
for i,model in enumerate(models[:int(nx/2)]):
	for j, model2 in enumerate(models):
		if model2 != model:
			plt.subplot(ny,ny,(i+1)*nx+j+1)
			plt.set_cmap('jet')
			if not 'schism' in str(type(model)):
				if not 'schism' in str(type(model2)):
					#exec('ph{:d}=plt.pcolormesh(model.lon,model.lat,model.T[0,:,:])'.format(i))
					#plt.colorbar()
					pass
				else:
					model.Mdiff[isin2d]=model2.slab[nn]
					mplot=np.ma.masked_invalid(model.Mdiff-model.slab)
					exec('ph{:d}=plt.pcolormesh(model.lon,model.lat,mplot,cmap=cm)'.format((i+1)*nx+j+1))
					#exec('ph{:d},ch=model2.plotAtelems(model2.slab-model.slab[model.nns[schism]],cmap=cm)'.format((i+1)*nx+j+1))
					eval('phs.append(ph{:d})'.format((i+1)*nx+j+1))
					plt.colorbar(extend='both')
			plt.plot(x[bdnodes],y[bdnodes],'k')
			plt.title(names[j]+' - '+names[i])
			plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
			plt.clim((-difflim,difflim))
			plt.gca().tick_params(labelbottom=False,labelleft=False)

###  update loop #######
ostia.Mdiff[isin2d]=s.slab[nn]
mplot=np.ma.masked_invalid(ostia.Mdiff-ostia.slab)
updata=[ostia.slab,s.slab,mplot]

# update labels         
plt.subplot(nx,ny,1)
ttl0=plt.title(names[0] +str(ostia.t))
plt.subplot(nx,ny,2)
ttl2=plt.title(names[1]+str(s.ti)[5:])

plt.subplot(nx,ny,4)
th1=plt.text(8.9,54.5,'bias: {:.2f} \n mae: {:.2f} '.format(np.nanmean(updata[2]),np.nanmean(np.abs(updata[2]))) )
ths=[th1]	#th1.remove()



plt.subplot(nx,ny,1)
plt.plot(lon,lat,'r+')
# update figures in loop


tdata=[]
ydata=[]
tdata2=[]
ydata2=[]
# plot timeseries
tdata.append(ostia.t)
tdata2.append(s.ti)
ydata.append(ostia.slab.flatten()[nn1])
ydata2.append(s.slab[nn2])	

plt.subplot(nx,ny,ny+1)
ax1 = plt.gca()
plt.ylim((vmin,vmax))

g1,=plt.plot(np.asarray(tdata),np.asarray(ydata))
g2,=plt.plot(np.asarray(tdata2),np.asarray(ydata2))

plt.grid()
plt.gcf().autofmt_xdate()
plt.ylabel(varname)
ax2 = ax1.twinx()
g3,=ax2.plot(np.asarray(tdata2),np.asarray(ydata2)-np.asarray(ydata),'k')
plt.ylim((-difflim,difflim))
plt.ylabel('difference')

lns = [g1,g2,g3]
#ax.legend(lns, labs, loc=0)
plt.legend(lns,['Ostia','SCHISM','diff'],loc='lower center',ncol=2,frameon=False)
xlim=ax1.get_xlim()
tdata=[]
ydata=[]
tdata2=[]
ydata2=[]

#################

for i,time in enumerate(times):
	t0=dt.datetime.now()
	print('doing plot'+str(i)+': '+str(time))

	#update data
	s.update(s,time)
	s.get_slab(s,time,varname)	
	ostia.get_slab(time=time)	

	t_nn=np.argmin(np.abs(s.time-time))
	idry=s.nc['wetdry_elem'][t_nn,:][s.nvplt2nvp]==1
	dryElem=(s.nc['wetdry_elem'][t_nn,:][s.nvplt2nvp][nn_element])==1
	slab=s.slab[nn]  # mask if node not member of wet elment
	slab[dryElem]=np.nan
	ostia.Mdiff[isin2d]=slab
	
	
	ostia.Mdiff[isin2d]=s.slab[nn]
	updata=[ostia.slab,s.slab,np.ma.masked_invalid(ostia.Mdiff-ostia.slab)]
	ostia.bias+=updata[-1]

	#update plot
	for phi,data,structi in zip(phs,updata,struct):
		if structi==1:
			phi.set_array(data[:-1,:-1].flatten())
		else:
			#phi.set_array(data[s.nvplt[:,:3]].mean(axis=1))
			phi.set_array(np.ma.masked_array(data[s.nvplt[:,:3]].mean(axis=1),mask=idry))
	for ii,thi in enumerate(ths):
		thi.set_text('bias: {:.2f} \n mae: {:.2f} '.format(updata[nx+ii].mean(),np.abs(updata[nx+ii]).mean()) )

	# update labels         
	plt.subplot(ny,nx,1)
	plt.title(names[0]+str(ostia.t))
	plt.subplot(ny,nx,2)
	plt.title(names[1]+str(s.ti)[5:])

	# plot timeseries
	tdata.append(ostia.t)
	tdata2.append(s.ti)
	ydata.append(ostia.slab.flatten()[nn1])
	ydata2.append(s.slab[nn2])	

	#plt.subplot(nx,ny,3)  subplot seems to kill it
	g1.set_xdata(np.asarray(tdata))
	g1.set_ydata(np.asarray(ydata))
	g2.set_xdata(np.asarray(tdata2))
	g2.set_ydata(np.asarray(ydata2))
	g3.set_xdata(np.asarray(tdata))
	deltay=np.asarray(ydata2)-np.asarray(ydata)
	g3.set_ydata(deltay)


	ax1.set_xlim(time-dt.timedelta(days=ndays),time+dt.timedelta(days=0.5))
	ax2.set_xlim(time-dt.timedelta(days=ndays),time+dt.timedelta(days=0.5))
	ax1.set_ylim( np.floor(np.min((np.min(ydata),np.min(ydata2)))), np.ceil(np.max((np.max(ydata),np.max(ydata2)))))
	dval=np.ceil(np.max(np.abs(deltay)))
	ax2.set_ylim((-dval,dval))

	#plt.xlim()
	#ax1.set_xlim((xlim[0],xlim[0]+(tdata[-1]-tdata[0]).total_seconds()/86400+0.5))
	ax1.set_xlim(tdata[0]-dt.timedelta(days=1),tdata[-1]+dt.timedelta(days=1))
	ax2.set_xlim(tdata[0]-dt.timedelta(days=1),tdata[-1]+dt.timedelta(days=1))
	print(str(tdata[-1]))
	#ax2.set_xlim((xlim[0],xlim[0]+(tdata[-1]-tdata[0]).total_seconds()/86400+0.5))
	#plt.gcf().autofmt_xdate()


	#fig.canvas.draw()
	plt.savefig(outdir2+'{:04d}_comp_ostia'.format(i)+'_'+varname,dpi=400)
	print('took '+str((dt.datetime.now()-t0).total_seconds())+' s')


ostia.bias/=len(times)
plt.close('all')
plt.pcolormesh(ostia.lon,ostia.lat,ostia.bias)
plt.plot(x[bdnodes],y[bdnodes],'k')
ch=plt.colorbar()
ch.set_label('bias')
plt.axis((x.min()-0.1,x.max()+0.1,y.min()-0.1,y.max()+0.1 ))
plt.title(names[2]+' - '+names[1]+'bias '+str(times[0])+'-'+str(times[-1]) )
plt.titght_layout()
plt.savefig(outdir2+'bias_'+varname,dpi=400)
