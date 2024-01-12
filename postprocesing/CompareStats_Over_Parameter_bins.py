"""
Compare Values of binned ranges
"""


from glob import glob
import os
import sys
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
from schism import* # import schism class to read grid structure
from matplotlib.path import Path   # polygon operations
plt.ion()
# parallel stuff
#%matplotlib inline
import xarray as xr
import matplotlib.pyplot as plt

from dask.diagnostics import ProgressBar
from dask_jobqueue import SLURMCluster
from distributed import Client, progress
from distributed.utils import tmpfile

import dask
import distributed
dask.__version__, distributed.__version__


setupdir=['/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/Veg_REF/',]
#####read grid_strtucture
os.chdir(setupdir[0])
s=schism_setup()

		
		
ds1=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/AnaMon10_quantile_storm/Ana_mveg_2_2017_10storm.nc')

ds2=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/AnaMon10_quantile_calm/Ana_mveg_2_2017_10calm.nc')



# netcdfdirs from experiments
ref='Veg_REF'
control='Veg_CNTRL'
experiment='Veg_max'
experiments=[ref,control,'Veg_max','Veg_HE','Veg_LE'] #,experiment]

expnames=['Ref','Blank','Veg$_{max}$','Veg$_{HE}$','Veg$_{LE}$']	
nfws=(8.03, 9.01, 54.32902335808929, 55.55471643942857)	
efws=(6.55, 8.45, 53.2, 53.83566262094643)
meadow=(6.64, 6.95,53.5, 53.62)
meadow_zoom=(6.842150912936795, 6.928917980826309,53.499620667175584, 53.55326546740791)	
	
domain_tag=['total','efws','nfws','meadow','meadow_zoom']
domains=domain_tag[1:]		

bins=np.arange(-4,10.5,.5)
lbin=bins[:-1]
ubin=bins[1:]
cbins=(lbin+ubin)/2
Dbins=bins

lon,lat=np.asarray(s.lon),np.asarray(s.lat)
D=np.asarray(s.depths)

for domain,axis_limit in zip(domains,[efws,nfws,meadow]):
	inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= axis_limit[3]) 
					
subsel=inreg
					
					
# by depths distributions	 # vergleiche stem density	
dirs=['Veg_REF','Veg_CNTRL','Veg_max','Veg_HE','Veg_LE']
setupdir=['/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/'+dir+'/' for dir in dirs]
cwd=os.getcwd()
gr3s=dict.fromkeys(experiments)
for exp,dir in zip(experiments,setupdir):
	os.chdir(dir)
	s.read_gr3('sav_N.gr3')
	gr3s[exp]=s.gr3['sav_N'].copy()
gr3s['Veg_CNTRL']*=0
					
					
# subdomain pre selected	
def bin_select(experiments,D,Dbins,data,gr3s,iref=0,inreg=None):
	""" bin data of different experiments,for aossosiated varaible D (e.g. depths bins) seperatly for cretarion gr3s. And optional preselector inreg (e.g. falste/true for region subsets) """

	lbins=Dbins[:-1]
	ubins=Dbins[1:]
	cbins=0.5*(lbins+ubins)

	if len(inreg)==1 and inreg== None:
		inreg=np.arange(len(D))
	D=D[inreg]
	truebins=np.zeros(len(cbins))
	falsebins=np.zeros(len(cbins))
	#truebinscntrl=np.zeros(len(cbins))
	truedeltas=np.zeros(len(cbins))
	falsedeltas=np.zeros(len(cbins))
	
	truedeltas_rel=np.zeros(len(cbins))
	falsedeltas_rel=np.zeros(len(cbins))

	#falsebinscntrl=np.zeros(len(cbins))	
	
	exptruebins=dict.fromkeys(experiments)
	expfalsebins=dict.fromkeys(experiments)
	
	exptruedeltas=dict.fromkeys(experiments)
	expfalsedeltas=dict.fromkeys(experiments)

	exptruedeltas_rel=dict.fromkeys(experiments)
	expfalsedeltas_rel=dict.fromkeys(experiments)
	
	
	refata=data[iref][inreg].copy()
	#from IPython import embed; embed()
	for ind in range(len(experiments)): #expinds:
	
		truebins=np.zeros(len(cbins))
		falsebins=np.zeros(len(cbins))
	
				
		indata=data[ind][inreg].copy()		
		if ind==1:
			ind=0
			itrue=gr3s[experiments[ind]][inreg]>0
			ifalse=gr3s[experiments[ind]][inreg]==0				
			ind=1
		else:	
			itrue=gr3s[experiments[ind]][inreg]>0
			ifalse=gr3s[experiments[ind]][inreg]==0				
	
	
			#truedata=delta[itrue]	
			
		print(ind)	
		for ibin in range(len(cbins)):
			din_false=(lbin[ibin] <= D[ifalse]) & (D[ifalse] < ubin[ibin])
			din_true=(lbin[ibin] <= D[itrue]) & (D[itrue] < ubin[ibin])
			if din_true.sum()>0:
				truebins[ibin]=np.nanmean(indata[itrue][din_true])
				truedeltas[ibin]=np.nanmean((indata[itrue]-refata[itrue])[din_true])
				truedeltas_rel[ibin]=truedeltas[ibin]/np.nanmean(refata[itrue][din_true])
				#truebinscntrl[ibin]=np.mean(Vcntrl[inreg][itrue][din_true])
				#np.nanmean(deltagrass[din_true])
			else:
				truebins[ibin]=np.nan
				#truebinscntrl[ibin]=np.nan
			if din_false.sum()>0:
				falsebins[ibin]=np.nanmean(indata[ifalse][din_false]) #deltafree[din_false]
				falsedeltas[ibin]=np.nanmean((indata[ifalse]-refata[ifalse])[din_false])
				falsedeltas_rel[ibin]=falsedeltas[ibin]/np.nanmean(refata[ifalse][din_false])
				#falsebinscntrl[ibin]=np.mean(Vcntrl[inreg][ifalse][din_false]) #deltafree[din_false]
			else:
				falsebins[ibin]=np.nan
				#falsebinscntrl[ibin]=np.nan	
		
		exptruebins[experiments[ind]]=truebins.copy()
		expfalsebins[experiments[ind]]=falsebins.copy()			
		exptruedeltas[experiments[ind]]=truedeltas.copy()   # keeps overwriting need copy
		expfalsedeltas[experiments[ind]]=falsedeltas.copy()			
		exptruedeltas_rel[experiments[ind]]=truedeltas_rel.copy()   # keeps overwriting need copy
		expfalsedeltas_rel[experiments[ind]]=falsedeltas_rel.copy()			
		
		#from IPython import embed; embed()
	return cbins,exptruebins,expfalsebins,exptruedeltas,expfalsedeltas,exptruedeltas_rel,expfalsedeltas_rel

def bin_scatter(experiments,D,Dbins,data,gr3s,iref=0,inreg=None):
	""" bin data of different experiments,for aossosiated varaible D (e.g. depths bins) seperatly for cretarion gr3s. And optional preselector inreg (e.g. falste/true for region subsets) """

	lbins=Dbins[:-1]
	ubins=Dbins[1:]
	cbins=0.5*(lbins+ubins)

	if len(inreg)==1 and inreg== None:
		inreg=np.arange(len(D))
	D=D[inreg]
	truebins=np.zeros(len(cbins))
	falsebins=np.zeros(len(cbins))
	#truebinscntrl=np.zeros(len(cbins))
	truedeltas=np.zeros(len(cbins))
	falsedeltas=np.zeros(len(cbins))
	
	truedeltas_rel=np.zeros(len(cbins))
	falsedeltas_rel=np.zeros(len(cbins))

	#falsebinscntrl=np.zeros(len(cbins))	
	
	exptruebins=dict.fromkeys(experiments)
	expfalsebins=dict.fromkeys(experiments)
	
	exptruedeltas=dict.fromkeys(experiments)
	expfalsedeltas=dict.fromkeys(experiments)

	exptruedeltas_rel=dict.fromkeys(experiments)
	expfalsedeltas_rel=dict.fromkeys(experiments)
	
	
	refata=data[iref][inreg].copy()
	#from IPython import embed; embed()
	for ind in range(len(experiments)): #expinds:
	
		truebins=np.zeros(len(cbins))
		falsebins=np.zeros(len(cbins))
	
				
		indata=data[ind][inreg].copy()		
		if ind==1:
			ind=0
			itrue=gr3s[experiments[ind]][inreg]>0
			ifalse=gr3s[experiments[ind]][inreg]==0				
			ind=1
		else:	
			itrue=gr3s[experiments[ind]][inreg]>0
			ifalse=gr3s[experiments[ind]][inreg]==0				
	
	
			#truedata=delta[itrue]	
			
		print(ind)	
		for ibin in range(len(cbins)):
			din_false=(lbin[ibin] <= D[ifalse]) & (D[ifalse] < ubin[ibin])
			din_true=(lbin[ibin] <= D[itrue]) & (D[itrue] < ubin[ibin])
			if din_true.sum()>0:
				truebins[ibin]=np.nanmean(indata[itrue][din_true])
				truedeltas[ibin]=np.nanmean((indata[itrue]-refata[itrue])[din_true])
				truedeltas_rel[ibin]=truedeltas[ibin]/np.nanmean(refata[itrue][din_true])
				#truebinscntrl[ibin]=np.mean(Vcntrl[inreg][itrue][din_true])
				#np.nanmean(deltagrass[din_true])
			else:
				truebins[ibin]=np.nan
				#truebinscntrl[ibin]=np.nan
			if din_false.sum()>0:
				falsebins[ibin]=np.nanmean(indata[ifalse][din_false]) #deltafree[din_false]
				falsedeltas[ibin]=np.nanmean((indata[ifalse]-refata[ifalse])[din_false])
				falsedeltas_rel[ibin]=falsedeltas[ibin]/np.nanmean(refata[ifalse][din_false])
				#falsebinscntrl[ibin]=np.mean(Vcntrl[inreg][ifalse][din_false]) #deltafree[din_false]
			else:
				falsebins[ibin]=np.nan
				#falsebinscntrl[ibin]=np.nan	
		
		exptruebins[experiments[ind]]=truebins.copy()
		expfalsebins[experiments[ind]]=falsebins.copy()			
		exptruedeltas[experiments[ind]]=truedeltas.copy()   # keeps overwriting need copy
		expfalsedeltas[experiments[ind]]=falsedeltas.copy()			
		exptruedeltas_rel[experiments[ind]]=truedeltas_rel.copy()   # keeps overwriting need copy
		expfalsedeltas_rel[experiments[ind]]=falsedeltas_rel.copy()			
		
		#from IPython import embed; embed()
	return cbins,exptruebins,expfalsebins,exptruedeltas,expfalsedeltas,exptruedeltas_rel,expfalsedeltas_rel


image_outdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/mon10_vs_REF/bar_plots_storm4/'




##
domain_tag=['total','efws','nfws','meadow','meadow_zoom']
level_tag=['_bottom','_surface']
varnames=['elevation','sigWaveHeight','turbulentKineticEner','depthAverageVel']
units=['m','m','W','m/s']
symbols=['$\zeta$','HS','TKE','$u_{<z>}$']

varnames+=['bottomStress','totalSuspendedLoad_surface','totalSuspendedLoad_bottom','turbulentKineticEner_surface','turbulentKineticEner_bottom','meanWavePeriod','horizontalVel_surface']
units+=['Pa','g/L','g/L','m$^2s^{-1}$','m$^2s^{-1}$','s','m/s']
symbols+=['$\\tau_{btm}$','$C_{totsurf}$','$C_{totbtm}$','$tke_{surf}$','$tke_{btm}$','T_{mean}','$V_{surf}$']		




varnames_use=['elevation','sigWaveHeight','bottomStress','totalSuspendedLoad_bottom','turbulentKineticEner_surface','turbulentKineticEner_bottom','depthAverageVel']



ref='Veg_REF'
control='Veg_CNTRL'
real=True
if real:
	#experiment='Veg_max'
	experiments=[ref,control,'Veg_max','Veg_HE','Veg_LE'] #,experiment]

	expnames=['Ref','Blank','Veg$_{max}$','Veg$_{HE}$','Veg$_{LE}$']	
else:
	experiments=[control,'Veg_REF','Veg_max','Veg_HE','Veg_LE'] #,experiment]
	expnames=['Blank','Ref','Veg$_{max}$','Veg$_{HE}$','Veg$_{LE}$']	


#varnames_use=varnames2
#
#
#fig,axs,cbar_ax,phs,qhs=compare_slabs_diff2(s,data1,clim0='auto',names=expnames,cmap=cmap,nrows=nrows,ncols=ncols,cblabel='$\Delta$'+label,axis_limit=axis_limit,clim=clim,figsize=figsize,cb_orientation=cb_orientation,relative=relative,nn=nn)
#
#plt.figure()
#ph,ch,ax=s.plotAtnodes(data1[2]-data1[0],cmap=plt.cm.jet)
#ph.set_clim((-.1,.1))
#
#cbins,exptruebins1,expfalsebins1=bin_select(experiments,D,Dbins,data1,gr3s,iref=0,inreg=inreg)	
#
#plt.figure()
#plt.clf()
#plt.scatter(D[inreg],(data1[1]-data1[0])[inreg])



def make_scatter(D,data,dataref,isegrass,inosegrass,varname,relative=False):
	if relative:
		delta=(data-dataref)/np.abs(dataref)*100
	else:
		delta=data-dataref
		
	inan=np.isnan(delta) | (delta >1e4)				
	delta=np.ma.masked_array(delta,mask=inan)
		
	deltafree=delta[inosegrass]
	deltagrass=delta[isegrass]
	plt.scatter(D[inosegrass],deltafree,marker='.',color='brown',s=0.4) 
	#isort=np.argsort(D[inosegrass])
	plt.scatter(D[isegrass],deltagrass,marker='.',color='green',s=0.4) 
	plt.legend(['free (avg:{:.4f})'.format( np.round(np.nanmean(deltafree),4)),'grass (avg:{:.4f})'.format( np.round(np.nanmean(deltagrass),4))],frameon=False)
	plt.ylabel('$\Delta$ ' +varname)
	if relative:
		plt.ylim((-100,100))
		if varname == 'elevation':
			plt.ylim((-30,30))	
	if irow==3:
		plt.xlabel('depth [m]')
	plt.grid()


def make_ratio_scatter(D,data,dataref,isegrass,inosegrass,varname,relative=False):
		
	ratio=np.abs(data/dataref)
	#inan=np.isnan(delta) | (delta >1e4)				
	#delta=np.ma.masked_array(delta,mask=inan)
	
	iinv=ratio<1 
	ratio[iinv]=-1/ratio[iinv]
	ratiofree=ratio[inosegrass]
	ratiograss=ratio[isegrass]
	plt.scatter(D[inosegrass],ratiofree,marker='.',color='brown',s=0.4) 
	#isort=np.argsort(D[inosegrass])
	plt.scatter(D[isegrass],ratiograss,marker='.',color='green',s=0.4) 
	plt.legend(['free (avg:{:.4f})'.format( np.round(np.nanmean(ratiofree),4)),'grass (avg:{:.4f})'.format( np.round(np.nanmean(ratiograss),4))],frameon=False)
	plt.ylabel('R ' +varname)
	if relative:
		plt.ylim((-100,100))
		if varname == 'elevation':
			plt.ylim((-30,30))	
	if irow==3:
		plt.xlabel('depth [m]')
	plt.grid()
	
#scatter plots


def set_label(ind):
	if ind+1< nexp:
		plt.xticks([])
	else:
		plt.xlabel('Depth [m]')
	plt.xlim((-4,6))
ncols=4

image_outdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/scatter_plots_storm4_refblank/'
if not os.path.isdir(image_outdir): os.mkdir(image_outdir)	
nexp=len(experiments)
iref=1
for i in range(len(varnames)):  #[:2]
	#if i < 2:
	#	continue
	vari=varnames[i]
	if vari not in varnames_use:
		print('skipping')
		continue
		
	print(vari)
	unit=units[i]
	for j,quanti in enumerate(['mean','quantile95']): #enumerate(quantaties):
		print(quanti)
		symbol=symbols[i]
		if j==[0]:
			symbol=symbol.join(('<','>'))
		else:	
			symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
		label0='{:s} [{:s}]'.format(symbol,unit)
		try: #eventually variable not saved
			data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
			data2=[ds2[expi+'_'+vari+'_'+quanti].values for expi in experiments]
		except:
			try: #eventually variable not saved
				try: #eventually variable not saved
					data1=[ds1[quanti+'_'+expi+'_'+vari].values for expi in experiments]
					data2=[ds2[quanti+'_'+expi+'_'+vari].values for expi in experiments]
				except:
					data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
					data2=[ds2[expi+'_'+vari+'_'+quanti].values for expi in experiments]
			except:
				try: #eventually variable not saved
					vari1,vari2=vari.split('_')
					data1=[ds1[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
					data2=[ds2[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
				except:
					continue

		data1=[np.ma.masked_array(datai,mask=drymask) for datai in data1]
		data2=[np.ma.masked_array(datai,mask=drymask) for datai in data2]
		label=label0

		
		
		for domain,axis_limit in zip(domains,[efws,nfws,meadow]):
			inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= axis_limit[3]) 
			plt.close('all')
			Dreg=D[inreg]
			dataref1=data1[iref][inreg]
			dataref2=data2[iref][inreg]
			plt.figure(figsize=(20,16))	
			count=-1
			
			for ind,expname in enumerate(experiments):
				if ind != iref:
					count+=1
					if ind==1:
						ind=0
						isegrass=gr3s[experiments[ind]][inreg]>0
						inosegrass=gr3s[experiments[ind]][inreg]==0				
						ind=1
					else:	
						isegrass=gr3s[experiments[ind]][inreg]>0
						inosegrass=gr3s[experiments[ind]][inreg]==0				
					
					#plt.subplot(4,2,(count+1)*2)
					plt.subplot(4,ncols,count*ncols+1)
					make_scatter(Dreg,data1[ind][inreg],dataref1,isegrass,inosegrass,varname=symbol+' '+unit,relative=False)
					set_label(ind)															
					plt.title('storm ' + expname)
					plt.subplot(4,ncols,count*ncols+2)
					make_scatter(Dreg,data1[ind][inreg],dataref1,isegrass,inosegrass,varname=symbol+' %',relative=True)
					set_label(ind)															
					plt.title('storm ' + expname)
					
					plt.subplot(4,ncols,count*ncols+3)
					make_scatter(Dreg,data2[ind][inreg],dataref2,isegrass,inosegrass,varname=symbol+' '+unit,relative=False)
					set_label(ind)
					plt.title('calm ' + expname)
					
					plt.subplot(4,ncols,(count+1)*ncols)
					make_scatter(Dreg,data2[ind][inreg],dataref2,isegrass,inosegrass,varname=symbol+' %',relative=True)
					set_label(ind)
					plt.title('calm ' + expname)
			plt.tight_layout()
			plt.suptitle(domain)
			plt.savefig(image_outdir+'_'.join((domain,vari,quanti))+'.png',dpi=300)


# ratios
for i in range(len(varnames)):  #[:2]
	#if i < 2:
	#	continue
	vari=varnames[i]
	if vari not in varnames_use:
		print('skipping')
		continue
		
	print(vari)
	unit=units[i]
	for j,quanti in enumerate(['mean','quantile95']): #enumerate(quantaties):
		print(quanti)
		symbol=symbols[i]
		if j==[0]:
			symbol=symbol.join(('<','>'))
		else:	
			symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
		label0='{:s} [{:s}]'.format(symbol,unit)
		try: #eventually variable not saved
			data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
			data2=[ds2[expi+'_'+vari+'_'+quanti].values for expi in experiments]
		except:
			try: #eventually variable not saved
				try: #eventually variable not saved
					data1=[ds1[quanti+'_'+expi+'_'+vari].values for expi in experiments]
					data2=[ds2[quanti+'_'+expi+'_'+vari].values for expi in experiments]
				except:
					data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
					data2=[ds2[expi+'_'+vari+'_'+quanti].values for expi in experiments]
			except:
				try: #eventually variable not saved
					vari1,vari2=vari.split('_')
					data1=[ds1[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
					data2=[ds2[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
				except:
					continue

		data1=[np.ma.masked_array(datai,mask=drymask) for datai in data1]
		data2=[np.ma.masked_array(datai,mask=drymask) for datai in data2]
		label=label0

		
		
		for domain,axis_limit in zip(domains,[efws,nfws,meadow]):
			inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= axis_limit[3]) 
			plt.close('all')
			Dreg=D[inreg]
			dataref1=data1[iref][inreg]
			dataref2=data2[iref][inreg]
			plt.figure(figsize=(20,16))	
			count=-1
			
			for ind,expname in enumerate(experiments):
				if ind != iref:
					count+=1
					if ind==1:
						ind=0
						isegrass=gr3s[experiments[ind]][inreg]>0
						inosegrass=gr3s[experiments[ind]][inreg]==0				
						ind=1
					else:	
						isegrass=gr3s[experiments[ind]][inreg]>0
						inosegrass=gr3s[experiments[ind]][inreg]==0				
					
					#plt.subplot(4,2,(count+1)*2)
					plt.subplot(4,ncols,count*ncols+1)
					make_ratio_scatter(Dreg,data1[ind][inreg],dataref1,isegrass,inosegrass,varname=symbol+' '+unit,relative=False)
					#plt.ylim((-10,10))
					set_label(ind)															
					plt.title('storm ' + expname)
					plt.subplot(4,ncols,count*ncols+2)
					make_ratio_scatter(Dreg,data2[ind][inreg],dataref2,isegrass,inosegrass,varname=symbol+' %',relative=True)
					set_label(ind)		
					#plt.ylim((-10,10))					
					plt.title('calm ' + expname)
			plt.tight_layout()
			plt.suptitle(domain)
			plt.savefig(image_outdir+'_'.join((domain,vari,quanti))+'.png',dpi=300)

			
# bar plots
#plt.figure(figsize=(20,16))	
if not os.path.isdir(image_outdir): os.mkdir(image_outdir)	
for mon in [10,]: #6]:	
	nametag='{:02d}'.format(mon)
	#ds=xr.open_dataset('newAna_2017_{:02d}.nc'.format(mon))
	#ds=xr.open_dataset('Ana_mveg_2017_{:02d}.nc'.format(mon))
	#ds=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/AnaMon10Restrict/Ana_mveg_drymasked_2017_{:02d}.nc'.format(mon))
	#ds=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/AnaMon10Restrict/Ana_mveg_2017_{:02d}.nc'.format(mon))
	#ds=xr.open_dataset(ncoutdir+'/Ana_mveg_2_2017_10.nc')
	#ds=xr.open_dataset(ncfile)
	
	for i in range(len(varnames)):  #[:2]
		#if i < 2:
		#	continue
		vari=varnames[i]
		if vari not in varnames_use:
			print('skipping')
			continue
			
		print(vari)
		unit=units[i]
		for j,quanti in enumerate(['mean','quantile95']): #enumerate(quantaties):
			print(quanti)
			symbol=symbols[i]
			if j==[0]:
				symbol=symbol.join(('<','>'))
			else:	
				symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
			label0='{:s} [{:s}]'.format(symbol,unit)
			try: #eventually variable not saved
				data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
				data2=[ds2[expi+'_'+vari+'_'+quanti].values for expi in experiments]
			except:
				try: #eventually variable not saved
					try: #eventually variable not saved
						data1=[ds1[quanti+'_'+expi+'_'+vari].values for expi in experiments]
						data2=[ds2[quanti+'_'+expi+'_'+vari].values for expi in experiments]
					except:
						data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
						data2=[ds2[expi+'_'+vari+'_'+quanti].values for expi in experiments]
				except:
					try: #eventually variable not saved
						vari1,vari2=vari.split('_')
						data1=[ds1[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
						data2=[ds2[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
					except:
						continue
			#break
		#break
	
			data1=[np.ma.masked_array(datai,mask=drymask) for datai in data1]
			data2=[np.ma.masked_array(datai,mask=drymask) for datai in data2]
			label=label0
			for domain,axis_limit in zip(domains,[efws,nfws,meadow]):
				# bins
				inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= axis_limit[3]) 
				# bins
				plt.clf()
			#for relative in [False,True]:
			#	plt.close('all')
			#	if relative:
			#		label=label0[:label0.rindex('[')]+'[%]'
			#		addtxt='_rel'
			#		#clim=((-25,25))
			#		clim='auto'
			#	else:
			#		label=label0
		    #
			#data=[results[expi][vari][quanti] for expi in experiments]
			# mask dry

				
				# bins
				cbins,exptruebins1,expfalsebins1,exptruedeltas1,expfalsedeltas1,exptruedeltas_rel1,expfalsedeltas_rel1=bin_select(experiments,D,Dbins,data1,gr3s,iref=0,inreg=inreg)	
				
				cbins,exptruebins2,expfalsebins2,exptruedeltas2,expfalsedeltas2,exptruedeltas_rel2,expfalsedeltas_rel2=bin_select(experiments,D,Dbins,data2,gr3s,iref=0,inreg=inreg)	
				
				##### plottings
				plt.close('all')
				for irow,expname in enumerate(experiments):
				
					# abs
					plt.figure(1,figsize=(20,16))
					#if irow==0:
					#	ax1=plt.subplot(5,2,irow*2+1)
					#else:	
					#	plt.subplot(5,2,irow*2+1,sharex=ax1,sharey=ax1)
					plt.subplot(5,2,irow*2+1)	
					plt.bar(cbins-0.1,expfalsebins1[expname],width=.25,color='brown')
					plt.bar(cbins+0.1,exptruebins1[expname],width=0.25,color='green')
					plt.ylabel(label)
					plt.xlim((lbin[0],ubin[-1]))
					if irow==0:
						plt.legend(('free','grass'))
						
					if irow==len(experiments):	
						plt.xlabel('Depth [m]')
					else:
						plt.xticks([])
					plt.title(expname + 'storm')	
					plt.subplot(5,2,irow*2+2)
					plt.bar(cbins-0.1,expfalsebins2[expname],width=.25,color='brown')
					plt.bar(cbins+0.1,exptruebins2[expname],width=0.25,color='green')
					plt.title(expname + 'calm')	
					plt.xlim((lbin[0],ubin[-1]))
					plt.xlabel('Depth')# edit
					plt.suptitle(domain)
					
					plt.figure(2,figsize=(20,16))
					#if irow==0:
					#	ax2=plt.subplot(5,2,irow*2+1)
					#else:	
					#	plt.subplot(5,2,irow*2+1,sharex=ax2,sharey=ax2)					
					plt.subplot(5,2,irow*2+1)
					plt.bar(cbins-0.1,expfalsedeltas1[expname],width=.25,color='brown')
					plt.bar(cbins+0.1,exptruedeltas1[expname],width=0.25,color='green')
					plt.ylabel('$\Delta$ '+label)
					plt.xlim((lbin[0],ubin[-1]))
					if irow==0:
						plt.legend(('free','grass'))
						
					if irow==len(experiments):	
						plt.xlabel('Depth [m]')
					else:
						plt.xticks([])
					plt.title(expname + 'storm')	
					plt.subplot(5,2,irow*2+2)
					plt.bar(cbins-0.1,expfalsedeltas2[expname],width=.25,color='brown')
					plt.bar(cbins+0.1,exptruedeltas2[expname],width=0.25,color='green')
					plt.title(expname + 'calm')	
					plt.xlim((lbin[0],ubin[-1]))
					plt.xlabel('Depth')# edit	
					plt.suptitle(domain)
				
				
					plt.figure(3,figsize=(20,16))
					if irow==0:
						ax3=plt.subplot(5,2,irow*2+1)
					else:	
						plt.subplot(5,2,irow*2+1,sharex=ax3,sharey=ax3)							
					#plt.subplot(5,2,irow*2+1)
					#plt.bar(cbins-0.1,expfalsedeltas1[expname]/expfalsebins1[experiments[iref]]*100,width=.25,color='brown')
					plt.bar(cbins-0.1,expfalsedeltas_rel1[expname]*100,width=.25,color='brown')
					plt.bar(cbins+0.1,exptruedeltas_rel1[expname]*100,width=0.25,color='green')
					#for ibin in range(len(cbins)):
					#	val=expfalsedeltas_rel1[expname][ibin]
					#	if val !=0:
					#		plt.text(cbins[ibin]-0.1,0,str(val*100),rotation=90*np.sign(val))
					#	val=exptruedeltas_rel1[expname][ibin]
					#	if val !=0:
					#		plt.text(cbins[ibin]+0.1,0,str(val*100),rotation=90*np.sign(val))
					plt.ylabel('$\Delta$ '+label[:label0.rindex('[')]+'[%]')
					plt.xlim((lbin[0],ubin[-1]))
					if irow==0:
						plt.legend(('free','grass'))
						
					if irow==len(experiments):	
						plt.xlabel('Depth [m]')
					else:
						plt.xticks([])
					plt.ylim((-100,100))	
					plt.title(expname + 'storm')	
					plt.subplot(5,2,irow*2+2)
					plt.bar(cbins-0.1,expfalsedeltas_rel2[expname]*100,width=.25,color='brown')
					plt.bar(cbins+0.1,exptruedeltas_rel2[expname]*100,width=0.25,color='green')
					#for ibin in range(len(cbins)):
					#	val=expfalsedeltas_rel2[expname][ibin]
					#	plt.text(cbins[ibin]-0.1,0,str(val*100),rotation=90*np.sign(val))
					#	val=exptruedeltas_rel2[expname][ibin]
					#	plt.text(cbins[ibin]+0.1,0,str(val*100),rotation=90*np.sign(val))
					
					plt.ylim((-100,100))
					plt.title(expname + 'calm')	
					plt.xlim((lbin[0],ubin[-1]))
					plt.xlabel('Depth')# edit	
					plt.suptitle(domain)	
				plt.figure(1)	
				plt.savefig(image_outdir+domain+'_'+vari+'_'+quanti+'abs.png',dpi=300)
				plt.figure(2)
				plt.savefig(image_outdir+domain+'_'+vari+'_'+quanti+'delta_abs.png',dpi=300)
				plt.figure(3)
				plt.savefig(image_outdir+domain+'_'+vari+'_'+quanti+'delta_rel.png',dpi=300)				
				
#				
#				###
#				
#		
#				## bar
#				#w=np.diff(cbins)[0]/(len(experiments)+1)
#				#dxs=np.arange((-np.diff(cbins)[0]+w)/2,(np.diff(cbins)[0]+w)/2,w)
#				#plt.figure()
#				#plt.clf()
#				#plt.subplot(2,2,1)
#				#plt.title('storm')
#				#clrs=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
#				#
#				
#				# abs diff abs diff rel storm calm
#				for irow,expname in enumerate(experiments):
#
#					# abs
#					plt.subplot(5,3,irow*3+1)
#					plt.bar(cbins-0.07,exptruebins1[expname],width=0.15)
#					plt.bar(cbins+0.07,exptruebins2[expname],width=0.15)
#					plt.title(expname)
#					#if irow<4:
#					#	plt.xticks([])
#					plt.ylabel('abs value')	
#					# diff	
#					plt.xlim((lbins[0],ubins[-1]))
#					if irow==0:
#						plt.legend(('storm','calm'))
#					plt.subplot(5,3,irow*3+2)
#					#plt.bar(cbins-0.07,exptruebins1[expname]-exptruebins1[experiments[iref]],width=0.15)
#					#plt.bar(cbins+0.07,exptruebins2[expname]-exptruebins2[experiments[iref]],width=0.15)
#					plt.bar(cbins-0.07,exptruedeltas1[expname],width=0.15)
#					plt.bar(cbins+0.07,exptruedeltas2[expname],width=0.15)
#					plt.title(expname)
#					plt.xlim((lbins[0],ubins[-1]))
#					#if irow<4:
#					#	plt.xticks([])
#					#else:	
#					#	plt.xlabel('Depth [m]')
#					plt.ylabel('$\Delta$ abs')	
#					# rel diff	
#					plt.subplot(5,3,irow*3+3)
#					plt.bar(cbins-0.07,(exptruebins1[expname]-exptruebins1[experiments[iref]])/exptruebins1[experiments[iref]]*100,width=0.15)
#					plt.bar(cbins+0.07,(exptruebins2[expname]-exptruebins2[experiments[iref]])/exptruebins2[experiments[iref]]*100,width=0.15)
#					plt.title(expname)
#					if irow<4:
#						plt.xticks([])
#					plt.ylabel('$\Delta$ rel [%]')	
#					plt.xlim((lbins[0],ubins[-1]))
#				plt.suptitle(domain + ' '+ label)	
#
#
#				plt.savefig(image_outdir+'{:s}_{:02d}_{:s}{:s}'.format(domain,mon,quanti,vari+addtxt),dpi=300)

	
	
# edit
# abs diff abs diff rel storm calm
iref=0

	
	
	
	
	
	# abs diff abs diff rel storm calm
iref=0
for irow,expname in enumerate(experiments):

	# abs
	barwidth=0.035
	plt.figure(1)
	plt.subplot(5,2,irow*2+1)
	plt.bar(cbins-0.1,expfalsebins1[expname],width=.25,color='brown')
	plt.bar(cbins+0.1,exptruebins1[expname],width=0.25,color='green')
	plt.ylabel(quanti +' '+vari)
	plt.xlim((lbin[0],ubin[-1]))
	if irow==0:
		plt.legend(('free','grass'))
		
	if irow==len(experiments):	
		plt.xlabel('Depth [m]')
	else:
		plt.xticks([])
	plt.title(expname + 'storm')	
	plt.subplot(5,2,irow*2+2)
	plt.bar(cbins-0.1,expfalsebins2[expname],width=.25,color='brown')
	plt.bar(cbins+0.1,exptruebins2[expname],width=0.25,color='green')
	plt.title(expname + 'calm')	
	plt.xlim((lbin[0],ubin[-1]))
	plt.xlabel('Depth')
	#if irow<4:
	#	plt.xticks([])


	
		
		
		
		
		
	plt.subplot(5,3,irow*3+2)
	#plt.bar(cbins-0.07,exptruebins1[expname]-exptruebins1[experiments[iref]],width=0.15)
	#plt.bar(cbins+0.07,exptruebins2[expname]-exptruebins2[experiments[iref]],width=0.15)
	#plt.bar(cbins-0.07,exptruedeltas1[expname],width=0.15)
	#plt.bar(cbins+0.07,exptruedeltas2[expname],width=0.15)
	plt.bar(cbins-2*barwidth,expfalsedeltas1[expname],width=barwidth)
	plt.bar(cbins-barwidth,exptruedeltas1[expname],width=barwidth)
	plt.bar(cbins+barwidth,expfalsedeltas2[expname],width=barwidth)
	plt.bar(cbins+2*barwidth,exptruedeltas2[expname],width=barwidth)

	plt.title(expname)
	plt.xlim((lbin[0],ubin[-1]))
	#if irow<4:
	#	plt.xticks([])
	#else:	
	#	plt.xlabel('Depth [m]')
	plt.ylabel('$\Delta$ abs')	
	# rel diff	
	plt.subplot(5,3,irow*3+3)
	plt.bar(cbins-2*barwidth,(expfalsebins1[expname]-expfalsebins1[experiments[iref]])/expfalsebins1[experiments[iref]]*100,width=0.15)
	plt.bar(cbins-barwidth,(exptruebins1[expname]-exptruebins1[experiments[iref]])/exptruebins1[experiments[iref]]*100,widt=hbarwidth)
	plt.bar(cbins+barwidth,(expfalsebins2[expname]-expfalsebins2[experiments[iref]])/expfalsebins2[experiments[iref]]*100,widthbarwidth)
	plt.bar(cbins+2*barwidth,(exptruebins2[expname]-exptruebins2[experiments[iref]])/exptruebins2[experiments[iref]]*100,width=barwidth)
	plt.title(expname)
	if irow<4:
		plt.xticks([])
	plt.ylabel('$\Delta$ rel [%]')	
	plt.xlim((lbin[0],ubin[-1]))
plt.suptitle(domain + ' '+ label)	



		
for stat in ['mean','q95']:
	for varname in varnames:
	
		for relative in [False,True]:
			plt.close('all')
			plt.figure(figsize=(20,16))
			count=0
			print(varname)
			for irow,ind in enumerate([1,2,3,4]):
				V=stats[experiments[ind]][varname][stat]
				Vcntrl=stats['Veg_REF'][varname][stat]
				
				for domain,axis_limit in zip(domains,[efws,nfws,meadow]):
					inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= 
									axis_limit[3]) 
									
					count+=1
					plt.subplot(4,3,count)	
					
					isegrass=gr3s[experiments[ind]][inreg]>0
					inosegrass=gr3s[experiments[ind]][inreg]==0				
					if relative:
						delta=(V[inreg]-Vcntrl[inreg])/Vcntrl[inreg]*100
					else:
						delta=V[inreg]-Vcntrl[inreg]
						
					inan=np.isnan(delta) | (delta >1e4)				
					delta=np.ma.masked_array(delta,mask=inan)
						
					deltafree=delta[inosegrass]
					deltagrass=delta[isegrass]
					plt.scatter(D[inreg][inosegrass],deltafree,marker='.',color='brown',s=0.4) 
					
					isort=np.argsort(D[inreg][inosegrass])
					
					plt.scatter(D[inreg][isegrass],deltagrass,marker='.',color='green',s=0.4) 
					plt.legend(['free (avg:{:.4f})'.format( np.round(np.nanmean(deltafree),4)),'seagrass (avg:{:.4f})'.format( np.round(np.nanmean(deltagrass),4))])
					plt.ylabel('$\Delta$ ' +varname +stat)
					if relative:
						plt.ylim((-100,100))
					if irow==3:
						plt.xlabel('depth [m]')
					plt.title(experiments[ind] + domain+' (avg:{:.4f})'.format( np.round(np.nanmean(delta),4)))
					plt.xlim((-4,5))
			plt.tight_layout() 
			plt.savefig(image_outdir+'scatter'+varname +stat+addtxt[relative]+'png',dpi=300)
			
			

dsstorm=