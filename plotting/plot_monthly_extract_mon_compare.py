""" Plot results of extract monthly 
colorbar positions need to be fidddeled around"""

from glob import glob
import os
import sys
#sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
import cartopy.crs as ccrs
import matplotlib
matplotlib.use('qtagg')
from matplotlib import pyplot as plt
from schism import* # import schism class to read grid structure

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
import cmocean
from matplotlib.path import Path   # polygon operations
dask.__version__, distributed.__version__


# 
file='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg_jul/Veg_REF/mon_Ana_mveg_jan_2017_07.nc' #stats file
rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg_jul/Veg_REF/'  # schism dircetory with grid

file1='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg_jul/Veg_REF/mon_Ana_mveg_jan_2017_07.nc'
file2='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg_jan/Veg_REF/mon_Ana_mveg_jan_2017_01.nc'

cwd=os.getcwd()



mon=1 # month of analyiss
# is store it in the fike abmes


#surface bottom
image_outdir='./mon_comp_period_stats/'
image_outdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg_jul/period_stats/'
if not os.path.exists(image_outdir):
	os.mkdir(image_outdir)

prefix='2017_07-_2017-08' # for output nc file
prefix='2017_01-_2017-02' # for output nc file

prefix='2017_Jan_vs_Jul'

## selected varnames
varnames=['elevation','sigWaveHeight','depthAverageVel','totalSuspendedLoad_bottom']

units={'elevation':'m','sigWaveHeight':'m','depthAverageVel':'m/s','totalSuspendedLoad_bottom':'g/L'}
symbols={'elevation':'$\zeta$','sigWaveHeight':'HS','depthAverageVel':'$u_{<z>}$','totalSuspendedLoad_bottom':'spm'}


# define plot regions --
#domain_tag=['total','efws','nfws']
#efws=(6.7722882425499, 8.341052998071572, 53.51960339871429, 53.83566262094643)
#nfws=(8.03, 9.01, 54.32902335808929, 55.55471643942857)
#efws=(6.55, 8.45, 53.2, 53.83566262094643)


# plot layout for sublot comparison of scnearios:
regions={
	'Domian':
		{'figsize':(8,6),
		'axis_limit':None,
		'nrows':2,
		'ncols':3,
		'cb_orientation':'horizontal',
		'cbarpos':(0.675,0.325,0.3,0.02), # set cbar position in subplot (try and error)
		'h_t_spacing':(None,None)        # hspace top space adjust # for subplotspacing ( hspace=-0.25,top=0.99)

		},
		
	'EFWS':
		{'figsize':(4.8,4.8),
		'axis_limit':(6.55, 8.45, 53.2, 53.83566262094643),
		'nrows':3,
		'ncols':2,
		'cb_orientation':'horizontal',
		'cbarpos':(0.565,0.2,0.35,0.02), # set cbar position in subplot (try and error)
		'h_t_spacing':(-0.25,0.99)        # hspace top space adjust # for subplotspacing ( hspace=-0.25,top=0.99) for  fig.subplots_adjust
		},

	'NFWS':
		{'figsize':(8,4.1),
		'axis_limit':(6.55, 8.45, 53.2, 53.83566262094643),
		'nrows':1,  #1,5
		'ncols':5,
		'cb_orientation':'horizontal',
		'cbarpos':(0.3,0.125,0.6,0.02), # set cbar position in subplot (try and error)
		'h_t_spacing':(-0.25,1.00)  
		}		
		
}



cmap=plt.cm.turbo
cmap0=plt.cm.viridis
cmap=cmocean.cm.balance





if not os.path.isdir(image_outdir):
	os.mkdir(image_outdir)
quantaties=['mean','quantile95','max']	
level_tag=['_bottom','_surface']



# plotting function for scenarios
def compare_slabs(s,data,names=None,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None):
	if nrows*ncols>=len(data):
		fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': ccrs.Mercator()},figsize=(11,8.5))
		axs = axs.flatten()
		alldata=np.hstack(data)		
		vmin,vmax=np.quantile(alldata,0.1),np.quantile(alldata,0.99)
		for i,datai in enumerate(data):
			ph,ch,ax=s.plotAtnodesGeo(datai,ax=axs[i])
			ch.remove()
			if names != None:
				ax.set_title(names[i])
			if axis_limit!=None:
				ax.set_extent(axis_limit)
			ph.set_clim((vmin,vmax))	
		# Add a colorbar axis at the bottom of the graph
		cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation='horizontal',extend='both',label=cblabel)
		plt.tight_layout()
	else:
		print('nrows*ncols < len(data)!')
	return cbar_ax
		
def compare_slabs_diff(s,data,names=None,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None):
	""" first entry of list data is reference substracete from experiments"""
	if nrows*ncols>=len(data)-1:
		fig, axs = plt.subplots(nrows=nrows,ncols=ncols,								subplot_kw={'projection': ccrs.Mercator()},figsize=(11,8.5))
		axs = axs.flatten()
		if nrows*ncols < 2:
			axs=[axs]
		for i,datai in enumerate(data[1:]):
			ph,ch,ax=s.plotAtnodesGeo(datai-data[0],ax=axs[i])
			ch.remove()
			if names != None:
				ax.set_title('-'.join((names[i+1],names[0])))
			if axis_limit!=None:
				ax.set_extent(axis_limit)
		# Add a colorbar axis at the bottom of the graph
		cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation='horizontal',extend='both',label=cblabel)
		plt.tight_layout()
	else:
		print('nrows*ncols < len(data)!')		
		
def get_quiver_locs(s,narrows=30,axis_limit=None):		
	if tuple(axis_limit)==None:
		axis_limit=(np.min(s.lon),np.max(s.lon),np.min(s.lat),np.max(s.lat))

	xlim=axis_limit[:2]
	ylim=axis_limit[2:]
	x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/narrows)
	y=np.arange(ylim[0],ylim[1],(ylim[1]-ylim[0])/narrows)
	X, Y = np.meshgrid(x,y)
	d,qloc=s.node_tree_latlon.query((np.vstack([X.ravel(), Y.ravel()])).transpose()) #quiver locations				
	#xref,yref=np.asarray(plt.axis())[[1,2]] +  np.diff(plt.axis())[[0,2]]*[- 0.2, 0.1]
	#if (self.shape==(self.nt,self.nnodes)) or (self.shape==(self.nt,self.nnodes,self.nz)): 
	#	vmax=1.5#np.percentile(np.sqrt(u[qloc]**2+v[qloc]**2),0.95)
	#else:
	#	vmax=np.double(self.maxfield.get())
	return qloc
	
		
def compare_slabs_diff2(s,data,names=None,clim0='auto',cmap0=plt.cm.viridis,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None,clim=None,figsize=(11,8.5),cb_orientation='horizontal',geo=True,relative=False):
	""" first entry of list data is reference substracete from experiments"""
	phs=[]
	qhs=[]
	def add_quiver(ivs,ax,x,y,u,v,qloc,axis_limit,scale=15,vmax=1,color='k'):
		axis_limit=np.asarray(axis_limit)
		xref,yref=axis_limit[[1,2]] +  np.diff(axis_limit)[[0,2]]*[-0.3, 0.85]
		if ivs:
			qh=ax.quiver(np.concatenate((x[qloc],[xref,])),np.concatenate((y[qloc],[yref,])),np.concatenate((u[qloc],[vmax,])),np.concatenate((v[qloc],[0,])),scale=scale,scale_units='inches',color=color)
			ax.text(xref,yref,'\n {:.2f} \n m/s'.format(vmax),fontsize='x-small')
		else:
			qh=[]	
		return qh
			
	if nrows*ncols>=len(data)-1:
		if geo:	
			import cartopy.crs as ccrs
			import cartopy.feature as cfeature
			from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

			fig, axs = plt.subplots(nrows=nrows,ncols=ncols,								subplot_kw={'projection': ccrs.Mercator()},figsize=figsize)
		else:
			fig, axs = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize)		
		axs = axs.flatten()
		
		if data[0].shape[0]==2:
			ivs=True #vector
			udata=[data[0][0,:]]
			vdata=[data[0][1,:]]
			for i,datai in enumerate(data[1:]):
				udata.append(data[i+1][0,:]-udata[0])
				vdata.append(data[i+1][1,:]-vdata[0])
			data=[(datai*2).sum(axis=0) for datai in data]	
			
			qloc=get_quiver_locs(s,narrows=40,axis_limit=axis_limit)
			proj=ccrs.Mercator()
			outproj=proj.transform_points(ccrs.Geodetic(),np.asarray(s.lon),np.asarray(s.lat))
			projx,projy=outproj[:,0],outproj[:,1]
			
			axis_limit2=np.asarray(axis_limit)
			outproj=proj.transform_points(ccrs.Geodetic(),axis_limit2[:2],axis_limit2[2:])
			projlim=tuple(np.hstack((outproj[:,0],outproj[:,1])))
		else:
			ivs=False
			projx,projy=np.asarray(s.lon),np.asarray(s.lat)

		if relative:
			diffs=[(datai-data[0])/data[0]*100 for datai in data[1:]]
		else:
			diffs=[datai-data[0] for datai in data[1:]]
		if clim=='auto': # limit to in regiond 5 and 95 percentile
			if axis_limit ==	None:
				inreg=np.arange(s.nnodes)
			else:
				lon,lat=np.asarray(s.lon),np.asarray(s.lat)
				inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= axis_limit[3])
			reg_diffs=np.hstack([np.abs(diffsi[inreg]) for diffsi in diffs])
			vmin,vmax=np.nanquantile(reg_diffs,0.1),np.nanquantile(reg_diffs,0.99)		
			clim=(-vmax*(1-(vmin==0)),vmax)
		if geo:	
			ph,ch,ax=s.plotAtnodesGeo(data[0],ax=axs[0],cmap=cmap0)
			if ivs:
				qh=add_quiver(ivs,ax,projx,projy,udata[0],vdata[0],qloc,projlim,scale=15,vmax=1,color='k')
				qhs.append(qh)		
		else:	
			ph,ch,ax=s.plotAtnodes(data[0],ax=axs[0],cmap=cmap0)
			qh=add_quiver(ivs,ax,np.asarray(s.lon),np.asarray(s.lat),udata[0],vdata[0],qloc,axis_limit,scale=1,vmax=1,color='w')
			qhs.append(qh)		
		if clim0=='auto': # limit to in regiond 5 and 95 percentile
			vmin0,vmax0=np.nanquantile(data[0],0.1),np.nanquantile(data[0],0.99)
			clim0=(-vmax0,vmax0)
		if clim != None:
				ph.set_clim(clim0)	
		phs.append(ph)		
		
		if names != None:
				ax.set_title(names[0])
		if axis_limit!=None:
				if geo:
					ax.set_extent(axis_limit)
				else:	
					ax.axis(axis_limit)
		if nrows*ncols < 2:
			axs=[axs]
		for i,datai in enumerate(data[1:]):
			if geo:	
				ph,ch,ax=s.plotAtnodesGeo(diffs[i],ax=axs[i+1],cmap=cmap)		
				if ivs:
					qh=add_quiver(ivs,ax,projx,projy,udata[i+1],vdata[i+1],qloc,projlim,scale=3,vmax=0.2,color='k')
					qhs.append(qh)	
			else:	
				ph,ch,ax=s.plotAtnodes(diffs[i],ax=axs[i+1],cmap=cmap)
				if ivs:
					qh=add_quiver(ivs,ax,np.asarray(s.lon),np.asarray(s.lat),udata[i],vdata[i],qloc,axis_limit,scale=1,vmax=0.5,color='w')
					qhs.append(qh)	
			ch.remove()
			if names != None:
				ax.set_title('-'.join((names[i+1],names[0])))
			if axis_limit!=None:
				if geo:
					ax.set_extent(axis_limit)
				else:	
					ax.axis(axis_limit)
			if clim != None:
				ph.set_clim(clim)
			phs.append(ph)	
		# Add a colorbar axis at the bottom of the graph
		cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation=cb_orientation,extend='both',label=cblabel)
		
		if geo==False:
			plt.tight_layout()
		
		for axi in axs[len(data):]:
			axi.remove()
		
		return fig,axs,cbar_ax,phs,qhs
	else:
		print('nrows*ncols < len(data)!')	


os.chdir(rundir)
s=schism_setup()		
s.init_node_tree(latlon=True)		
		
nametag='{:02d}'.format(mon)


ds1=xr.open_dataset(file1)
ds2=xr.open_dataset(file2)

expkeys=[key for key in list(ds1.keys()) if  'dryFlagNode_mean' in key  ]
experiments=[key.split('_dryFlagNode_mean')[0] for key in expkeys]
expnames=['Ref','Blank','Veg$_{max}$','Veg$_{HE}$','Veg$_{LE}$']


for i in range(len(varnames)):

	vari=varnames[i]
	
	print(vari)
	unit=units[vari]
	
	for j,quanti in enumerate(quantaties):
		print(quanti)

		symbol=symbols[vari]
		if j==[0]:
			symbol=symbol.join(('<','>'))
		else:	
			symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
		label0='{:s} [{:s}]'.format(symbol,unit)
		try: #eventually variable not saved
			data1=[ds1[expi+'_'+vari+'_'+quanti].values for expi in experiments]
			data2 = [ds2[expi + '_' + vari + '_' + quanti].values for expi in experiments]
		except:
			try: #eventually variable not saved
				data1=[ds1[quanti+'_'+expi+'_'+vari].values for expi in experiments]
				data2=[ds2[quanti+'_'+expi+'_'+vari].values for expi in experiments]

			except:
				try: #eventually variable not saved
					vari1,vari2=vari.split('_')
					data1=[ds1[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
					data2=[ds2[expi+'_'+vari1+'_'+quanti+'_'+vari2].values for expi in experiments]
				except:
					continue

		#data=[results[expi][vari][quanti] for expi in experiments]
		
		for reg in regions:

			# set plot properties
			figsize=regions[reg]['figsize']
			axis_limit=regions[reg]['axis_limit']
			nrows=regions[reg]['nrows']
			ncols=regions[reg]['ncols']
			cb_orientation=regions[reg]['cb_orientation']
			cbarpos=regions[reg]['cbarpos']
			hspace,top=regions[reg]['h_t_spacing']		
			
			for relative in [False,True]:
				plt.close('all')
				if relative:
					label=label0[:label0.rindex('[')]+'[%]'
					addtxt='_rel'
					clim=((-25,25))
				else:
					label=label0
					clim='auto'
					addtxt=''						
				if len(data[0].shape)==2: # surface and bottom
					for ilevel in [0,-1]:
						level_data=[datai[:,ilevel] for datai in data]
						fig,axs,cbar_ax,phs,qhs=compare_slabs_diff2(s,level_data,names=expnames,cmap=cmap,nrows=nrows,ncols=ncols,cblabel='$\Delta$'+label,axis_limit=axis_limit,clim=clim,figsize=figsize,cb_orientation=cb_orientation,relative=relative)#clim=(-0.01,0.01))
						plt.suptitle( '{:s} {:s}_{:s} \n {:s} - {:s}'.format(quanti,vari,level_tag[ilevel],date0,date1)  )
						
						
						if hspace!=None:
							fig.subplots_adjust( hspace=hspace,top=top)					
						cbar_ax.set_position(cbarpos)
						plt.savefig(image_outdir+nametag+prefix+'_{:s}_{:02d}_{:s}{:s}_{:s}.png'.format(reg,mon,quanti,vari+addtxt,level_tag[ilevel]),dpi=300)
						
				else: #2d
					fig,axs,cbar_ax,phs,qhs=compare_slabs_diff2(s,data,names=expnames,cmap=cmap,nrows=nrows,ncols=ncols,cblabel='$\Delta$'+label,axis_limit=axis_limit,clim=clim,figsize=figsize,cb_orientation=cb_orientation,relative=relative)	
										
					# adjust this part to set colorbar toi appropriate postion
					if reg=='Domain':
						xt=axs[0].get_xticks()
						for axi in axs[:len(data)]:
							axi.set_xticks(xt[::2])
							
					if hspace!=None:
						fig.subplots_adjust( hspace=hspace,top=top)					
						
					cbar_ax.set_position(cbarpos)
					plt.tight_layout()
					plt.savefig(image_outdir+prefix+'_{:s}_{:02d}_{:s}{:s}.png'.format(reg,mon,quanti,vari+addtxt),dpi=300)
ds.close()