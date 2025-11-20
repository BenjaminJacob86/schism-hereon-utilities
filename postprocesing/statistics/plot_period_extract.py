""" compare statistics from different simulations """
from glob import glob
import os
import sys
import pandas as pd
sys.path.insert(0, '/work/gg0028/SCHISM/schism-hzg-utilities/')
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
from schism import *  # import schism class to read grid structure

# parallel stuff
# %matplotlib inline
import xarray as xr

from dask.diagnostics import ProgressBar
from dask_jobqueue import SLURMCluster
from distributed import Client, progress
from distributed.utils import tmpfile

import dask
import distributed

plt.ion()




# compare storm max levels
rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir1='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/stats_extract_Herwart/Veg_CNTRL/'
ncdir2='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/stats_extract2/Veg_CNTRL/'

ncdir1b='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/stats_extract_Herwart/Veg_max/'
ncdir2b='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/stats_extract2/Veg_max/'


#names=['Herwart','Xaver']
names=['1997','2090','2090-1997']
ncdir1='/work/gg0028/g260114/RUNS/Rest-Coast/analysis/1997_stats/19971001_19971101/Veg_CNTRL/'
ncdir2='/work/gg0028/g260114/RUNS/Rest-Coast/analysis/2090_stats/20901001_20901101/Veg_CNTRL/'

ncdir1b='/work/gg0028/g260114/RUNS/Rest-Coast/analysis/1997_stats/19971001_19971101/Veg_max/'
ncdir2b='/work/gg0028/g260114/RUNS/Rest-Coast/analysis/2090_stats/20901001_20901101/Veg_max/'



image_dir='/work/gg0028/g260114/RUNS/Rest-Coast/analysis/pics/'

#/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/stats_extract_2017_JJA/
#/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/stats_extract_2017_JJA/

#a=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/stats_extract_2017_JJA/Veg_CNTRL/bottomStressX_quantile95_1.nc')
#plt.figure()
#s.plotAtnodes(a['bottomStressX_quantile95'].values)


os.chdir(rundir)
s=schism_setup()

#quad
# get element areas
A=[]
for i in range(s.nvplt.shape[0]):
	nodes=s.nvplt[i,:]+1
	A.append(s.proj_area(nodes))
A=np.asarray(A)
#  s.A=s.compute_element_areas()

## element centers
lon=np.asarray(s.lon)
lat=np.asarray(s.lat)
faces=s.nvplt
cx=np.mean(lon[faces],axis=1)
cy=np.mean(lat[faces],axis=1)

elcoords=[[cx[i],cy[i]] for i in range(len(cx))] # pooint pairs of nodes
elem_nn_tree = cKDTree(elcoords) # next neighbour search tree	     
#
# get element indices for basins
eleminds={}
#for ifile,tag in enumerate(areas):
#	areaPoly=Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
#	eleminds[tag]=areaPoly.contains_points(elcoords)	



 

varnames=['elevation','sigWaveHeight','depthAverageVel','totalSuspendedLoad_bottom','bottomStress']
units={'elevation':'m','sigWaveHeight':'m','depthAverageVel':'m/s','totalSuspendedLoad_bottom':'g/L','bottomStress':'Pa'}
symbols={'elevation':'$\zeta$','sigWaveHeight':'HS','depthAverageVel':'$u_{<z>}$','totalSuspendedLoad_bottom':'spm','bottomStress':'$\\tau_{btm}$'}

varname='elevation'
quantaties=['mean','max','quantile95']


# plot layout for sublot comparison of scnearios:
regions={
	'EFWS':
		{'figsize':(4.8,4.8),
		'axis_limit':(6.55, 8.45, 53.2, 53.83566262094643),
		'nrows':3,
		'ncols':2,
		'cb_orientation':'horizontal',
		'cbarpos':(0.565,0.2,0.35,0.02), # set cbar position in subplot (try and error)
		'h_t_spacing':(-0.25,0.99)        # hspace top space adjust # for subplotspacing ( hspace=-0.25,top=0.99) for  fig.subplots_adjust
		},
}

if not os.path.exists(image_dir):
    os.mkdir(image_dir)



axis_limit=regions['EFWS']['axis_limit']


inbox_elems= (cx >axis_limit[0]) &  (cx < axis_limit[1]) & (cy >axis_limit[2]) &  (cy < axis_limit[3])
inbox_weights=A[inbox_elems]/A[inbox_elems].sum()


# plotting function for scenarios
def compare_slabs(s,data,names=None,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None,figsize=(11,8.5)):
	if nrows*ncols>=len(data):
		fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': ccrs.Mercator()},figsize=figsize)
		axs = axs.flatten()
		alldata=np.hstack(data)		
		vmin,vmax=np.nanquantile(alldata,0.1),np.nanquantile(alldata,0.99)
		for i,datai in enumerate(data):
			ph,ch,ax=s.plotAtnodesGeo(datai,ax=axs[i])
			ch.remove()
			if names != None:
				ax.set_title(names[i])
			if axis_limit!=None:
				ax.set_extent(axis_limit)
			ph.set_clim((vmin,vmax))	
		# Add a colorbar axis at the bottom of the graph


        # per default centered filling hakf width    
		x0=1/(2*ncols)
		w=1/ncols
		cbar_ax = fig.add_axes([x0, 0.2, w, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation='horizontal',extend='both',label=cblabel)
		plt.tight_layout()
	else:
		print('nrows*ncols < len(data)!')
	return fig,axs,cbar_ax



# plotting function for scenarios
def compare_slabs_with_diff(s, data, names=None, cmap=plt.cm.turbo, cblabel='', cblabel2='', axis_limit=None,
                            figsize=(11, 8.5)):
    """ plot Two scenarios in abs values and their differnce, each row difference betwwen data set """
    nrows = int(np.floor(np.sqrt(len(data))))
    ncols = int(np.ceil(np.sqrt(len(data)))) + 1


    
    show_reg_avg=True # plot region average in title
    addtxt=''
    
    if nrows * ncols >= len(data):
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, subplot_kw={'projection': ccrs.Mercator()}, figsize=figsize)
        axs = axs.flatten()
        alldata = np.hstack(data)
        vmin, vmax = np.nanquantile(alldata, 0.1), np.nanquantile(alldata, 0.99)
        for i, datai in enumerate(data):
        
            if show_reg_avg:            
                elem_vals=datai[s.nvplt].mean(axis=1)
                regmean=np.nansum(elem_vals[inbox_elems]*inbox_weights)
                addtxt=' (avg={:.2f})'.format(regmean)
            
            ph, ch, ax = s.plotAtnodesGeo(datai, ax=axs[i])
            ch.remove()
            if names != None:
                ax.set_title(names[i]+addtxt)
            if axis_limit != None:
                ax.set_extent(axis_limit)
            ph.set_clim((vmin, vmax))
        # Add a colorbar axis at the bottom of the graph

        # per default centered filling hakf width
        x0 = 1 / (2 * ncols)
        w = 1 / ncols
        cbar_ax = fig.add_axes([x0, 0.2, w, 0.02])

        # Draw the colorbar
        cbar = fig.colorbar(ph, cax=cbar_ax, orientation='horizontal', extend='both', label=cblabel)
        plt.tight_layout()

        # add difference plot
        Delta = data[1] - data[0]
        
        if show_reg_avg:            
            elem_vals=Delta[s.nvplt].mean(axis=1)
            regmean=np.nansum(elem_vals[inbox_elems]*inbox_weights)
            addtxt=' (avg={:.2f})'.format(regmean)        
        
        ph, ch, ax = s.plotAtnodesGeo(Delta, ax=axs[-1])
        ch.remove()
        if names != None:
            ax.set_title(names[-1]+addtxt)
        if axis_limit != None:
            ax.set_extent(axis_limit)
        ph.set_cmap('RdBu_r')
        vmax = np.nanquantile(np.abs(Delta), 0.95)
        ph.set_clim((-vmax, vmax))

        # per default centered filling hakf width
        plt.pause(0.0001)  # need to pause otherwise wrong values
        temp = axs[2].get_position()  # strange results in one go
        temp = temp.get_points()[0][1]
        # wa_axis=axs[2].get_position().get_points()[0][1] # subplot axis width
        wa_axis = temp
        w = wa_axis * 0.75  # colorbar
        x0 = 2 * 1 / (ncols) + (wa_axis - w) / 2

        cbar_ax2 = fig.add_axes([x0, 0.2, w, 0.02])
        cbar2 = fig.colorbar(ph, cax=cbar_ax2, orientation='horizontal', extend='both', label=cblabel2)
        plt.pause(0.0001) 
    else:
        print('nrows*ncols < len(data)!')
    return fig, axs, cbar_ax, cbar_ax2


# side by side comparison
if False:
    for varname in varnames:
        vari=varname
        unit=units[vari]

        for quanti in quantaties:

            print('doing '+ ' '.join((vari,quanti)))
            break
            # build uint
            symbol=symbols[vari]
            if quanti=='mean':
                symbol=symbol.join(('<','>'))
            else:	
                symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
            label0='{:s} [{:s}]'.format(symbol,unit)

            
            #change name order in creation file
            varname_temp=varname.replace('_bottom','')
            
            if 'bottom' in varname:
                files=glob.glob(ncdir1+'*'+varname_temp+'*'+quanti+'*')
                filename=files[0].split(ncdir1)[1]
                ds1=xr.open_dataset(ncdir1+filename)
                variable=list(ds1.variables.keys())[0]
            elif 'surface' in varname:            
                files=glob.glob(ncdir1+'*'+varname_temp+'*'+quanti+'*')
                filename=files[1].split(ncdir1)[1]
                ds1=xr.open_dataset(ncdir1+filename)
                variable=list(ds1.variables.keys())[0]
            else:
                filename=varname+'_'+quanti+'_1.nc'
                variable=varname+'_'+quanti
            
            ds1=xr.open_dataset(ncdir1+filename)
            variable=list(ds1.variables.keys())[0]
            data1=ds1[variable].values

            ds2=xr.open_dataset(ncdir2+filename)
            variable=list(ds2.variables.keys())[0]
            data2=ds2[variable].values


            ds1b=xr.open_dataset(ncdir1b+filename)
            variable=list(ds1b.variables.keys())[0]
            data1b=ds1b[variable].values

            ds2b=xr.open_dataset(ncdir2b+filename)
            variable=list(ds2b.variables.keys())[0]
            data2b=ds2b[variable].values

            # compare no Seagrasss periods absolute

            data=[data1,data2]
            
            
            fig,axs,cbar_ax=compare_slabs(s,data,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=label0,axis_limit=axis_limit,figsize=(8,3.6))
            fig.subplots_adjust( hspace=0.1,top=1.2)		
            plt.savefig(image_dir+varname+'_'+quanti+'_SLR_Veg_CNTRL.png',dpi=300)
            plt.close()

            # compare reduction by  Seagrasss periods absolute
            Ddata=[data1b-data1,data2b-data2]
            clabel='$\Delta$'+label0
            fig,axs,cbar_ax=compare_slabs(s,Ddata,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=clabel,axis_limit=axis_limit,figsize=(8,3.6))
            fig.subplots_adjust( hspace=0.1,top=1.2)		
            plt.savefig(image_dir+varname+'_'+quanti+'_SLR_Veg_max_minus_CNTRL.png',dpi=300)
            plt.close()

            # compare reduction by  Seagrasss periods relative
            Ddata=[(data1b-data1)/data1*100,(data2b-data2)/data2*100]
            
            # replace /0 nans by nan
            for item in Ddata:
                iinf=np.isinf(item)
                item[iinf]=np.nan
            np.isinf(Ddata).max()
            
            clabel='$\Delta$'+label0
            clabel=clabel[:-4]+'[%]'
            fig,axs,cbar_ax=compare_slabs(s,Ddata,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=clabel,axis_limit=axis_limit,figsize=(8,3.6))
            fig.subplots_adjust( hspace=0.1,top=1.2)		
            plt.savefig(image_dir+varname+'_'+quanti+'_SLR_Veg_max_minus_CNTRL_rel.png',dpi=300)
            plt.close()

# difference of seagrass during storms

## element centers
#cx=np.mean(x[faces],axis=1)
#cy=np.mean(y[faces],axis=1)
#elcoords=[[cx[i],cy[i]] for i in range(len(cx))] # pooint pairs of nodes
#elem_nn_tree = cKDTree(elcoords) # next neighbour search tree	     
##
## get element indices for basins
#eleminds={}
##for ifile,tag in enumerate(areas):
##	areaPoly=Path(list(zip(seas[tag][:,0],seas[tag][:,1])))
#	eleminds[tag]=areaPoly.contains_points(elcoords)	


for varname in varnames:
    vari=varname
    unit=units[vari]

    for quanti in quantaties:

        print('doing '+ ' '.join((vari,quanti)))

        # build uint
        symbol=symbols[vari]
        if quanti=='mean':
            symbol=symbol.join(('<','>'))
        else:	
            symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
        label0='{:s} [{:s}]'.format(symbol,unit)

        
        #change name order in creation file
        varname_temp=varname.replace('_bottom','')
        
        if ('bottom' in varname) & (varname!='bottomStress'):
            files=glob.glob(ncdir1+'*'+varname_temp+'*'+quanti+'*')
            filename=files[0].split(ncdir1)[1]
            ds1=xr.open_dataset(ncdir1+filename)
            variable=list(ds1.variables.keys())[0]
        elif 'surface' in varname:            
            files=glob.glob(ncdir1+'*'+varname_temp+'*'+quanti+'*')
            filename=files[1].split(ncdir1)[1]
            ds1=xr.open_dataset(ncdir1+filename)
            variable=list(ds1.variables.keys())[0]
        else:
            filename=varname+'_'+quanti+'_1.nc'
            variable=varname+'_'+quanti
        
        ds1=xr.open_dataset(ncdir1+filename)
        variable=list(ds1.variables.keys())[0]
        data1=ds1[variable].values

        ds2=xr.open_dataset(ncdir2+filename)
        variable=list(ds2.variables.keys())[0]
        data2=ds2[variable].values


        ds1b=xr.open_dataset(ncdir1b+filename)
        variable=list(ds1b.variables.keys())[0]
        data1b=ds1b[variable].values

        ds2b=xr.open_dataset(ncdir2b+filename)
        variable=list(ds2b.variables.keys())[0]
        data2b=ds2b[variable].values

        # compare no Seagrasss periods absolute

        data=[data1,data2]
        
        
        fig,axs,cbar_ax,cbar_ax2=compare_slabs_with_diff(s,data,names=names,cmap=plt.cm.turbo,cblabel=label0,cblabel2='$\Delta$ '+label0,axis_limit=axis_limit,figsize=(10,3.6))
        fig.subplots_adjust( hspace=0.1,top=1.2)		
        plt.savefig(image_dir+varname+'_'+quanti+'_Veg_CNTRL.png',dpi=300)
        plt.close()
        

        # compare reduction by  Seagrasss periods absolute
        Ddata=[data1b-data1,data2b-data2]
        clabel='$\Delta$'+label0
        fig,axs,cbar_ax,cbar_ax2=compare_slabs_with_diff(s,Ddata,names=names,cmap=plt.cm.turbo,cblabel=clabel,cblabel2='$\Delta$ '+clabel,axis_limit=axis_limit,figsize=(10,3.6))
        fig.subplots_adjust( hspace=0.1,top=1.2)		
        plt.savefig(image_dir+varname+'_'+quanti+'_storms_Veg_max_minus_CNTRL.png',dpi=300)
        plt.close()


        
        # compare reduction by  Seagrasss periods relative
        Ddata=[(data1b-data1)/data1*100,(data2b-data2)/data2*100]
        clabel='$\Delta$'+label0
        clabel=clabel[:-4]+'[%]'

        # replace /0 nans by nan
        for item in Ddata:
            iinf=np.isinf(item)
            item[iinf]=np.nan
        #np.isinf(Ddata).max()
        # fix extremse
        for item in Ddata:        
            inan=np.abs(item)>250
            item[inan]=np.nan        
        #fig,axs,cbar_ax=compare_slabs(s,Ddata,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=clabel,axis_limit=axis_limit,figsize=(8,3.6))
        fig,axs,cbar_ax,cbar_ax2=compare_slabs_with_diff(s,Ddata,names=names,cmap=plt.cm.turbo,cblabel=clabel,cblabel2='$\Delta$ '+clabel,axis_limit=axis_limit,figsize=(10,3.6))        
        fig.subplots_adjust( hspace=0.1,top=1.2)		
        plt.savefig(image_dir+varname+'_'+quanti+'_storms_Veg_max_minus_CNTRL_rel.png',dpi=300)
        plt.close()


### scatter events 
scatter_dir=image_dir+'scatter/'
if not os.path.exists(scatter_dir):
    os.mkdir(scatter_dir)

for varname in varnames:
    vari=varname
    unit=units[vari]

    for quanti in quantaties:

        print('doing '+ ' '.join((vari,quanti)))

        # build uint
        symbol=symbols[vari]
        if quanti=='mean':
            symbol=symbol.join(('<','>'))
        else:	
            symbol=symbol.join((quanti.replace('quantile','prc')+'(',')'))
        label0='{:s} [{:s}]'.format(symbol,unit)
        clabel='$\Delta$'+label0
        
        #change name order in creation file
        varname_temp=varname.replace('_bottom','')
        
        if ('bottom' in varname) & (varname!='bottomStress'):
            files=glob.glob(ncdir1+'*'+varname_temp+'*'+quanti+'*')
            filename=files[0].split(ncdir1)[1]
            ds1=xr.open_dataset(ncdir1+filename)
            variable=list(ds1.variables.keys())[0]
        elif 'surface' in varname:            
            files=glob.glob(ncdir1+'*'+varname_temp+'*'+quanti+'*')
            filename=files[1].split(ncdir1)[1]
            ds1=xr.open_dataset(ncdir1+filename)
            variable=list(ds1.variables.keys())[0]
        else:
            filename=varname+'_'+quanti+'_1.nc'
            variable=varname+'_'+quanti
        
        ds1=xr.open_dataset(ncdir1+filename)
        variable=list(ds1.variables.keys())[0]
        data1=ds1[variable].values

        ds2=xr.open_dataset(ncdir2+filename)
        variable=list(ds2.variables.keys())[0]
        data2=ds2[variable].values


        ds1b=xr.open_dataset(ncdir1b+filename)
        variable=list(ds1b.variables.keys())[0]
        data1b=ds1b[variable].values

        ds2b=xr.open_dataset(ncdir2b+filename)
        variable=list(ds2b.variables.keys())[0]
        data2b=ds2b[variable].values

        # compare no Seagrasss periods absolute

        data=[data1,data2]
        # compare reduction by  Seagrasss periods absolute
        Ddata=[data1b-data1,data2b-data2]


        # make scatter plots on reduction
                
        Ddata=[data1b-data1,data2b-data2]        

        V1=Ddata[0][s.nvplt].mean(axis=1)[inbox_elems]
        D1=data1[s.nvplt].mean(axis=1)[inbox_elems]

        V2=Ddata[1][s.nvplt].mean(axis=1)[inbox_elems]
        D2=data2[s.nvplt].mean(axis=1)[inbox_elems]


        Dep=np.asarray(s.depths)[s.nvplt].mean(axis=1)[inbox_elems]

        isveg=(Dep<4) & (Dep>-4)




        plt.figure(figsize=(12,7))
        plt.clf()
        ax1=plt.subplot(2,2,1)
        plt.plot(Dep,V1,'.')
        plt.xlabel('depth [m]')
        plt.ylabel(clabel)
        plt.title('2017')
        plt.xlim((-5,15))
        plt.grid()
        m1=np.mean(V1[isveg])
        q1=np.quantile(V1[isveg],0.05)
        plt.hlines(m1,-4,4,color='k')
        plt.hlines(q1,-4,4,color='r')
        plt.vlines(-4,-2,1.5,color='g')
        plt.vlines(4,-2,1.5,color='g')
        plt.legend(['Data','Mean: {:.3f}'.format(m1),'Q95: {:.3f}'.format(q1),'Veg'],frameon=False)

        ax3=plt.subplot(2,2,3)
        #plt.vlines(4,-2,1.5,color='g')
        plt.plot(D1[isveg],V1[isveg],'.',color='g')
        plt.plot(D1[~isveg],V1[~isveg],'.',color='gray')
        m1, b1 = np.polyfit(D1[isveg],V1[isveg], 1)
        m2, b2 = np.polyfit(D1[~isveg],V1[~isveg], 1)
        xlim=plt.xlim()
        xq=np.linspace(0,xlim[1])
        plt.plot(xq,m1*xq+b1,color='k',linewidth=2)
        plt.plot(xq,m2*xq+b2,color='b',linewidth=2)
        plt.legend(['Veg','VegFree','{:.3f} x + {:.3f}'.format(m1,b1),'{:.3f} x + {:.3f}'.format(m2,b2)],frameon=False)

        #plt.xlabel('magnitude')
        #plt.ylabel('reduction')
        plt.xlabel(label0)
        plt.ylabel(clabel)

        plt.grid()
        ax2=plt.subplot(2,2,2,sharey=ax1)
        plt.plot(Dep,V2,'.')
        plt.xlim((-5,15))
        m1=np.mean(V2[isveg])
        q1=np.quantile(V2[isveg],0.05)
        plt.hlines(m1,-4,4,color='k')
        plt.hlines(q1,-4,4,color='r')
        plt.vlines(-4,-2,1.5,color='g')
        plt.vlines(4,-2,1.5,color='g')
        plt.legend(['Data','Mean: {:.3f}'.format(m1),'Q95: {:.3f}'.format(q1),'Veg'],frameon=False)
        plt.xlabel('depth [m]')
        plt.ylabel('reduction')
        plt.title('2090')
        plt.xlim((-5,15))
        plt.grid()
        ax4=plt.subplot(2,2,4,sharey=ax3,sharex=ax3)
        plt.plot(D2[isveg],V2[isveg],'.',color='g')
        plt.plot(D2[~isveg],V2[~isveg],'.',color='gray')
        m1, b1 = np.polyfit(D2[isveg],V2[isveg], 1)
        m2, b2 = np.polyfit(D2[~isveg],V2[~isveg], 1)
        xlim=plt.xlim()
        xq=np.linspace(0,xlim[1])
        plt.plot(xq,m1*xq+b1,color='k',linewidth=2)
        plt.plot(xq,m2*xq+b2,color='b',linewidth=2)
        plt.legend(['Veg','VegFree','{:.3f} x + {:.3f}'.format(m1,b1),'{:.3f} x + {:.3f}'.format(m2,b2)],frameon=False)
        plt.xlabel(label0)
        plt.ylabel(clabel)
        plt.grid()
        plt.suptitle(' '.join((vari,quanti)))
        plt.tight_layout()
        plt.savefig(scatter_dir+varname+'_'+quanti+'_Veg_max_minus_CNTRL_rel.png',dpi=300)
        plt.close()