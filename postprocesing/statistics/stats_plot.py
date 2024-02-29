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
rundir1='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir1=rundir1+'outputs_all2/'

os.chdir(rundir1)
s=schism_setup()


# compare storm max levels
rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir1='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/stats_extract_Herwart/Veg_CNTRL/'
ncdir2='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/stats_extract2/Veg_CNTRL/'

ncdir1b='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/stats_extract_Herwart/Veg_max/'
ncdir2b='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/stats_extract2/Veg_max/'



#seasonal
rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir1='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/extract_JF/Veg_CNTRL/'
ncdir2='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/extract_JA/Veg_CNTRL/'

ncdir1b='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/extract_JF/Veg_CNTRL/Veg_max2/'
ncdir2b='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/extract_JA/Veg_max2/'




 
    
 
#get data

## statistics to compute: # potential quantaties (percentiles, eg.   'quantile95','quantile5', ...
quantaties = ['mean', 'min', 'max', 'std', 'quantile95']  # ['mean','min','max','std','quantile95']

## variables to analyse
# vector vars are analysed  by magnitude computed from compoent files
varnames = ['dryFlagNode', 'elevation', 'sigWaveHeight', 'totalSuspendedLoad']
vector_vars = ['depthAverageVel', 'bottomStress']  # add variables with x and y component

# reduce varlist selection to speed up loading
varlist = ['out2d', 'turbulentKineticEner', 'totalSuspendedLoad', 'horizontalVelX', 'horizontalVelY']

ds=schism_outputs_by_variable(ncdir1, max_stack=-1, varlist=varlist)
ds2=schism_outputs_by_variable(ncdir2, max_stack=-1, varlist=varlist)



###### get time in datetime	and time subsets
def get_dates(setupdir,access):
    """ get dates from SCHISM seconds since start based on param.nml """
    import datetime as dt
    p = param(setupdir + '/')
    reftime = np.asarray(dt.datetime(int(p.get_parameter('start_year')),
                                     int(p.get_parameter('start_month')),
                                     int(p.get_parameter('start_day')),
                                     int(p.get_parameter('start_hour')), 0, 0), np.datetime64)
    a = access
    t = a.get('elevation').time.values
    startdate = reftime + 3600 * np.timedelta64(1, 's')
    dates = reftime + np.asarray(t, np.timedelta64(1, 's'))
    #years = dates.astype('datetime64[Y]').astype(int) + 1970
    #months = dates.astype('datetime64[M]').astype(int) % 12 + 1
    #days = (dates - dates.astype('datetime64[M]') + 1) / np.timedelta64(1, 'D')
    return dates

dates1=get_dates(rundir1,ds)
dates2=get_dates(rundir2,ds2)


os.chdir('/work/gg0028/g260114/RUNS/Rest-Coast/')
s.read_reg('efws_s')

# make Path from region coordinates for point inside checkings	
plt.ion()
s.create_region()
from matplotlib.path import Path		
xy=s.reg['efws_s']
RegPoly=Path(list(zip(xy[:,0],xy[:,1])))
grid_coords=list(zip(s.lon,s.lat))
pt_in_poly=RegPoly.contains_points(grid_coords)	
# seagrass <=4

e=ds.get('elevation')
e2=ds2.get('elevation')
efws_m1=np.hstack([np.mean(e[ti,:].values[pt_in_poly]) for ti in range(len(dates1)) ]) #
efws_m2=np.hstack([np.mean(e2[ti,:].values[pt_in_poly]) for ti in range(len(dates2)) ]) #


# Just by points cleand would be by area

# 2d
e=ds.get('elevation')
e2=ds2.get('elevation')
efws_m1=np.hstack([np.mean(e[ti,:].values[pt_in_poly]) for ti in range(len(dates1)) ]) #
efws_m2=np.hstack([np.mean(e2[ti,:].values[pt_in_poly]) for ti in range(len(dates2)) ]) #

#2d vector
vals1=[(ds.ds['bottomStress']['bottomStress']**2).sum(axis=0)[ti].values[pt_in_poly].mean() for ti in range(len(dates1))]
vals2=[(ds2.ds['bottomStress']['bottomStress']**2).sum(axis=0)[ti].values[pt_in_poly].mean() for ti in range(len(dates2))]

#3d vector

# stats curver

q5_1=np.quantile(efws_m1,.05)
m1=np.mean(efws_m1)
q95_1=np.quantile(efws_m1,.95)

q5_2=np.quantile(efws_m2,.05)
m2=np.mean(efws_m2)
q95_2=np.quantile(efws_m2,.95)



varnames=['elevation','sigWaveHeight','depthAverageVel','totalSuspendedLoad_bottom','bottomStress']
units={'elevation':'m','sigWaveHeight':'m','depthAverageVel':'m/s','totalSuspendedLoad_bottom':'g/L','bottomStress':'Pa'}
symbols={'elevation':'$\zeta$','sigWaveHeight':'HS','depthAverageVel':'$u_{<z>}$','totalSuspendedLoad_bottom':'spm','bottomStress':'$\\tau_{btm}$'}

varname='elevation'
ylab=symbols[varname]+' '+units[varname]
data1,data2=box_compare(ds,ds2,varname=varname,ylab=ylab)


# save time series varname

prefix='Veg_CNTRL'
for varname in varnames:

    print(varname)
    
    e=ds.get(varname)
    e2=ds2.get(varname)
    data1=np.hstack([np.mean(e[ti,:].values[pt_in_poly]) for ti in range(len(dates1)) ]) #
    data2=np.hstack([np.mean(e2[ti,:].values[pt_in_poly]) for ti in range(len(dates2)) ]) #

    ds_efws1=xr.Dataset(
    data_vars=dict(
    efws_bstress=(["time"],np.asarray(data1)),
   

    coords=dict(time=dates1),
    attrs=dict(description="Wadden Sea spatial average")
    )

    ds_efws2=xr.Dataset(
    data_vars=dict(
    efws_bstress=(["time"],np.asarray(data2)),
    ),
    coords=dict(time=dates2),
    attrs=dict(description="Wadden Sea spatial average")
    )
    
    ds_efws1.to_netcdf(prefix+'efws'+year+varname+'.nc')
    ds_efws1.to_netcdf(prefix+'efws'++varname+'.nc')



def box_compare(ds,ds2,varname='elevation',ylab='ylab'):


    e=ds.get(varname)
    e2=ds2.get(varname)
    data1=np.hstack([np.mean(e[ti,:].values[pt_in_poly]) for ti in range(len(dates1)) ]) #
    data2=np.hstack([np.mean(e2[ti,:].values[pt_in_poly]) for ti in range(len(dates2)) ]) #


    ds_efws1=xr.Dataset(
    data_vars=dict(
    efws_bstress=(["time"],np.asarray(data1)),
    ),
    coords=dict(time=dates1),
    attrs=dict(description="Wadden Sea spatial average")
    )


    ds_efws2=xr.Dataset(
    data_vars=dict(
    efws_bstress=(["time"],np.asarray(data2)),
    ),
    coords=dict(time=dates2),
    attrs=dict(description="Wadden Sea spatial average")
    )


    g1=ds_efws1.groupby('time.month')
    g2=ds_efws2.groupby('time.month')


    plt.figure()
    ax1=plt.subplot(2,1,1)
    plt.boxplot(monthly,patch_artist = True,
               boxprops = dict(facecolor = "lightblue"))           
    plt.ylabel(ylab)           
    plt.grid()           
    ax2=plt.subplot(2,1,2,sharey=ax1)
    plt.boxplot(monthly2,patch_artist = True,
               boxprops = dict(facecolor = "lightgreen"))           
    plt.grid()
    plt.ylabel(ylab)

    return data1,data2





[g1[1]

def moments(vals):
    return np.quantile(vals,.05),np.mean(vals),np.quantile(vals,.95)


q5_1,m1,q95_1=moments(vals1)
q5_2,m2,q95_2=moments(vals2)


monthly=[g1[i]['efws_bstress'] for i in range(1,13)]
monthly2=[g2[i]['efws_bstress'] for i in range(1,13)]


plt.subplot(2,1,1)
plt.boxplot(x=np.arange(1,12)-0.5,data=monthly,patch_artist = True,
           boxprops = dict(facecolor = "lightblue"))









           
plt.boxplot(monthly2,color='r')

# 11


def ts_moment_plot(dates1,dates2,vals1,vals2):
    q5_1,m1,q95_1=moments(vals1)
    q5_2,m2,q95_2=moments(vals2)
    
    # plot ts    
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(dates1,vals1)
    plt.hlines(q5_1,dates1[0],dates1[-1],'r')
    plt.hlines(m1,dates1[0],dates1[-1],'r')
    plt.hlines(q95_1,dates1[0],dates1[-1],'r')
    plt.legend(['ts','q05','mean','q95'],ncol=4,frameon=False)
    plt.subplot(2,1,2)
    plt.plot(dates2,vals2)
    plt.hlines(q5_2,dates2[0],dates2[-1],'r')
    plt.hlines(m2,dates2[0],dates2[-1],'r')
    plt.hlines(q95_2,dates2[0],dates2[-1],'r')
    
    #bar
    plt.figure()
    plt.clf()
    plt.bar(-.125+np.arange(3),[q5_1,m1,q95_1],width=0.25)
    plt.bar(0.125+np.arange(3),[q5_2,m2,q95_2],width=0.25)
    plt.xticks([0,1,2,])
    plt.gca().set_xticklabels(['Q05','mean','Q95'])
    plt.grid()

plt.close('all')
ts_moment_plot(dates1,dates2,vals1,vals2)


plt.figure()
plt.boxplot([vals1,vals2])
plt.grid()

#Quartille
#25 ,25
75

q5_1=np.quantile(efws_m1,.05)
m1=np.mean(efws_m1)
q95_1=np.quantile(efws_m1,.95)

q5_2=np.quantile(efws_m2,.05)
m2=np.mean(efws_m2)
q95_2=np.quantile(efws_m2,.95)



plt.figure()
plt.subplot(2,1,1)
plt.plot(dates1,efws_m1)
plt.hlines(q5_1,dates1[0],dates1[-1],'r')
plt.hlines(m1,dates1[0],dates1[-1],'r')
plt.hlines(q95_1,dates1[0],dates1[-1],'r')
plt.subplot(2,1,2)
plt.plot(dates2,efws_m2)
plt.hlines(q5_2,dates2[0],dates2[-1],'r')
plt.hlines(m2,dates2[0],dates2[-1],'r')
plt.hlines(q95_2,dates2[0],dates2[-1],'r')

plt.figure()
plt.clf()
plt.bar(-.125+np.arange(3),[q5_1,m1,q95_1],width=0.25)
plt.bar(0.125+np.arange(3),[q5_2,m2,q95_2],width=0.25)
plt.xticks([0,1,2,])
plt.gca().set_xticklabels(['Q05','mean','Q95'])
plt.grid()



plt.figure()
plt.subplot(2,1,1)
#plt.cla()
plt.plot(dates1,vals1)
#plt.hlines(q5_1,dates1[0],dates1[-1],'r')
#plt.hlines(m1,dates1[0],dates1[-1],'r')
#plt.hlines(q95_1,dates1[0],dates1[-1],'r')
plt.subplot(2,1,2)
plt.plot(dates2,vals2)
plt.title()






access = [ for ncdiri in ncdirs]



ncdir2='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/stats_extract2/Veg_CNTRL/'

ncdir1b='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/stats_extract_Herwart/Veg_max/'
ncdir2b='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/stats_extract2/Veg_max/'








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



#
image_dir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/herwart_vs_xaver/'
if not os.path.exists(image_dir):
    os.mkdir(image_dir)


names=['Herwart','Xaver']
axis_limit=regions['EFWS']['axis_limit']



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
        plt.savefig(image_dir+varname+'_'+quanti+'_storms_Veg_CNTRL.png',dpi=300)
        plt.close()

        # compare reduction by  Seagrasss periods absolute
        Ddata=[data1b-data1,data2b-data2]
        clabel='$\Delta$'+label0
        fig,axs,cbar_ax=compare_slabs(s,Ddata,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=clabel,axis_limit=axis_limit,figsize=(8,3.6))
        fig.subplots_adjust( hspace=0.1,top=1.2)		
        plt.savefig(image_dir+varname+'_'+quanti+'_storms_Veg_max_minus_CNTRL.png',dpi=300)
        plt.close()

        # compare reduction by  Seagrasss periods relative
        Ddata=[(data1b-data1)/data1*100,(data2b-data2)/data2*100]
        clabel='$\Delta$'+label0
        clabel=clabel[:-4]+'[%]'
        fig,axs,cbar_ax=compare_slabs(s,Ddata,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=clabel,axis_limit=axis_limit,figsize=(8,3.6))
        fig.subplots_adjust( hspace=0.1,top=1.2)		
        plt.savefig(image_dir+varname+'_'+quanti+'_storms_Veg_max_minus_CNTRL_rel.png',dpi=300)
        plt.close()


#ncols=2
#x0=1/(2*ncols)
#w=1/ncols
#cbar_ax.set_position((x0,0.1,w,0.02))
#plt.tight_layout()

# difference of seagrass during storms
fig,axs,cbar_ax=compare_slabs(s,data,names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=clabel,axis_limit=axis_limit,figsize=(8,3.6))
fig.subplots_adjust( hspace=0.1,top=1.2)		

plt.savefig(image_dir+'Zeta_max_storms_Veg_CNTRL.png',dpi=300)

### differences:


fig,axs,cbar_ax=compare_slabs(s,[data1,data1b],names=names,cmap=plt.cm.turbo,nrows=1,ncols=2,cblabel=label0,axis_limit=axis_limit,figsize=(8,3.6))
fig.subplots_adjust( hspace=0.1,top=1.2)		