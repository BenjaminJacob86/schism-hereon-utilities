"""
Compare schism against monthly mean SPM from CMEMS sattelite product:
https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description
"""
import os
import sys
import pandas as pd
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
from scipy.spatial import cKDTree


rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir=rundir+'outputs_all/'

rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/SPM_tune/'
ncdir=rundir+'outputs_all_temp/'

spmdir='/work/gg0028/g260114/DATA/CMEMS/SPM/'       # CMEMS spm data monthly
varlist=['out2d','totalSuspendedLoad']

varname='totalSuspendedLoad' # compare to data
year=2017


# units schism : g / L # schism
# units data:  g/m³   monthly mean
factor=1000            # factor schis#-> *1000     
label='SPM  [g/m^3]'

limL=0
lumU=150


################# Program Start #######################################

os.chdir(rundir)
s=schism_setup()
s.init_node_tree(latlon=True)
s.ds=schism_outputs_by_variable(ncdir,max_stack=-1,varlist=varlist).ds



## get time information
p=param()
reftime = np.asarray(dt.datetime(int(p.get_parameter('start_year')),
                                 int(p.get_parameter('start_month')),
                                 int(p.get_parameter('start_day')),
                                 int(p.get_parameter('start_hour')), 0, 0), np.datetime64)
t=s.ds['out2d'].time.values
#convert in seconds
dates = reftime + np.asarray(t, np.timedelta64(1, 's'))
years = dates.astype('datetime64[Y]').astype(int) + 1970
months = dates.astype('datetime64[M]').astype(int) % 12 + 1
days = (dates - dates.astype('datetime64[M]') + 1) / np.timedelta64(1, 'D')



# subset
def get_timesteps_seconds(date0,date1):
    """ select time values in seconds since start corresponding to dates date0 and date1 """
    #date0 = np.asarray(pd.to_datetime(date0, format='%Y%m%d'), np.datetime64)
    #date1 = np.asarray(pd.to_datetime(date1, format='%Y%m%d'), np.datetime64)

    # select time slices
    # date0=dt.datetime(startdate)

    index = (dates >= date0) & (dates <= date1)
    levels = [len(s.vgrid[1]) - 1, 0]
    label = ['surface', 'bottom']
    # schism time indices for selectio in seconds since starts (depends on schism version)
    t0sec = t[index][0]
    t1sec = t[index][-10]
    return t0sec,t1sec
######################


## load spm data
spm_files=np.sort(glob.glob(spmdir+'*'+str(year)+'*.nc'))
dspm=xr.open_mfdataset(spm_files)


# limit data to modeldomain
lonmin=np.min(s.lon)-1
lonmax=np.max(s.lon)+1
latmin=np.min(s.lat)-1
latmax=np.max(s.lat)+1

dspm=dspm.sel(lon=slice(lonmin,lonmax),lat=slice(latmax,latmin)) # !!inverse lat lon order
lon,lat=dspm.lon.values,dspm.lat.values
LON,LAT=np.meshgrid(lon,lat)
##

## nearest neighbour indices
nntree = cKDTree(list(zip(LON.flatten(),LAT.flatten())))
nns=nntree.query(list(zip(np.asarray(s.lon),np.asarray(s.lat))))[1]


#
##means=[]
#for mon in range(1,12+1):
#    date0=np.datetime64('2017-{:02d}-01'.format(mon))
#    date1=np.datetime64('{:d}-{:02d}-01'.format(year+(mon==12),mon+1))
#    dspmi=dspm.sel(time=slice(date0,date1-np.timedelta64(1,'D'))) #select data in range
#    t0sec,t1sec=get_timesteps_seconds(date0,date1) #select data in range
#    #monthly mean
#    dssi=s.ds[varname].sel(time=slice(t0sec,t1sec),nSCHISM_vgrid_layers=-1)
#    means.append(dssi.mean(dim='time')[varname].values)
#

month_avail=np.unique(months)
keep=[days[np.where(months==mon)[0][-1]]>29 for mon in month_avail]
month_avail=month_avail[keep]


#month loop
means=[]
for mon in month_avail:
    print(mon)
 
    date0=np.datetime64('2017-{:02d}-01'.format(mon))
    date1=np.datetime64('{:d}-{:02d}-01'.format(year+(mon==12),mon+1))


    dspmi=dspm.sel(time=slice(date0,date1-np.timedelta64(1,'D'))) #select data in range
    t0sec,t1sec=get_timesteps_seconds(date0,date1) #select data in range

    #monthly mean
    dssi=s.ds[varname].sel(time=slice(t0sec,t1sec),nSCHISM_vgrid_layers=-1)
    mean_schism=dssi.mean(dim='time')[varname].values
    means.append(mean_schism)

    # meshgrid for binning
    plt.clf()
    plt.subplot(2,2,1)
    ph0=dspmi.SPM.plot()
    ph0.set_clim((limL,lumU))
    plt.subplot(2,2,2)
    ph,ch,ax=s.plotAtnodes(mean_schism*factor)
    ax.set_label(label)
    ph.set_clim((limL,lumU))
    ch.set_label(label)
    plt.subplot(2,2,3)
    ph,ch,ax=s.plotAtnodes(mean_schism*factor- dspmi.SPM.values[0,:].flatten()[nns],cmap=plt.cm.jet)
    ph.set_clim((-50,50))
    ch.set_label('\Delta '+label)

    plt.subplot(2,2,4)
    ph,ch,ax=s.plotAtnodes(dspmi.SPM.values[0,:].flatten()[nns])
    ch.set_label('interp_data')
    plt.suptitle('Monthly mean '+str(year)+'\{:02d}'.format(mon))
    plt.tight_layout()
    
    plt.savefig('Monthly_mean'+str(year)+'{:02d}.png'.format(mon))

     
