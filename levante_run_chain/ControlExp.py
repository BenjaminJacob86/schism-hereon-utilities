# example illustrating grid and ploting hovmöller diagramm
#export OMP_NUM_THREADS=1 
from glob import glob
import os
import sys
sys.path.insert(0,'/home/g/g260114/schism-hzg-utilities/')
import xarray as xr
from matplotlib import pyplot as plt
from schism import* # import schism class to read grid structure
s=schism_setup() #read grid_strtucture
# hgrid.ll has grid information in lon lat or colocation
import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
plt.ion()



######################## Settings ############################
basedir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/'
image_outdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Analyse_storm/'
if not os.path.isdir(image_outdir):
	os.mkdir(image_outdir)
#glob.glob('Veg_*')
#ref='Veg_LE', 'Veg_CNTRL', 'Veg_CNTRL_no_wave', 'Veg_REF', 'Veg_HE', 'Veg_max'


os.chdir(basedir)
s=schism_setup()


ref='Veg_REF'
control='Veg_CNTRL'
experiment='Veg_max'


#expnames=['Ref','Blank','Veg$_{max}$','Veg$_{LE}$','Veg$_{HE}$']	
expnames=['Ref','Blank','Veg$_{max}$','Veg$_{HE}$','Veg$_{LE}$']	
experiments=[ref,control,'Veg_max','Veg_HE','Veg_LE'] #,experiment]
#ncdirs=[basedir+exp+'/outputs01/'  for exp in experiments]
ncdirs=[basedir+exp+'/outputs_merged/'  for exp in experiments]

plt.ion()



a=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/outputs01/out2d_27.nc')
b=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/outputs/out2d_27.nc')




a=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_max_control/outputs01/out2d_23.nc')
b=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_max_control/outputs_regular_hotstart/out2d_23.nc')
c=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_max_control/outputs_hotsart_sel2/out2d_23.nc')
d=xr.open_dataset('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_max_control/outputs/out2d_23.nc')


varname='sigWaveHeight'

ti=0
data1=a[varname][ti,:].values
data2=b[varname][ti,:].values
data3=c[varname][ti,:].values
data4=d[varname][ti,:].values

#sieht nict aus wie initial condition

# hotstart not equal
plt.close('all')
fig, axs = plt.subplots(nrows=1,ncols=3)
ph,ax,cax=s.plotAtnodes(data1,ax=axs[0])
ph.set_clim((0,0.75))
ph,ax,cax=s.plotAtnodes(data2,ax=axs[1])
ph.set_clim((0,0.75))


data1=a[varname][ti,:].values
data2=b[varname][ti,:].values

#sieht nict aus wie initial condition

# hotstart not equal
varname='elevation'

# nach 24 wieder

ti=0
data1=a[varname][ti,:].values # normal run
data2=b[varname][ti,:].values # regular restart
data3=c[varname][ti,:].values # restart select hot2
data4=d[varname][ti,:].values # restat select hot1
plt.close('all')
fig, axs = plt.subplots(nrows=1,ncols=4)
ph,ax,cax=s.plotAtnodes(data1,ax=axs[0],cmap=plt.cm.turbo)
axs[0].set_title('run continous')
ph2,ax2,cax2=s.plotAtnodes(data2,ax=axs[1],cmap=plt.cm.turbo)
axs[1].set_title('2 ts hotfile')
ph3,ax3,cax3=s.plotAtnodes(data3,ax=axs[2],cmap=plt.cm.turbo)
axs[2].set_title('hot 2st extracted')
ph4,ax4,cax4=s.plotAtnodes(data4,ax=axs[3],cmap=plt.cm.turbo)
axs[3].set_title('hot 1st extracted')

ph.set_clim((0,0.75))
ph2.set_clim((0,0.75))

# it looks for hotstart file
# selection correspondends to time extrat
# it seemling does not check time


# written hotstart is wrnong?
# ocean time

# initialisation must be wrong
#removed second entry


#start time set correct
# boundary forcing red xcorrect


# Falscher hotstart?


#hotstart entspricht nicht run


plt.subplot(2,1,2)

ph.set_clim((0,0.75))
