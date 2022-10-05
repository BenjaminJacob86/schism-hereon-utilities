#
import os
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
import numpy as np
from schism import *
from matplotlib import pyplot as plt


def get_parameter(dir,param='dt'):
	f=open(dir+'param.nml')
	for line in f.readlines():
		if param+' =' in line:
			param= line.split('=')[1].split('!')[0]
			try:
				param=float(param)
			except:
				pass
			break		
	return param

plt.ion()

dir1='/gpfs/work/jacobb/data/RUNS/IDEAL/BSidealSlopeSmoothBasinTriCoarsen1/' #'/gpfs/work/jacobb/data/RUNS/IDEAL/BSidealTransition/'
dir2='/gpfs/work/jacobb/data/RUNS/IDEAL/BSidealSlopeSmoothBasinTriCoarsen2/'
 #'/gpfs/work/jacobb/data/RUNS/IDEAL/BSidealTransition2xRes/'

 
name1='BasinTriCoarsen1'
name2='BasinTriCoarsen2' 

# variable and time step to compare
ti=24*4
varname='salt'

 
#dir3='/gpfs/work/jacobb/data/RUNS/IDEAL/BSidealTransition2xRes_onlyPave/'


	
	
dt1=get_parameter(dir1,param='dt')
dt2=get_parameter(dir2,param='dt')


os.chdir(dir1)
ncdir1=dir1+'/v0/combined/'
s1=schism_setup()
os.chdir(dir2)
s2=schism_setup()
ncdir2=dir2+'combined/'

#ncdir3=dir3+'combined/'
#os.chdir(dir3)
#s3=schism_setup()


r1=np.asarray(list(s1.resolution_by_element.values()))
r2=np.asarray(list(s2.resolution_by_element.values()))
#r3=np.asarray(list(s3.resolution_by_element.values()))

combined=np.concatenate((r1,r2))
vmin=np.min(combined)
vmax=np.max(combined)

plt.clf()
plt.subplot(2,1,1)
ph,ch=s1.plotAtelems(r1[s1.nvplt2nvp],latlon=False)
plt.title(name2)
ch.set_label('Resolution [m]')
plt.clim(vmin,vmax)
plt.subplot(2,1,2)
ph,ch=s2.plotAtelems(r2[s2.nvplt2nvp],latlon=False)
plt.title(name2)
ch.set_label('Resolution [m]')
plt.clim(vmin,vmax)
plt.savefig('resolution.png',dpi=300)

plt.clf()
cfl1=s1.compute_cfl(dt=dt1)
cfl2=s2.compute_cfl(dt=dt2)
combined=np.concatenate((cfl1,cfl2))
vmin=np.min(combined)
vmax=np.max(combined)
plt.subplot(2,1,1)
ph,ch=s1.plotAtelems(cfl1[s1.nvplt2nvp],latlon=False)
plt.title(name1)
ch.set_label('cfl')
plt.clim((vmin,vmax))
plt.subplot(2,1,2)
ph,ch=s2.plotAtelems(cfl2[s2.nvplt2nvp],latlon=False)
plt.title(name2)
ch.set_label('cfl')
plt.clim(vmin,vmax)
plt.clim((vmin,vmax))
plt.savefig('cfl.png',dpi=300)

plt.clf()
plt.subplot(2,1,1)
s1.plot_domain_boundaries(latlon=False,append=True)
plt.title(name1)
s1.plot_mesh(latlon=False)
plt.axis('equal')
plt.subplot(2,1,2)
s2.plot_domain_boundaries(latlon=False,append=True)
plt.title(name2)
s2.plot_mesh(latlon=False)
plt.axis('equal')
plt.savefig('mesh.png',dpi=300)


nz1=np.asarray([s1.vgrid[key].count() for key in s1.vgrid.keys()])
nz2=np.asarray([s2.vgrid[key].count() for key in s2.vgrid.keys()])

plt.clf()
plt.subplot(2,1,1)
ph,ch=s1.plotAtnodes(nz1,latlon=False)
#s1.plot_mesh(latlon=False)
ch.set_label('# vert')
plt.title(name1)
plt.subplot(2,1,2)
ph,ch=s2.plotAtnodes(nz2,latlon=False)
#s2.plot_mesh(latlon=False)
ch.set_label('# vert')
#plt.clim(vmin,vmax)
plt.title(name2)
plt.savefig('nvert.png',dpi=300)



s1.nc=schism_output2(ncdir1).nc
s2.nc=schism_output2(ncdir2).nc


val1=s1.nc[varname][ti,:,-1].values
val2=s2.nc[varname][ti,:,-1].values

plt.clf()
plt.subplot(2,1,1)
ph,ch=s1.plotAtnodes(val1,latlon=False)
plt.clim((vmin,vmax))
ch.set_label(varname)
plt.axis('equal')
plt.title(name1)
plt.axis(ax)
plt.subplot(2,1,2)
ph,ch=s2.plotAtnodes(val2,latlon=False)
plt.clim((vmin,vmax))
ch.set_label(varname)
plt.title(name2)
plt.suptitle(str(s1.nc['time'][ti].values)[:19])
plt.axis('equal')
plt.axis(ax)
plt.savefig(varname+str(ti)+'png',dpi=300)







#s3.nc=schism_output2(ncdir3).nc


plt.figure()

varname='salt'


ti=0

val1=s1.nc[varname][ti,:,-1].values
val2=s2.nc[varname][ti,:,-1].values

vmin=np.min(np.concatenate((val1,val2)))
vmax=np.max(np.concatenate((val1,val2)))+2

name1='R4kmB100mC'
name2='R2kmB50mC'
name3='R4kmB100mCTri'

ti=24*15
val1=s1.nc[varname][ti,:,-1].values
val2=s2.nc[varname][ti,:,-1].values
val3=s3.nc[varname][ti,:,-1].values

vmin=np.min(np.concatenate((val1,val2)))
vmax=np.max(np.concatenate((val1,val2)))


ti=24*2
val1=s1.nc['hvel'][ti,:,-1,0].values
val2=s2.nc['hvel'][ti,:,-1,0].values


varname='Ux'
plt.clf()
plt.subplot(2,1,1)
ph,ch=s1.plotAtnodes(val1,latlon=False)
plt.clim((vmin,vmax))
ch.set_label(varname)
plt.axis('equal')
plt.title(name1)
plt.axis(ax)
plt.subplot(2,1,2)
ph,ch=s2.plotAtnodes(val2,latlon=False)
plt.clim((vmin,vmax))
ch.set_label(varname)
plt.title(name2)
plt.suptitle(str(s1.nc['time'][ti].values)[:19])
plt.axis('equal')
plt.axis(ax)


Hsurf1=(s1.nc['zcor'][ti,:,-1]-s1.nc['zcor'][ti,:,-2]).values
Hsurf2=(s2.nc['zcor'][ti,:,-1]-s2.nc['zcor'][ti,:,-2]).values

vmin=np.min(np.concatenate((Hsurf1,Hsurf2)))
vmax=np.max(np.concatenate((Hsurf1,Hsurf2)))

plt.clf()
plt.subplot(2,1,1)
ph,ch=s1.plotAtnodes(Hsurf1,latlon=False)
plt.clim((vmin,vmax))
ch.set_label('H surf layer')
plt.axis('equal')
plt.title(name1)
#plt.axis(ax)
plt.axis(ax)
plt.subplot(2,1,2)
ph,ch=s2.plotAtnodes(Hsurf2,latlon=False)
plt.clim((vmin,vmax))
ch.set_label('H surf layer')
plt.title(name2)
plt.suptitle(str(s2.nc['time'][ti].values)[:19])
plt.axis('equal')
plt.axis(ax)


s1.init_node_tree(latlon=False)
dd,nn=s1.node_tree_xy.query(list(zip(s2.x,s2.y)))

plt.figure()
s2.plotAtnodes(val2-val1[nn],latlon=False)





