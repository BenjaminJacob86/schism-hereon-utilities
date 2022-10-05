import numpy as np
from matplotlib import pyplot as plt
import sys
import os
sys.path.insert(0,'//gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *

setupdir=os.getcwd()+'/' #'/gpfs/work/jacobb/data/RUNS/IDEAL/BSideal/'
outdir=setupdir+'combined/'
s=schism_setup()

#trans=bp_file('transect.bp')
s.nc=schism_output2(outdir).nc

# get transect
#x=np.hstack((np.arange(-50000,-15000,5000),np.arange(-15000,15000,1000),np.arange(15000,70000,5000)))
#x=np.hstack((np.arange(-50000,-15000,4000),np.arange(-15000,15000,1000),np.arange(15000,70000,5000)))
#x=np.arange(-50000,70000,100)#np.hstack((,np.arange(-15000,15000,1000),np.arange(15000,70000,5000)))


x=np.arange(-100000,120000,100)#np.hstack((,np.arange(-15000,15000,1000),np.arange(15000,70000,5000)))


x=np.arange(-20144.12093154148, -14393.984333674438,100)#np.hstac


x=np.arange(-30000, -14393.984333674438,100)#np.hstac

x=np.arange(-25000, -14000.984333674438,100)#np.hstac


x=np.arange(14393.984333674438,30000,100)#np.hstac
y=0*x
# save transect
m=np.zeros((len(x),4))
m[:,0]=np.arange(1,len(x)+1)
m[:,1]=x
m[:,2]=y
np.savetxt('transect.bp',m,fmt='%d %f %f %f')



s.trans=bp_transect(s, s.nc, x, y, latlon=False)

s.init_node_tree(latlon=False)
#dd,nn=s.node_tree_xy.query(list(zip(x,y)))
#plt.clf()
#s.plotAtnodes(s.depths,latlon=False)
#s.plot_mesh(latlon=False)
#plt.plot(np.asarray(s.x)[nn],np.asarray(s.y)[nn],'w.')

#d=np.asarray(s.depths)[nn]
#plt.figure()
#plt.clf()
#plt.plot(s.trans.l[:,0],-d,'.-')

fixlim=True

varname='temp'
varname='zcor'
varname='salt'
for	ti in range(0,24*24,6):
	print(ti)

	Z=s.trans.ds['zcor'][ti,:].values
	V=s.trans.ds[varname][ti,:].values
	time=str(s.trans.ds['time'][ti].values)[:19]
	VF=s.nc[varname][ti,:,-1].values

	plt.clf()
	plt.subplot(2,1,1)
	s.plotAtnodes(VF,latlon=False)
	#s.plot_mesh(latlon=False)
	#s.plot_mesh(latlon=False)
	#plt.colorbar(latlon=False)
	L=x.max()-x.min()
	if x.max()>0:
		plt.xlim((x.min()*1.1,x.max()*1))
	else:	
		plt.xlim((x.min()*0.9,x.max()*0.9))
	plt.ylim((y.mean()-L/4,y.mean()+L/4))
	plt.plot(x,y,'.-',color='w',linewidth=1,markersize=1)
	#plt.axis('equal')
	plt.title(time)
	
	if ti==0:
		vmin0, vmax0 = plt.gci().get_clim()
		
	if fixlim==True:
		vmin0, vmax0 = 18,25
		vmin0, vmax0 = 15,24
	
	plt.clim((vmin0,vmax0))
	plt.subplot(2,1,2)
	#ph=plt.pcolor(s.trans.l/1000,Z,V,shading='flat',cmap=plt.cm.jet,edgecolors='w',linewidth=0.4)
	ph=plt.pcolor(s.trans.l/1000,Z,V,shading='flat',cmap=plt.cm.jet)
	#plt.ylim((-150,0))
	plt.ylim((-30,0))
	#plt.ylim((-150,150))
	ch=plt.colorbar()
	#ch.set_label('salinity [-]')
	ch.set_label(varname)
	#plt.pause(0.01)
	if ti==0:
		vmin, vmax = plt.gci().get_clim()
	plt.clim((vmin,vmax))
	plt.savefig('{:04d}_{:s}.png'.format(ti,varname),dpi=130)

###

varname='zcor'
varname='hvel'
for	ti in range(0,24*24,6):
	print(ti)

	Z=s.trans.ds['zcor'][ti,:].values
	V=s.trans.ds[varname][ti,:].values[:,:,0]
	time=str(s.trans.ds['time'][ti].values)[:19]
	VF=s.nc[varname][ti,:,-1,0].values

	plt.clf()
	plt.subplot(2,1,1)
	s.plotAtnodes(VF,latlon=False)
	s.plot_mesh(latlon=False)
	#plt.colorbar(latlon=False)
	L=x.max()-x.min()
	plt.xlim((x.min()*1,x.max()*1))
	plt.ylim((y.mean()-L/4,y.mean()+L/4))
	plt.plot(x,y,'.-',color='w',linewidth=1)
	#plt.axis('equal')
	plt.title(time)
	
	if ti==0:
		vmin0, vmax0 = plt.gci().get_clim()
	plt.clim((-0.3,0.3))
	plt.subplot(2,1,2)
	#ph=plt.pcolor(s.trans.l,Z,V,shading='flat',cmap=plt.cm.jet,edgecolors='w',linewidth=0.4)
	ph=plt.pcolor(s.trans.l,Z,V,shading='flat',cmap=plt.cm.jet)
	#plt.ylim((-150,0))
	plt.ylim((-100,0))
	ch=plt.colorbar()
	ch.set_label('Eastern Velocity [m/s]')
	#plt.pause(0.01)
	if ti==0:
		vmin, vmax = plt.gci().get_clim()
	plt.clim((vmin,vmax))
	plt.clim((-0.3,0.3))
	plt.savefig('{:04d}_{:s}.png'.format(ti,varname),dpi=130)