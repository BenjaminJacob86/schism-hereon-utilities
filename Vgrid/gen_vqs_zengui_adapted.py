#!/usr/bin/env
#from pylib import *
"""
Adapted from
https://github.com/wzhengui/pylibs/tree/master/scripts/gen_vqs.py
for my workflow
"""


close('all')

import numpy as np
from matplotlib import pyplot as plt

s=schism_setup()
s.depths=np.asarray(s.depths)
#------------------------------------------------------------------------------
#input
#------------------------------------------------------------------------------
grd='./hgrid.gr3'        #hgrid name (grid.npz, or hgrid.gr3)

#hsm: depth for each master grid;  nv: number for layers for each depth
#hsm=np.array([0.5,1.0,1.5, *np.arange(2,11),*[11.5,14,17,21],*np.arange(25,51,5),*np.arange(60,81,10),95,125])
#hsm=np.array([0.5,1.0,1.5, *np.arange(2,11),*[11.5,14,17,21],*np.arange(25,51,5),*np.arange(60,81,10),95,125])


# hack  -to reuse forcing set open boundary nodes to amx depths


##


hsm=np.array([0.5,1.0,1.5, *np.arange(2,11),*[11.5,14,17,21],*np.arange(25,51,5)]) #,*np.arange(60,81,10),95,125])


nhm=len(hsm); nv=2+np.arange(nhm)

# restict
nv=np.minimum(21,nv) # to keep aligned with german bight

# hack for german bight whcih was sigma  keep the maximum layes at the open bd to the default of 21 to reuse bd forcing
bdind=np.asarray(s.bdy_segments[0])-1
s.depths[bdind]=hsm[-1]

#check transect
bname='./transect.bp'   #transect bpfile
#------------------------------------------------------------------------------
#compute master grid
#------------------------------------------------------------------------------
nvrt=nv[-1]; z_mas=np.ones([nhm,nvrt])*np.nan; eta=0.0
for m, [hsmi,nvi] in enumerate(zip(hsm,nv)):
    #strethcing funciton
    theta_b=0; theta_f=2.5; hc=10.5
    hc=np.minimum(hsmi,hc)

    for k in np.arange(nvi):
        sigma= k/(1-nvi)  #zi=-sigma #original sigma coordiante

        #compute zcoordinate
        cs=(1-theta_b)*np.sinh(theta_f*sigma)/np.sinh(theta_f)+theta_b*(np.tanh(theta_f*(sigma+0.5))-np.tanh(theta_f*0.5))/2/np.tanh(theta_f*0.5)
        z_mas[m,k]=eta*(1+sigma)+hc*sigma+(hsmi-hc)*cs

    #normalize z_mas
    z_mas[m]=-(z_mas[m]-z_mas[m,0])*hsmi/(z_mas[m,nvi-1]-z_mas[m,0])
s_mas=np.array([z_mas[i]/hsm[i] for i in np.arange(nhm)])

#check master grid
for i in np.arange(nhm-1):
    if min(z_mas[i,:nv[i]]-z_mas[i+1,:nv[i]])<0: print('check: master grid layer={}, hsm={}, nv={}'.format(i+1,hsm[i+1],nv[i+1]))

#plot master grid
plt.figure(figsize=[10,5])
for i in np.arange(nhm): plt.plot(i*np.ones(nvrt),z_mas[i],'k-',lw=0.3)
for k in np.arange(nvrt): plt.plot(np.arange(nhm),z_mas.T[k],'k-',lw=0.3)

ax = plt.gca()
ax.set_xlim([-0.5, nhm-0.5])
ax.set_ylim([-hsm[-1], 0.5])

plt.gcf().tight_layout()
#move_figure(plt.gcf(), 0, 0)
plt.savefig('Master_Grid', dpi=200)

# sys.exit()

#------------------------------------------------------------------------------
#compute vgrid
#------------------------------------------------------------------------------
#read hgrid
#gd=loadz(grd).hgrid if grd.endswith('.npz') else read_schism_hgrid(grd)
# using my schism_setup here
fpz=s.depths<hsm[0]; s.depths[fpz]=hsm[0]

#find hsm index for all points
rat=np.ones(s.nnodes)*np.nan; nlayer=np.zeros(s.nnodes).astype('int');
ind1=np.zeros(s.nnodes).astype('int'); ind2=np.zeros(s.nnodes).astype('int')
for m, hsmi in enumerate(hsm):
    if m==0:
        fp=s.depths<=hsm[m];
        ind1[fp]=0; ind2[fp]=0; rat[fp]=0; nlayer[fp]=nv[0]
    else:
        fp=(s.depths>hsm[m-1])*(s.depths<=hsm[m])
        ind1[fp]=m-1; ind2[fp]=m
        rat[fp]=(s.depths[fp]-hsm[m-1])/(hsm[m]-hsm[m-1]); nlayer[fp]=nv[m]

znd=z_mas[ind1]*(1-rat[:,None])+z_mas[ind2]*rat[:,None]; #z coordinate
for i in np.arange(s.nnodes): znd[i,nlayer[i]-1]=-s.depths[i]
snd=znd/s.depths[:,None]; #sigma coordinate

#check vgrid
for i in np.arange(s.nnodes):
    for k in np.arange(nvrt-1):
        if znd[i,k]<=znd[i,k+1]:
            sys.exit('wrong vertical layers')

#write vgrid.in
fid=open('vgrid.in','w+')
fid.write('1  !average # of layers={:0.2f}\n{}\n'.format(np.mean(nlayer),nvrt))
for i in np.arange(s.nnodes):
    nlayeri=nlayer[i]; si=np.flipud(snd[i,:nlayeri])
    fstr='{:6}    {:2}   '+'{:10.6f}  '*nlayeri+'\n'
    fid.write(fstr.format(i+1,nvrt-nlayeri+1,*si))
fid.close()
print('Average number of layers is: {:0.2f}'.format(np.mean(nlayer)))

#------------------------------------------------------------------------------
#plot transect
#------------------------------------------------------------------------------
if os.path.exists(str(bname)):
    bp=read_schism_bpfile(str(bname))

    #compute dist
    dist=[0,]
    for i in arange(bp.nsta-1):
        disti=abs((bp.x[i+1]-bp.x[i])+1j*(bp.y[i+1]-bp.y[i]))+dist[i]
        dist.append(disti)
    dist=array(dist)

    #compute zcor
    sindp=near_pts(c_[bp.x,bp.y],c_[gd.x,gd.y]); zi=znd[sindp]
    for i in arange(bp.nsta): fpn=isnan(zi[i]); zi[i][fpn]=min(zi[i])


    #plot
    figure(figsize=[10,5])
    for k in arange(nvrt): plot(dist,zi[:,k],'k',lw=0.5)
    for i in arange(bp.nsta): plot(ones(nvrt)*dist[i],zi[i],'k',lw=0.5)
    setp(gca(),ylim=[zi.min()-1,0.5],xlim=[0,dist.max()])
    gcf().tight_layout()
    # move_figure(gcf(),0,0)