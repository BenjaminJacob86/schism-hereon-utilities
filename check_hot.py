import os

plt.ion()


os.chdir('/gpfs/work/kochw/schism-routine-BS/RUN24d/')
s0=schism_setup()
ds0=xr.open_dataset('/gpfs/work/kochw/schism-routine-BS/RUN24d/hotstart.nc0')

zi=-1
plt.figure()
plt.subplot(1,2,1)
a=ds['tr_nd'][:,zi,0].values
s.plotAtnodes(a)
ax=plt.axis()
b=ds0['tr_nd'][:,zi,0].values
plt.subplot(1,2,2)
s0.plotAtnodes(b)
plt.axis(ax)

plt.figure()
plt.subplot(1,2,1)
a=ds['tr_nd'][:,0,0].values
s.plotAtnodes(s.depths)
ax=plt.axis()
b=ds0['tr_nd'][:,0,0].values
plt.subplot(1,2,2)
s0.plotAtnodes(s0.depths)
plt.axis(ax)
plt.clim(-40,120)
s.plot_domain_boundaries(append=True)