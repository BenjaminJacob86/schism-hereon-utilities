# cut Black Sea
# cut out grid of setup within domain of second setup
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *
from matplotlib import path
import os

# reference setup to take grid from
os.chdir('/gpfs/work/jacobb/data/RUNS/BS-routine') #/gpfs/work/kochw/schism-routine-BS/RUN24d/')
s=schism_setup()

trim_setup=False
if trim_setup:
	# setup with ne domain to reduce to
	os.chdir('/gpfs/work/jacobb/data/SETUPS/NWBlackSea/setup3')
	s1=schism_setup() # sammalle setup to cut iff
	bds=[]
	for i in range(len(s1.bdy_segments)):
		bds+=s1.bdy_segments[i]+s1.land_segments[i][1:]
	bds=np.asarray(bds)
	xpoly=np.asarray(s1.lon)[bds-1]
	ypoly=np.asarray(s1.lat)[bds-1]
else: # trim to regile
	regfile='basinSelect.reg' # extract grid within region to stand alone grid
	m=np.loadtxt(regfile,skiprows=3)
	xpoly,ypoly=m[:,0],m[:,1]



#os.chdir('/gpfs/work/jacobb/data/SETUPS/NWBlackSea/bscut/')



#s1.plot_domain_boundaries()
#nvrt=np.asarray([s.vgrid[key].count()    for key in s.vgrid.keys()])
#plt.figure()
#s.plotAtnodes(nvrt)


p=path.Path(list(zip(xpoly,ypoly)))

x=np.asarray(s.x)
y=np.asarray(s.y)
x2=np.asarray(s.lon)
y2=np.asarray(s.lat)

if ypoly.max()> 90:
	cx=x[s.nvplt].mean(axis=1)
	cy=y[s.nvplt].mean(axis=1)
else:	#latlon
	cx=x2[s.nvplt].mean(axis=1)
	cy=y2[s.nvplt].mean(axis=1)
coords=list(zip(cx,cy))
#
#ikeep=cx<31.5 # keep elements with center west of easter boundary
ikeep=p.contains_points(coords) # keep elements within polygon
triscut=s.nvplt[ikeep,:].copy()

uniqueNr=np.unique(triscut)
newNoderNr=np.asarray([i for i in range(len(uniqueNr))])
switcher=dict(zip(uniqueNr, zip(newNoderNr))) # relate new number to old numbers
vgrid2={i+1:[] for i in range(1,len(uniqueNr))}

d0=np.asarray(s.depths)
d=np.zeros(len(uniqueNr))

xout=np.zeros(len(uniqueNr)) # lon 
yout=np.zeros(len(uniqueNr))  #  lat

xoutproj=np.zeros(len(uniqueNr)) # lon 
youtproj=np.zeros(len(uniqueNr)) #   lat


# subsest elementable, depth and vgrid 
for oldIndex in uniqueNr:
	newIndex=switcher[oldIndex][0]
	triscut[triscut==oldIndex]=newIndex
	d[newIndex]=d0[oldIndex]
	xout[newIndex]=x2[oldIndex]
	yout[newIndex]=y2[oldIndex]
	vgrid2[newIndex+1]=s.vgrid[oldIndex+1]	
	xoutproj[newIndex]=x[oldIndex]
	youtproj[newIndex]=y[oldIndex]

# sortierung der list
import collections
vgrid3 = collections.OrderedDict(sorted(vgrid2.items()))

#nvrt=np.asarray([s.vgrid[key].count()    for key in s.vgrid.keys()])
#nvrt2=np.asarray([vgrid2[key].count()    for key in vgrid2.keys()])
#nvrt3=np.asarray([vgrid3[key].count()    for key in vgrid3.keys()])
#
#plt.figure()
#plt.subplot(2,1,1)
#s.plotAtnodes(nvrt)
#plt.subplot(2,1,2)
#s2.plotAtnodes(nvrt3)
#
#
#plt.figure()
#plt.subplot(2,1,1)
#s.plotAtnodes(s.depths)
#plt.subplot(2,1,2)
#s2.plotAtnodes(s2.depths)
#
	
	
	
# save reduced grid	
s.x=xoutproj
s.y=youtproj
s.lon=xout
s.lat=yout
s.depths=d
s.vgrid=vgrid3
s.nnodes=len(xout)
s.nelements=triscut.shape[0]
s.nv=triscut+1

	
s.dump_hgridgr3('hgrid_cut.gr3')	
s.dump_hgridll('hgrid_cut.ll')	
s.dump_vgrid('vgrid_cut.in')



# test
s2=schism_setup(hgrid_file='hgrid_cut.gr3',ll_file='hgrid_cut.ll',vgrid_file='vgrid_cut.in')
nvrt2=np.asarray([s2.vgrid[key].count()    for key in s2.vgrid.keys()])
s2.plotAtnodes(nvrt2)

triscut=s2.nvplt.copy()




np.asarray((triscut==inode).sum(axis=1),bool)

for inode in range(s2.nnodes):
	ilem=np.asarray((triscut==inode).sum(axis=1),bool)
	np.where(triscut==inode)
	break
	
# keep only none outspiking elements i.e. remove elements of nodes occuring only in one element	
keep_elems=np.asarray([np.where((triscut==inode).sum(axis=1)>0)[0] for inode in range(s2.nnodes) if len(np.where((triscut==inode).sum(axis=1)>0)[0]) >1] )
keep_elems=np.unique(np.hstack(keep_elems))


plt.figure()
s2.plot_domain_boundaries()
s2.plot_mesh(latlon=True)
plt.triplot(s2.lon,s2.lat,triscut[keep_elems,:],color='k')



