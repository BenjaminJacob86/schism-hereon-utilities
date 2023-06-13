
import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
s=schism_setup()
plt.ion()


ds=xr.open_dataset('ww3.2090_incJade_spec.nc') # open wave süectra file

#names=[''.join(ds.station_name[i,:].values) for i in range(int(ds.station[-1].values))]

lon=ds.longitude.values[0,:]
lat=ds.latitude.values[0,:]


#import matplotlib
#from matplotlib.collections import PolyCollection as PolyCollection
from ww3 import WW3_mesh
ww3=WW3_mesh('NBSext_bl.msh')

ww3.plot_mesh()

plt.clf()
ww3.plotAtnodes(ww3.d)
s.plot_domain_boundaries(append=True)
plt.plot(lon,lat,'k.')
plt.axis((7.861384840637267, 8.66835755783823, 53.27858432360115, 53.897847246893726))
plt.savefig('WaveForcingLocations.png',dpi=300)