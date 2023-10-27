import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
plt.ion()
s=schism_setup()
s.read_gr3('TEM_nudge.gr3')
include=s.gr3['TEM_nudge']>0
field=np.zeros(s.nnodes)


ds=xr.open_dataset('TEM_nu.nc')
field[include]=ds.tracer_concentration[0,:,-1,0].values
include.sum()