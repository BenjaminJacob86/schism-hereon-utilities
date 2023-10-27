# make LockExchange .py
import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *

# import path fore in region check
from matplotlib.path import Path



s=schism_setup()



				
				
				
#write temperature
s.dump_gr3('temp.ic',const=4,comment='homogenousTempForLockExchange')

# initialize Salinity for one half
S=np.ones(s.nnodes)*18

# fille other part by define region
s.read_reg('LockExchangeLeft.reg')
reg=s.reg['LockExchangeLeft']

Poly=Path(list(zip(reg[:,0],reg[:,1])))
coords=list(zip(s.x,s.y))
inarea_nodes=np.asarray(np.where(Poly.contains_points(coords)))
S[inarea_nodes]=36
s.dump_gr3_spat_var('salt.ic',S,comment='18/36 LockExchange')

# switch off wind
s.dump_gr3('windfactor.gr3',const=0,comment='ZeroWInd')
