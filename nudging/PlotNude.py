import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hereon-utilities/')
from schism import *

plt.ion()
s=schism_setup()


s.read_gr3('include.gr3')
inc=np.asarray(s.gr3['include'],bool)

T,S=np.zeros(s.nnodes),np.zeros(s.nnodes)


dsT=xr.open_dataset('TEM_nu.nc')
dsS=xr.open_dataset('SAL_nu.nc')

SST=dsT.tracer_concentration[0,:,-1,0].values
SSS=dsS.tracer_concentration[0,:,-1,0].values

S[inc]=SSS

ph,ch,ax=s.plotAtnodes(S)