
import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')

from schism import *

s=schism_setup()
bdnodes=np.hstack(s.bdy_segments[:-1])-1   # zeros baesed numbering
wwmbnd=np.zeros(s.nnodes)
wwmbnd[bdnodes]=2

#s.plotAtnodes(wwmbnd)
s.depths=wwmbnd
s.dump_hgridll('wwmbnd.gr3')
