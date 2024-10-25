# crate homogenous and depth threshold based  gr3 files
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *

# .gr3 name without filetype: vale  - spatially homogenous
#homogenous={'albedo':0.0,'diffmax':1.0,'diffmin':1.0000000e-06,'watertype':7.0000000e+00}

homogenous={'diffmax':1.0,'diffmin':1.0000000e-06,'watertype':7.0000000e+00,'rough':5e-4} #'shapiro':5e-4}

# depth dpendend parameterisation
#depth_dependent={'drag':{2:[5e-3,1.5e-3]}} # gr3 name:{Threshold:[ValueWhereShallower,ValueWhereDeeper]}

depth_dependent={}
s=schism_setup()
depths=np.asarray(s.depths).copy()	

values=np.zeros(s.nnodes)
for key in homogenous.keys():
	values[:]=homogenous[key]
	s.depths=values[:]
	s.dump_hgridgr3(filename=key+'.gr3',comment=key+'spat. hom. ='+str(homogenous[key]))
	
for key in depth_dependent.keys():
	threshold=list(depth_dependent[key].keys())[0]
	values[depths<=threshold]=depth_dependent[key][threshold][0]
	values[depths>threshold]=depth_dependent[key][threshold][1]
	s.depths=values
	s.dump_hgridgr3(filename=key+'.gr3', comment=key+ 'depth dependend shaollower/deeper {:f} use {:f}/{:f}'.format(threshold,depth_dependent[key][threshold][0],depth_dependent[key][threshold][1]))

s.dump_tvd_prop()

s.lon=np.ones(len(s.x))*34
s.lat=np.ones(len(s.y))*43.3

s.dump_hgridll('hgrid.ll')

