"""
Create Points list for WW3 to write
SCHISM wave boundary forcing for WWM.
And Append Wave bouy stations also.
"""

import os
import netCDF4
import sys
import matplotlib
import datetime as dt
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')
from schism import * # import schism functions

######### SETTINGS ###########################

#bdlistname='GB_bd_points.list'
bdlistname='SNS_bd_points.list'
rundir='/gpfs/work/jacobb/data/SETUPS/GermanBight/GB2018/'  # schism rundir
schism_waveforce_bds=[0,1]    # list of schism open boundaries to appliy wave forcing at [0,1,3,...]
use_each_bd_node=1          # sub stepping for bd nodes to be used 

# lon lat  #name
stations={'Elbe':[8.1105,53.9970 ],
'Helgoland':[7.8190,54.2190], 
'Westerland':[8.2210,54.9097 ],
'Fino-1':[6.5838,54.0143 ],
'Fino-3':[7.1500,55.2000 ]
}
prefix='SNS'
############################################################################

# I assume interpolations to points is done
#    Lon			'LAt     # Noder Nr Schism
# -4.698938      48.573900  'N192146'

pwd=os.getcwd() 
os.chdir(rundir)
s=schism_setup()
os.chdir(pwd)

lon,lat=np.asarray(s.lon),np.asarray(s.lat)
s.plot_domain_boundaries(append=True,nr=True)
lines=''
for obd in schism_waveforce_bds:
	bdnodes=np.asarray(s.bdy_segments[obd])[::use_each_bd_node]
	ph=plt.plot(lon[bdnodes-1],lat[bdnodes-1],'k.')
	plt.legend(ph,['WaveForcePoints'])
	for node in bdnodes:
		lines+="{:.6f} {:.6f} '{:s}{:d}'\n".format(lon[node-1],lat[node-1],prefix,node)	
		
plt.savefig('WaveForcingBoundaries.png',dpi=200)
plt.close()

for key in stations.keys():
	lines+="{:.6f} {:.6f} '{:s}'\n".format(stations[key][0],stations[key][1],key)

with open(bdlistname,'w') as f:
	f.writelines(lines)
