"""
gen station out based on CMEMS tide gauge and marnet data
within tolerance distance to model grid nodes
"""


__copyright__   = "Copyright 2021 - Helmholtz-Zentrum hereon GmbH"
__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hereon.de"
__status__ = "Development"

import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import *
from netCDF4 import Dataset

########### Settings ####################################################################
rundir='/gpfs/work/jacobb/data/SETUPS/NorthSea/NorthSea200mHelgoland/RERUN2/'
tgdir='/gpfs/work/ksddata/observation/insitu/CMEMS/NorthWestShelf/TG/201801/'
mardir='/gpfs/work/ksddata/observation/insitu/Marnet_Sebastian/' 					  # myocean
dtol=0.05      # distance tolerance in degree lon/lat
year=2017      # for which year (Marnet data caries year in name)
f=open('station.in','w')
f.write('1 0 0 0 1 1 1 1 1 !on (1)|off(0) flags for elev, air pressure, windx, windy, T, S, u, v, w\n')

cartesian=False  # convert lon lat to cartesian coordinates (nearaset neighbour) for nws = 2
##########################################################################################

######### identify Gauge Data ##########################################
os.chdir(rundir)
s=schism_setup()
s.init_node_tree()
files=glob(tgdir+'*.nc')
TGcoords,TGnames=[],[]
for file in files:
	#if 'Newhaven' in file:
	#	break
	nc=Dataset(file)
	if 'SLEV' in nc.variables.keys():
		try:
			coord=((nc['LONGITUDE'][0],nc['LATITUDE'][0]))
		except:
			coord=((float(nc.geospatial_lon_min),float(nc.geospatial_lat_min)))
		if s.node_tree_latlon.query(coord)[0]  <= dtol:
			if cartesian:
				nn=s.find_nearest_node(coord[0],coord[1])-1
				TGcoords.append((s.x[nn],s.y[nn])+(0,))
				#coord
				#s.lon[nn],s.lat[nn]
			else:
				TGcoords.append(coord+(0,))
			TGnames.append(file[file.rfind('/')+1:].split('_')[-2])
	nc.close()
#for name in TGnames:
#	if 'Newhaven' in name:
#		break
	
### marnet
marnet=glob(mardir+'*_{:d}.nc'.format(year))	
Mcoords=[]
Mnames=[]
nns=[]
for file in marnet:
	nc=Dataset(file)
	coord=(nc['ppos'][1],nc['ppos'][0])
	if s.node_tree_latlon.query(coord)[0]  <= dtol:
		# get zlevels of WT  WS
		depths=[]
		for key in nc.variables.keys():
			if ('WT' in key) | ('WS' in key):
				d=np.float(key.split('m')[0][2:])
				depths.append(d)
				if cartesian:
					nn=s.find_nearest_node(coord[0],coord[1])-1 # minus vergessen vorher
					#nn=s.find_nearest_node(coord[0],coord[1]) #  faöse minus vergessen vorher
					Mcoords.append((s.x[nn],s.y[nn])+(d,))
					nns.append(nn)
				else:
					Mcoords.append(coord+(d,))
					nn=s.find_nearest_node(coord[0],coord[1])-1 # minus vergessen vorher
					nns.append(nn)
				Mnames.append(file[file.rfind('/')+1:].split('_')[-2]+ ' z'+str(d))
		nc.close()
#################	

lon=np.asarray(s.lon)
lat=np.asarray(s.lat)
#plt.clf()
#s.plot_domain_boundaries()
#plt.plot(np.asarray(Mcoords)[:,0],np.asarray(Mcoords)[:,1],'ko')
#plt.plot(lon[nns],lat[nns],'r+')

f.write('{:d} !# of stations\n'.format(len(TGcoords)+len(Mcoords)))
f.write('1 {:f} {:f} {:f} ! {:s}  Format: station #,x,y,z; if ics=2, x,y are degr in lon/lat. z is z-coord (not distance from surface!). For 3D variables, code will extrapolate above surface/below bottom if necessary\n'.format(TGcoords[0][0],TGcoords[0][1],TGcoords[0][2],TGnames[0]))

for i, packed in enumerate(zip(TGcoords[1:]+Mcoords,TGnames[1:]+Mnames)):
	coord, name = packed
	f.write('{:d} {:f} {:f} {:f} ! {:s}\n'.format(i+2,coord[0],coord[1],coord[2],name))
f.close()

# plot
plt.ion()
if cartesian:
	s.lon=s.x
	s.lat=s.y
s.plot_domain_boundaries(append=True)
for i, packed in enumerate(zip(TGcoords,TGnames)):
	coord, name = packed
	ph1=plt.plot(coord[0],coord[1],'ko')
	plt.text(coord[0],coord[1],name[:5])
for i, packed in enumerate(zip(Mcoords,Mnames)):
	coord, name = packed
	ph2=plt.plot(coord[0],coord[1],'ro')
	plt.text(coord[0],coord[1],name[:5])
plt.legend((ph1,ph2),['TG','Marnet'])
plt.savefig('station_in_map.png',dpi=300)	
plt.legend([ph1[0],ph2[0]],['TG','Marnet'])