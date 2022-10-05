"""
load station output variables base don station input
"""

import numpy as np
import sys
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
from schism import *
import datetime as dt
class schism_station_output:
	def __init__(self,fname='station.in',output='outputs',comments='!'):
		with open('station.in') as f:
			#line=f.readline().replace('\n','')
			lines=f.readlines()
			active_flag,var_names=lines[0].replace('\n','').split('for')
			active_flag=np.asarray([np.int(val) for val in active_flag.split('!')[0].split(' ')[:9]])
			var_names=[name.replace(' ','') for name in var_names.split(',')]
			nstations=np.int(lines[1].split('!')[0])
			
			line=lines[2]
			num,txt=line.split('!')[:2]
			self.coords=np.asarray([np.float(i) for i in num.split(' ')[1:4]])
			self.stations=[txt.split(' ')[1]]
			for line in lines[3:]:
				num,txt=line.replace('\n','').split('!')
				self.coords=np.vstack((self.coords,np.asarray([np.float(i) for i in num.split(' ')[1:4]])))
				self.stations.append(txt[1:])
			
			# coords to lon lat based on nearset neighbour
			if (self.coords[:,1]>90).sum() > 0:
				s=schism_setup()
				s.init_node_tree(latlon=False)
				nn=s.node_tree_xy.query(list(zip(self.coords[:,0],self.coords[:,1])))[1]
				lon,lat=np.asarray(s.lon),np.asarray(s.lat)	
				self.coords[:,0]=lon[nn]
				self.coords[:,1]=lat[nn]
				
			self.station_out={name:0 for name in var_names}
			for nr in np.where(active_flag)[0]:
				values=np.loadtxt(output+'/staout_{:d}'.format(nr+1))
				self.station_out[var_names[nr]]=values[:,1:]
			station_out['time']=values[:,0]
			
			with open('param.nml') as f:
				for line in f.readlines():
					if 'start_year' in line:
						year=np.int(line.split('=')[1].split('!')[0])
					elif 'start_month' in line:
						month=np.int(line.split('=')[1].split('!')[0])
					elif 'start_day' in line:
						day=np.int(line.split('=')[1].split('!')[0])
					elif 'start_hour' in line:
						hour=np.float(line.split('=')[1].split('!')[0])
						minute=np.int((hour-np.floor(hour))*60)
						hour=np.int(np.floor(hour))
						break
			refdate=dt.datetime(year,month,day,hour,minute)
			self.time=refdate+dt.timedelta(seconds=station_out['time'][0])*np.arange(1,len(station_out['time'])+1)
			
	def plot_station(self,istat=0):		
		if type(istat)==int:
			i=0				
		else: # get from name
			for i,name in enumerate(stations):
				if istat.lower() in name.lower():
					break
			
		plt.plot(self.time,
		self.station_out['elev'][:,i])
		plt.title(self.stations[i])
		plt.grid()
		plt.show()
		plt.gcf().autofmt_xdate()
		plt.tight_layout()

