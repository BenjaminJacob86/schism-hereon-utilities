


class dwd_stations(object):
	"""
	read station data files from https://opendata.dwd.de/climate_environment/CDC/observations_germany/
	"""
	import pandas as pd
	import numpy as np
	from glob import glob
	
	def __init__(self,dwddir='/work/gg0028/g260114/DATA/dwd/obs'):
		self.dwddir=dwddir
		with open(self.dwddir+'/'+'stations_dwd.txt') as f:
			header=f.readline().split()
			nheader=len(header)
			d={h:[] for h in header}

			line=f.readline()

			for line in f.readlines():
				vals=line.split()
				if len(vals) > nheader: # space in names, merge
					vals[-3]=' '.join((vals[-3],vals[-2]))
					vals.remove(vals[-2])
				#else:
				for key,val in zip(header,vals):
					d[key]+=[val,]
			
			#convert strings
			for i,key in enumerate(header):
				if i==0:
					d[key]=np.asarray(d[key],int)
				elif i <3:
					d[key]=np.asarray(pd.to_datetime(d[key]),np.datetime64)
				elif i < 6:			
					d[key]=np.asarray(d[key],float)
			
			self.stations=pd.DataFrame.from_dict(d)
			self.stationdata={}

	def get_id_from_name(self,name='Helgoland'):
		return self.stations.loc[self.stations.Stationsname==name].Stations_id

	def get_data_for_id(self,id):
		self.stationdata[id]={}
		sel=self.stations.loc[self.stations.Stations_id==id]
		self.stationdata[id]['name']=sel.Stationsname.values[0]
		self.stationdata[id]['lat']=sel.geoBreite.values[0]
		self.stationdata[id]['lon']=sel.geoLaenge.values[0]

		files=glob(self.dwddir+'/*{:d}*.txt'.format(id))
		for file in files:
			print('loading ' +file)
			varname=file.split('_')[3]
			dsadd=pd.read_csv(file,delimiter=';',parse_dates=['MESS_DATUM'])
			for key in dsadd.keys(): #mask dummy -999
				dsadd[key]=dsadd[key].mask(dsadd[key]==-999)
			if varname in self.stationdata[id].keys():
				pd.merge(self.stationdata[id][varname],dsadd)
			else:
				self.stationdata[id][varname]=dsadd#pd.read_csv(file,delimiter=';',parse_dates=['MESS_DATUM'])
		return self.stationdata
