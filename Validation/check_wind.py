import sys
import os
# check wwind
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities//Lib/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/Lib/')

from schism import *
from glob import glob


#schism_dir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'

schism_dir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/'
nc_dir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs01/'
os.chdir(schism_dir)
s=schism_setup()
s.ds=schism_outputs_by_variable(nc_dir,varlist='out2d').ds
s.init_node_tree()
s.init_node_tree(latlon=False)

reftime=np.asarray(dt.datetime(np.int(p.get_parameter('start_year')),\
	np.int(p.get_parameter('start_month')),\
	np.int(p.get_parameter('start_day')),\
	np.int(p.get_parameter('start_hour')),0,0),np.datetime64)
tschism=reftime+s.ds['out2d'].time.values *np.timedelta64(1,'s')

swind=s.ds.get('windSpeed')
shs=s.ds['out2d']# s.ds.get('sigWaveHeight')


# helgo post
XOUTS=422968.911
YOUTS=6008763.77
nn2=s.node_tree_xy.query((XOUTS,YOUTS))[1]

plt.figure()
plt.subplot(2,1,1)
s.plot_domain_boundaries(latlon=False,append=True)
plt.plot(XOUTS,YOUTS,'ko')
plt.subplot(2,1,2)
s.plot_domain_boundaries(latlon=True,append=True)
plt.plot(data[2115]['lon'],data[2115]['lat'],'ko')


dwd=dwd_stations()		
data=dwd.get_data_for_id(dwd.get_id_from_name('Helgoland').values[0])

t=data[2115]['ff'].MESS_DATUM.values
vabs=data[2115]['ff'].FF_10.values
dwd.stationdata[2115]

#nearest neighbour
nn=s.node_tree_latlon.query((data[2115]['lon'],data[2115]['lat']))[1]
vabs_schism=np.sqrt((swind['windSpeed'][:,:,nn]**2).sum(axis=0).values)


vabs_schism=np.sqrt((swind['windSpeed'][:,:,nn]**2).sum(axis=0).values)
hs_schism=shs['sigWaveHeight'][:,nn].values
hs_schism2=shs['sigWaveHeight'][:,nn2].values


# sflux
dsair=xr.open_mfdataset(np.sort(glob(schism_dir+'sflux/*air*')))
#
nny,nnx=np.unravel_index(np.argmin((dsair['lon'].values-data[2115]['lon'])**2 + (dsair['lat'].values-data[2115]['lat'])**2),dsair['lon'].shape)
dsair=dsair.sel(ny_grid=[nny,],nx_grid=[nnx,]) #, method="nearest"
sflux_abs=np.sqrt(dsair['uwind']**2+dsair['vwind']**2).values.flatten()
sflux_time=dsair['time'].values

#wave

#def load_wwm_site(file):
#	""" line by line reading only first to columns to circumvent issues """
#	with open(file) as f:
#		line=f.readline()
#		header=line.split()[:-1][:2]
#		data=dict.fromkeys(header)
#		t=[]
#		hs=[]
#		tm01=[]
#		tm02=[]
#		#tm10=[]
#		for line in f.readlines():
#			ti,hsi,tm1i,tm2i=line.split()[:4]
#			t.append(ti)
#			hs.append(hsi)
#			tm01.append(tm1i)
#			tm02.append(tm2i)
#		data['TIME']=np.asarray(t,float)	
#		data['HS']=np.asarray(hs,float)
#		data['tm01']=np.asarray(tm01,float)
#		data['tm02']=np.asarray(tm02,float)
#		min=(data['TIME'])*100
#		min=np.asarray(np.fix(100*(min-np.fix(min))),int)
#		hour=(data['TIME'])
#		hour=np.asarray(np.round(100*(hour-np.fix(hour))),int)
#		day=(data['TIME']/100)
#		day=np.asarray(np.fix(100*(day-np.fix(day))),int)
#		mon=(data['TIME']/10000)
#		mon=np.asarray(np.fix(100*(mon-np.fix(mon))),int)
#		year=np.asarray(np.fix(data['TIME']/10000),int)
#		data['dates']=np.asarray(dt.datetime(year[0],mon[0],day[0],hour[0],0,0)+np.arange(len(year))*dt.timedelta(hours=(np.diff(min[:2])/60+np.diff(hour[:2]))[0]),np.datetime64)
#
#	return data




file='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/wwm_out01/HELGON.site'

file='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_REF/wwm_out01/HELGON.site'

hel=load_wwm_site(file)


def load_wwm_site_multiple_files(files,write_joint_file=False):
	""" line by line reading only first to columns to circumvent issues """

	
	for i,file in enumerate(files):
		with open(file) as f:
			line=f.readline()
			
			if i==0:
				header=line.split()
				#data=dict.fromkeys(header)
				data={h:[] for h in header}#dict.fromkeys(header)
				nheader=len(header)
			
			for line in f.readlines():
				line=line.replace('***************',' nan')
				ar=line.split()[:nheader]
				for key,val in zip(header,ar):
					data[key].append(val)


    #uniqute time steps (remove double times due to restarts)						
	unqindex=np.unique(data['TIME'],return_index=True)[1]
	for key in header:
		data[key]=np.asarray(data[key],float)[unqindex]


	if write_joint_file:
		outname=files[0].split('/')[-1].replace('.','_joint.')		
		#header='\t TIME	\t HS \t TM01 \t TM02'
		#m=np.vstack((data['TIME'],data['HS'],data['tm01'],data['tm02'])).T
		m=np.vstack((data[key] for key in data.keys() )).T
		np.savetxt(outname,m,fmt='%.6f \t '*m.shape[1],header=header)

	import pandas as pd
	data['dates']=pd.to_datetime(hel['TIME'], format='%Y%m%d.%H%M%S')	
	
	return data	
	

	
# wind schism
plt.clf()
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,vabs,'.-',label='DWDobs')
plt.plot(sflux_time,sflux_abs,'.-',label='DWDmod')
plt.plot(tschism,vabs_schism,'.-',label='SCHISM')
wwm_abs=np.sqrt(hel['WIND-X']**2+hel['WIND-Y']**2)
plt.plot(hel['dates'],wwm_abs,'.-',label='WWMout')
plt.xlim(tschism[[0,-1]])
plt.ylabel('wind abs [m/s]')
plt.legend()
#plt.gcf().autofmt_xdate()
plt.grid()
plt.subplot(2,1,2)
plt.plot(hel['dates'],hel['HS'],'.-',label='WWMsite')
plt.plot(tschism,hs_schism,'.-',label='SCHISMnn')
plt.plot(tschism,hs_schism2,'r--',label='SCHISMnn2')
plt.ylabel('Hs [m]')
plt.xlim(tschism[[0,-1]])
plt.grid()
plt.legend()
plt.gcf().autofmt_xdate()