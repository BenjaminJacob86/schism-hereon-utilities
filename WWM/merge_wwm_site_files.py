"""
mv wwms *.site files from within different folders
"""

from glob import glob
import numpy as np

sitenames=[name.split('/')[1] for name in glob('wwm_out01/*.site')]
sites=dict.fromkeys(sitenames)
for name in sitenames:
	files=np.sort([ os.getcwd() + '/' + file for file in glob('wwm_out*/'+name)	 ])
	sites[name]=load_wwm_site_multiple_files(files,write_joint_file=True)

def load_wwm_site_multiple_files(files,write_joint_file=False):
	""" line by line reading only first to columns to circumvent issues """
	t=[]
	hs=[]
	tm01=[]
	tm02=[]
	
	for file in files:
		with open(file) as f:
			line=f.readline()
			header=line.split()[:-1][:2]
			data=dict.fromkeys(header)
			for line in f.readlines():
				ti,hsi,tm1i,tm2i=line.split()[:4]
				t.append(ti)
				hs.append(hsi)
				tm01.append(tm1i)
				tm02.append(tm2i)

    #uniqute time steps (remove double times due to restarts)						
	unqindex=np.unique(t,return_index=True)[1]

	data['TIME']=np.asarray(t,float)[unqindex]	
	data['HS']=np.asarray(hs,float)[unqindex]
	data['tm01']=np.asarray(tm01,float)[unqindex]
	data['tm02']=np.asarray(tm02,float)[unqindex]
	min=(data['TIME'])*100
	min=np.asarray(np.fix(100*(min-np.fix(min))),int)
	hour=(data['TIME'])
	hour=np.asarray(np.round(100*(hour-np.fix(hour))),int)
	day=(data['TIME']/100)
	day=np.asarray(np.fix(100*(day-np.fix(day))),int)
	mon=(data['TIME']/10000)
	mon=np.asarray(np.fix(100*(mon-np.fix(mon))),int)
	year=np.asarray(np.fix(data['TIME']/10000),int)
	data['dates']=np.asarray(dt.datetime(year[0],mon[0],day[0],hour[0],0,0)+np.arange(len(year))*dt.timedelta(hours=(np.diff(min[:2])/60+np.diff(hour[:2]))[0]),np.datetime64)
	
	if write_joint_file:
		outname=files[0].split('/')[-1].replace('.','_joint.')		
		header='\t TIME	\t HS \t TM01 \t TM02'
		m=np.vstack((data['TIME'],data['HS'],data['tm01'],data['tm02'])).T
		np.savetxt(outname,m,fmt='%.6f \t %.6f \t %.6f \t %.6f \t',header=header)
	return data	
