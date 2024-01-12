"""
mv wwms *.site files from within different folders
"""
import os
from glob import glob
import numpy as np
import datetime as dt


# 8 digits after point in wwm output but not always
# some outputs onöly ****


def load_wwm_site_multiple_files(files,varlist=False,write_joint_file=False):
	""" line by line reading only first to columns to circumvent issues """
	
	for file in files:
		with open(file) as f:
			line=f.readline()
			header=line.split()#[:-1][:2]
			if file==files[0]:
				data=dict.fromkeys(header)
				for key in data.keys():
					data[key]=[]
			for line in f.readlines():
			


				# fix stars destroing format
				isstar=np.asarray([char=='*' for char in line])
				ireplace=np.where(~isstar[:-1] & isstar[1:])[0]+1
				for index in ireplace:
					line=line[:index-1]+' '+line[index+1:]

					line=line.replace('*','9')
				#also with digits after comma it as issues not eperating numbers
				items=line.split()
				if len(items) != len(header):
					items2=[]
					for itemi in items:
						if itemi.count('.')>1: # connectect numbers:
							isplit=itemi.index('.')+8 # 8 digets after points
							items2.append(itemi[:isplit+1])
							items2.append(itemi[isplit+1:])
						else:
							items2.append(itemi)
					items=items2
				#if len(line.split()) != len(header):
				#	print(line)
				#	line2=line
				#	from IPython import embed; embed()									
					
					
				#for d,val in zip(list(data.keys()),line.split()):
				for d,val in zip(list(data.keys()),items):
					data[d].append(val.replace('*','9'))

					
    #uniqute time steps (remove double times due to restarts)						
	unqindex=np.unique(data['TIME'],return_index=True)[1]
	try:
		for d in (data.keys()):
			data[d]=np.asarray(data[d],float)[unqindex]
	except:
		print('error assigning data')			
		from IPython import embed; embed()	
	min=(data['TIME'])*100
	min=np.asarray(np.fix(100*(min-np.fix(min))),int)
	hour=(data['TIME'])
	hour=np.asarray(np.round(100*(hour-np.fix(hour))),int)
	day=(data['TIME']/100)
	day=np.asarray(np.fix(100*(day-np.fix(day))),int)
	
	mon=(data['TIME']/10000)
	mon=np.asarray(np.fix(100*(mon-np.fix(mon))),int)
	year=np.asarray(np.fix(data['TIME']/10000),int)
	
	#day=np.asarray(np.fix((100*day-np.fix(100*day))),int)
	
	day=np.asarray(np.fix(data['TIME']-year*10000)-mon*100,int)
	
	data['dates']=np.asarray(dt.datetime(year[0],mon[0],day[0],hour[0],0,0)+np.arange(len(year))*dt.timedelta(hours=(np.diff(min[:2])/60+np.diff(hour[:2]))[0]),np.datetime64)
	#seconds=
	#data['time'])(data['dates']-data['dates'][0])/np.timedelta64(1,'s')
	
	if varlist!=False:
		keys=list(data.keys())
		for key in keys:
			if key not in ['TIME']+varlist+['dates']:
				header.remove(key)
				data.pop(key)
		
		
				
	if write_joint_file:
		outname=files[0].split('/')[-1].replace('.','_joint.')	
		keys=list(data.keys())
		if 'dates' in keys: #header:
			#header.remove('dates')		
			data.pop('dates')		
			keys.remove('dates')
		#header='%\t TIME' ' \t '.join(varlist) # 
		header=' \t'.join(data.keys())
		try:
			m=np.vstack((data[key] for key in keys)).T # data.keys()
		except:
			print('error stacking data')		
			from IPython import embed; embed()				
		fmt=('%.6f \t'*len(data.keys()))[:-2]
		np.savetxt(outname,m,fmt=fmt,header=header)

		#np.savetxt(outname,m,fmt='%.6f \t %.6f \t %.6f \t %.6f \t',header=header)
	return data	

sitenames=[name.split('/')[1] for name in glob('wwm_out01/*.site')]
sites=dict.fromkeys(sitenames)

varlist=['HS','DM','DSPR','TPP','ELEVATION']
#varlist=False
for name in sitenames:
	files=np.sort([ os.getcwd() + '/' + file for file in glob('wwm_out*/'+name)	 ])
	sites[name]=load_wwm_site_multiple_files(files,varlist=varlist,write_joint_file=True)