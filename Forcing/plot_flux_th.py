"""
	Plot flux.th rivers and label based on comments
	within bctides.in
"""
import numpy as np
import datetime as dt

with open('param.nml') as f:
	for line in f.readlines():
		line=line.split('!')[0]
		if 'start_year' in line:
			year=np.int(line.split('=')[1])
		if 'start_month' in line:
			month=np.int(line.split('=')[1])
		if 'start_day' in line:
			day=np.int(line.split('=')[1])
		if 'start_hour' in line:
			hour=np.float(line.split('=')[1])
		if 'utc_start' in line:
			utc=np.float(line.split('=')[1])
			break
			
minute=np.int((hour - np.floor(hour))*60)
hour=np.int(hour)   
refdate=dt.datetime(year,month,day,hour,minute)   
	

with open('bctides.in') as f:
	lines=f.readlines()
	nlines=len(lines)
	names=[]
	for i in range(5,nlines):
		line=lines[i]
		if '!' in line:
			vars,comment=line.split('!')
			vars=[np.float(v) for v in (vars.split(' ')) if len(v)>0]
			if len(vars)>3:
				if vars[2]==1:
					names.append(comment.split('\n')[0])

Q=np.loadtxt('flux.th')					
t=Q[:,0]
Q=Q[:,1:]