"""
Set WW3 start and end dates
based ond dates assumed as start and end in 
ww3_shel.nml, update with start and enddates in DWD windfile name (ICONglob*)
and call replace script.
"""

from glob import glob
import datetime as dt


# load current dates from nml
with open('ww3_shel.nml') as f:
	for line in f.readlines():
		if 'DOMAIN%START' in line and '!' not in line:
			startdate0=line[line.index("'")+1:line.index("'")+9]	
		if 'DOMAIN%STOP' in line and '!' not in line:
			enddate0=line[line.index("'")+1:line.index("'")+9]	
			break

# load desired dates from windfile name
windfile=glob('ICONglob*')[0]		
str1=windfile[windfile.index('_')+1:]
startdate1=str1[:str1.index('_')]
enddate1=str1[str1.index('_')+1:str1.index('.')]

#update dates in replace script
script='replaceDate.sh'
backup='./'+script+'_'+str(dt.date.today())
import shutil
shutil.copy('./'+script, './'+script+'_'+str(dt.date.today())) #backup

with open(backup) as fin:
	with open(script,'w') as fout:
		for line in fin.readlines():
			print(line)
			if 'startdate=' in line:
				fout.write('startdate='+startdate0+'\n')
			elif 'startdatenew=' in line:
				fout.write('startdatenew='+startdate1+'\n')
			elif 'enddate=' in line:
				fout.write('enddate='+enddate0+'\n')
			elif 'enddatenew=' in line:
				fout.write('enddatenew='+enddate1+'\n')
			else:
				fout.write(line)
import subprocess
subprocess.check_output('bash replaceDate.sh', shell=True)# call replace script