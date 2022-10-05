# calculate speed up of schism run based on latest model output step

import os
import datetime

testfile = "mirror.out"

# read start time
f=open(testfile)
line=f.readline()
ind1=line.index('20') # rund starts 20 something
ind2=line.index(',')
ind3=line.index('\n')

startDate=line[ind1:ind2]
startTime=line[ind2+2:ind3]

year=int(startDate[0:4])
mon=int(startDate[4:6])
day=int(startDate[6:])

hour=int(startTime[0:2])
minute=int(startTime[2:4])
second=int(float(startTime[4:]))


starttime=datetime.datetime(year,mon,day,hour,minute,second)

mtime=os.path.getmtime(testfile)
lastEditDate=datetime.datetime.fromtimestamp(mtime)




content=f.readlines()
# search firtst time output
for line in content:
    if 'TIME=' in line:
        break
ind1=line.index('TIME=')
ind3=line.index('\n')
line=line[ind1+5:ind3]
simtime0=float(line)
#f.close()

#f=open(testfile)
# search last time output

finished=0
for line in reversed(content):
    if 'Run completed' in line:
	    ind1=line.index('20') # rund starts 20 something
	    ind2=line.index(',')
	    ind3=line.index('\n')
	    startDate=line[ind1:ind2]
	    startTime=line[ind2+2:ind3]
	    year=int(startDate[0:4])
	    mon=int(startDate[4:6])
	    day=int(startDate[6:])
	    hour=int(startTime[0:2])
	    minute=int(startTime[2:4])
	    second=int(float(startTime[4:]))
	    endtime=datetime.datetime(year,mon,day,hour,minute,second)
	    finished=1	
    break

print('run completeted: '+ str(finished==1))

if finished==0:
	endtime=lastEditDate


runtime=endtime-starttime

for line in reversed(content):
	if 'TIME=' in line:
		break
ind1=line.index('TIME=')
ind3=line.index('\n')
line=line[ind1+5:ind3]
simtime1=float(line)
f.close()
simtime=simtime1-simtime0

speedup=(simtime)/float(runtime.seconds)
print('current model day is ' + str(simtime1/86400))
print('speed up is ' + str(speedup))

# output
outname='speedup' + str(datetime.datetime.now())  + '.txt'
outname=outname.replace(' ','_')
f2=open(outname,'w')
f2.write('model start at \t' +  str(starttime) +'\n')
if finished==0:
	f2.write('latest model step at ' +  str(endtime) +'\n')
else:
	f2.write('run completeted at ' +  str(endtime) +'\n')
f2.write('total runtime ' +  str(runtime) +'\n')

f2.write('start time in model days ' +  str(simtime0/86400) +'\n')
f2.write('end time in model days ' +  str(simtime1/86400) +'\n')
f2.write('simulated time ' +  str(simtime/3600) + ' h =  ' + str(simtime/86400)  +'days \n')
f2.write('speed up against realtime  ' +  str(speedup) +'\n')
f2.close();
