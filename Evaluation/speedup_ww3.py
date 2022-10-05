import datetime as dt

f=open('ww3run.out')
lines=f.readlines()

for i,line in enumerate(lines):
	if 'WAVEWATCH III calculating' in line:
		line0=line
		break	
		
for i,line in enumerate(lines[::-1]):
	if 'WAVEWATCH III calculating' in line:
		line1=line
		break
		
simtime0=line0.split('for')[1].split('UTC')[0][1:-1]
simtime0=dt.datetime.strptime(simtime0,'%Y/%m/%d %H:%M:%S')
starttime=dt.datetime.strptime(line0.split(' at ')[1][:-1],'%H:%M:%S')


simtime1=line1.split('for')[1].split('UTC')[0][1:-1]
simtime1=dt.datetime.strptime(simtime1,'%Y/%m/%d %H:%M:%S')
endtime=dt.datetime.strptime(line1.split(' at ')[1][:-1],'%H:%M:%S')

speedup=(simtime1-simtime0).total_seconds()/(endtime-starttime).total_seconds()
runtime=(endtime-starttime).total_seconds()/3600

# output
outname='speedup' + str(dt.datetime.now())  + '.txt'
outname=outname.replace(' ','_')
f2=open(outname,'w')
f2.write('model start at \t' +  str(starttime) +'\n')
#if finished==0:
#	f2.write('latest model step at ' +  str(endtime) +'\n')
#else:
#	f2.write('run completeted at ' +  str(endtime) +'\n')
f2.write('total runtime ' +  str(runtime) +'h \n')

f2.write('start time in model days ' +  str(simtime0) +'\n')
f2.write('end time in model days ' +  str(simtime1) +'\n')
f2.write('simulated time ' +  str((simtime1-simtime0).total_seconds()/3600) + ' h =  ' + str((simtime1-simtime0).total_seconds()/86400)  +'days \n')
f2.write('speed up against realtime  ' +  str(speedup) +'\n')
f2.close();

print('start time in model days ' +  str(simtime0) +'\n')
print('end time in model days ' +  str(simtime1) +'\n')
print('total runtime ' +  str(runtime) +'h \n')
print('simulated time ' +  str((simtime1-simtime0).total_seconds()/3600) + ' h =  ' + str((simtime1-simtime0).total_seconds()/86400)  +'days \n')
print('speed up against realtime  ' +  str(speedup) +'\n')