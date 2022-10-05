import sys
sys.path.insert(0,'/pf/g/g260114/git/hzg/schism-hzg-utilities')
import datetime as dt
from schism import *
from matplotlib import pyplot as plt

# append river data repeating years before
t0=dt.datetime(2012,6,1,0,0,0)
t1=dt.datetime(2020,1,1,0,0,0)
ndays=t1-dates[-1]


s=schism_setup()
river=sources(s)



dates=t0+ dt.timedelta(seconds=np.float(np.diff(river.vtime)[0]))*np.arange(len(river.vtime))
months=np.asarray([dates[i].month for i in range(len(dates))])
days=np.asarray([dates[i].day for i in range(len(dates))])

dates2=dates.copy()
msource2=river.msource.copy()
vsource2=river.vsource.copy()
for date2 in dates[-1]+np.arange(1,np.int(ndays.days))*dt.timedelta(seconds=np.float(np.diff(river.vtime)[0])):
	ind=np.where((date2.month==months) & 	(date2.day==days))[0]
	if len(ind)==0:
		ind=np.where((date2.month==months) & 	(date2.day-1==days))[0] # schalt jahr nutze tag zuve
	ind=ind[0]
	
	dates2=np.hstack((dates2,date2))
	msource2=np.vstack((msource2,msource2[ind,:]))	
	vsource2=np.vstack((vsource2,vsource2[ind,:]))	
	
	
#plt.plot(dates2,msource2[:,0])
seconds=np.asarray([(date-dates2[0]).total_seconds() for date in dates2])
msource2=np.hstack((seconds.reshape(len(seconds),1),msource2))	
vsource2=np.hstack((seconds.reshape(len(seconds),1),vsource2))	

np.savetxt('msource.th_added21_km_a'+str(str(t0).replace('-','')[:8]) + '_' +str(str(t1).replace('-','')[:8]),msource2)
np.savetxt('vsource.th_added21_km_a'+str(str(t0).replace('-','')[:8]) + '_' +str(str(t1).replace('-','')[:8]),vsource2)
	
	
	
	
	
	
