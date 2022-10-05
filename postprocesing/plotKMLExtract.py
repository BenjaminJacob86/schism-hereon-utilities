##############
fname='EuropeBasinExtracts.nc'
plt.ion()
nc=Dataset(fname)   
t0=dt.datetime.strptime(nc['time'].units[11:],'%Y-%m-%d %H:%M:%S')
dates=t0+np.arange(1,len(nc['time'])+1)*dt.timedelta(hours=6)
year=np.asarray([ti.year for ti in dates])
#names=[key[:key.index('_')]  for key in nc.variables.keys() if key != 'time' ]
names=[key[:key.rindex('_')]  for key in nc.variables.keys() if key != 'time' ]

meanref=dates[3::4]
daymeans=dict.fromkeys(list(nc.variables.keys())[1:])
monmeans=dict.fromkeys(list(nc.variables.keys())[1:])
yearmeans=dict.fromkeys(list(nc.variables.keys())[1:])
years=np.asarray([date.year for date in dates])		
months=np.asarray([date.month for date in dates])		
for varname in nc.variables.keys():
	if varname != 'time':
		daymeans[varname]=(nc[varname][3:-3]+nc[varname][4:-2]+nc[varname][5:-1]+nc[varname][6:])[::4]/4
		
		ts=nc[varname][:]
		tsmean=[]

		if ts.shape[-1]==2:
			yearmeans[varname]=np.asarray([ts[years==year].mean(axis=0) for year in range(2012,2017)])
			for year in range(2012,2017):
				for mon in range(1,13):
					tsmean.append(ts[(years==year) & (months==mon), :].mean(axis=0))
		else:
			yearmeans[varname]=np.asarray([ts[years==year].mean(axis=0) for year in range(2012,2017)])
			for year in range(2012,2017):
				for mon in range(1,13):
					tsmean.append(ts[(years==year) & (months==mon) ].mean(axis=0))
		monmeans[varname]=np.asarray(tsmean)

monref=[]		
for year in range(2012,2017):
	for mon in range(1,13):
		monref.append(dt.datetime(year,mon,15))
monref=np.asarray(monref)

yearref=np.asarray([dt.datetime(year,6,15) for year in range(2012,2017)])
ivalid=~np.isnan(monmeans['AlboranSea_mean_elev'])
ivalid=~np.isnan(monmeans[list(monmeans.keys())[0]])
fit=np.polyfit(np.arange(len(monref))[ivalid]+1,monmeans[varname][ivalid],1)




# daily means
plt.clf()
count=0
for i,varname in enumerate(list(nc.variables.keys())[1:]):
	#if 'elev' in varname:
	if len( daymeans[varname].shape)==1:
		plt.clf()
		plt.plot(meanref[:-1],daymeans[varname],linewidth=1)
		plt.plot(monref,monmeans[varname],linewidth=3)
		plt.plot(yearref,yearmeans[varname],'k',linewidth=2)
		fit=np.polyfit(np.arange(len(monref))[ivalid]+1,monmeans[varname][ivalid],1)
		plt.plot(monref[ivalid],np.polyval(fit,np.arange(len(monref))+1)[ivalid],'r')
		plt.text(monref[ivalid][-15],plt.ylim()[0]+np.diff(plt.ylim())*0.9,'y= {:f} x + {:f}'.format(fit[0],fit[1]),color='r')
		plt.title(names[i])
		plt.grid()
		#plt.ylim((-0.8,0))
		if 'elev' in varname: 
			plt.ylabel('ssh [m]')
		elif 'salt' in varname: 
			plt.ylabel('salinity [psu]')
		elif 'temp' in varname: 
			plt.ylabel('tempertature [degC]')
		elif 'pressure' in varname: 
			plt.ylabel('pressure [Pa]')
		plt.legend(['daily mean','monthly mean','annual mean','monthly trend'],ncol=2,loc='upper left')
		#plt.tight_layout()
		plt.savefig(varname+'.png',dpi=300)


!mkdir salt temp elev wind pressure
!mv *temp.png temp/
!mv *salt.png salt/
!mv *elev.png elev/
!mv *pressure.png pressure/
!mv *speed.png wind/

		
## daily means
#plt.clf()
#count=0
#for i,varname in enumerate(list(nc.variables.keys())[1:]):
#	if 'elev' in varname:
#		count+=1
#		plt.subplot(3,2,count)
#		plt.plot(meanref[:-1],daymeans[varname],linewidth=1)
#		plt.plot(monref,monmeans[varname],linewidth=3)
#		plt.plot(yearref,yearmeans[varname],'k',linewidth=2)
#		fit=np.polyfit(np.arange(len(monref))[ivalid]+1,monmeans[varname][ivalid],1)
#		plt.plot(monref[ivalid],np.polyval(fit,np.arange(len(monref))+1)[ivalid],'r')
#		plt.text(monref[ivalid][-15],plt.ylim()[1]*0.9,'y= {:f} x + {:f}'.format(fit[0],fit[1]),color='r')
#		plt.title(names[i])
#		plt.grid()
#		plt.ylim((-0.8,0))
#		if i%2==0:
#			plt.ylabel('ssh [m]')
#plt.legend(['daily mean','monthly mean','annual mean','monthly trend'],ncol=4)
#plt.tight_layout()