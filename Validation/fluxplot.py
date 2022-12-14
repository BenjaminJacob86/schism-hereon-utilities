# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 11:40:17 2019

@author: JacobB
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 17:13:37 2018

@author: JacobB
"""
#execute with /gpfs/home/jacobb/anaconda3/bin/python
#			/gpfs/home/jacobb/anaconda3/bin/ipython

import numpy as np
from matplotlib import pyplot as plt
import sys
#sys.path.insert(0,'/p/home/jusers/jacob4/juwels/shared')
sys.path.insert(0,'/gpfs/home/jacobb/code/python')
sys.path.insert(0,'/mnt/lustre01/pf/g/g260114/Programs/python/scripts')

sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')

from schism import *
import os
import datetime as dt



##### settings

setupdir='/gpfs/work/jacobb/data/discharge/'
fluxfile='/gpfs/work/jacobb/data/discharge/flux.out'


# look for repeating timestemp

#  0)totflux  1) pos flix 2) neg flux 3) pos salt flux 4) neg salt flux
#  0)totflux  1) pos flix 2) neg flux 3) pos salt flux 4) neg salt flux

# looks:												  
colors={'outflow':'r','inflow':'g'}
#MS=8											#markersize										  
# in: positive is in basin out: positive is out basin (flaxes are from larger to smaller number in fluxflag)

prefix='BS_full_'					# prefix for diagaram file names
##################################

# temp flux
os.chdir(setupdir)
s=schism_setup()

fluxflag=np.loadtxt('fluxflag.prop')[:,1]
#s.plotAtelems(fluxflag)
s.plotAtelems(fluxflag[s.nvplt2nvp])

p=param(setupdir+'/param.nml')
reftime=dt.datetime(np.int(p.get_parameter('start_year')),
np.int(p.get_parameter('start_month')),
np.int(p.get_parameter('start_day')),
np.int(p.get_parameter('start_hour')),0,0)

# rows are types  regions rage columns
m=np.loadtxt(fluxfile)
#n=734657

MS=8

nregs=len(m[0,:])-1
ntypes=np.where(np.diff(m[:100,0])>0)[0][0]+1
t=m[::ntypes,0]  # time in days
flux={'vol':{'tot':m[::ntypes,1:],'pos':m[1::ntypes,1:],'neg':m[2::ntypes,1:]},'salt':{'pos':m[3::ntypes,1:],'neg':m[4::ntypes,1:]}}
endtime=reftime+dt.timedelta(days=t[-1])

#cross_section_names=['Asov','Bosporus N','Bosporus C','Bosporus S','Dardanelle N','Dardanelle C','Dardanelle S','Bosporus narrow','Gibraltar','little Belt','Great Belt','Sound']
# in: positive is in basin out: positive is out basin (flaxes are from larger to smaller number in fluxflag)

# positive in elbe

direction=['in']*nregs #['out','out','out','out','out','out','out','in','in','in','in','in']

#transient
transient_days=0  # use data after 30 days
#startdate+=dt.timedelta(days=transient_days)
flux2=dict.fromkeys(flux.keys())
for key in flux.keys():
	flux2[key]=dict.fromkeys(flux[key].keys())
	for key2 in flux[key].keys():
		flux2[key][key2]=flux[key][key2][t>transient_days,:]


cross_section_names=['Elbe '+str(i) for i in range(nregs)]		
#cross_section_names=['Asov','Bosporus N','Bosporus C','Bosporus S','Dardanelle N','Dardanelle C','Dardanelle S']
#locnames=['Asov','Bosporus N','Bosporus C','Bosporus S','Dardanelle N','Dardanelle C','Dardanelle S','Bosporus','Gibraltar','little Belt','Great Belt','Sound']
#direction=['out','out','out','out','out','out','out']


# kerch asov (f)
limits=[]
for i,name in enumerate(cross_section_names):
	if 'Bosporus' in name:
		limits.append(5*10**4)
	elif 'Dardanelle' in name:
		if 'C' in name:
			limits.append(8*10**4)			
		else:
			limits.append(2.5*10**4)

	elif i > 9:
			limits.append(2.2*10**5)
		
	else:	
		limits.append(0)
#limits[cross_section_names.index('little Belt')]=0.8*10**5
#limits[cross_section_names.index('Asov')]=2.5*10**4
#limits[cross_section_names.index('Dardanelle C')]=8*10**4

# singular images - volume
# outflow is positive
# inflow is negative

#i=cross_section_names.index('Gibraltar')

for i in range(nregs):
		if direction[i]=='out':
			outflow=flux2['vol']['pos'][:,i]
			inflow=flux2['vol']['neg'][:,i]
			totflow=flux2['vol']['tot'][:,i]
		else:
			outflow=-flux2['vol']['neg'][:,i]
			inflow=-flux2['vol']['pos'][:,i]
			totflow=-flux2['vol']['tot'][:,i]

		# make plot	
		plt.clf()
		plt.clf()
		plt.plot(totflow,-inflow,','+colors['inflow'])#,Markersize=MS)              
		plt.plot(totflow,outflow,','+colors['outflow'])#,',Markersize=MS)              
		plt.xlabel('total flux [m^3/s]' )
		plt.ylabel('partial flux [m^3/s]')
		leg=plt.legend(['inflow(-)','outflow(+)'],frameon=False)
		for line, text in zip(leg.get_lines(), leg.get_texts()):
			text.set_color(line.get_color())
		plt.grid()
		plt.title(cross_section_names[i])
		if limits[i] != 0:
			plt.xlim((-limits[i],limits[i]))
			plt.ylim((0,limits[i]))
		plt.tight_layout()	
		plt.savefig(cross_section_names[i].replace(' ','')+'.png',dpi=300)
		
		
		
cmap=plt.cm.Accent(range(8))		
plt.figure()
plt.clf()
plt.subplot(2,2,1)
s.plotAtelems(fluxflag[s.nvplt2nvp],cmap=plt.cm.Accent)
plt.axis((35.82223142092876, 37.490492529540724, 44.98029282307836, 45.61701014946948))
plt.subplot(2,2,2)
s.plotAtelems(fluxflag[s.nvplt2nvp],cmap=plt.cm.Accent)
plt.axis((28.89103133419248, 29.273887022808957, 40.90079848359462, 41.30672171545255))
plt.subplot(2,2,3)
s.plotAtelems(fluxflag[s.nvplt2nvp],cmap=plt.cm.Accent)
plt.axis((25.96472477188497, 27.217658207780175, 39.76084334865623, 40.87766676909033))
plt.subplot(2,2,4)
outflows=[]
inflows=[]
for i in range(nregs):
	outflows+=[flux2['vol']['pos'][:,i].mean()]
	inflows+=[flux2['vol']['neg'][:,i].mean()]
outflows=np.asarray(outflows)
inflows=np.asarray(inflows)
	
shortnames=np.asarray([name[:4]+name[-1] for name in cross_section_names])
#plt.clf()
plt.bar(shortnames,np.abs(outflows/inflows),color=cmap[1:8]) 
plt.xticks(shortnames, shortnames, rotation=45)
plt.grid()
plt.ylabel('Ratio  <outflow>/<inflow> / 1')
for i in range(len(outflows)):
	plt.text(i-0.125,0.2,'{:.2f}'.format(np.abs(outflows/inflows)[i]),rotation=90) #cross_section_names[i]
plt.title('Ratio outflow:inflow ')
plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
plt.tight_layout()
plt.savefig('outflow_inflow_ratios.png',dpi=300)	
	
m3pers_to_km_3pera=86400*365/1000**3
plt.figure()
plt.clf()
plt.subplot(3,1,1)
plt.bar(shortnames,np.abs(inflows)*m3pers_to_km_3pera,color=cmap[1:8]) 
plt.grid()
plt.title('inflow [km^3/a]')
plt.subplot(3,1,2)
plt.bar(shortnames,np.abs(outflows)*m3pers_to_km_3pera,color=cmap[1:8]) 
plt.grid()
plt.title('outflow [km^3/a]')
netout=(np.abs(outflows)-np.abs(inflows))*m3pers_to_km_3pera
plt.subplot(3,1,3)
plt.bar(shortnames,netout,color=cmap[1:8]) 
plt.grid()
plt.title('outflow - inflow [km^3/a]')

for i in range(len(outflows)):
	plt.subplot(3,1,1)
	plt.text(i-0.125,50,'{:.2f}'.format(np.abs(inflows[i])*m3pers_to_km_3pera),rotation=90)	
	plt.subplot(3,1,2)
	plt.text(i-0.125,50,'{:.2f}'.format(np.abs(outflows[i])*m3pers_to_km_3pera),rotation=90)
	plt.subplot(3,1,3)
	plt.text(i-0.125,50,'{:.2f}'.format(netout[i]),rotation=90)
plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
plt.tight_layout()	
plt.savefig('outflow_inflow.png',dpi=300)

	
	
	
	
	
	
	
	

plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
fig.subplots_adjust(top=0.9)
plt.savefig('flux_diag.png',dpi=300)
plt.close()			


# out of bosporus is positive
plt.rc('font',**{'family':'normal','size':12})
i=3
plt.figure()
plt.plot(flux2['vol']['tot'][:,i],-flux2['vol']['neg'][:,i],'r*',Markersize=MS)              
plt.plot(flux2['vol']['tot'][:,i],flux2['vol']['pos'][:,i],'g*',Markersize=MS)   
           
plt.xlim((-2.2*1e5,2.2*1e5))
plt.ylim((0*1e4,2.2*1e5))
plt.title('central ' + cross_section_names[i] + ' EU tuned')


flux['vol']['tot'].mean()
              
plt.ion()

nrows=int(np.ceil(np.sqrt(nregs)))
ncols=int(np.round(np.sqrt(nregs)))


MS=3 # Markersize

# Black Sea 5 x 10^4
# Baltic Sea 2.5 10^5



plt.close('all')
fig=plt.figure()
for i in range(nregs):
    plt.subplot(nrows,ncols,i+1)
    tot=flux['salt']['pos'][:,i]+flux['salt']['neg'][:,i]
    plt.plot(t, np.cumsum(tot[:n]),'k')              
    plt.plot(t, np.cumsum(flux['salt']['pos'][:,i][:n]),'--')              
    plt.plot(t, np.cumsum(flux['salt']['neg'][:,i][:n]),'--')              
    plt.ylabel('vol salt [kg m^3/s]')
    plt.title(locnames[i])
    plt.legend(['tot','pos','neg'],frameon=False)
plt.tight_layout()
plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
fig.subplots_adjust(top=0.9)
plt.savefig('flux_salt_cumts.png',dpi=300)
plt.close()


fig.clf()
for i in range(nregs):
    plt.subplot(nrows,ncols,i+1)
    plt.plot(t, np.cumsum(flux['vol']['tot'][:,i][:n]),'k')              
    plt.plot(t, np.cumsum(flux['vol']['pos'][:,i][:n]),'--')              
    plt.plot(t, np.cumsum(flux['vol']['neg'][:,i][:n]),'--')              
    plt.ylabel('vol [m^3/s]')
    plt.title(locnames[i])
    plt.legend(['tot','pos','neg'],frameon=False)
plt.tight_layout()
plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
fig.subplots_adjust(top=0.9)
plt.savefig('flux_cumts.png',dpi=300)
plt.close()


# plt.figure()
# for i in range(nregs):
    # plt.subplot(nrows,ncols,i+1)
    # plt.plot(t, np.cumsum(flux['salt']['pos'][:,i]),'--')              
    # plt.plot(t, np.cumsum(flux['salt']['neg'][:,i]),'--')              
    # plt.plot(t, np.cumsum(flux['salt']['pos'][:,i]-np.cumsum(flux['vol']['neg'][:,i])),'k')              
    # plt.ylabel('vol salt [kg m^3/s]')
    # plt.title(locnames[i])
    # plt.legend(['pos','neg'],frameon=False)
# plt.tight_layout()
# plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
# plt.savefig('flux_saltts.png',dpi=300)
# plt.close()

fig.clf()
for i in range(nregs):
    plt.subplot(nrows,ncols,i+1)
    plt.plot(flux['vol']['tot'][:,i],flux['vol']['pos'][:,i][:n],',',Markersize=MS)              
    plt.plot(flux['vol']['tot'][:,i],-flux['vol']['neg'][:,i][:n],',',Markersize=MS)              
    plt.xlabel('total flux [m^3/s]' )
    plt.ylabel('partial flux [m^3/s]')
    plt.legend(['pos','neg'],frameon=False)
    plt.title(locnames[i])
plt.tight_layout()
plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
fig.subplots_adjust(top=0.9)
plt.savefig('flux_diag.png',dpi=300)
plt.close()

fig.clf()
for i in range(nregs):
    plt.subplot(nrows,ncols,i+1)
    tot=flux['salt']['pos'][:,i]-flux['salt']['neg'][:,i]
    plt.plot(tot,flux['salt']['pos'][:,i][:n],',',Markersize=MS)              
    plt.plot(tot,-flux['salt']['neg'][:,i][:n],',',Markersize=MS)              
    plt.xlabel('total salt \n flux [kg m^3/s]' )
    plt.ylabel('partial salt \n flux [kg m^3/s]')
    plt.legend(['pos','neg'],frameon=False)
    plt.title(locnames[i])
plt.tight_layout()
plt.suptitle(str(reftime)[:10] + ' - ' + str(endtime) )
fig.subplots_adjust(top=0.9)
plt.savefig('flux_diag_salt.png',dpi=300)
plt.close()



