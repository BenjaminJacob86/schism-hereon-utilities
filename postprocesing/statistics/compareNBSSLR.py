#/gpfs/work/jacobb/data/shared/schism-hzg-utilities/
#calc_stats(rundir='./',ncdir='./outputs',varnames=['hvel','wind_speed','zcor','salt','temp'],lvl=-1,periods=['year']):

import sys
sys.path.insert(0,'/gpfs/home/jacobb/code/python/')   #strand
sys.path.insert(0,'/mnt/lustre01/pf/g/g260114/Programs/python/scripts/')  #mistral
sys.path.insert(0,'/sciclone/home20/bjacob01/git/schism-hzg-utilities/')   #vims
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/')   #strand
sys.path.insert(0,'/gpfs/work/jacobb/data/shared/schism-hzg-utilities/postprocesing/statistics')   #strand

from schism import *
from schism_statistics import calc_stats


# baroline ref


dir1='/gpfs/work/jacobb/data/SETUPS/NorthBalticSea/NBS/'
dir2='/gpfs/work/jacobb/data/SETUPS/NorthBalticSea/NBSSLR/'

# NBS combined is varocline
# NBS combined outputs is barotrope

# SLR combined ist barotrop
# SLR combined 1 ist baroklin

# baroclinic
dir1='/gpfs/work/jacobb/data/SETUPS/NorthBalticSea/NBS/outputs'
dir2='/gpfs/work/jacobb/data/SETUPS/NorthBalticSea/NBSSLR/combined1'



# org from outputs
#Ref=calc_stats(rundir=dir1, ncdir=dir1+'combined/sub/', varnames=['elev'], lvl=-1, periods=['year'])
Ref=calc_stats(rundir=dir1, ncdir=dir1+'combined_barotrop/', varnames=['elev'], lvl=-1, periods=['year'])
SLR=calc_stats(rundir=dir2, ncdir=dir2+'combined/', varnames=['elev'], lvl=-1, periods=['year'])

statsRef=Ref[0]
dryRef=Ref[1]
s1=statsRef['setup']

statsSlr=SLR[0]
drySlr=SLR[1]
s2=statsSlr['setup']

s2.plotAtnodes(np.asarray(s2.depths)-np.asarray(s1.depths))
plt.clim((1.99,2.01))
plt.savefig('Delta H')


plt.figure()
plt.subplot(2,2,1)
pdata=(np.sqrt(statsRef['elev']['var'][0]))
ph,ch=s1.plotAtnodes(pdata,mask=dryRef)
plt.xlim([-6,10])
plt.ylim([48,63])
plt.clim((0,1.8))
ch.set_label('std(ssh) [m]')
plt.title('Reference')

plt.subplot(2,2,2)
ph,ch=s2.plotAtnodes(np.sqrt(statsSlr['elev']['var'][0]),mask=drySlr)
plt.xlim([-6,10])
plt.ylim([48,63])
plt.clim((0,1.8))
plt.title('2 m SLR')
ch.set_label('std(ssh) [m]')


#baroclinic



plt.subplot(2,2,3)
plt.figure()
pdata=((np.sqrt(statsSlr['elev']['var'][0])-np.sqrt(statsRef['elev']['var'][0])))
#s=np.sign(pdata)
#pdata2=np.log10(np.abs(pdata))
#pdata2[s>0]=pdata2[s>0].min()+7
#pdata2[s<0]=-(pdata2[s<0].min()+7)
# shift
#inver
#pdata2[pdata2>0].min()+1
#pdata2[pdata2>0]
ph,ch=s1.plotAtnodes(pdata,mask=drySlr)
plt.clim((-0.1,0.1))
plt.xlim([-6,14])
plt.ylim([48,63])
plt.clim((-0.07,0.07))
ch.set_label('Delta std(ssh) [m]')
plt.title('2 m SLR -Ref')
plt.set_cmap('RdBu_r')

plt.subplot(2,2,4)
ph,ch=s2.plotAtnodes(((np.sqrt(statsSlr['elev']['var'][0])-np.sqrt(statsRef['elev']['var'][0]))),mask=drySlr)
plt.xlim([6,10])
plt.ylim([53,55.5])
plt.clim((-0.02,.1))
ch.set_label('Delta std(ssh) [m]')
plt.title('2 m SLR')


plt.figure()
plt.set_cmap('RdBu_r')
pdata=((np.sqrt(statsSlr['elev']['var'][0])-np.sqrt(statsRef['elev']['var'][0])))
ph,ch=s2.plotAtnodes(pdata,mask=drySlr)
#plt.ylim([48,63])
plt.clim([-0.1,0.1])
plt.title('2 m SLR')
ch.set_label('std(ssh) [m]')


