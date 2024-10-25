import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hereon-utilities/')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
s=schism_setup()

# control plots for
istart=0
iend=0
#plt.figure()
plt.clf()
i=0
iend+=len(s.bdy_segments[i])
inds=np.asarray(range(istart,iend))
istart=iend
frcbdnodes=s.bdy_segments[i]
bdnodes=np.asarray(frcbdnodes)[inds]#-1	
bddepths=np.asarray([s.vgrid[inode].filled(-1)*s.depthsdict[inode] for inode in bdnodes])
xcoord=np.tile(inds,[len(bddepths[0]),1]).T


schism_step=1261
tib4=int(np.floor(1261*100/3600))
tinext=int(np.ceil(1261*100/3600))

nancheck=[]
for file in ['TEM_3D.th.nc','SAL_3D.th.nc','uv3D.th.nc']:
    ds=xr.open_dataset('TEM_3D.th.nc')
    temp=np.isnan(ds.time_series[tib4:tinext+1,:,:,0])
    nancheck.append(temp)

    plt.figure()
    plt.subplot(1,2,1)
    plt.pcolor(xcoord,bddepths,temp[0,:])
    plt.colorbar()
    plt.title(tib4)
    plt.subplot(1,2,2)
    plt.pcolor(xcoord,bddepths,temp[1,:])
    plt.colorbar()
    plt.title(tinext)
    plt.suptitle(file + ' nan check')
    try:
        plt.savefig('nan_check_step'+str(schism_step)+file[:file.index('.')]+'.png',dpi=300)
    except:
        pass