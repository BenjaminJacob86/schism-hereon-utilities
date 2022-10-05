"""
Link Temperature and salinity files for nudging script
gen_nudge_from_AMM15_check_varids.f90
"""

__author__ = "Benjamin Jacob"
__copyright__   = "Copyright 2018 - 03\2021 Helmholtz-Zentrum Geesthacht"
__copyright__   = "Copyright 04\2021 - Helmholtz-Zentrum Hereon GmbH"
__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"
from glob import glob
import numpy as np
import os

indir='/gpfs/work/jacobb/data/SETUPS/NorthSea/NorthSea200mHelgoland/RERUN2/nudging/neu/tsfiles0/' #'/work/gg0028/g260114/Data/Forcing/cmems_dayly/use_for_nudge_2007_2014/'
linkdir='tsfiles/'


input_is_integer=True   # different sorting necessary for T_1.nc T_n.nc ... then for   *20120601.nc  

startnr=353  # to append unfinished file

# assume to be same number of files for T and S
infilesT=list(np.sort(glob(indir+'*T_*.nc')))[startnr-1:]
infilesS=list(np.sort(glob(indir+'*S_*.nc')))[startnr-1:]

if input_is_integer:  # assumes number at end between '_' and '.nc'
	nrs=[np.int(file.split('_')[-1].split('.')[0]) for file in infilesT]
	isort=np.argsort(nrs)
	infilesT=np.asarray(infilesT)[isort]
	infilesS=np.asarray(infilesS)[isort]

n=0
for i,file in enumerate(infilesS):
        n=i+1
        cmd='ln -s {:s} {:s}S_{:d}.nc'.format(file,linkdir,n)
        os.system(cmd)
print('linked ' +str(n) +' files')

n=0
for i,file in enumerate(infilesT):
        n=i+1
        cmd='ln -s {:s} {:s}T_{:d}.nc'.format(file,linkdir,n)
        os.system(cmd)
print('linked ' +str(n) +' files')
