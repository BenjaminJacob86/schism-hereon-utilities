"""
merge Era5 wind files
"""

__author__ = "Benjamin Jacob"
__copyright__   = "Copyright 2018 - 03\2021 Helmholtz-Zentrum Geesthacht"
__copyright__   = "Copyright 04\2021 - Helmholtz-Zentrum Hereon GmbH"
__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"

import os
from netCDF4 import Dataset,MFDataset
from glob import glob
import numpy as np
import shutil
import datetime as dt
from cftime import utime


datadir='/h/ksddata06/routine-ksd/dwd1/'
dates=np.loadtxt('run_period')
date0=dt.datetime.strptime(str(int(dates[0])),'%Y%m%d')
date1=dt.datetime.strptime(str(int(dates[-1])),'%Y%m%d')

files=[datadir+'dwdna'+ dt.datetime.strftime(date0+dt.timedelta(days=day),'%Y%m%d')+'.nc' for day in range((date1-date0).days)]

mon0='{:02d}'.format(date0.month)
mon1='{:02d}'.format(date1.month)
year='{:02d}'.format(date0.year)
fnameout=''.join(('dwdwnd_',year,mon0,'to',mon1,'.nc'))

os.remove('dwdwind.nc')
os.symlink(fnameout,'dwdwind.nc')

#shutil.copy('./'+files[0], './'+fnameout)
shutil.copy(files[0], './'+fnameout)
files=files[1:]


nca=Dataset(fnameout, "a")
utref=utime(nca['time'].units)
for file in files:
	print('appending '+ file)
	ncb=Dataset(file, "r")
	utadd=utime(ncb['time'].units)
	lena=len(nca['time']) 
	lenb=len(ncb['time']) 
	
	nca['time'][lena:lena+lenb]=utref.date2num(utadd.num2date(ncb['time'])) #ncb['time'][:] # zero start
	for var in 'u10','v10':
		append_data=ncb[var][:]
		growdata=nca[var]
		growdata[lena:lena+lenb,:]=append_data
	ncb.close()
nca.sync()	
nca.close()	