"""
merge schism nudginging files created for adajcent periods
by appending into the first file.
Both files are required to have the same time step and to be one time step appart
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

files=np.sort(glob('ERA5wnd_2017??.nc'))[:4]

mon0=files[0]
mon0=mon0[mon0.index('.')-2:mon0.index('.')]
mon1=files[-1]
mon1=mon1[mon1.index('.')-2:mon1.index('.')]
year=files[0]
year=year[year.index('.')-6:year.index('.')-2]

fnameout=''.join(('ERA5wnd_',year,mon0,'to',mon1,'.nc'))


shutil.copy('./'+files[0], './'+fnameout)
files=files[1:]


nca=Dataset(fnameout, "a")
for file in files:
	print('appending '+ file)
	ncb=Dataset(file, "r")
	lena=len(nca['time']) 
	lenb=len(ncb['time']) 
	
	nca['time'][lena:lena+lenb]=ncb['time'][:] # zero start
	for var in 'u10','v10':
		append_data=ncb[var][:]
		growdata=nca[var]
		growdata[lena:lena+lenb,:]=append_data
	ncb.close()
nca.sync()	
nca.close()	
