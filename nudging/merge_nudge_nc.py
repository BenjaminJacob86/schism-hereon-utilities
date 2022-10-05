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

import netCDF4
import numpy as np

#files1=['SAL_nu.nc_to_352','TEM_nu.nc_to_366'] # files to be appendend
#files2=['SAL_nu.nc_353_to','TEM_nu.nc_353_to'] # files to append from


files1=['SAL_nu.nc_bis715','TEM_nu.nc_bis715'] # files to be appendend
files2=['SAL_nu.nc_from_716','TEM_nu.nc_from_716'] # files to append from


for ifile in range(len(files1)):
	appendfile=files1[ifile]
	addfile=files2[ifile]

	nca=netCDF4.Dataset(appendfile, "a")
	ncb=netCDF4.Dataset(addfile, "r")

	if ifile== 0:
		lena=len(nca['time']) 
		lenb=len(ncb['time']) 
		
	nca['time'][lena:lena+lenb]=nca['time'][lena-1]+ncb['time'][0]+ncb['time'] # zero start
		
	# apppend
	for i in range(len(ncb['time'])):
		nca['tracer_concentration'][lena+i,:]=ncb['tracer_concentration'][i,:]
	nca.sync()	
	nca.close()	
	ncb.close()	
