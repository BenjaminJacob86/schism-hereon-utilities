""" 
Add one line at beginning to shift th.nc
"""

import xarray as xr
import numpy as np
import os


from glob import glob

pattern='*416_days_nn.th.nc'

files=glob(pattern)

for file in files:
    print(file)
    
    ds=xr.open_dataset(file)
    
    time=ds.time.values
    timenew=np.hstack((time[0],time[1]+time))  # beginns with zerro
    
    data=ds.time_series.values
    datanew=np.concatenate((data[:1,:],data),axis=0)  # beginns with zerro
    
    dsnew = xr.Dataset(
		data_vars=dict(
			time_step=(["one"], ds.time_step.values),
			time_series=(["time", "nOpenBndNodes", "nLevels", "nComponents"], datanew),
		),
		coords=dict(
			time=timenew,   # copy time vector from prc
		),
		attrs=dict(description="Time lined added nc."),
	)
    os.rename(file,'org_'+file)
    dsnew.to_netcdf(file)
    
    ds.close()
    dsnew.close()