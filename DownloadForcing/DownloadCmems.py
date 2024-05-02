#!/gpfs/home/jacobb/anaconda3/envs/cmt_1.0/bin/python

import datetime as dt
import copernicusmarine as cm

# time 
t00=dt.datetime(2023,1,1,0,0,0,0) #start time download 
t11=dt.datetime(2023,1,3,23,0,0,0) #endtime download time download 
step=dt.timedelta(hours=24)  		 # step witdth e.g. each 24h a new file
endhour=step-dt.timedelta(hours=1)

#folder
outdir='./NWS/'

#name
#default name currently
#download_NWS_ibi.py
# elevation
# Only from AUgust 2023 for 
dmax=300

for i in range((t11-t00).days):
    t0=t00+step*i
    t1=t0+endhour
    print(t0)
    print(t1)
    
    # instantanours elev
    cm.subset(
      dataset_id="cmems_mod_nws_phy_anfc_0.027deg-2D_PT15M-i",
      #dataset_version="202309",
      variables=["zos"],
      minimum_longitude=-15.996014595031738,
      maximum_longitude=9.977008819580078,
      minimum_latitude=46.003639221191406,
      maximum_latitude=61.28188705444336,
      start_datetime=t0.strftime("%Y-%m-%dT%H:%M:%S"),
      end_datetime=t1.strftime("%Y-%m-%dT%H:%M:%S") ,
      #start_datetime="2023-02-01T00:00:00",
      #end_datetime="2023-02-01T01:00:00",
      force_download=True,
      output_directory=outdir
       #--force-download
    )

    #3D vars
    cm.subset(
      dataset_id="cmems_mod_nws_phy_anfc_0.027deg-3D_PT1H-m",
      #dataset_version="202309",
      variables=["uo","vo","so","thetao"],
      #variables=["uo"],
      minimum_longitude=-15.996014595031738,
      maximum_longitude=9.977008819580078,
      minimum_latitude=46.003639221191406,
      maximum_latitude=61.28188705444336,
      start_datetime=t0.strftime("%Y-%m-%dT%H:%M:%S"),
      end_datetime=t1.strftime("%Y-%m-%dT%H:%M:%S") ,
      minimum_depth=0.4940253794193268,
      maximum_depth=dmax,
      output_directory=outdir
    )
