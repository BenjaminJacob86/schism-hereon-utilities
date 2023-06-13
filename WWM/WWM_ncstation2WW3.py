The frequencies of the model are put in an exponential grid FR(1:MSC)
with SC)=FRHIGH,FR(1)=FRLOW, FR(M and FR(i+1) = XFR FR(i). The
value of XFR should be around 1.1


f=FR(1)=FRLOW, FR(MSC)=FRHIGH, and FR(i+1) = XFR FR(i)

import numpy as np

msc=ds.dims['nfreq']
f=np.zeros(msc)
f[0]=ds.frlow
for i in range(1,len(f)):
	f[i]=f[i-1]*1.1


# what is frequency 1 and frequency 2 in ww3?
! - frequency : unit Hz, center band frequency - linear log scale (XFR
! - frequency1 : unit Hz, lower band frequency
! - frequency2 : unit Hz, upper band frequency

ds.spsig

# code from carolina about computation of f in WAM
# Is it what WWM is doing? I dont get the formulation in wwm
f_min=ds.frlow.values
f_max=ds.frhigh.values
log_freq = np.logspace(np.log10(f_min), np.log10(f_max), num=msc)


dsin=xr.open_mfdataset('wwm_station_0001.nc')

#
msc=ds.dims['nfreq']
f_min=dsin.frlow.values
f_max=dsin.frhigh.values
log_freq = np.logspace(np.log10(f_min), np.log10(f_max), num=msc)
frequency=log_freq

msd=ds.dims['ndir']
direction=np.linspace(0,360,msd)

dsin.ocean_time


(ocean_time: 16, ndir: 24, nfreq: 36, nbstation: 33)
efth=dsin.AC.values
efth.swapaxes(3,1)

 * time       (time) datetime64[ns] 2013-11-01 ... 2014-01-01
  * station    (station) float64 1.0 2.0 3.0 4.0 5.0 ... 954.0 955.0 956.0 957.0
  * frequency  (frequency) float32 0.04 0.044 0.0484 ... 0.929 1.022 1.124
  * direction  (direction) float32 90.0 75.0 60.0 45.0 ... 135.0 120.0 105.0

  
# get coordindates 

def split_line():
	varname,vals=line.split('=')
	if '!' in vals:
		vals,comment=vals.split('!'	)
		comment='!'+comment
	else:
		comment='\n'
	return varname,vals,comment
	
with open('wwminput.nml') as fin:
	for line in fin.readlines():
		if ('XOUTS' in line):
			varname,vals,comment= split_line()
			xouts=vals.strip()
		elif ('YOUTS' in line):
			varname,vals,comment= split_line()
			youts=vals.strip()

xouts=[np.float(xi) for xi in xouts.split(',')]
youts=[np.float(yi) for xi in youts.split(',')]			

lon
lat


# write netcdf  
dstime = xr.DataArray(name="time",data=seconds,dims=["time"], coords=dict(
						time=(["time"],seconds),
					),					attrs=attrs,
					)			  
  
