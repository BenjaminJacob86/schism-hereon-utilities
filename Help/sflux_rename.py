from glob import glob
import numpy as np
import os

files=glob('*1.*.nc')

n=0
for i,file in enumerate(files):
		prefix,nr,suffix=file.split('.')
		nr='{:04d}'.format(int(nr))
		outname='.'.join((prefix,nr,suffix))
		
		os.rename(file,outname)


 /work/gg0028/g260114/RUNS/GermanBight/GB_timeslices/GB1995_wave/
 outputs_0