import xarray as xr
import datetime as dt
import numpy as np
import os

class param:		
	"""	functions for param.in for reading and editing. Operates in local directory """
	import os
	def __init__(self,fname='param.nml',comments='!'):
		#self.param_in=np.asarray(np.loadtxt(fname,comments=comments)[1:-1],int)-1		
		
		if '/' in fname:
			islash=fname.rindex('/')
			self.dir=fname[:islash]		
		else:
			self.dir='./'
			
		f=open(self.dir+'param.nml')	
		self.lines=f.readlines()
		f.close()
		
	def get_parameter(self,param='dt'):
		""" read parameter from param.in"""
		
		for line in self.lines:
			if param+' =' in line:
				param= line.split('=')[1].split('!')[0]
				try:
					param=float(param)
				except:
					param=str(param)
				break
		return param

	def set_parameter(self,params,values,outname='param.nml',outdir='./'):
		"""set_parameter(self,params,values,outname='param.nml',outdir='./') change parameters in param.in """
		if outname=='param.nml':
			try:
				os.rename('param.nml','param.nml.bkp')
			except:
				pass
		
		if type(params) == str:
			params=[params,]
			values=[values,]
		fout=open(outdir+outname,'w') 
		for line in self.lines:
			for param,value in zip(params,values):
				if param+' =' in line:
					line=' {:s} = {:.0f} !'.format(param,value)+line.split('!')[1]+'\n'
					values.remove(value)
					params.remove(param)
			fout.write(line)		

		fout.close()		
		# load updated param.nml
		print('updated param.nml has been loaded and will be accessed by get_parameters')	
		f=open(outdir+outname)	
		self.lines=f.readlines()
		f.close()
			
	def set_time_step(self,dt,outname='param.nml',outdir='./',rnday=None):
		""" set_time_step(self,dt,outname='param.nml',outdir='./',rnday=None) updates time step (dt) and all time related parameters (nspool,ihfskip,nhot_write,nspool_sta) to maintain output timings. rnday != None changes simulate length to specified Nr of days. dt new time step """
		params=['dt','nspool','ihfskip','nhot_write','nspool_sta','wtiminc']
		values=[self.get_parameter(param=param) for param in params]
		values=[dt]+list(np.asarray(values[1:])*values[0]/dt)
		if rnday != None:
				params.append('rnday')
				values.append(rnday)
		self.set_parameter(params=params,values=values,outname=outname,outdir='./')		



# schism
p=param()
reftime=np.asarray(dt.datetime(np.int(p.get_parameter('start_year')),
np.int(p.get_parameter('start_month')),
np.int(p.get_parameter('start_day')),
np.int(p.get_parameter('start_hour')),0,0),np.datetime64)

ds=xr.open_dataset('hotstart.nc')
ds['time']
hottime_schism=reftime+ds['time'].values*np.timedelta64(1,'s')


with open('runnr') as f:
	count=np.int(f.readline()[:2])

wwmhotname='wwm_hotfile_out_{:d}.nc'.format(count)
ds2=xr.open_dataset(wwmhotname)
hottime_wwm=ds2['ocean_time'].values

date=np.asarray(hottime_wwm[-1],dt.datetime)
strdate=str(date)

strdate=strdate.replace('-','')[:8]

with open('wwminput.nml-template') as f:
	lines=f.readlines()
	linenrs=[]
	for iline,line in enumerate(lines):
		if ('BEGTC' in line) and  ('Begin time of the wave boundary file' not in line) and (' PROC%BEGTC' not in line):
			lines[iline]=line.replace(line[line.index("'")+1:line.rindex(".0")],strdate)
		elif 	'FILEHOT_IN' in line: # update hotstart filenake
			lines[iline]=line.replace(line[line.index("'")+1:line.rindex("'")],wwmhotname)
		elif line[:6]==' LHOTR':
			lines[iline]=line.replace('= F','= T')
			break
	
with open('wwminput_new.nml','w') as fout:	
	fout.writelines(lines)		


exit()

