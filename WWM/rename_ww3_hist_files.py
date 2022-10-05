# after each restart wwm3 starts counting from 1
# rename after timestep
from glob import glob
folders=glob('wwmhist*')
import xarray as xr
import os
import numpy as np
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
		params=['dt','nspool','ihfskip','nhot_write','nspool_sta']
		values=[self.get_parameter(param=param) for param in params]
		values=[dt]+list(np.asarray(values[1:])*values[0]/dt)
		if rnday != None:
				params.append('rnday')
				values.append(rnday)
		self.set_parameter(params=params,values=values,outname=outname,outdir='./')		


p=param()				

tref=np.datetime64('{:4.0f}-{:02.0f}-{:02.0f}T00:00:00.000000000'.format(p.get_parameter('start_year'),p.get_parameter('start_month'),p.get_parameter('start_day')))+np.timedelta64(0,'D')
#2017-01-26T00:00:00.000000000

# assume daily		
for folder in folders:
		os.chdir(folder)
		files=np.sort(glob('wwm_hist_*.nc'))
		for file in files:
			file
			dsi=xr.open_dataset(file)
			day=(dsi['ocean_time'][0].values-tref)/np.timedelta64(1,'D')+1
			dsi.close()
			os.rename(file,'wwm_hist_{:04.0f}.nc'.format(day))
		os.chdir('../')
		
