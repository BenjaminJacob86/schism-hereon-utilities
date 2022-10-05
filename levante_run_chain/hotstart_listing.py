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
reftime=np.asarray(dt.datetime(int(p.get_parameter('start_year')),
int(p.get_parameter('start_month')),
int(p.get_parameter('start_day')),
int(p.get_parameter('start_hour')),0,0),np.datetime64)


from glob import glob

folders=glob('Veg_*')

for folder in folders:
	os.chdir(folder)
	
	schismhots=glob('hotstart.nc_*')
	wwmhots=glob('wwm_hotfile_out_*_org.nc')

	
	with open('hottime_list','w') as f:
		for wwmhot in np.sort(wwmhots):
			try:
				ds2=xr.open_dataset(wwmhot)
				hottime_wwm=ds2['ocean_time'].values
				ds2.close()
				f.write(wwmhot +' ' + str(hottime_wwm) + '\n')
			except:	
				f.write(wwmhot +' invalid \n')
		for schismhot in np.sort(schismhots):
			try:
				ds=xr.open_dataset(schismhot)
				hottime_schism=reftime+ds['time'].values*np.timedelta64(1,'s')
				f.write(schismhot +' ' + str(hottime_schism) + '\n')
				ds.close()
			except:	
				f.write(schismhot +' invalid \n')
	os.chdir('../')