""" Counduct Experiment with different TimeSteps and Compare """
import os
import numpy as np

### time experiments
dts=[20,40,60]


def get_parameter(dir,param='dt'):
	f=open(dir+'param.nml')
	for line in f.readlines():
		if param+' =' in line:
			param= line.split('=')[1].split('!')[0]
			try:
				param=float(param)
			except:
				pass
			break		
	return param		

def set_parameters(dir,dir2,params,values):
	f=open(dir+'param.nml')
	fout=open(dir2+'param.nml','w')
	for line in f.readlines():
		for param,value in zip(params,values):
			if param+' =' in line:
				line=' {:s} = {:.0f} !'.format(param,value)+line.split('!')[1]+'\n'
				values.remove(value)
				params.remove(param)
		fout.write(line)
	#os.rename(dir+'param.nml',dir+'param.nml0')	
	#os.rename(dir+'param_new.nml',dir+'param.nml')	

	
def set_time(dir1,dir2,dt,rnday=15):
	params=['dt','nspool','ihfskip','nhot_write','nspool_sta']	
	values=[get_parameter(dir1,param=param) for param in params]
	values=[dt]+list(np.asarray(values[1:])*values[0]/dt)
	params.append('rnday')
	values.append(rnday)
	set_parameters(dir1,dir2,params=params,values=values)
	
maindir=os.getcwd()+'/'
for dt in dts:	
	
	if get_parameter(maindir,param='dt') != dt:
		expdir='run_dt_{:d}'.format(dt)
		os.system('mkdir {:s}'.format(expdir))		
		os.chdir(expdir)
		os.system('mkdir outputs combined hotstarts')		
		set_time(maindir,os.getcwd()+'/',dt,rnday=15)
		os.system('ln -s ../* .')		
		os.chdir(maindir)
	

# launch Jobs	
#RES=$(sbatch --dependency=afterok:$runJobID subroutines/createLastHot.sh $rundir $hotcombine $nproc $ntracer )



import numpy as np
import glob
ncdir='./combined/'
nstaskMax=48 #
step=4       # stacks to compile per tasks

def get_parameter(dir,param='dt'):
	f=open(dir+'param.nml')
	for line in f.readlines():
		if param+' =' in line:
			param= line.split('=')[1].split('!')[0]
			try:
				param=float(param)
			except:
				pass
			break		
	return param

ihfskip=get_parameter('./','ihfskip')
dt=get_parameter('./','dt')
	
schismfiles=[] 
for iorder in range(8): # check for schout_nc files until 99999
	schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
max_stack_combined=list(np.asarray(nrs)[np.argsort(nrs)])[-1]

with open('outputs/mirror.out') as f:
	for line in reversed(f.readlines()):
		if 'TIME=' in line:
			break
	ind1=line.index('TIME=')
	ind3=line.index('\n')
	line=line[ind1+5:ind3]
	simtime1=float(line)
max_stack_computed=np.floor(simtime1/(ihfskip*dt))

stacksToCombine=max_stack_computed-max_stack_combined

tasks=np.floor(stacksToCombine/step)

if tasks/nstaskMax <=1 1:
	nstaskMax=tasks
	iterations=1
else:
	iterations=np.ceil(stacksToCombine/nstaskMax)
	
with open('CombineMultipleNCs.batch') as filein
	with open('CombineMultipleNCs.batch_new','w') as filout:
		for line in filein.readlines():
			if '#SBATCH --ntasks-per-node=' in line:
				line='#SBATCH --ntasks-per-node='+str(np.int(ntasks)))
			if 'stack0=' in line:
				line='stack0='+str(np.int(max_stack_combined+1)))
			if 'stack1=' in line:
				line='stack0='+str(np.int(max_stack_computed+1)))
				
				
			fileout.writeline(line)

import fileinput
with fileinput.FileInput('CombineMultipleNCs.batch', inplace=True, backup='.bak') as file:
	ntasks=12
	for line in file:
		if '#SBATCH --ntasks-per-node=' in line:
			print(line.replace('#SBATCH --ntasks-per-node=', '#SBATCH --ntasks-per-node='+str(np.int(ntasks)))+'!', end='\n')
			#print(line.replace(line, '#SBATCH --ntasks-per-node='+str(np.int(ntasks))), end='')
			

	
import glob	

schout_10.nc  schout_11.nc  schout_12.nc  schout_1.nc  schout_2.nc  schout_3.nc  schout_4.nc  schout_5.nc  schout_6.nc	schout_7.nc  schout_8.nc  schout_9.nc