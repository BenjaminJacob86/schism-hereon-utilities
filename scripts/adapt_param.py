# take over parameters.py
import numpy as np

a=open('param.nml.jp','r')
b=open('param.nml-code-template','r')

linesA=a.readlines()
linesB=b.readlines()

params=[]
vals=[]
for line in linesA:
	if (line[0] != '!') & (len(line)>1) & ('=' in line): 
		param=line.split(' ')
		while '' in param:
			param.remove('')
		paramout=[]   # deal with missing psace after =
		

		
		for p in param:
			if ('=' in p) and (len(p)>1):
				p.split('=')
				paramout+=[p.split('=')]
			else:
				paramout+=[p]
				#paramout.append()
		param=np.hstack(paramout)				
			
		params.append(param[0])
		vals.append(param[2])


nomatch=[]
c=open('param_adapted_new.nml','w')
for line in linesB:

	#if 'lev_tr' in line:
	#	print('da')
	#	from IPython import embed; embed()
	#	break	

	# take over values
	if (line[0] != '!') & (len(line)>1) & ('=' in line): 
		line=line.replace('\n','')
		param=line.split(' ')
		#param=list(np.hstack([p.split('=') for p in param]))
		while '' in param:
			param.remove('')
		#param=[p.split('=') for p in param]# if no space after
		paramout=[]   # deal with missing psace after =
		for p in param:
			if ('=' in p) and (len(p)>1):
				p.split('=')
				paramout+=[p.split('=')]
			else:
				paramout+=[p]
				#paramout.append()
		param=np.hstack(paramout)	
		paramsi=param[0]
		vali=param[2]
		if '!rnday' in line:
			param[100]

		#if 'lev_tr_source' in line:
		#	from IPython import embed; embed()
		if paramsi in params:
			line=line.replace(vali,vals[ params.index(paramsi)],1)+'\n'
			
		else:
			line+='\n'
			# find parameter
			nomatch.append(paramsi)
	c.write(line)	
c.close()		

d=open('param_nwe__not_found_parameters_from_source.nml','w')
for line in nomatch:
	d.write(line +'\n')
d.close()

