import numpy as np
import warnings
xy=np.loadtxt('points_dissa.txt',comments='%')[:,:2]
 

# append in
np=xy.shape[0]

PID='P' # point identifier string to be numerate
off=1   # offset 1: start coutnign  from one 



names=","+",".join(["'P"+"{:d}'".format(off+i) for i in range(np)])
XOUTS=','+','.join(['{:f}'.format(x) for x in xy[:,0]])
YOUTS=','+','.join(['{:f}'.format(x) for x in xy[:,1]])
CUTOFF=(','+'0.0,'*np)[:-1]

def split_line():
	varname,vals=line.split('=')
	if '!' in vals:
		vals,comment=vals.split('!'	)
		comment='!'+comment
	else:
		comment='\n'
	return varname,vals,comment


with open('wwminput.nml-template0') as fin, open('wwminput.nml-template','w') as fout:
	for line in fin.readlines():
		if ('ILOUTS' in line) or ('IOUTS' in line):
			varname,vals,comment= split_line()
			vals='{:d}\t'.format(int(vals)+np)
			if int(vals) > 50:
				warnings.warn('Number of specified stations succeds wwm limit of 50')	
			line=varname+'='+ ' ' + vals + comment
		elif ('NLOUTS' in line):
			varname,vals,comment= split_line()
			vals=vals.strip()
			vals+=names
			line=varname+'='+ ' ' + vals + comment
		elif ('NOUTS' in line):
			varname,vals,comment= split_line()
			vals=vals.strip()
			vals+=names
			line=varname+'='+ ' ' + vals + comment
		elif ('XOUTS' in line):
			varname,vals,comment= split_line()
			vals=vals.strip()
			vals+=XOUTS
			line=varname+'='+ ' ' + vals + comment
		elif ('YOUTS' in line):
			varname,vals,comment= split_line()
			vals=vals.strip()
			vals+=YOUTS
			line=varname+'='+ ' ' + vals + comment
		elif ('CUTOFF' in line):
			varname,vals,comment= split_line()
			vals=vals.strip()
			vals+=CUTOFF
			line=varname+'='+ ' ' + vals + comment
		fout.write(line)	