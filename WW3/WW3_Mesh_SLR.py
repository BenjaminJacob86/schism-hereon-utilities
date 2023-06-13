# ww3 modify deptsh
import numpy as np

SLR=1
input_grid='NBSext_bl.msh'
output_grid=input_grid[:input_grid.index('.')]+'{:d}mSLR.msh'.format(SLR)

fin=open(input_grid,'r')
fout=open(output_grid,'w')
fin=open(input_grid,'r')

with open(input_grid,'r')as fin, open(output_grid,'w') as fout:
	for i,line in enumerate(fin.readlines()):
		if i==4:
			counter=np.float(line)
		if (i>4) & ((counter-1)>=0):
			counter-=1
			sep=line.split()
			sep[-1]=str(np.float(sep[-1])+SLR) #add SLR
			line=' '.join(sep)+ ' \n'
		fout.write(line)
