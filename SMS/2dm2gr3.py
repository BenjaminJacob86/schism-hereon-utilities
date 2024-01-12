
from glob import glob

nodes=[]
elems=[]

#fname='new15-upgraded_Bosporus100mNeu.2dm'

fname=glob('*.2dm')[0]
f=open(fname)
for line in f.readlines():
	
	if 'E3T' in line:
		splitted=line.split()
		elems.append(splitted[1] + ' 3 ' + ' '.join(splitted[2:-1])+'\n')
	elif 'E4Q' in line:	
		splitted=line.split()
		elems.append(splitted[1] + ' 4 ' + ' '.join(splitted[2:-1])+'\n')
	elif 'ND' in line:		
		nodes.append(line[3:])
gr3name=fname.split('.')[0]+'.gr3'

with open(gr3name,'w') as f:
	f.write(gr3name+'\n'+'{:d} {:d}\n'.format(len(elems),len(nodes)))		
	for node in nodes:
		f.write(node)
	for elem in elems:
		f.write(elem)