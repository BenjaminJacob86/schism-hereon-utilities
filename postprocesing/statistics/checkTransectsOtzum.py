from glob import glob
import nu

s.init_node_tree(latlon=False)

trans=[]
for file in glob('*.bp'):
	if not '_ll' in file:
		trans.append(np.loadtxt(file,skiprows=2)[:,1:3])



nns=[s.node_tree_xy.query(list(zip(trans[i][:,0],trans[i][:,1])))[1] for i in range(len(trans))]

s.depths=np.asarray(s.depths)
depths=[s.depths[nns[i]] for i in range(len(trans))]
L=[ np.hstack((0,np.cumsum(np.sqrt((np.diff(trans[i],axis=0)**2).sum(axis=1))))) for i in range(len(trans)) ]

dists=[

trans
s.depths[nns[i]] for i in range(len(trans))]

for i in range(len(trans)):
	plt.subplot(3,2,i+1)
	plt.plot(L[i],-depths[i])