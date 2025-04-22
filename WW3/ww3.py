""" read ww3 grid """

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import xarray as xr
plt.ion()

class WW3_mesh():
	import numpy as np
	from matplotlib import pyplot as plt
	from matplotlib.collections import PolyCollection
	import xarray as xr
	def __init__(self,file):
		with open(file) as f:
			
			lines=f.readlines()

			nodestart=lines[lines.index('$Nodes\n')+1]
			nodestart=lines.index(nodestart)
			nnodes=int(lines[lines.index('$Nodes\n')+1].split(' ')[0])

			elemstart=lines[lines.index('$Elements\n')+1]
			elemstart=lines.index(elemstart)
			nelem=int(lines[lines.index('$Elements\n')+1].split(' ')[0])

			self.nodes=np.loadtxt('NBSext_bl.msh',skiprows=nodestart+1,max_rows=nnodes)[:,1:]
			self.elems=np.loadtxt('NBSext_bl.msh',skiprows=elemstart+1,max_rows=nelem)[:,6:]-1
			self.elems=np.asarray(self.elems,int)
			self.x,self.y,self.d=self.nodes[:,0],self.nodes[:,1],self.nodes[:,2]
			
		  # plot functions 
	def plotAtnodes(self,nodevalues,cmap=plt.cm.jet,mask=None):
		"""
		visualisation routine plotting triangles at nodes (quads are splitted)
		"""
		ph=plt.tripcolor(self.x,self.y,self.elems,nodevalues,shading='flat',mask=mask,cmap=cmap)	  
		ch=plt.colorbar()

		return ph,ch		
		
	def plot_mesh(self,tri_color='k',quad_color='m',linewidth=0.2):	  
		xy=np.c_[self.x,self.y]
		tripc = PolyCollection(xy[self.elems,:3],facecolors='none',edgecolors=tri_color,linewidth=linewidth) #, **kwargs)
		plt.gca().add_collection(tripc)
		
	def load_points(self,fname,sep=' '):	  
		with open(fname) as f:
			x,y,names=[],[],[]
			for line in f.readlines():
				line=line.split(' ')
				xi,yi,name=[li for li in line if li != '']
				x.append(np.float(xi))
				y.append(np.float(yi))
				names.append(name.split('\n')[0])
			self.px=np.asarray(x)
			self.py=np.asarray(y)
			self.pnames=np.asarray(names)
	def plot_points(self,showname=False):	  
		plt.plot(self.px,self.py,'ko')
		if showname:
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,' '+ name,color='w')
		else:	
			count=0
			for xi,yi,name in zip(self.px,self.py,self.pnames):
				plt.text(xi,yi,str(count))		
				count+=1

				
ww3=WW3_mesh('NBSext_bl.msh')
ww3.load_points('GB_bd_points.list')
#ww3.load_points('../points.list',sep=' ')
plt.axis((4, 9.628056711758582, 53.054607859066536, 55.8))
#ww3.load_points('/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd_2017/points.list')
	
				

ds=xr.open_dataset('ww3.201701_spec.nc')
x=ds['longitude'][0,:].values
y=ds['latitude'][0,:].values






plt.figure()
ww3.plotAtnodes(ww3.d,cmap=plt.cm.jet,mask=None)
plt.plot(ww3.px[:321:2],ww3.py[:321:2],'ko')

xreduced=np.hstack((ww3.px[:321:2],ww3.px[362:]))
yreduced=np.hstack((ww3.py[:321:2],ww3.py[362:]))
nreduced=np.hstack((ww3.pnames[:321:2],ww3.pnames[362:]))
f=open('points_reduced.list','w')
for xi,yi,ni in zip(xreduced,yreduced,nreduced):
			f.write('   {:f}      {:f}  {:s}\n'.format(xi,yi,ni))
f.close()

ww3.load_points('points_reduced.list')

plt.figure()
ww3.plotAtnodes(ww3.d,cmap=plt.cm.jet,mask=None)
#ww3.plot_mesh()
ww3.plot_points()
plt.plot(x,y,"2",color='r')
plt.clim((0,50))
plt.show()

