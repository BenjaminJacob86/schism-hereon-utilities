

#  adapted  from:
# https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids 
# 

# interpolation
#import scipy.interpolate as spint
import numpy as np
import scipy.spatial.qhull as qhull
from scipy.spatial import cKDTree



class unstructured_interpolant:
	""" Interpolation Object From Unstructrured Data or triangular grid nodes to xy coordinates. Call unstructured_interpolant( s,xun,yun,xq,yq,eff0=1e-12). Checks for Parents for Barycentric Interpolation and Next neighbours for nn interpolation as altornative for points outside triangles."""

		
	def __init__(self, s,xun,yun,xq,yq,eff0=1e-12):

		# points grid from
		self.xy=[[xun[i],yun[i]] for i in range(len(xun))] # pooint pairs of nodes

		# points grid to
		self.xyq=[[xq[i],yq[i]] for i in range(len(xq))] # pooint pairs of nodes



		def interp_weights(self):
		    tri = qhull.Delaunay(self.xy)
			tri.simplices=s.nvplt
		    simplex = tri.find_simplex(self.xyq) # -1 when not found
			# can i overwrite simlices?
			vertices = np.take(tri.simplices, simplex, axis=0)
		    temp = np.take(tri.transform, simplex, axis=0) # transform to barycentric coordinates
		    delta = self.xyq - temp[:, 2]

		    bary = np.einsum('njk,nk->nj', temp[:, :2, :], delta) # einstein summation
		    return tri,simplex,vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


		self.delauny,self.parent_tris,self.parent_nodes,self.parent_weights=interp_weights(self)


		self.no_parents=np.min(self.parent_weights,axis=1)<0


		# build next neighbour map
		self.xy_nn_tree = cKDTree(self.xy) # next neighbour search tree	
		self.nn_nodes=[]
		self.nn_dist=[]
		for ind in range(len(xq)):
			d,wo=self.xy_nn_tree.query((xq[ind],yq[ind]))
			self.nn_dist.append(d)
			self.nn_nodes.append(wo)
		self.nn_dist=np.asarray(self.nn_dist)
		self.nn_nodes=np.asarray(self.nn_nodes)


		# effective zero
		self.parent_weights[np.abs(self.parent_weights)<eff0]=0

	def interp_bary(self,nodevalues,fill_value=np.nan):
	    ret = np.einsum('nj,nj->n', np.take(nodevalues, self.parent_nodes), self.parent_weights)
	    ret[np.any(self.parent_weights < 0, axis=1)] = fill_value
	    return ret

	
	def interp_nn(self,nodevalues,nn_inds,tol_dist=300,fill_value=np.nan):
		ret=np.asarray(nodevalues)[self.nn_nodes[nn_inds]]
		ret[self.nn_dist[nn_inds]>tol_dist]=fill_value
		return ret





#interpolant=unstructured_interpolant(sin.x,sin.y,sout.x,sout.y)
##


# barycentric interpolation
#xout=interpolant.interp_bary(s.x)
#yout=interpolant.interp_bary(s.y)

# next neighbour interpolation for points outside triangle
#xout[interpolant.no_parents]=interpolant.interp_nn(np.asarray(s.x),interpolant.no_parents,tol_dist=643)
#yout[interpolant.no_parents]=interpolant.interp_nn(np.asarray(s.y),interpolant.no_parents,tol_dist=643)



