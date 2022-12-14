# quick and dirty setup conversion to other grid
# input: 
#	old setup
#	new grid
#
# from old setup translate:
#	bathymetry
#	vgrid
#	gr3 files
#	boundaries
#	bnd forcing
#
# requirements: 
#	grid is subdomain
#	and shares land and ocean boundaries with other setup
#	yet no quads
#	new grid should have only one land boundary created with xmgredit starting where the first ocean boundary shall begin.

# libraries
import os
import netCDF4
import sys
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import netCDF4
import pickle 

#own libraries
sys.path.insert(0,'/gpfs/home/jacobb/code/python/')
sys.path.insert(0,'/gpfs/work/jacobb/data/PreProc/Setup2Other/')
from schism import *
from schism_plots import *
from schism_interpolations import *



################# S E T T I N G S ##################################
# source_setup
#source_setup='/gpfs/work/jacobb/data/SETUPS/Europe/HighRes/'# '/gpfs/work/jacobb/data/SETUPS/Europe/HighRes/'
#dest_setup='/gpfs/work/jacobb/data/RUNS/HighResCut/' #/gpfs/work/jacobb/data/SETUPS/Europe/ImprovedDanishStraits/'
#source_setup='/gpfs/work/jacobb/data/RUNS/HighResCut/hot/'

source_setup='/gpfs/work/jacobb/data/SETUPS/NWBlackSea/setup/'
dest_setup='/gpfs/work/jacobb/data/SETUPS/NWBlackSea/setup2/'

if source_setup[-1]!='/':
	source_setup+='/'

if dest_setup[-1]!='/':
	dest_setup+='/'


## what inputs to clone from different setup
do_bnd=1        # transfer boundaries
do_vgrid=0      # transfer vgrid
do_gr3=1        # transfer gr3 files
do_ic=1        # transfer ic files
do_hot=0        # transfer hotsart
do_vm_source=0  # transfer sources 
do_bath=1       # transfer bathymetries

hot_outname='hotstart_transfered.nc'
###################################################################





######### transfer open bounadries ############################################

def transfer_bnd(sin,sout):

	# just one initial boundary to have concave hull
	# ensure that boundary order starts with open boundary 1 corrsponding to  'from'  grid
	if len(sout.land_segments)==1:
		## nn search try for bd nodes
		bd_nn_tree = cKDTree(list(zip(np.asarray(sout.x)[np.asarray(sout.land_segments[0])-1],\
		np.asarray(sout.y)[np.asarray(sout.land_segments[0])-1]))) # next neighbour search tree	
		
		open_bd = sin.bdy_segments[0]

		d0,i0=bd_nn_tree.query([sin.x[open_bd[0]-1],sin.y[open_bd[0]-1]])
		d1,i1=bd_nn_tree.query([sin.x[open_bd[-1]-1],sin.y[open_bd[-1]-1]])

		# first boundary if i0 > i1  boundary start is shifteted
		# >> reorder for more intuitive start
		if i0 > i1:
			bounding=sout.land_segments[0][i0:]
			bounding=bounding+sout.land_segments[0][:i1+1]
			bounding=bounding+sout.land_segments[0][i1+1:i0]
		else:
			bounding=sout.land_segments[0]

		# create new tree
		bd_nn_tree = cKDTree(list(zip(np.asarray(sout.x)[np.asarray(bounding)-1],\
		np.asarray(sout.y)[np.asarray(bounding)-1]))) # next neighbour search tree	

		# build boundaries adapted from old mesh
		new_land_segments=[]  # use indexing starting from one to be consistent with richards stuff
		new_bd_segments=[]
		bd_limits=[]

		f=open('river_merge.txt','w')

		# land starts after ocean boundary
		for i,open_bd in enumerate(sin.bdy_segments): #[:10]
			d0,i0=bd_nn_tree.query([sin.x[open_bd[0]-1],sin.y[open_bd[0]-1]])
			d1,i1=bd_nn_tree.query([sin.x[open_bd[-1]-1],sin.y[open_bd[-1]-1]])

			# len only 1 -> needs min 2 nodes
			if i0==i1:
				i1+=1
			elif i0>i1:
				i0,i1 = i1,i0
	

			# distinct rivers not resolved >> merge
			#if i>0 and new_bd_segments[-1]==bounding[i0:i1+1]:
			#	print('merge rivers ' +str(i)  + ' ' + str(i+1))
			#	f.writelines('merge rivers ' +str(i)  + ' ' + str(i+1) +'\n')
			#else:	
			new_bd_segments.append(bounding[i0:i1+1])
			bd_limits.append([i0,i1])
	
			if i > 0:  # shifted into elese ; ennough ?
				new_land_segments.append(bounding[bd_limits[i-1][1]:bd_limits[i][0]+1])
				#if len(new_land_segments[-1])==1:
				#	print('land seg ' +str(i)+ ' only 1')		

		if  bd_limits[-1][-1]!=len(bounding):  # does not end with land boundary
			#new_land_segments.append(bounding[bd_limits[i-1][1]:]) # here append
			new_land_segments.append(bounding[bd_limits[-1][-1]:])


		# corrections 
		# IF DUE TO COARSEN RIVERS OVERLAP Delete boundary and merge fluxes
		idble=[]
		for i in range(1,len(new_bd_segments)):
			if new_bd_segments[i]==new_bd_segments[i-1]:
				idble.append(i)

		for i in reversed(idble):
			del new_bd_segments[i-1]
			del new_land_segments[i-1]
			del bd_limits[i-1]
			f.writelines('merge rivers of open boundaries ' +str(i)  + ' ' + str(i+1) +'\n')
		f.close()


		#len land bd== 1 >> rivers adjecent remove land bd 


		# update boundary information in grid object
		sout.num_bdy_segments=len(new_bd_segments)
		sout.bdy_segments=new_bd_segments
		sout.bdy_nodes=[]
		for i in range(sout.num_bdy_segments):
			sout.bdy_nodes += sout.bdy_segments[i][:] 
		sout.num_bdy_nodes=len(sout.bdy_nodes)

		sout.num_land_segments=len(new_land_segments)+len(sout.island_segments)
		sout.land_segments=new_land_segments
		sout.land_nodes=[]
		for i in range(len(sout.land_segments)):
			sout.land_nodes += sout.land_segments[i][:] 
		for i in range(len(sout.island_segments)):
			sout.land_nodes += sout.island_segments[i][:] 
		sout.num_land_nodes=len(sout.land_nodes)

		return sout


	else:
		print('There are already several boundaries. For the boundary transfer from source to destination grid the latter has to have only one initial Land boundary which matches the domain boundary')



###### program start #####

# load source setup
os.chdir(source_setup)
sin=schism_setup(hgrid_file='hgrid.gr3',ll_file='hgrid.ll',vgrid_file='vgrid.in')

# load dest_setup
os.chdir(dest_setup)
sout=schism_setup(hgrid_file='hgrid.gr3')
#sout=schism_setup(hgrid_file='hgrid.gr3',ll_file='hgrid.ll',vgrid_file='vgrid.in') # reload for testing


# create unstructured interpolant object from  sin to sout
interpolant=unstructured_interpolant(sin.x,sin.y,sout.x,sout.y)

if sin.x[0]!=sin.lon[0]:
	in_cartesian=1
else:
	in_cartesian=0

if sout.x[0]!=sout.lon[0]:
	out_cartesian=1
else:
	out_cartesian=0


########### interp bathymetry
if do_bath:
	sout.depths=interpolant.interp_bary(sin.depths) # barycentric interpolation
	sout.depths[interpolant.no_parents]=interpolant.interp_nn(np.asarray(sin.depths),interpolant.no_parents,tol_dist=643) # next neighbour interpolation for  points outside triangle

	# calc at elements
#	for  key in sout.nvdict.keys()#sout.element_tree_ids:
#		sout.element_depth[key]=np.mean(sout.depths[np.asarray(sout.nvdict[key])-1])







# quads
do_plot=0
if do_plot:


	plt.subplot(2,1,1)
	schism_plotAtnodes(sin,sin.depths)
	plt.subplot(2,1,2)
	schism_plotAtnodes(sout,sout.depths)
	plt.tight_layout()
	plt.savefig('depth_change.png',dpi=300)
	#plt.show()
	plt.close()


	plt.subplot(2,1,1)
	if in_cartesian:
		schism_plotAtnodes(sin,sin.resolution_by_nodes.values())
	else:
		schism_plotelems(sin,sin.cpp_resolution_by_element.values())
	plt.subplot(2,1,2)
	if out_cartesian:
		schism_plotAtnodes(sout,sout.resolution_by_nodes.values())
	else:
		schism_plotelems(sout,np.asarray(sout.cpp_resolution_by_element.values()))
	#plt.tight_layout()
	plt.savefig('eu3.0_coarse_change.png',dpi=600)
	#plt.show()
	plt.close()
	
#quit()
####################################################


############## transfer boundaries ##########################
if do_bnd:
	os.chdir(dest_setup)

	sout=transfer_bnd(sin,sout)
if do_bnd or do_bath:
	sout.dump_hgridgr3('hgrid_bnd_transefered.gr3') # write new grid
	#sout=schism_setup(hgrid_file='hgrid_bnd_transefered.gr3',ll_file='hgrid.ll',vgrid_file='vgrid.in') # reload for testing
	print('done transeferring boundaries')
	os.system("mv hgrid.gr3 hgrid.gr30")
	os.system("ln -s hgrid_bnd_transefered.gr3 hgrid.gr3")
##################################################




##### copy vgrid ##############################################################################
# sin .vgrid i
# masked_array(data = [-1.       -0.984539 -0.972944 -0.961348 -0.953618 -0.945888 -0.938157
# -0.933133 -0.928881 -0.925402 -0.922697 -0.920764 -0.919605 -0.918832
# -0.918252 -0.917865 -0.917479 -0.917286 -0.904648 -0.890678 -0.875257
# -0.858265 -0.839589 -0.819129 -0.796802 -0.772547 -0.746337 -0.718181
# -0.688134 -0.656298 -0.622828 -0.587933 -0.551869 -0.514936 -0.477468
# -0.439818 -0.402349 -0.365416 -0.329352 -0.294457 -0.260988 -0.229152
# -0.199104 -0.170948 -0.144738 -0.120484 -0.098156 -0.077696 -0.059021
# -0.042029 -0.026607 -0.012638  0.      ],
#             mask = False,
#       fill_value = -9999.0)


## quick next neighbour transfer
#  maybe  check 3 closest neighbours adn chose vertical layers of that with most similar depth
# for now just take from next neighbour - one should check how much depth deviates

# copy next neighbour vgrid profile from old setup dictionary
if do_vgrid:
	print('transfering vgrid')
	sout.znum=sin.znum
	sout.vgrid={}
	sout.bidx={} # counter of masked arrays
	for key in sout.xdict.keys():
		dists,inn=interpolant.xy_nn_tree.query((sout.xdict[key],sout.ydict[key]),k=1)
		sout.vgrid.update({key:sin.vgrid[inn+1]})
		sout.bidx.update({key:sout.znum-sin.vgrid[inn+1].count()+1})

	# write vgrid
	f=open('vgrid_out.in','w')
	f.write('%8i\n'%(1))
	f.write('%8i\n'%sout.znum)
	for key in sout.xdict.keys():
		f.write('%8i %8i'%(key,sout.bidx[key]))
		for val in sout.vgrid[key][sout.bidx[key]-1:]:
			f.write('%12.6f'%val) 
		f.write('\n') 
	f.close()
	os.system("ln -s vgrid_out.in vgrid.in")


	nlayers_old=[sin.vgrid[key].count() for key in sin.xdict.keys()]
	nlayers_new=[sout.vgrid[key].count() for key in sout.xdict.keys()]


	plt.subplot(2,1,1)
	ph,ch=sin.plotAtnodes(nlayers_old)
	plt.title('old')
	ch.set_label('# vert layers')
	plt.subplot(2,1,2)
	ph2,ch2=sout.plotAtnodes(nlayers_new)
	plt.title('coarse')
	ch2.set_label('# vert layers')
	plt.tight_layout()
	plt.savefig('vert_layers_transefered',dpi=300)
	#plt.show()

############################################## Done Vgrid #################################################################






################################### copy gr3 files  ##############################################################################
class gr3file():  # for tris - add quads
	def __init__(self,hgrid_file='hgrid.gr3'):

		self.hgrid_file=hgrid_file

		# parse hgrid file
		f = open(hgrid_file)
		self.description = f.readline().rstrip()
		dat = f.readline().split()
		self.nelements = int(dat[0])
		self.nnodes = int(dat[1])

		n=[]
		x=[]
		y=[]
		d=[]
		self.depthsdict={}
		self.xdict={}
		self.ydict={}
		for nn in range(self.nnodes):
			dat=f.readline().split()
			n.append(int(dat[0]))
			x.append(float(dat[1]))
			y.append(float(dat[2]))
			d.append(float(dat[3]))
			self.depthsdict[n[-1]] = d[-1]
			self.ydict[n[-1]] = y[-1]
			self.xdict[n[-1]] = x[-1]
		self.inodes = n
		self.x = x
		self.y = y
		self.depths = d

		# gr3 file might lack element table
		n=[]
		nv = []
		nvdict = {}
		for nn in range(self.nelements):
			dat=f.readline().split()
			if len(dat)>0:
				n.append(int(dat[0]))
				nvnum = int(dat[1])
				nv.append([ int(ii) for ii in dat[2:2+nvnum]])
				nvdict[n[-1]] = nv[-1]
		self.ielement = n
		self.nv = nv
		self.nvdict = nvdict

if do_gr3:
	gr3files=glob(source_setup+'*gr3')
	gr3files=np.asarray(gr3files)
	delinds=[]
	for i in range(len(gr3files)):
		if 'hgrid' in gr3files[i]:
			delinds.append(i)
	
	gr3files=np.delete(gr3files,np.asarray(delinds,int))

	for filei in gr3files:
		print('interpolating file'+filei)
		try:
			gr3=gr3file(filei)
			nodevalues=interpolant.interp_bary(gr3.depths) # barycentric interpolation
			nodevalues[interpolant.no_parents]=interpolant.interp_nn(np.asarray(gr3.depths),interpolant.no_parents,tol_dist=643) # next neighbour interpolation for
			name=filei[filei.rfind('/')+1:]
			sout.dump_gr3_spat_var(name,nodevalues)
		except:
			print('error with file '+filei)
			pass

# interpolate .ic files			
if do_ic:
	gr3files=glob(source_setup+'*.ic')
	gr3files=np.asarray(gr3files)
	delinds=[]
	
	gr3files=np.delete(gr3files,np.asarray(delinds,int))

	for filei in gr3files:
		print('interpolating file'+filei)
		try:
			gr3=gr3file(filei)
			nodevalues=interpolant.interp_bary(gr3.depths) # barycentric interpolation
			nodevalues[interpolant.no_parents]=interpolant.interp_nn(np.asarray(gr3.depths),interpolant.no_parents,tol_dist=643) # next neighbour interpolation for
			name=filei[filei.rfind('/')+1:]
			sout.dump_gr3_spat_var(name,nodevalues)
		except:
			print('error with file '+filei)
			pass			
			
sout.dump_tvd_prop() 
##################################################









#copy hostart ####################################
if do_hot:

	class schismhot():

	  def __init__(self,file,hotsetup):
	    tnc = netCDF4.Dataset(file)
	    tv = tnc.variables
	    #snc = netCDF4.Dataset(file)
	    #sv = snc.variables
	 

	    self.lon = np.asarray(hotsetup.lon)
	    self.lat = np.asarray(hotsetup.lat)
	    self.vgrid = hotsetup.vgrid
	    self.d = hotsetup.depths
	    self.time = tv['time'][:]
	    self.tidx = 0
	    self.s = tv['tr_nd0'][:,:,1] # node vert tracer
	    self.t = tv['tr_nd0'][:,:,0] # node vert tracer
	  
	    tnc.close()	

	  def interpolate(self,depths,nodelon,nodelat,bidx=1):
	    # start
	    t = np.zeros((len(depths),))
	    s = np.zeros((len(depths),))

	    inn=np.argmin( (nodelon-self.lon)**2 + (nodelat-self.lat)**2) # horizontal next neighbour
	    d=self.vgrid[inn+1]*self.d[inn]

		
	    #vertical linear 
	    for ik,dep in enumerate(depths[bidx-1:]):
	      if np.sum(d<=dep):
	        ibelow=np.nonzero(d<=dep)[-1][-1] #  find(d<=dep)[-1]
	        iabove=ibelow+1 #- (ibelow+1==len(d))
	        dists=d[ibelow:iabove+1]-dep
	        if dists[0]==0:
	          w=np.float64(np.logical_not(dists))
	        else:
	          w=1/abs(dists)/sum(1/abs(dists))
	        s[ik]=np.sum(w*self.s[inn,ibelow:iabove+1])
	        t[ik]=np.sum(w*self.t[inn,ibelow:iabove+1])
	      else:
	        s[ik]=np.nan
	        t[ik]=np.nan


	    return (t,s)

	  # much faster using np	
	  def interpolate2(self,depths,nodelon,nodelat,bidx=1):
	    
		# start
	    t = np.zeros((len(depths),))
	    s = np.zeros((len(depths),))

	    inn=np.argmin( (nodelon-self.lon)**2 + (nodelat-self.lat)**2) # horizontal next neighbour
	    d=self.vgrid[inn+1]*self.d[inn]
		
	    bidx0=d.nonzero()[0][0]
	    s[bidx:]=np.interp(depths[bidx:],d[bidx0:],self.s[inn,bidx0:])
	    t[bidx:]=np.interp(depths[bidx:],d[bidx0:],self.t[inn,bidx0:])
	    s[:bidx],t[:bidx]=s[bidx],t[bidx]
	    #bidx=1
	    #vertical linear 
	    #for ik,dep in enumerate(depths[bidx-1:]):
		#  #ikeff=ik+bidx	
		#  if np.sum(d<=dep):
		#	ibelow=np.nonzero(d<=dep)[-1][-1] #  find(d<=dep)[-1]
		#	iabove=ibelow+1 #- (ibelow+1==len(d))
		#	dists=d[ibelow:iabove+1]-dep
		#	if dists[0]==0:
		#		w=np.float64(np.logical_not(dists))
		#	else:
		#		w=1/abs(dists)/sum(1/abs(dists))
		#	s[ik]=np.sum(w*self.s[inn,ibelow:iabove+1])
		#	t[ik]=np.sum(w*self.t[inn,ibelow:iabove+1])
		# # else:
		#	s[ik]=np.nan
		#	t[ik]=np.nan
	    return (t,s)


	#file='/work/gg0028/g260114/RUNS/Europe3.0/hotstarts/hotstart_cmems_jannsen_nbs_bs.nc'
	file=source_setup+'hotstart.nc'
	hotstart_in=schismhot(file,sin)
	
	
	#import time
	#s = {}
	#t = {}
	#t0=time.time()
	#for i in range(100):
	#	t[i],s[i] = hotstart_in.interpolate(depths,nodelon,nodelat,bidx=1)
	#t1=time.time()
	#dt1=t1-t0
	#print(dt1)
	
	#s2 = {}
	#t2 = {}
	#t0=time.time()
	#for i in range(100):
	#	t2[i],s2[i] = hotstart_in.interpolate2(depths,nodelon,nodelat,bidx=bidx-1)
	#t1=time.time()
	#dt2=t1-t0
	#print(dt2)

	
	# interpolate hotstart fields  - horizontal next neighbour vertical linear
	# write t,s on nodes
	s = {}
	t = {}
	# create t,s fields:
	for i,nodelon,nodelat,d in zip(sout.inodes,sout.lon,sout.lat,sout.depths):
		if (i%10000) == 0:
			print('  interpolate i = %d'%i)
	
		depths = sout.vgrid[i].filled(-1)*d
		bidx =sout.bidx[i]

		# interpolate from large domain
		#t[i],s[i] = hotstart_in.interpolate(depths,nodelon,nodelat,bidx=1)
		t[i],s[i] = hotstart_in.interpolate2(depths,nodelon,nodelat,bidx=bidx-1)

	#write pickle
	f = open('ts.pickle','wb')
	pickle.dump((t,s),f)
	f.close()

	#t2,s2=pickle.load(open('ts.pickle','rb'))
	#t2=np.asarray(t2.values())
	#s2=np.asarray(s2.values())



	t=np.asarray(list(t.values()))
	s=np.asarray(list(s.values()))

	# not necessary with new interpolation ?
	# replace nan with value above
	#for i in range(len(s)):
	#	naninds=np.isnan(t[i])
	#	ibtm=np.sum(naninds)
	#	t[i][naninds]=t[i][ibtm]
	#	s[i][naninds]=s[i][ibtm]

	# write hotstart
	tr_nd=np.stack((t,s),axis=2) 
	sout.write_hotstart(tr_nd,filename=hot_outname,time=0.0,iths=0,ifile=0,elev=0.0)




	# tr_nd.shape has to be (nelements,nvrt,ntracers)
	# interpolate from node to elements via avaraging
	#tris=np.asarray(sin.nvdict.values())-1
	#t_elem=np.mean(hotstart_in.t[tris],axis=1)
	#s_elem=np.mean(hotstart_in.s[tris],axis=1)
	#tr_nd=np.stack((t_elem,s_elem),axis=2) 

	#tr_nd=np.stack((hotstart_in.t,hotstart_in.s),axis=2) 
	#sin.write_hotstart(tr_nd,filename='hotstart_test.nc',time=0.0,iths=0,ifile=0,elev=0.0,uvel=0.0,vvel=0.0)

	#sin.write_hotstart(tr_nd,filename='hotstart_test.nc',time=0.0,iths=0,ifile=0,elev=0.0)


##########################################################


# source sink elments
if do_vm_source:
	os.chdir(dest_setup)
	sout=schism_setup()

	fin=open(source_setup+'source_sink.in','r')
	fout=open('source_sink.in','w')
	line=fin.readline()
	fout.write(line)
	
	# copy sources
	nelements=int(line.split()[0])
	for i in range(nelements):
		line=fin.readline()
		elem_in=int(line.split()[0])

		xs=[]
		ys=[]
		for inode in sin.nvdict[elem_in]:
			#xs.append(sin.x[inode-1])
			#ys.append(sin.y[inode-1])
			xs.append(sin.lon[inode-1]) #force lon lat
			ys.append(sin.lat[inode-1])

			
		elem_in_check=sin.find_nearest_element(np.mean(xs),np.mean(ys),latlon=True)
		elem_out=sout.find_nearest_element(np.mean(xs),np.mean(ys),latlon=True)
		fout.write('%i ! global element # for source %i \n'%(elem_out,i+1))

	while not 'sinks' in line:
		line=fin.readline()
		fout.write(line)
	
	# copy sinks
	nelements=int(line.split()[0])
	for i in range(nelements):
		line=fin.readline()
		elem_in=int(line.split()[0])

		xs=[]
		ys=[]
		for inode in sin.nvdict[elem_in]:
			xs.append(sin.x[inode-1])
			ys.append(sin.y[inode-1])

		elem_in_check=sin.find_nearest_element(np.mean(xs),np.mean(ys),latlon=True)
		elem_out=sout.find_nearest_element(np.mean(xs),np.mean(ys),latlon=True)
		fout.write('%i ! global element # for sink %i \n'%(elem_out,i+1))

	# add bottom
	while line!='':
		line=fin.readline()
		fout.write(line)

	fin.close()
	fout.close()



# check
#	plt.triplot(sout.x,sout.y,triangles=np.asarray(sout.nvdict.values())-1,color='k')
#	plt.plot(np.mean(xs),np.mean(ys),'rx')
#	plt.plot([sin.x[i-1] for i in sin.nvdict[elem_in] ],[sin.y[i-1] for i in sin.nvdict[elem_in] ],color='b')
#	plt.plot([sout.x[i-1] for i in sout.nvdict[elem_out] ],[sout.y[i-1] for i in sout.nvdict[elem_out] ],color='r')




