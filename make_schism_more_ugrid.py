
from glob import glob 
import xarray as xr
import os 
import numpy as np
import sys
import datetime as dt


## what to attach for ugrid 

infiles=['2_combined_variables/schism_raw_1.nc',]  # input files to convert
fillvalue=-9999   # fill value to overwrite in data
cf_role=xr.DataArray(int(2),name='schism_mesh',attrs={'cf_role':"mesh_topology",'long_name' : "Topology data of Mesh2","topology_dimension":"2","node_coordinates":"SCHISM_hgrid_node_x SCHISM_hgrid_node_y","face_node_connectivity":"SCHISM_hgrid_face_nodes"})
outdir='./3_ugrid_conversion/'

class hgrid(object):
    """ read a schsim grid file """
    def __init__(self,hgrid_file):
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
      
      n=[]
      nv = []
      nvdict = {}
      
      nvtridict = {}
      for nn in range(self.nelements):
        dat=f.readline().split()
        n.append(int(dat[0]))
        nvnum = int(dat[1])
        nv.append([ int(ii) for ii in dat[2:2+nvnum]])
        nvdict[n[-1]] = nv[-1]
      self.ielement = n
      self.nv = nv
      self.nvdict = nvdict
      
      nvtridict={}
      k=0
      for elem in nvdict.values():
          k+=1
          if len(elem)==4:
              nvtridict[k]=elem[0:3]
              k+=1
              nvtridict[k]=[elem[i] for i in [0,1,3]]
          else:
              nvtridict[k]=elem
      self.nvtridict = nvtridict
      
      # compute sides
      # first get sequence array
      self.nx = {}
      self.nx[3] = [[1,2],[2,0],[0,1]]
      self.nx[4] = [[1,2],[2,3],[3,0],[0,1]]

      # get inverse nv - neighbouring elements
      self.node_neighbour_elements = { i:[] for i in self.inodes }
      self.n_node_neighbours = { i:0 for i in self.inodes }
      for i,nv in zip(self.ielement,self.nv):
        for nnode in nv:
          self.n_node_neighbours[nnode] += 1
          self.node_neighbour_elements[nnode].append(i)
      # find neighbouring elements around element (ic3)
      self.element_sides={}
      self.side_nodes={}
      for i,nv in zip(self.ielement,self.nv):
        isides = []
        # loop around element for existing sides
        for iloop,(ind1,ind2) in enumerate(self.nx[len(nv)]):
          iside = 0
          nd1,nd2 = nv[ind1],nv[ind2]
          for checkelement in self.node_neighbour_elements[nd1]:
            if (checkelement != i) and (nd2 in self.nvdict[checkelement]):
              iside = checkelement
              break
          isides.append(iside)
        self.element_sides[i] = isides
      # count sides
      self.nsides = 0
      element_ids = self.element_sides.keys()
      #element_ids.sort()
      #element_ids=np.sort(element_ids)
      for i in element_ids:
        for ii,iside in enumerate(self.element_sides[i]):
          if iside==0 or i<iside:
            self.nsides += 1
            iinds = self.nx[len(self.element_sides[i])][ii]
            self.side_nodes[self.nsides] = [self.nvdict[i][iinds[0]],self.nvdict[i][iinds[1]]]


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
		
		

## Load longitude latitude coordinates from gridfile
## to overwrite the cartesian coordinates in the netcdf		

ll=hgrid('hgrid.ll')
lon=np.asarray(ll.x)
lat=np.asarray(ll.y)

# load start time from namelist to overwrite time format in netcdf (original it is just counting seconds sinc emodel start)
p=param()
reftime=dt.datetime(int(p.get_parameter('start_year')),
int(p.get_parameter('start_month')),
int(p.get_parameter('start_day')),
int(p.get_parameter('start_hour')),0,0)		



# overwrite FIll values:

# variables #to exckude for fill value overwerting
exclude=['time','SCHISM_hgrid', 'SCHISM_hgrid_face_nodes', 'SCHISM_hgrid_edge_nodes', 'SCHISM_hgrid_node_x',
         'SCHISM_hgrid_node_y', 'bottom_index_node', 'SCHISM_hgrid_face_x', 'SCHISM_hgrid_face_y', 
         'ele_bottom_index', 'SCHISM_hgrid_edge_x', 'SCHISM_hgrid_edge_y', 'edge_bottom_index',
         'sigma', 'dry_value_flag', 'coordinate_system_flag', 'minimum_depth', 'sigma_h_c', 'sigma_theta_b', 
         'sigma_theta_f', 'sigma_maxdepth', 'Cs', 'dryFlagElement', 'dryFlagSide', 'dryFlagNode','depth','zCoordinates'] # exclude for plot selection


## modify 
outnames=[]

for file in infiles:
	print(file)
	ds=xr.open_dataset(file)
	# change xcoord to lon lat
	ds['SCHISM_hgrid_node_x'].attrs["standard_name"] = "longitude"
	ds['SCHISM_hgrid_node_x'].attrs["units"] = "degrees_east"
	ds['SCHISM_hgrid_node_x'][:]=lon
	
	ds['SCHISM_hgrid_node_y'].attrs["standard_name"] = "latitude"
	ds['SCHISM_hgrid_node_y'].attrs["units"] = "degrees_north"
	ds['SCHISM_hgrid_node_y'][:]=lat
	
	# change time
	ds['time'].attrs['units']=reftime.strftime('seconds since %Y-%m-%d %H:%M:%S')
	
	try:
		outname= file[file.rindex('/')+1:]
	except:
		outname=file
	outname='conv_'+outname	
	
	#mask variables
	drynodes=np.asarray(ds.dryFlagNode.values,bool)
	for var in ds.variables:
		if var not in exclude:
			print(var)
			values=ds[var].values
			if len(values.shape)==2: #2d
				for ti in range(values.shape[0]):
					values[ti,drynodes[ti,:]]=fillvalue
			else:		
				for ti in range(values.shape[0]):
					values[ti,drynodes[ti,:],:]=fillvalue

			ds[var][:]=values
			ds[var].attrs['_FillValue']=fillvalue
	
	ds.to_netcdf(outdir+outname)
	cf_role.to_netcdf(outdir+outname,mode='a')
	ds.close()

	
# rename files withd ate in name	to date 
## rename
#outnames=[]
#for file in infiles:
#	ds=xr.open_dataset(file)
#	if 'SCHISM' in str(ds.dims.keys()):
#		t0=str(ds['time'][0].values)[:19] #[:19]  #[:10]: just date  [:19] date and time
#		t1=str(ds['time'][-1].values)[:19]#[:19]
#		outstring=outdir+'schism_'+t0+'_'+t1+'.nc'
#		outnames.append(outstring.replace('-','').replace(':',''))
#	ds.close()		 