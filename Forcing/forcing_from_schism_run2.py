import sys
import os
""" Code to expand schism downscaling for child model boundaries that liew within dry follaing areay of mother model
so that the schism interpolate_vraibale utilities where used with limitted boudanry.
In mother dry areay child model takes over vertical mean T and S values of nn profiles in mothder model both if wet and and dry
the simulated elevation if wet and sea level 0 if dry. For uv the vertical mean of nn profile if wet and zero if dry
 """
from time import time
sys.path.insert(0,'/home/g/g260114/git/schism-hereon-utilities/')
from schism import* # import schism class to read grid structure
import glob
# gen forcing using nudign scripts for reduced boundary
# overwrite neareas neighbour for missing boundaries


for scenario in ['r219','r217']:
	print(scenario)
	#ncdir='/work/bg1186/g260094/SNS/SNSE3D_01a_MSLRnsob_prc50_r212/'
	ncdir='/work/bg1186/g260094/SNS/SNSE3D_01a_MSLRnsob_prc50_'+scenario+'/'

	max_stack=367

	# schout style:
	schismfiles=[] 
	for iorder in range(8): # check for schout_nc files until 99999
		schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
	nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
	schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
	nrs=list(np.asarray(nrs)[np.argsort(nrs)])


	s0dir='/work/gg0028/g260094/SNS/SNSE3D_MSLRnsob_prc50_r212/'
	s1dir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/ClimateProj/'




	os.chdir(s0dir)
	s0=schism_setup()#vgrid_file='vgrid.in.old')
	os.chdir(s1dir)
	s1=schism_setup()#vgrid_file='vgrid.in.old')#
	os.chdir('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/ClimateProj/forcingFromSNS/')
	#s0=schism_setup(hgrid_file='bg.gr3',vgrid_file='vigrid.bg')

	#s1=schism_setup(hgrid_file='fg.gr3',vgrid_file='vigrid.fg')
	#os.chdir('/work/gg002F8/g260114/RUNS/GermanBight/GB_2017_wave_sed/ClimateProj/forcingFromSNS/')
	nvertmax=21 #max number of vertical layers in destination grid

	#acc=schism_output2('/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/ClimateProj/forcingFromSNS/',max_stack=max_stack).nc # to much space
	acc=schism_output2(ncdir,max_stack=max_stack).nc

	############## load schism info ###############
	s=s1
	lon=np.asarray(s.lon)
	lat=np.asarray(s.lat)
	# detect at which boundaries forcing has to be applied
	with open('bctides.in') as f:
		openbd_segs=[]
		bd=0
		for line in f.readlines()[4:]:
			print(line)
			splitted=line.split()
			bd+=np.sum([val.isdigit() for val in splitted[:5]])==5 #is openboundary 
			if (splitted[1:5]==['4', '4', '4', '4']) | (splitted[1:5]==['5', '5', '5', '5']):
				openbd_segs.append(bd-1)
	print('create forcing for open boundaries '+str(openbd_segs))	
	ibd=np.hstack([np.asarray(s.bdy_segments[bdseg])-1 for bdseg in openbd_segs])		

	d=-np.asarray(s.depths)[ibd] # ant positive (here negative later multiplired with sigma)) values to allow mor effective interpolation
	bdlon=lon[ibd]
	bdlat=lat[ibd]
	depths=np.asarray(s.depths)
	bdcoords=np.asarray(list(zip(bdlon,bdlat)))

	#check if same projection before
	bdx=np.asarray(s1.x)[ibd]
	bdy=np.asarray(s1.y)[ibd]
	bdcoords=np.asarray(list(zip(bdx,bdy)))
	nbd=len(bdcoords)
	#######################################


	#appendbds=np.asarray([ibd-1 for ibd in s1org.bdy_segments[0] if ibd not in s1.bdy_segments[0]])
	#appendbds_indx=[ibd-1 for ibd in s1org.bdy_segments[0] if ibd not in s1.bdy_segments[0]]

	#xbd=np.asarray(s1.x)[appendbds]
	#ybd=np.asarray(s1.y)[appendbds]
	#bdcoords=list(zip(xbd,ybd))

	#insert_index=[np.where(s1org.bdy_segments[0]==ibd+1)[0][0] for ibd in appendbds]


	# find parent elements in mother grid
	#parents,ndeweights=s0.find_parent_tri(bdlon,bdlat,dThresh=0.4,latlon=True)
	parents,ndeweights=s0.find_parent_tri(bdx,bdy,dThresh=5000,latlon=False)

	plt.ion()

	# find nearest element centers for dry assessment
	# element centers
	x0,y0=np.asarray(s0.x),np.asarray(s0.y)
	cx=np.mean(x0[s0.nvplt],axis=1)
	cy=np.mean(y0[s0.nvplt],axis=1)
	elcoords=[[cx[i],cy[i]] for i in range(len(cx))] # pooint pairs of nodes
	elem_nn_tree = cKDTree(elcoords) # next neighbour search tree	     
	gets_dry=acc.wetdry_elem.max(axis=0).values
	gets_dry_ind=np.where(gets_dry)[0]

	# maniupluate coordinates to avoid selection of dry replacement nn elements
	cx_mnpl,cy_mnpl  = cx.copy(),cy.copy()
	cx_mnpl[gets_dry_ind]=np.inf
	elcoords_wet=[[cx_mnpl[i],cy_mnpl[i]] for i in range(len(cx))] # pooint pairs of nodes
	elem_nn_tree_wet = cKDTree(elcoords_wet) # next neighbour search tree	    
	#cy_mnpl

	force_elems=np.zeros(s0.nelements)
	force_elems[parents]=1
	dry_force_elem=(gets_dry>0) & (force_elems > 0)

	ireplace=np.where(dry_force_elem)[0]
	ireplacement=elem_nn_tree_wet.query(list(zip(cx[ireplace],cy[ireplace])))[1]

	elems=s0.nvplt

	## re check this part
	#plt.figure()
	##s0.plotAtelems(dry_force_elem,latlon=False)
	#ph,ch=s0.plotAtelems(gets_dry,latlon=False)
	#ch.set_label('Falls Dry')
	#s0.plot_domain_boundaries(latlon=False,append=True)
	#plt.triplot(s0.x,s0.y,elems[ireplace],color='r')
	#plt.triplot(s0.x,s0.y,elems[ireplacement],color='g')
	#s0.plot_mesh(latlon=False)
	#plt.triplot(s0.x,s0.y,elems[ireplace],color='r')
	#plt.triplot(s0.x,s0.y,elems[ireplacement],color='g')


	# replace  element with closest non dry element pseudo interpoaltion replacement
	parents_wet=parents.copy()
	for ind_replace, ind_replacment in zip(ireplace,ireplacement):
		ind_replace, ind_replacment
		parents_wet[ind_replace==ind_replace]=ind_replacment


	nt=acc.elev.shape[0]
	t=((acc.time-acc.time[0])/np.timedelta64(1,'s')).values

	# output arrays
	nz=21
	elev=np.zeros((nt,nbd,1,1))
	T=np.zeros((nt,nbd,nz,1))
	S=np.zeros((nt,nbd,nz,1))
	uv=np.zeros((nt,nbd,nz,2))


	if len(openbd_segs)>0:
		frcbdnodes=[]
		for seg in openbd_segs:
			frcbdnodes+=s.bdy_segments[seg]
			bdyvgrid = np.asarray([s.vgrid[ii].filled(-1.) for ii in frcbdnodes ])
	else:
		frcbdnodes=s.bdy_nodes


	stack_size=acc.elev.chunksizes['time'][0]
	nstacks=int(acc.time.shape[0]/stack_size)
		
		

	method='nn' # 'along_sigma'
		
	if method=='nn':	
		t0=time()
		# nearest
		s0.init_node_tree(latlon=False)
		nnpoints=s0.node_tree_xy.query(bdcoords)[1]
		nn=np.asarray(s0.node_tree_xy.query(list(zip(bdx,bdy)))[1])
		acc_sel=acc.sel(nSCHISM_hgrid_node=nn)


		t_total=0
		dt=0
		for stack in range(nstacks):
			print('doing stack ' + str(stack)+'/'+str(nstacks))
			t0=time()
			#load stacks
			istart=stack*stack_size
			iend=(stack+1)*stack_size
			elevin=acc.elev[istart:iend].values
			saltin=acc.salt[istart:iend].values
			tempin=acc.temp[istart:iend].values
			uvin=acc.hvel[istart:iend].values

			elev[istart:iend,:,0,0]=acc_sel.elev[istart:iend,:].values
			T[istart:iend,:,:,0]=acc_sel.temp[istart:iend,:].values
			S[istart:iend,:,:,0]=acc_sel.salt[istart:iend,:].values
			uv[istart:iend,:]=acc_sel.hvel[istart:iend,:].values
			# time estimate	
			dt+=time()-t0
			dtavg=dt/(stack+1)
			print('remaining time estimate ' + str(round((nstacks-(stack+1))*dtavg/60,2)) + ' min' )
		
		# maybe if parallel otherise loop
		#elev[:,:,0,0]=acc_sel.elev.values
		#T[:,:,:,0]=acc_sel.temp.values
		#S[:,:,:,0]=acc_sel.salt.values
		#uv=acc_sel.hvel.values
		#print('nn interp took ' +  str((time()-t0)/60) + ' min')
	elif method=='along_sigma':
		
		nz_in=21
		ndeweights_depth=np.tile(ndeweights,(nz_in,1,1)).swapaxes(0,1).swapaxes(1,2)
		ndeweights_hvel=np.tile(ndeweights_depth,(2,1,1,1)).swapaxes(0,1).swapaxes(1,2).swapaxes(2,3)


		# alternative nearset node


		# interpolate from SCHISM along sigma layers
		
		t_total=0
		dt=0
		for stack in range(2):#(nstacks):

			t0=time()
			#load stacks
			istart=stack*stack_size
			iend=(stack+1)*stack_size
			elevin=acc.elev[istart:iend].values
			saltin=acc.salt[istart:iend].values
			tempin=acc.temp[istart:iend].values
			uvin=acc.hvel[istart:iend].values

			# stupid along sigma layers
			print('doing stack ' + str(stack))
			for ti in range(24):
				elev[t_total,:,0,0]=(elevin[ti][nodes]*ndeweights).sum(axis=1)
				S[t_total,:,:,0]=(saltin[ti][nodes,:]*ndeweights_depth).sum(axis=1)
				T[t_total,:,:,0]=(tempin[ti][nodes,:]*ndeweights_depth).sum(axis=1)
				uv[t_total,:]=(uvin[ti][nodes,:]*ndeweights_hvel).sum(axis=1)
				t_total+=1	
				
			# time estimate	
			dt+=time()-t0
			dtavg=dt/(stack+1)
			print('remaining time estimate ' + str(round((nstacks-(stack+1))*dtavg/60,2)) + ' min' )

			
	print('eliminate nans')
	inan_profiles=np.where(np.isnan(S[0,:,:,:].sum(axis=1)))[0] #where are nans in porfile
	inan=[ np.where(np.isnan(S[0,inode,:,0]))[0] for inode in inan_profiles] # where in vertical
	i1stvalid=[column[-1]+1 for column in inan]																		
	for inode,ireplace,ireplacement in zip(inan_profiles,inan,i1stvalid):
		for isub in range(len(ireplace)):
			isub 
			S[:,inode,ireplace[isub],:]=S[:,inode,ireplacement,:]
			T[:,inode,ireplace[isub],:]=T[:,inode,ireplacement,:]
			uv[:,inode,ireplace[isub],:]=uv[:,inode,ireplacement,:]
	print('done eliminating nans')



	ndays=int(round(t[-1]/86400,1))
	filename='elev_{:s}_{:d}_days_nn.th.nc'.format(scenario,ndays)
	s1.write_bdy_netcdf(filename,t,elev,frcbdnodes)
	filename='TEM_{:s}_3D_{:d}_days_nn_th.nc'.format(scenario,ndays)
	s1.write_bdy_netcdf(filename,t,T,frcbdnodes)
	filename='SAL_{:s}_3D_{:d}_days_nn.th.nc'.format(scenario,ndays)
	s1.write_bdy_netcdf(filename,t,S,frcbdnodes)
	filename='uv3D_{:s}_{:d}_days_nn.th.nc'.format(scenario,ndays)
	s1.write_bdy_netcdf(filename,t,uv,frcbdnodes)

	#write


	# improve interp
	#bddepths0=np.asarray([s0.vgrid[inode].filled(-1)*s0.depthsdict[inode] for inode in nn+1])
	#T0=acc.temp[0,:].values
	bdprofiles2={'uv':uv,'salt':S,'temp':T,'elev':elev}


	name_ssh='elev'
	name_salt='salt'
	name_temp='temp'

	plot_bd=True
	if plot_bd:
		plt.figure()
		names={'ssh':name_ssh,'salt':name_salt,'temp':name_temp}
		for varname in 'ssh','salt','temp':
			name=names[varname]
			ti=np.int(nt/2)
			#t=t0+ti*deltaT
			#tt=t.strftime("%Y-%m-%dT%H:%M:%S")
			datain=acc[name][ti].load().values
			plt.clf()
			profil=bdprofiles2[name][ti,:,:]
			if len(datain.shape)==2:
				datain=datain[:,-1]
				profil=profil[:,-1,:]
			vmin=np.nanmin(datain)
			vmax=np.nanmax(datain)
			#vmin=bdprofiles2['ssh'][ti,:,:].min()
			#vmax=bdprofiles2['ssh'][ti,:,:].max()
			#plt.pcolormesh(x2d,y2d,datain,vmin=vmin,vmax=vmax)
			
			#plt.scatter(bdx,bdy,s=15,c=profil,vmin=vmin,vmax=vmax)#,edgecolors=None,linewidth=0.01)
			ph,ch,ax=s0.plotAtnodes(datain,latlon=False)
			ph1=plt.scatter(bdx,bdy,s=8,c=profil,vmin=vmin,vmax=vmax,edgecolors='k',linewidth=0.1)
			#ch=plt.colorbar()
			ch.set_label(varname)
			plt.suptitle(str(acc.time[ti].values))
			plt.title('surface '+varname)
			plt.legend([ph1],['interpolated bd points'])
			plt.tight_layout()
			plt.savefig(scenario+'halftime_surface_'+varname+'.png',dpi=300)

		# control plots for
		istart=0
		iend=0
		#plt.figure()
		plt.clf()
		for i in openbd_segs: #range(len(s.bdy_segments)):
			iend+=len(s.bdy_segments[i])
			inds=np.asarray(range(istart,iend))
			istart=iend
			bdnodes=np.asarray(frcbdnodes)[inds]#-1	

			bddepths=np.asarray([s.vgrid[inode].filled(-1)*s.depthsdict[inode] for inode in bdnodes])
			xcoord=np.tile(inds,[len(bddepths[0]),1]).T

			label=['start','end']
			for ti in range(0,-2,-1):
				plt.clf()
				plt.subplot(2,2,1)
				plt.pcolor( xcoord,bddepths,bdprofiles2['uv'][ti,inds,:,0])
				ch=plt.colorbar()
				plt.xlabel('bdy node')
				plt.ylabel('bddepths')
				ch.set_label('u [m/s]')
			
				plt.subplot(2,2,2)
				plt.pcolor( xcoord,bddepths,bdprofiles2['uv'][ti,inds,:,1])
				ch=plt.colorbar()
				plt.xlabel('bdy node')
				plt.ylabel('bddepths')
				ch.set_label('v [m/s]')

				plt.subplot(2,2,3)
				plt.pcolor( xcoord,bddepths,bdprofiles2['temp'][ti,inds,:,0])
				ch=plt.colorbar()
				plt.xlabel('bdy node')
				plt.ylabel('bddepths')
				ch.set_label('t [deg C]')
				
				
				plt.subplot(2,2,4)
				plt.pcolor( xcoord,bddepths,bdprofiles2['salt'][ti,inds,:,0])
				ch=plt.colorbar()
				plt.xlabel('bdy node')
				plt.ylabel('bddepths')
				ch.set_label('s [g/kg]')
				#plt.suptitle(str(t0+timesq[ti]/timesq[1]*deltaT))	
				#plt.suptitle(str(t0+timesq[ti]/timesq[1]*deltaT))	
				plt.suptitle(str(acc.time[ti].values))
				plt.tight_layout()
				
				plt.savefig(scenario+'bd_%i_%s'%(i,label[ti]),dpi=400)
				#plt.show()
				plt.close()
