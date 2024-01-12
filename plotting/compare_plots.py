def compare_slabs(s,data,names=None,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None):
	if nrows*ncols>=len(data):
		fig, axs = plt.subplots(nrows=nrows,ncols=ncols,subplot_kw={'projection': ccrs.Mercator()},figsize=(11,8.5))
		axs = axs.flatten()
		alldata=np.hstack(data)		
		vmin,vmax=np.quantile(alldata,0.1),np.quantile(alldata,0.99)
		for i,datai in enumerate(data):
			ph,ch,ax=s.plotAtnodesGeo(datai,ax=axs[i])
			ch.remove()
			if names != None:
				ax.set_title(names[i])
			if axis_limit!=None:
				ax.set_extent(axis_limit)
			ph.set_clim((vmin,vmax))	
		# Add a colorbar axis at the bottom of the graph
		cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation='horizontal',extend='both',label=cblabel)
		plt.tight_layout()
	else:
		print('nrows*ncols < len(data)!')
	return cbar_ax
		
def compare_slabs_diff(s,data,names=None,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None):
	""" first entry of list data is reference substracete from experiments"""
	if nrows*ncols>=len(data)-1:
		fig, axs = plt.subplots(nrows=nrows,ncols=ncols,								subplot_kw={'projection': ccrs.Mercator()},figsize=(11,8.5))
		axs = axs.flatten()
		if nrows*ncols < 2:
			axs=[axs]
		for i,datai in enumerate(data[1:]):
			ph,ch,ax=s.plotAtnodesGeo(datai-data[0],ax=axs[i])
			ch.remove()
			if names != None:
				ax.set_title('-'.join((names[i+1],names[0])))
			if axis_limit!=None:
				ax.set_extent(axis_limit)
		# Add a colorbar axis at the bottom of the graph
		cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation='horizontal',extend='both',label=cblabel)
		plt.tight_layout()
	else:
		print('nrows*ncols < len(data)!')		
		
		
s.init_node_tree(latlon=True)
def get_quiver_locs(s,narrows=30,axis_limit=None):		
	if tuple(axis_limit)==None:
		axis_limit=(np.min(s.lon),np.max(s.lon),np.min(s.lat),np.max(s.lat))

	xlim=axis_limit[:2]
	ylim=axis_limit[2:]
	x=np.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/narrows)
	y=np.arange(ylim[0],ylim[1],(ylim[1]-ylim[0])/narrows)
	X, Y = np.meshgrid(x,y)
	d,qloc=s.node_tree_latlon.query((np.vstack([X.ravel(), Y.ravel()])).transpose()) #quiver locations				
	#xref,yref=np.asarray(plt.axis())[[1,2]] +  np.diff(plt.axis())[[0,2]]*[- 0.2, 0.1]
	#if (self.shape==(self.nt,self.nnodes)) or (self.shape==(self.nt,self.nnodes,self.nz)): 
	#	vmax=1.5#np.percentile(np.sqrt(u[qloc]**2+v[qloc]**2),0.95)
	#else:
	#	vmax=np.double(self.maxfield.get())
	return qloc
		
		
def compare_slabs_diff2(s,data,names=None,clim0='auto',cmap0=plt.cm.viridis,cmap=plt.cm.turbo,nrows=1,ncols=1,cblabel='',axis_limit=None,clim=None,figsize=(11,8.5),cb_orientation='horizontal',geo=True,relative=False):
	""" first entry of list data is reference substracete from experiments"""
	phs=[]
	qhs=[]
	def add_quiver(ivs,ax,x,y,u,v,qloc,axis_limit,scale=15,vmax=1,color='k'):
		axis_limit=np.asarray(axis_limit)
		xref,yref=axis_limit[[1,2]] +  np.diff(axis_limit)[[0,2]]*[-0.3, 0.85]
		if ivs:
			qh=ax.quiver(np.concatenate((x[qloc],[xref,])),np.concatenate((y[qloc],[yref,])),np.concatenate((u[qloc],[vmax,])),np.concatenate((v[qloc],[0,])),scale=scale,scale_units='inches',color=color)
			ax.text(xref,yref,'\n {:.2f} \n m/s'.format(vmax),fontsize='x-small')
		else:
			qh=[]	
		return qh
			
	if nrows*ncols>=len(data)-1:
		if geo:	
			import cartopy.crs as ccrs
			import cartopy.feature as cfeature
			from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

			fig, axs = plt.subplots(nrows=nrows,ncols=ncols,								subplot_kw={'projection': ccrs.Mercator()},figsize=figsize)
		else:
			fig, axs = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize)		
		axs = axs.flatten()
		
		if data[0].shape[0]==2:
			ivs=True #vector
			udata=[data[0][0,:]]
			vdata=[data[0][1,:]]
			for i,datai in enumerate(data[1:]):
				udata.append(data[i+1][0,:]-udata[0])
				vdata.append(data[i+1][1,:]-vdata[0])
			data=[(datai*2).sum(axis=0) for datai in data]	
			
			qloc=get_quiver_locs(s,narrows=40,axis_limit=axis_limit)
			proj=ccrs.Mercator()
			outproj=proj.transform_points(ccrs.Geodetic(),np.asarray(s.lon),np.asarray(s.lat))
			projx,projy=outproj[:,0],outproj[:,1]
			
			axis_limit2=np.asarray(axis_limit)
			outproj=proj.transform_points(ccrs.Geodetic(),axis_limit2[:2],axis_limit2[2:])
			projlim=tuple(np.hstack((outproj[:,0],outproj[:,1])))
		else:
			ivs=False
			projx,projy=np.asarray(s.lon),np.asarray(s.lat)

		if relative:
			diffs=[(datai-data[0])/data[0]*100 for datai in data[1:]]
		else:
			diffs=[datai-data[0] for datai in data[1:]]
		if clim=='auto': # limit to in regiond 5 and 95 percentile
			if axis_limit ==	None:
				inreg=np.arange(s.nnodes)
			else:
				lon,lat=np.asarray(s.lon),np.asarray(s.lat)
				inreg= (axis_limit[0] <= lon) & (lon <= axis_limit[1]) & (axis_limit[2] <= lat) & (lat <= axis_limit[3])
			reg_diffs=np.hstack([np.abs(diffsi[inreg]) for diffsi in diffs])
			vmin,vmax=np.nanquantile(reg_diffs,0.1),np.nanquantile(reg_diffs,0.99)		
			clim=(-vmax,vmax)
		if geo:	
			ph,ch,ax=s.plotAtnodesGeo(data[0],ax=axs[0],cmap=cmap0)
			if ivs:
				qh=add_quiver(ivs,ax,projx,projy,udata[0],vdata[0],qloc,projlim,scale=15,vmax=1,color='k')
				qhs.append(qh)		
		else:	
			ph,ch,ax=s.plotAtnodes(data[0],ax=axs[0],cmap=cmap0)
			qh=add_quiver(ivs,ax,np.asarray(s.lon),np.asarray(s.lat),udata[0],vdata[0],qloc,axis_limit,scale=1,vmax=1,color='w')
			qhs.append(qh)		
		if clim0=='auto': # limit to in regiond 5 and 95 percentile
			vmin0,vmax0=np.nanquantile(data[0],0.1),np.nanquantile(data[0],0.99)
			clim0=(-vmax0,vmax0)
		if clim != None:
				ph.set_clim(clim0)	
		phs.append(ph)		
		
		if names != None:
				ax.set_title(names[0])
		if axis_limit!=None:
				if geo:
					ax.set_extent(axis_limit)
				else:	
					ax.axis(axis_limit)
		if nrows*ncols < 2:
			axs=[axs]
		for i,datai in enumerate(data[1:]):
			if geo:	
				ph,ch,ax=s.plotAtnodesGeo(diffs[i],ax=axs[i+1],cmap=cmap)		
				if ivs:
					qh=add_quiver(ivs,ax,projx,projy,udata[i+1],vdata[i+1],qloc,projlim,scale=3,vmax=0.2,color='k')
					qhs.append(qh)	
			else:	
				ph,ch,ax=s.plotAtnodes(diffs[i],ax=axs[i+1],cmap=cmap)
				if ivs:
					qh=add_quiver(ivs,ax,np.asarray(s.lon),np.asarray(s.lat),udata[i],vdata[i],qloc,axis_limit,scale=1,vmax=0.5,color='w')
					qhs.append(qh)	
			ch.remove()
			if names != None:
				ax.set_title('-'.join((names[i+1],names[0])))
			if axis_limit!=None:
				if geo:
					ax.set_extent(axis_limit)
				else:	
					ax.axis(axis_limit)
			if clim != None:
				ph.set_clim(clim)
			phs.append(ph)	
		# Add a colorbar axis at the bottom of the graph
		cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
		
		# Draw the colorbar
		cbar=fig.colorbar(ph, cax=cbar_ax,orientation=cb_orientation,extend='both',label=cblabel)
		
		if geo==False:
			plt.tight_layout()
		
		for axi in axs[len(data):]:
			axi.remove()
		
		return fig,axs,cbar_ax,phs,qhs
	else:
		print('nrows*ncols < len(data)!')	