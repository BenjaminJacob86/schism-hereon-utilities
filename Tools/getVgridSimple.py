np.max(np.log(np.asarray(s.depths))*10)
np.min(np.log(np.asarray(s.depths))*10)


n=np.asarray(np.round((np.sqrt(np.asarray(s.depths)))*3),int)
nmax=n.max()
vgrid={i:np.ma.masked_array(s.depths[i]*np.concatenate((-np.ones(nmax-n[i]),np.linspace(-1,0,n[i]))),mask=np.concatenate((np.zeros(nmax-n[i],bool),np.ones(n[i],bool)))) for i n range(1,1+s.nnodes)}

