s=schism_setup()
maxDepth=np.zeros(s.nelements)
for ielement in s.ielement:
        maxDepth[ielement-1]=max([s.depths[idx-1] for idx in s.nv[ielement-1]])
dx=np.asarray(list(s.resolution_by_element.values()))
cfl=0.5
g=9.81


p=param()
dt=p.get_parameter('dt')
cfl=s.compute_cfl(dt=dt
s.plotAtelems(cfl)

dt=dx*cfl/np.sqrt(g*maxDepth)
plt.ion()
s.plotAtelems(dt)

cfl_by_elements=np.sqrt(g*maxDepth)*dt/dxs


