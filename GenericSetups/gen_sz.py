import numpy as np
i0=25
z0=-30
z1=-60
dz=2

D=np.arange(-30,z1+dz,dz)
nr=np.arange(i0,i0+len(D))

np.savetxt('sz_section.ascii',np.vstack((nr,D)).T,fmt='%d %d')

S=np.linspace(-1,0,8)
nr2=np.arange(1,len(S)+1)

np.savetxt('s_section.ascii',np.vstack((nr2,S)).T,fmt='%d %f')
