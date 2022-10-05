  def parse_vgrid(self,vgrid_file='vgrid.in'):
      #import numpy as np
      f = open(vgrid_file)
      line=f.readline()
      if '!' in line:
        line=line.split('!')[0]	
      first = int(line)
      if first == 1: # unsructured vertical
        znum = int(f.readline())
        self.znum = znum
        a = {}
        self.bidx = {}
        for line in f.readlines():
          sigma1d = -9999.*np.ones((znum,))
          data = line.split()
          self.bidx[int(data[0])] = int(data[1])
          sigma1d[int(data[1])-1:] = np.asarray([float(ii) for ii in data[2:]])
          a[int(data[0])] = np.ma.masked_equal(sigma1d,-9999.)
        f.close()
        self.vgrid = a
      else: # sigma - z assume z
        line=	  
        self.znum,self.n_zlevels,self.h_s=f.readline().split(' ')[:3]
        self.znum,self.n_zlevels,self.h_s=np.int(self.znum),np.int(self.n_zlevels),np.int(self.h_s)
        f.readline()
		# zlevels
        nz=np.int(f.readline().split(' ')[0])
        f.close()
		
		self.zlevels=np.loadtxt(vgrid_file,skiprows=3,comments='!',max_rows=nz)[:,1]
		self.slevels=np.loadtxt(vgrid_file,skiprows=5+nz,comments='!')[:,1]
		self.vgrid={}
        if nz > 1: #s_z
          comb=np.hstack((zlevels,h_s*slevels))
          for i,d in enumerate(s.depths): #pure slevel when shallower?
            if d < self.h_s: # pure sigma ?		   
              combi=np.hstack((self.zlevels,s.depths[i]*self.slevels))
              self.vgrid[i+1]=np.ma.masked_array(combi,mask=bomb<-d)
			else:
              self.vgrid[i+1]=np.ma.masked_array(comb,mask=bomb<-d)
        else: # pure sigma
          for i,d in enumerate(s.depths):
            self.vgrid[i+1]=s.depths[i]*self.slevels