""" Donwload CMEMS data for a grid  """

#/opt/miniforge/envs/forcing/bin/ipython3

import sys
import os
import numpy as np

#boto list
import boto3
from botocore import UNSIGNED
import s3fs   # To access files on user dir directly using AWS
import datetime as dt
from pathlib import Path
import shutil

grid_selection = prefix = os.getenv('GRID_FOLDER') # FOlger in user my files to determine the grid
forcing_selection = os.getenv('FORCING_FOLDER') # FOlger in user my files to determine the grid
userdir=os.environ.get("EDITO_INFRA_OUTPUT")    # where to write the output

os.chdir('/app/setup/')

if userdir==None:
    local_test=True
    userdir='./outputs/'
else:
    local_test=False

if local_test:
    bucket = "oidc-jacobb"
    endpoint_url = 'https://minio.dive.edito.eu' # new one
    
    s3 = boto3.client("s3",endpoint_url = 'https://'+'minio.dive.edito.eu',
                    aws_access_key_id= 'V2WOPE9NZOJD9SFT0T7T', 
                    aws_secret_access_key= 'K9RemsfC39k3k+IBScXiPdBZ3rSjj44qtgg4xjf0', 
                    aws_session_token = 'eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9.eyJhY2Nlc3NLZXkiOiJWMldPUEU5TlpPSkQ5U0ZUMFQ3VCIsImFjciI6IjAiLCJhbGxvd2VkLW9yaWdpbnMiOlsiKiJdLCJhdWQiOlsibWluaW8iLCJhY2NvdW50Il0sImF1dGhfdGltZSI6MTc2NTgxODI4NSwiYXpwIjoib255eGlhLW1pbmlvIiwiZW1haWwiOiJiZW5qYW1pbi5qYWNvYkBoZXJlb24uZGUiLCJlbWFpbF92ZXJpZmllZCI6dHJ1ZSwiZXhwIjoxNzY1OTczMjMwLCJmYW1pbHlfbmFtZSI6IkphY29iIiwiZ2l2ZW5fbmFtZSI6IkJlbmphbWluIiwiZ3JvdXBzIjpbIkVESVRPX1VTRVIiLCJmb2NjdXMiLCJvbWkiXSwiaWF0IjoxNzY1ODg2ODMwLCJpc3MiOiJodHRwczovL2F1dGguZGl2ZS5lZGl0by5ldS9hdXRoL3JlYWxtcy9kYXRhbGFiIiwianRpIjoiYmVmYmRmZjYtODhhNi00NmFmLTk4MDUtNjEwYThlOGRhY2NjIiwibmFtZSI6IkJlbmphbWluIEphY29iIiwicG9saWN5Ijoic3Rzb25seSIsInByZWZlcnJlZF91c2VybmFtZSI6ImphY29iYiIsInJlYWxtX2FjY2VzcyI6eyJyb2xlcyI6WyJkZWZhdWx0LXJvbGVzLWRhdGFsYWIiLCJvZmZsaW5lX2FjY2VzcyIsInVtYV9hdXRob3JpemF0aW9uIl19LCJyZXNvdXJjZV9hY2Nlc3MiOnsiYWNjb3VudCI6eyJyb2xlcyI6WyJtYW5hZ2UtYWNjb3VudCIsIm1hbmFnZS1hY2NvdW50LWxpbmtzIiwidmlldy1wcm9maWxlIl19LCJtaW5pbyI6eyJyb2xlcyI6WyJzdHNvbmx5Il19fSwic2NvcGUiOiJvcGVuaWQgZW1haWwgcHJvZmlsZSIsInNlc3Npb25fc3RhdGUiOiI3MzFlZjc3ZS04MzY0LTQzYjYtOGJkNi01ODNlZGRjNDBlZjEiLCJzaWQiOiI3MzFlZjc3ZS04MzY0LTQzYjYtOGJkNi01ODNlZGRjNDBlZjEiLCJzdWIiOiI5ODY1MzkzOC01YTY5LTQyYzEtYWJkZi0wYjI5YWYzMWNhNWIiLCJ0eXAiOiJCZWFyZXIifQ.3eHIh2yl2agQQ1mkhyCuRqSyWLgbowv8-i1pVWlWpVL1ZKWlFKWLOGZ1ex4R-iAP-Mk77CXzNkSrwZP_PGnwSA')

    grid_selection = prefix = 'NWS_NEW'  #os.getenv('GRID_FOLDER') # FOlger in user my files to determine the grid
    forcing_selection = 'NWS_NEW/forcing_output/'
    os.makedirs('outputs', exist_ok=True)
else:

    grid_selection = prefix = os.getenv('GRID_FOLDER') # FOlger in user my files to determine the grid
    forcing_selection = prefix = os.getenv('FORCING_FOLDER')
    userdir=os.environ.get("EDITO_INFRA_OUTPUT")    # where to write the output

    endpoint_url = 'https://'+os.environ.get("AWS_S3_ENDPOINT" )  # get Url form env variable
    AWS_S3_ENDPOINT=os.environ.get("AWS_S3_ENDPOINT")  
    AWS_ACCESS_KEY_ID=os.environ.get("AWS_ACCESS_KEY_ID")  
    AWS_SECRET_ACCESS_KEY=os.environ.get("AWS_SECRET_ACCESS_KEY")  
    AWS_SESSION_TOKEN=os.environ.get("AWS_SESSION_TOKEN")  

    #  access user my Files directory
    userMyFilesdir=os.environ.get("User_Name").replace('user','oidc')  # get user name from environment variabl user-jacobb
    bucket=os.environ.get("User_Name").replace('user','oidc')

    # unvertified:
    s3 = boto3.client("s3",endpoint_url = endpoint_url,
                  aws_access_key_id= AWS_ACCESS_KEY_ID, 
                  aws_secret_access_key= AWS_SECRET_ACCESS_KEY, 
                  aws_session_token = AWS_SESSION_TOKEN)
                  

    # link the outputs to Edito Folder  -this is actually seemling super slow because it is an S3 bucket
    #from pathlib import Path
    os.makedirs('outputs', exist_ok=True)

# Helper function to download files preserving directory structure
def download_with_structure(s3_client, bucket, prefix, base_dir='.'):
    """
    Download all files from S3 prefix, preserving directory structure.
    
    Args:
        s3_client: boto3 S3 client
        bucket: S3 bucket name
        prefix: S3 prefix (folder path) to download from
        base_dir: Local base directory to download to (default: current directory)
    """
    # Normalize prefix: remove trailing slash for consistent matching
    prefix = prefix.rstrip('/')
    
    paginator = s3_client.get_paginator('list_objects_v2')
    pages = paginator.paginate(Bucket=bucket, Prefix=prefix)
    
    for page in pages:
        if 'Contents' not in page:
            continue
        for obj in page['Contents']:
            key = obj['Key']
            # Skip if it's a directory marker (ends with /)
            if key.endswith('/'):
                continue
            
            # Get relative path after the prefix
            if key.startswith(prefix):
                # Remove prefix and any leading slashes
                rel_path = key[len(prefix):].lstrip('/')
            else:
                # Fallback: use just the filename
                rel_path = key.split('/')[-1]
            
            # Construct local file path
            local_path = os.path.join(base_dir, rel_path)
            
            # Create directory if needed
            local_dir = os.path.dirname(local_path)
            if local_dir:
                os.makedirs(local_dir, exist_ok=True)
            
            # Download file
            s3_client.download_file(Bucket=bucket, Key=key, Filename=local_path)
            print(f"Downloaded: {key} -> {local_path}")

# Download setup files (grid files)
print(f"Downloading grid files from: {grid_selection}")
download_with_structure(s3, bucket, grid_selection)

# Download forcing files (including sflux subdirectory if present)
print(f"Downloading forcing files from: {forcing_selection}")
download_with_structure(s3, bucket, forcing_selection)

print(' '.join(('grid_selection : ',grid_selection)))
print(' '.join(('forcing_selection : ',forcing_selection)))


print(os.listdir())
print('----------------------------------------------------')    
print(os.environ)    
# make sure to link outputs to output folder on edito
#cp /SCHISM/exec.sh .

#with open(fname,'w') as f:
#    f.write('{:s}\n{:s}'.format(start_date,end_date))
   
# should be produced in foring gehnerator
with open('simulation_period.txt','r') as f:    
    start_date=f.readline().split()[0]
    end_date=f.readline().split()[0]
t00=dt.datetime.strptime(start_date,'%Y-%m-%d') #start time download 
t11=dt.datetime.strptime(end_date,'%Y-%m-%d') #start time download 
ndays=(t11-t00).days+1 # check schism forcing

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
        
        
p=param()

values=[int(value) for value in  start_date.split('-') ]+[ndays,]
p.set_parameter('start_year',values[0])
p.set_parameter(params=['start_year','start_month','start_day','rnday'],values=[int(value) for value in  start_date.split('-') ])        
p.set_parameter('rnday',ndays)

# Check if sflux folder exists (atmospheric forcing)
use_sflux = os.path.exists('sflux') and os.path.isdir('sflux')
if use_sflux:
    print('sflux available - enabling atmospheric forcing (nws=2)')
    p.set_parameter('nws', 2)  # default is 0
    # Export environment variable for shell scripts
    os.environ['USE_SFLUX'] = '1'
else:
    print('sflux not available - using default nws=0 (no atmospheric forcing)')
    os.environ['USE_SFLUX'] = '0'

# repair if tvd misses
if not os.path.isfile('tvd.prop'):
    with open('hgrid.gr3') as f:
        line=f.readline()
        line=f.readline()
    nelem=int(line.split()[0])
    tvd=np.zeros((nelem,2),int)
    tvd[:,0]=1+np.arange(nelem)
    tvd[:,1]=1
    np.savetxt('tvd.prop',tvd,fmt='%d %d')


# Source and destination paths
with open(os.path.join(userdir,'description.txt'),'w') as f:
    f.write('run_description:\n grid_dir: {:s}\n forcing_dir: {:s}\n start_date: {:s}\n end_date: {:s}'.format(grid_selection,forcing_selection,start_date,end_date))

# copy files to output to have all in one places for analysis
cp_files=['hgrid.gr3','hgrid.ll','vgrid.in']
for file in cp_files:
    shutil.copy(file, os.path.join(userdir,file))
    
# create execscript - select binary based on sflux availability
# The exec.sh script checks for sflux directory directly (more reliable than env var)
exec_script = """#!/bin/bash
# Select SCHISM binary based on atmospheric forcing availability
# Check if sflux directory exists in current directory
if [ -d "sflux" ] && [ "$(ls -A sflux 2>/dev/null)" ]; then
    echo "sflux directory found - Using SCHISM with atmospheric forcing (PREC_EVAP)"
    mpirun -n 8 /SCHISM/GBconfig/bin/pschism_PREC_EVAP_TVD-SB 5
else
    echo "sflux directory not found - Using SCHISM without atmospheric forcing"
    mpirun -n 8 /SCHISM/GBconfig/bin/pschism_TVD-SB 5
fi
"""
with open('exec.sh','w') as f:
    f.write(exec_script)
