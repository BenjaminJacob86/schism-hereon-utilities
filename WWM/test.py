###### BSH ########
class bsh_spec:
	""" Read BSH Wave spectra *.spec files with format
	Line 1:  Buoy-Type Station Datum+Time
	Line 2:  up to 12 Parameters
    01 Peakfrequency (Hz),
	02 Energy Maximum (m*m*s),
	03 Dir at Peak (degN), 90 degrees is from east
	04 Spread at Peak (deg),
	05 Hs (m),
	06 Tp (s),
	07 Tm1 (s),
	08 Tm2 (s),
	09 Tm-1 (s),
	10 Nspc ,
	11 Hmax (m),
	12 T(Hmax) (s)
	% Line 3 - Nspc+3: f (Hz), e (m*m*s), dir (degN), spr (deg)
	
	Input:
	filename of spectral files as str
	or chronologically ordered list of files of bouy data:
	(bsh_spec(file=list(np.sort(glob.glob('<path>/*<station>*.spec')))   )
	
	Output class with fields:
	.name := station name as specified in file e.g. FN3
	.B.type :=  Buoy-Type - wave rider etc.
	.dates := list of dates of measurements
	.integral_header := header of integral parameters sotred in respective according numpy array (.integral_values)
	.integral_values := time times header item numpy array of integra parameters e.g. 
						bouy.integral_values[:,bouy.integral_header.index('Hs')] containing HS time series
						according to bouy.dates	
	.spectal_header := ['f', 'e', 'dir', 'spr'] header of spectral parmeters within .spectra
	.spectra := list with spectra occoring to .dates each item being an numpy array
				with spectral freuqencies in rows anbd colums according to spectal_header
	"""
	# not always all frequencies in spectra spectra
	def __init__(self,file):
	
		if type(file) == list:
			files=file
			file=files[0]
			nfiles=len(files)
		else:
			nfiles=1
			
		with open(file) as f:
			names=['Peakfrequency','Energy_Maximum','Dir_at_Peak','Spread_at_Peak',
			'Hs','Tp','Tm1','Tm2','Tm-1','Nspc','Hmax','T_Hmax']
			integral_header=names
			spectal_header=[ 'f', 'e' , 'dir', 'spr']
			lines=f.readlines()
			date,type,name=[tmp for tmp in lines[0].rstrip('\n').split(' ') if tmp != '']
			dates=[dt.datetime.strptime(date,'%Y%m%d%H%M')]
			i_nspec=names.index('Nspc')
			vals=[ float(val) for val in lines[1].rstrip('\n').split()]
			integral_values=[vals]
			count=2
			nspec=int(vals[i_nspec]) 
			
			spectrum=[[ float(item) for item in lines[count].rstrip('\n').split()]]
			for i in range(count+1,count+nspec):
				spectrum.append([ float(item) for item in lines[i].rstrip('\n').split()])
			spectra=[np.asarray(spectrum)]
			
			count+=nspec
			
			append_data(count,lines)

		if nfiles > 1:
			for file in files[1:]:
				#print(file)
				count=0
				try:
					with open(file) as f:
						append_data(count,f.readlines())
				except:
					print('error loading file '+ file)
					pass
		
		# convert lists to numpy arrays	
		integral_values=np.asarray(integral_values)
		
		
	def append_data(self,count,lines):
		while count < len(lines):
			date=lines[count].rstrip('\n').split(' ')[-1]
			dates.append(dt.datetime.strptime(date,'%Y%m%d%H%M'))
			
			count+=1 # read integral parameters
			integral_values.append([ float(val) for val in lines[count].rstrip('\n').split()])
			nspec=int(integral_values[-1][i_nspec]) 
			
			count+=1 # read spectra
			spectrum=[[ float(item) for item in lines[count].rstrip('\n').split()]]
			for i in range(count+1,count+nspec):
				spectrum.append([ float(item) for item in lines[i].rstrip('\n').split()])
			spectra.append(np.asarray([spectrum]))
			count+=nspec
		
	def get_parameter(self,parameter='Hs'):
		""" get buoy data """
		return integral_values[:,integral_header.index(parameter)]