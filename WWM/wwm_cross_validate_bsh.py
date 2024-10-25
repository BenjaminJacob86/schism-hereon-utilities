""""
SCHISM Wave modeling

Perform intercomparison of bsh bouys against
station output from WW3 and WWM model runs.
Making plots of time series qqplot and taylordiagram
exportet as images and pdf.

"""

import os
from glob import glob
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
plt.ion()
import datetime as dt
import sys
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/Lib/')
# own libraries 
from data_and_model_classes import bsh_spec, WW3_mesh
from validation_statistics_and_plots import interp_to_data_time,QQplot, TaylorDiagram,plotTaylor

########### SETTINGS #########


skip_n_first_days=0 #skip first days of WWM where model is building up
#skip_n_first_days=np.timedelta64(skip_n_first_days,'D')
skip_n_first_days=dt.timedelta(days=skip_n_first_days)

bouydir='/gpfs/work/jacobb/data/Validation/Waves/bouy/bouys/' #bsh buoy directory containing subfolders with the bouys and in those a file with name <subfolder>_BSH.spec

# WWM run directories containing model station output
WWMdirs=['/gpfs/work/jacobb/data/Validation/Waves/REF_seagrass/']

#'/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_from_mistral/FRCbdCHeck/',
#2nd '/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_900s_onpiece/'
#2nd 'WWM_Bmax_1.6'
wwm_names=['REF']  # names of scenario for display
#,': WWM run with betamax inreased to 1.6'
wwm_descriptions=['']
comments='' #sentences following wwm and ww3 descriptions in pdf
display_names=wwm_names

# names of stations
keysOBS=['ELB','FNO','FN3','HEL']
keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland']
keysWWM=['prep_ELBE', 'prep_FINO1', 'prep_FINO3', 'prep_HELGON']

#keysOBS=['ELB','FNO','FN3','HEL']
#keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland']
#keysWWM=['prep_ELBE', 'prep_FINO1', 'prep_FINO3', 'prep_HELGON']

## WWW3 ###
addWW3=False
ww3dirs=['']
ww3mesh='NBSext_bl.msh'
ww3PointList='points.list'
#ww3GridOutput==np.sort(glob(ww3dir+'mix_select_prev_post_run/'+'ww3.??????.nc'))# np.sort(glob(ww3dir+'defect_ncs/'+'ww3.????01.nc')) # ww3.201701.nc
#ww3PointOutput=[ww3dirs[0]+'ww3.2017_cat_01_03.nc',ww3dirs[1]+'ww3.2017cat0103.spec.nc',ww3dirs[2]+'ww3.2017cat0103.spec.nc',ww3dirs[3]+'ww3.2017cat0103.spec.nc',ww3dirs[4]+'ww3.2017cat0103.spec.nc'] #'ww3.2017_spec.nc'  # ww3 spectral output
ww3_names=['']
ww3_descriptions=['']


pdfname="neu_WWM_Validation_2017.pdf" # name of image pdf
#####

image_output_dir='/gpfs/work/jacobb/data/Validation/Waves/REF_seagrass/'

varnameWWM='HS' ## WWM parameter to validate
varnameBSH='Hs'
unit='[m]'
##### appearance

def set_FS(Fincrease=1.4):
	#SMALL_SIZE = 10*Fincrease
	SMALL_SIZE = 8*Fincrease
	MEDIUM_SIZE = 10*Fincrease
	BIGGER_SIZE = 12*Fincrease
	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure
	###########################
set_FS(Fincrease=1.4)


#################################################



		
# make image output dir
if not os.path.exists(image_output_dir):
	os.mkdir(image_output_dir)		
		
# load buoys
subdirs=glob(bouydir+'*')
subdirs=np.sort([sdir[sdir.rindex('/')+1:] for sdir in subdirs])
bouys=dict.fromkeys(subdirs)
for subdir in subdirs:
	#bouys[subdir]=bsh_spec(bouydir+'{:s}/{:s}_BSH.spec'.format(subdir,subdir))
	bouys[subdir]=bsh_spec(bouydir+'{:s}/{:s}_BSH.spec'.format(subdir,subdir))
#B=bouys[subdirs[0]] B.get_parameter('Hs')
################
#ELB_BSH.spec



#######  Wave Watch 3   ######
if addWW3:
	parameters=['hs']		
	ww3mesh='NBSext_bl.msh'
	ww3PointList='points.list'

	ww3Bouys_specInt=dict.fromkeys(ww3_names)
	ww3_stations_intp=dict.fromkeys(ww3_names)
	for name,diri,Pout in zip (ww3_names,ww3dirs,ww3PointOutput):	
		ww3Bouys_specInt[name]=dict.fromkeys(keysWW3)
		ww3_stations_intp[name]=dict.fromkeys(keysWW3)
		ww3=WW3_mesh(diri+ww3mesh)
		ww3.load_points(diri+ww3PointList)
		ww3.open_point_output(Pout)
		nrs=[ww3.pnames.index(key) for key in keysWW3]
		#build dictionary
		for key,nr in zip(keysWW3,nrs):
			ww3Bouys_specInt[name][key]=dict.fromkeys(['dates']+parameters)
			ww3_stations_intp[name][key]=dict.fromkeys(['dates']+parameters)#dict.fromkeys(ww3Bouys_specInt.keys())	
			ww3Bouys_specInt[name][key]['hs']=ww3.Hs[:,nr]
			ww3Bouys_specInt[name][key]['dates']=ww3.tspec

			
#a=load_wwm_site('/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_from_mistral/BetaMax2.4/prep_ELBE.site')	
#b=load_wwm_site('/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/wave_from_mistral/FRCWW31.8/prep_ELBE.site')	
#plt.figure
#plt.plot(a['HS'])
#plt.plot(b['HS'])

def load_wwm_site(file):
	with open(file) as f:
		lines=f.readline()
	header=lines.split()[:-1]
	try: 
		m=np.loadtxt(file,skiprows=1)
	except:
		m=np.loadtxt(file,skiprows=2)
	data=dict.fromkeys(header)
	for i,key in enumerate(header):
		data[key]=m[:,i]

	min=(data['TIME'])*100
	min=np.asarray(np.fix(100*(min-np.fix(min))),int)
	hour=(data['TIME'])
	hour=np.asarray(np.round(100*(hour-np.fix(hour))),int)
	day=(data['TIME']/100)
	day=np.asarray(np.fix(100*(day-np.fix(day))),int)
	mon=(data['TIME']/10000)
	mon=np.asarray(np.fix(100*(mon-np.fix(mon))),int)
	year=np.asarray(np.fix(data['TIME']/10000),int)
	data['dates']=dt.datetime(year[0],mon[0],day[0],hour[0],0,0)+np.arange(len(year))*dt.timedelta(hours=(np.diff(min[:2])/60+np.diff(hour[:2]))[0])

	return data

	
wwm_stations=dict.fromkeys(display_names)	
for dir,name in zip(WWMdirs,display_names):

	os.chdir(dir)	
	print(dir)
	##### work around station output formatting issues
	# error * in files.
	files=glob('*.site')
	sites={file.split('.')[0]:None for file in files}
	for file in files:
		with open(file) as file_in:
			if not os.path.exists('prep_'+file) and (not 'prep' in file):
				with open('prep_'+file,'w') as file_out:
					for nr,line in enumerate(file_in.readlines()):
						#if nr==99:
						#	break
						#line=line.replace('*','0').replace('Infinity','inf').replace('NaN            NaN','NaN')#.replace('NaN\tNaN','NaN') #.replace('NaN','-99')
						line=line.replace('*','0').replace('Infinity','').replace('NaN            NaN','NaN')#.replace('NaN\tNaN','NaN') #.replace('NaN','-99')
						file_out.write(line)

	files=np.sort(glob('prep*.site'))
	sites={file.split('.')[0]:None for file in files}
	for site,file in zip(sites.keys(),files):
		sites[site]=load_wwm_site(file)
	wwm_stations[name]=sites
	os.chdir('../')							
#####################


n,m=np.int(np.ceil(np.sqrt(len(wwm_stations)))),np.int(np.floor(np.sqrt(len(wwm_stations))))
wwm_stations_intp=dict.fromkeys(display_names)
for i,name in enumerate(display_names):
	wwm_stations_intp[name]=dict.fromkeys(wwm_stations[name].keys())
	
	
	
################## Plot ##############################################
scatter=[]
scatter_ww3=[]
ts=[]
taylor=[]

for key1,key2 in zip(keysOBS,keysWWM):
	key1,key2
		
	B=bouys[key1]
	dataB=B.get_parameter(varnameBSH)
	date0=B.dates[0]
	date1=B.dates[-1]
	for i,name in enumerate(display_names):
		WWM=wwm_stations[name][key2]
		date0=np.maximum(date0,WWM['dates'][0])
		date1=np.minimum(date1,WWM['dates'][-1])
	
	date0=np.maximum(date0,WWM['dates'][0]+skip_n_first_days) #ofset days
	#date1=np.maximum(date0,date1)
	
	if addWW3:
		key3=keysWW3[keysOBS.index(key1)]	
		WW3s=[]
		for name in ww3_names:
			name
			WW3s.append(ww3Bouys_specInt[name][key3])
		
	for i,name in enumerate(display_names):
		
		WWM=wwm_stations[name][key2]
		wwm_stations_intp[name][key2]=dict.fromkeys(['dates',varnameWWM,'obs'+varnameWWM])
		
		tdata,param_intp,dataWWMintp=interp_to_data_time(B.dates,dataB,WWM['dates'],WWM[varnameWWM],min_maxdate=[date0,date1])
		wwm_stations_intp[name][key2]={'dates':tdata,varnameWWM:dataWWMintp,'obs'+varnameWWM:param_intp}
		
		#plt.subplot(m,n,i+1)
		plt.clf()
		if len(dataWWMintp)>0:
			QQplot(param_intp,dataWWMintp,stationname=key1,obsname='bouy',modname=name,parameter='Hs',unit='[m]')
		plt.suptitle(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])	
		plt.tight_layout()
		fname=image_output_dir+'wwm_'+key1+'scatter_'+name+'.png'
		plt.savefig(fname,dpi=300)
		scatter.append(fname)
		
	if addWW3:		
		for WW3,name in zip(WW3s,ww3_names):
			tdata,param_intp,dataWW3intp=interp_to_data_time(B.dates,dataB,WW3['dates'],WW3['hs'],min_maxdate=[date0,date1])
			ww3_stations_intp[name][key3]={'dates':tdata,'hs':dataWW3intp,'obs'+varnameWWM:param_intp}
			plt.clf()	
			if len(dataWW3intp)>0:
				QQplot(param_intp,dataWW3intp,stationname=key1,obsname='bouy',modname=name,parameter='Hs',unit='[m]')
			plt.suptitle(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])	
			plt.tight_layout()
			fname=image_output_dir+key1+'scatter_'+name+'.png'
			scatter_ww3.append(fname)
			plt.savefig(fname,dpi=300)
	
	
	# time series
	plt.clf()
	plt.plot(B.dates,dataB,'b.',  markersize=4)
	for i,name in enumerate(display_names):
		WWM=wwm_stations[name][key2]
		plt.plot(WWM['dates'],WWM[varnameWWM],'-')
	plt.xlim((date0,date1))
	plt.ylabel(varnameWWM+unit)
	plt.gcf().autofmt_xdate()
	plt.title(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])
	if addWW3:
		for WW3,name in zip(WW3s,ww3_names):
			plt.plot(WW3['dates'],WW3['hs'],'--')
		plt.legend(['buoy']+display_names+ww3_names,loc='best',frameon=False)		
	else:	
		plt.legend(['buoy']+display_names,loc='best',frameon=False)			
	plt.grid()
	plt.tight_layout()
	fname=image_output_dir+'wwm_'+key1+'timeseries_'+'_'.join(display_names)+'.png'
	ts.append(fname)
	plt.savefig(fname,dpi=300)
	


	
## Collect samplews per run for taylor diag	
samples=[]	
wwm_bias=[]
wwm_rmse=[]
wwm_corr=[]
for i,name in enumerate(display_names):	
	corr=np.asarray([np.corrcoef(wwm_stations_intp[name][key][varnameWWM],wwm_stations_intp[name][key]['obs'+varnameWWM])[0,1] for key in wwm_stations_intp[name].keys()])
	std_rel=np.asarray([np.std(wwm_stations_intp[name][key][varnameWWM])/np.std(wwm_stations_intp[name][key]['obs'+varnameWWM]) for key in wwm_stations_intp[name].keys()])
	rmse=np.asarray([np.sqrt(np.mean((wwm_stations_intp[name][key][varnameWWM]-wwm_stations_intp[name][key]['obs'+varnameWWM])**2)) for key in wwm_stations_intp[name].keys()])
	samples.append(list(zip(std_rel,corr,keysOBS)))
	
	bias=np.asarray([np.mean((wwm_stations_intp[name][key][varnameWWM]-wwm_stations_intp[name][key]['obs'+varnameWWM])) for key in wwm_stations_intp[name].keys()])
	wwm_bias.append(bias)
	wwm_rmse.append(rmse)
	wwm_corr.append(corr)

if addWW3:
	ww3_bias=[]
	ww3_rmse=[]
	ww3_corr=[]
	for i,name in enumerate(ww3_names):	
		corr=np.asarray([np.corrcoef(ww3_stations_intp[name][key]['hs'],ww3_stations_intp[name][key]['obs'+varnameWWM])[0,1] for key in ww3_stations_intp[name].keys()])
		std_rel=np.asarray([np.std(ww3_stations_intp[name][key]['hs'])/np.std(ww3_stations_intp[name][key]['obs'+varnameWWM]) for key in ww3_stations_intp[name].keys()])
		rmse=np.asarray([np.sqrt(np.mean((ww3_stations_intp[name][key]['hs']-ww3_stations_intp[name][key]['obs'+varnameWWM])**2)) for key in ww3_stations_intp[name].keys()])
		samples.append(list(zip(std_rel,corr,keysOBS)))

		bias=np.asarray([np.mean((ww3_stations_intp[name][key]['hs']-ww3_stations_intp[name][key]['obs'+varnameWWM])) for key in ww3_stations_intp[name].keys()])
		ww3_bias.append(bias)
		ww3_rmse.append(rmse)
		ww3_corr.append(corr)
		
display_names=wwm_names+ww3_names
	
from matplotlib.text import TextPath
label = TextPath((0,0), str(1))#, linewidth=1)
plt.close('all')
dia=plotTaylor(samples[0],stdref=1,extend=False) #negative
colors=plt.cm.tab10(range(len(samples)))  #['b','r']
#Add models to Taylor diagram
phs=[plt.plot(np.nan,np.nan,color='k',marker=label,linestyle='',markersize=14)[0]]
for nr,sample in enumerate(samples[1:]):
	for i,(stddev, corrcoef, name) in enumerate(sample):
		i,stddev, corrcoef, name
		dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
	phs+=plt.plot(np.nan,np.nan,marker=label,color=colors[nr,:],linestyle='',markersize=14)
plt.legend(phs,display_names)
#plt.tight_layout()
fname=image_output_dir+'wwm_'+key1+'taylour_'+'_'.join(display_names)+'.png'
taylor.append(fname)
plt.savefig(fname,dpi=300)
plt.close()


if addWW3:
	plt.close('all')
	dia=plotTaylor(samples[0],stdref=1,extend=False) #negative
	colors=plt.cm.tab10(range(len(samples)))  #['b','r']
	#Add models to Taylor diagram
	phs=[plt.plot(np.nan,np.nan,color='k',marker=label,linestyle='',markersize=14)[0]]
	for nr,sample in enumerate(samples[1:len(wwm_names)]):
		for i,(stddev, corrcoef, name) in enumerate(sample):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
		phs+=plt.plot(np.nan,np.nan,marker=label,color=colors[nr,:],linestyle='',markersize=14)
	plt.legend(phs,wwm_names)
	#plt.tight_layout()
	fname=image_output_dir+'wwm_only_'+key1+'taylour_'+'_'.join(display_names)+'.png'
	taylor.append(fname)
	plt.savefig(fname,dpi=300)
	plt.close()
	
	plt.close('all')
	dia=plotTaylor(samples[len(wwm_names)],stdref=1,extend=False) #negative
	#Add models to Taylor diagram
	phs=[plt.plot(np.nan,np.nan,color='k',marker=label,linestyle='',markersize=14)[0]]
	for nr,sample in enumerate(samples[len(wwm_names)+1:]):
		for i,(stddev, corrcoef, name) in enumerate(sample):
			i,stddev, corrcoef, name
			dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='',mfc=colors[nr,:], mec=colors[nr,:],label=name)
		phs+=plt.plot(np.nan,np.nan,marker=label,color=colors[nr,:],linestyle='',markersize=14)
	plt.legend(phs,ww3_names)
	#plt.tight_layout()
	fname=image_output_dir+'ww3_only_'+key1+'taylour_'+'_'.join(display_names)+'.png'
	taylor.append(fname)
	plt.savefig(fname,dpi=300)
	plt.close()
	
####### Report #############################


###### make report #######		
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
import time

		
chart_style = TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                          ('VALIGN', (0, 0), (-1, -1), 'CENTER')])

chart_style = TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTER')])

						  
doc = SimpleDocTemplate(image_output_dir+pdfname,pagesize=letter,
                        rightMargin=72,leftMargin=72,
                        topMargin=72,bottomMargin=18,title='Blubb')

doc.title="WWM Validation"
					
time.ctime()

title = 'Validation of WWM run against BSH bouys'							
formatted_time = time.ctime()
full_name = "Benjamin Jacob"


styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

Story=[]
Story.append(Paragraph(title, styles["title"]))
Story.append(Spacer(1, 12))
Story.append(Spacer(1, 12))

####
ptext = '%s' % formatted_time
Story.append(Paragraph(ptext, styles["Normal"]))
# Create return address
ptext = '%s' % full_name
Story.append(Paragraph(ptext, styles["Normal"]))    
Story.append(Spacer(4, 12))

ptext='Validation of WWM '
if len(ww3_names) >0:
	ptext+='and WW3 '
ptext+='run(s) for period ' +  str(date0)[:10] + ' - '+ str(date1)[:11] +'.'
Story.append(Paragraph(ptext, styles["Normal"]))    
if len(wwm_names) >0:
	ptext='WWM scenarios are: ' 
	Story.append(Paragraph(ptext, styles["Normal"]))    
	for name,desc in zip(wwm_names,wwm_descriptions):
		ptext=' '+name + desc +', '
		Story.append(Spacer(0.5, 5))
		Story.append(Paragraph(ptext, styles["Normal"]))    
#ptext=ptext[:-1]	

if len(ww3_names) >0:
	Story.append(Spacer(1, 12))
	ptext='WW3 scenarios are ' 
	Story.append(Paragraph(ptext, styles["Normal"]))    
	for name,desc in zip(ww3_names,ww3_descriptions):
		ptext=' '+name + desc
		Story.append(Spacer(0.5, 5))
		Story.append(Paragraph(ptext, styles["Normal"]))    
	
#ptext=ptext[:-1]	
ptext=comments
Story.append(Paragraph(ptext, styles["Normal"]))    
	
counter=1
Story.append(Spacer(4, 12))
Story.append(Paragraph('Overview', styles["h1"]))
im = Image(taylor[0],4.5*inch,4*inch)
Story.append(im)
Story.append(Paragraph('Fig {:d}: taylor diagram (all runs) of significant wave height performance against BSH buoys for period {:s} - {:s} -All samples.'.format(counter,str(date0),str(date1)), styles["Normal"]))
Story.append(Spacer(2, 12))

if addWW3:
	counter+=1
	Story.append(Spacer(4, 12))
	im = Image(taylor[1],4.5*inch,4*inch)
	Story.append(im)
	Story.append(Paragraph('Fig {:d}: taylor diagram (WWM only) of significant wave height performance against BSH buoys for period {:s} - {:s}.'.format(counter,str(date0),str(date1)), styles["Normal"]))
	Story.append(Spacer(2, 12))

	counter+=1
	Story.append(Spacer(4, 12))
	im = Image(taylor[2],4.5*inch,4*inch)
	Story.append(im)
	Story.append(Paragraph('Fig {:d}: taylor diagram (WW3 only) of significant wave height performance against BSH buoys for period {:s} - {:s}.'.format(counter,str(date0),str(date1)), styles["Normal"]))
	Story.append(Spacer(2, 12))


ptext='The station averaged values (each station weighted equally independent of nr of observtions) correspond to simulation: bias/rmse/correlation '
Story.append(Paragraph(ptext, styles["Normal"]))    
Story.append(Spacer(1, 12))
for i,name in enumerate(wwm_names):
	Story.append(Spacer(1, 12))
	ptext='\n{:s}: {:f}/{:f}/{:f}, '.format(name,wwm_bias[i].mean(),wwm_rmse[i].mean(),wwm_corr[i].mean())
	Story.append(Paragraph(ptext, styles["Normal"]))    
if len(ww3_names) >0:	
	for i,name in enumerate(ww3_names):
		Story.append(Spacer(1, 12))
		ptext='\n{:s}: {:f}/{:f}/{:f}, '.format(name,ww3_bias[i].mean(),ww3_rmse[i].mean(),ww3_corr[i].mean())
		Story.append(Paragraph(ptext, styles["Normal"]))    	
#ptext=ptext[:-1]+'.'
#Story.append(Paragraph(ptext, styles["Normal"]))    
Story.append(Spacer(4, 12))
Story.append(Paragraph('time series', styles["h1"]))



for i,fig in enumerate(ts):
	Story.append(Spacer(4, 12))
	im = Image(ts[i],4*inch,3*inch)
	Story.append(im)
	counter+=1
	Story.append(Paragraph('Fig {:d}: Time Series of significant wave height at station {:s}.'.format(counter,str(keysWW3[i])), styles["Normal"]))
	

Story.append(Spacer(4, 12))
Story.append(Paragraph('scatter plots', styles["h1"]))


# 1 row per station
nwwm=len(wwm_names)
nww3=len(ww3_names)

nruns=len(display_names)
sqrt=np.sqrt(nruns)	
ncols=np.minimum(nruns,4)
nrows=np.int(np.ceil(nruns/ncols))

inchW=8.2/ncols*1.1
inchH=6/ncols*1.1
for i in range(len(keysWW3)):
	Story.append(Spacer(4, 12))
	
	counter+=1	
	ims=[]
	for k in range(i*(nwwm),(i+1)*nwwm):
		ims.append(Image(scatter[k],inchW*inch,inchH*inch))
	
	if addWW3:
		for k in range((i)*nww3,(i+1)*nww3):
			ims.append(Image(scatter_ww3[k],inchW*inch,inchH*inch))
		
	for irow in range(nrows):
		i0=irow*ncols
		i1=i0+ncols	
		i1=np.minimum(i1,len(ims))
		Story.append(Table([ims[i0:i1]],
			colWidths=[inchW * inch, inchW * inch],
			rowHeights=[inchH * inch], style=chart_style))				
	Story.append(Spacer(2, 0))			
			
	Story.append(Paragraph('Fig {:d}: Scatter plot of significant wave height from top left to bottom right (counting along rows) simulated by '.format(counter) +' '.join(display_names) +' against BSH buoys at station {:s}.'.format(str(keysWW3[i])), styles["Normal"]))
doc.build(Story)				
		
