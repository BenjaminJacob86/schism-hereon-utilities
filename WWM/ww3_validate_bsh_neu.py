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
from matplotlib.collections import PolyCollection
########### SETTINGS #########
bouydir='/gpfs/work/jacobb/data/SETUPS/WW3/bouy/bouys/' #bsh buoy directory containing subfolders with the bouys and in those a file with name <subfolder>_BSH.spec

comments='' #sentences following wwm and ww3 descriptions in pdf

# names of stations
keysOBS=['ELB','FN1','FN3','HEL','WES']
keysWW3=['Elbe', 'Fino-1','Fino-3','Helgoland','Westerland']

## WWW3 ###
ww3dirs=['/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd_routine/',]#['/gpfs/work/jacobb/data/SETUPS/WW3/WW4NBSbnd/',]
ww3mesh='NBSext_bl.msh'
ww3PointList='points.list'
#ww3GridOutput==np.sort(glob(ww3dir+'mix_select_prev_post_run/'+'ww3.??????.nc'))# np.sort(glob(ww3dir+'defect_ncs/'+'ww3.????01.nc')) # ww3.201701.nc
ww3PointOutput=[ww3dirs[0]+ww3.202108_09_spec.nc',]#'ww3.2021_01to03_spec.nc',] #ww3.202108_spec.nc',]
ww3_names=['WW3_DWD']
ww3_descriptions=[': WW3 run with betamax of 1.5 and Global DWD NRT wind',]


pdfname="WW3_Validation_DWD_NRT_2021_2.pdf" # name of image pdf
#####

image_output_dir=ww3dirs[0]
varnameBSH='Hs'
unit='[m]'
##### appearance

def set_FS(Fincrease=1.4):
	#SMALL_SIZE = 10*Fincrease
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
class bsh_ncbouy:
	def __init__(self,files):
		self.ds=xr.open_mfdataset(files)
		self.dates=self.ds['time'].values
	def	get_parameter(self,varname=''):
		return self.ds[varname.lower()].values

	
# make image output dir
if not os.path.exists(image_output_dir):
	os.mkdir(image_output_dir)		
		
# load buoys
subdirs=glob(bouydir+'*')
subdirs=np.sort([sdir[sdir.rindex('/')+1:] for sdir in subdirs])
bouys=dict.fromkeys(subdirs)
for subdir in subdirs:
	files=list(np.sort(glob(bouydir+'{:s}/*.spt'.format(subdir))))
	if len(files)>0:
		bouys[subdir]=bsh_spec(files=files)
		ncbouys=False
	else: #checheck netcdf 	
		files=list(np.sort(glob(bouydir+'{:s}/*.nc'.format(subdir))))
		bouys[subdir]=bsh_ncbouy(files=files)
			
		ncbouys=True
		
#np.diff(bouys['ELB']['time'])

#tdata,param_intp,dataWW3intp=interp_to_data_time(B.dates,dataB,WW3['dates'],WW3['hs'],min_maxdate=[date0,date1])
		
#B=bouys[subdirs[0]] B.get_parameter('Hs')
################

ww3Bouys_specInt=dict.fromkeys(ww3_names)
ww3_stations_intp=dict.fromkeys(ww3_names)
ww3s=dict.fromkeys(ww3_names)
for name,diri,Pout in zip (ww3_names,ww3dirs,ww3PointOutput):	
	ww3Bouys_specInt[name]=dict.fromkeys(keysWW3)
	ww3_stations_intp[name]=dict.fromkeys(keysWW3)
	ww3=WW3_mesh(diri+ww3mesh)
	ww3s[name]=ww3
	


	
	
#######  Wave Watch 3   ######
#nrs=[ww3.pnames.index(key) for key in keysWW3]
#plt.plot(ww3.specds['time']ww3.specds['wnd'][:,nr]
#plt.figure()
#ww3.specds['wnd'][:].plot()


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

	
	
################## Plot ##############################################
scatter_ww3=[]
ts=[]
taylor=[]
from cycler import cycler
import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = cycler(color='rcmykbg')

for key1,key2 in zip(keysOBS,keysWW3):
	
		
	B=bouys[key1]
	dataB=B.get_parameter(varnameBSH)
	date0=B.dates[0]
	date1=B.dates[-1]
	
	key3=keysWW3[keysOBS.index(key1)]	
	WW3s=[]
	for name in ww3_names:
		WW3s.append(ww3Bouys_specInt[name][key2])
		WW3=ww3Bouys_specInt[name][key2]
		date0=np.maximum(date0,WW3['dates'][0])
		date1=np.minimum(date1,WW3['dates'][-1])
	
	
	for WW3,name in zip(WW3s,ww3_names):
		tdata,param_intp,dataWW3intp=interp_to_data_time(B.dates,dataB,WW3['dates'],WW3['hs'],min_maxdate=[date0,date1])
		ww3_stations_intp[name][key3]={'dates':tdata,'hs':dataWW3intp,'obs'+'Hs':param_intp}
		plt.clf()	
		QQplot(param_intp,dataWW3intp,stationname=key1,obsname='bouy',modname=name,parameter='Hs',unit='[m]')
		plt.suptitle(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])	
		plt.tight_layout()
		fname=image_output_dir+key1+'scatter_'+name+'.png'
		scatter_ww3.append(fname)
		plt.savefig(fname,dpi=300)
	
	
	# time series
	plt.clf()
	plt.plot(B.dates,dataB,'b.',  markersize=4)
	for WW3,name in zip(WW3s,ww3_names):
		plt.plot(WW3['dates'],WW3['hs'],'--')
	plt.legend(['buoy']+ww3_names,loc='best',frameon=False)		
	plt.xlim((date0,np.maximum(date1,WW3['dates'][-1])))
	plt.ylabel('Hs [m]')
	plt.gcf().autofmt_xdate()
	plt.title(key1 + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])
	plt.grid()
	plt.tight_layout()
	fname=image_output_dir+'ww3_'+key1+'timeseries_'+'_'.join(ww3_names)+'.png'
	ts.append(fname)
	plt.savefig(fname,dpi=300)
	


	
## Collect samplews per run for taylor diag	
samples=[]	
ww3_bias=[]
ww3_rmse=[]
ww3_corr=[]
for i,name in enumerate(ww3_names):	
	corr=np.asarray([np.corrcoef(ww3_stations_intp[name][key]['hs'],ww3_stations_intp[name][key]['obs'+'Hs'])[0,1] for key in ww3_stations_intp[name].keys()])
	std_rel=np.asarray([np.std(ww3_stations_intp[name][key]['hs'])/np.std(ww3_stations_intp[name][key]['obs'+'Hs']) for key in ww3_stations_intp[name].keys()])
	rmse=np.asarray([np.sqrt(np.mean((ww3_stations_intp[name][key]['hs']-ww3_stations_intp[name][key]['obs'+'Hs'])**2)) for key in ww3_stations_intp[name].keys()])
	samples.append(list(zip(std_rel,corr,keysOBS)))

	bias=np.asarray([np.mean((ww3_stations_intp[name][key]['hs']-ww3_stations_intp[name][key]['obs'+'Hs'])) for key in ww3_stations_intp[name].keys()])
	ww3_bias.append(bias)
	ww3_rmse.append(rmse)
	ww3_corr.append(corr)
ww3_names=ww3_names
	
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
plt.legend(phs,ww3_names)
#plt.tight_layout()
fname=image_output_dir+'ww3_'+key1+'taylour_'+'_'.join(ww3_names)+'.png'
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

title = 'Validation of WW3 run against BSH bouys'							
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

ptext='Validation of WW3 run(s) '
ptext+='run(s) for period ' +  str(date0)[:10] + ' - '+ str(date1)[:11] +'.'
if len(ww3_names) >0:
	Story.append(Spacer(1, 12))
	ptext='WW3 scenarios are ' 
	Story.append(Paragraph(ptext, styles["Normal"]))    
	for name,desc in zip(ww3_names,ww3_descriptions):
		ptext=' '+name + desc
		Story.append(Spacer(0.5, 5))
		Story.append(Paragraph(ptext, styles["Normal"]))    

ptext=comments
Story.append(Paragraph(ptext, styles["Normal"]))    
	
counter=1
Story.append(Spacer(4, 12))
Story.append(Paragraph('Overview', styles["h1"]))
im = Image(taylor[0],4.5*inch,4*inch)
Story.append(im)
Story.append(Paragraph('Fig {:d}: taylor diagram (all runs) of significant wave height performance against BSH buoys for period {:s} - {:s} -All samples.'.format(counter,str(date0),str(date1)), styles["Normal"]))
Story.append(Spacer(2, 12))


ptext='The station averaged values (each station weighted equally independent of nr of observtions) correspond to simulation: bias/rmse/correlation '
Story.append(Paragraph(ptext, styles["Normal"]))    
Story.append(Spacer(1, 12))
for i,name in enumerate(ww3_names):
	Story.append(Spacer(1, 12))
	ptext='\n{:s}: {:f}/{:f}/{:f}, '.format(name,ww3_bias[i].mean(),ww3_rmse[i].mean(),ww3_corr[i].mean())
	Story.append(Paragraph(ptext, styles["Normal"]))    	
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
nww3=len(ww3_names)

nruns=len(ww3_names)
sqrt=np.sqrt(nruns)	
ncols=np.minimum(nruns,4)
nrows=np.int(np.ceil(nruns/ncols))

inchW=8.2/np.max((ncols,2))*1.1
inchH=6/np.max((ncols,2))*1.1
for i in range(len(keysWW3)):
	Story.append(Spacer(4, 12))
	
	counter+=1	
	ims=[]
	

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
			
	Story.append(Paragraph('Fig {:d}: Scatter plot of significant wave height from top left to bottom right (counting along rows) simulated by '.format(counter) +' '.join(ww3_names) +' against BSH buoys at station {:s}.'.format(str(keysWW3[i])), styles["Normal"]))
doc.build(Story)				
		
