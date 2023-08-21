"""
Make Initial conditions (hvar_*.ic) for Water quality model ICM 
eitehr as constant values or taken from CMEMS
"""

import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
from scipy.spatial import cKDTree

#!---------------------------------------------------------------------------------
#!---------------------------state variables in ICM--------------------------------
#!---------------------------------------------------------------------------------
#! 1  ZB1   :  1st zooplankton                            g/m^3
#! 2  ZB2   :  2nd zooplankton                            g/m^3
#! 3  PB1   :  Diatom                                     g/m^3
#! 4  PB2   :  Green Algae                                g/m^3
#! 5  PB3   :  Cyanobacteria                              g/m^3
#! 6  RPOC  :  Refractory Particulate Organic Carbon      g/m^3
#! 7  LPOC  :  Labile Particulate Organic Carbon          g/m^3
#! 8  DOC   :  Dissolved Orgnaic Carbon                   g/m^3
#! 9  RPON  :  Refractory Particulate Organic Nitrogen    g/m^3
#! 10 LPON  :  Labile Particulate Organic Nitrogen        g/m^3
#! 11 DON   :  Dissolved Orgnaic Nitrogen                 g/m^3
#! 12 NH4   :  Ammonium Nitrogen                          g/m^3
#! 13 NO3   :  Nitrate Nitrogen                           g/m^3
#! 14 RPOP  :  Refractory Particulate Organic Phosphorus  g/m^3
#! 15 LPOP  :  Labile Particulate Organic Phosphorus      g/m^3
#! 16 DOP   :  Dissolved Orgnaic Phosphorus               g/m^3
#! 17 PO4t  :  Total Phosphate                            g/m^3
#! 18 SU    :  Particulate Biogenic Silica                g/m^3
#! 19 SAt   :  Available Silica                           g/m^3
#! 20 COD   :  Chemical Oxygen Demand                     g/m^3
#! 21 DOX   :  Dissolved Oxygen                           g/m^3

#----------- PhD module -----------------------------------------
#! 22 TIC   :  Total Inorganic Carbon                     g/m^3
#! 23 ALK   :  Alkalinity                                 g[CaCO3]/m^3
#! 24 CA    :  Dissolved Calcium                          g[CaCO3]/m^3
#! 25 CACO3 :  Calcium Carbonate                          g[CaCO3]/m^3
#!---------------------------------------------------------------------------------

######################### motu download of CMEMS files #######################################################
do_download=False
if do_download:
	import motuclient
	import subprocess

	date='2022-01-01'
	user='bjacob'
	passwd='FreFru86'

	pids={'cmems_mod_nws_bgc-no3_anfc_7km-3D_P1D-m':'no3',
	'cmems_mod_nws_bgc-chl_anfc_7km-3D_P1D-m':'chl',
	'cmems_mod_nws_bgc-phyc_anfc_7km-3D_P1D-m':'phyc',
	'cmems_mod_nws_bgc-o2_anfc_7km-3D_P1D-m':'o2',
	'cmems_mod_nws_bgc-po4_anfc_7km-3D_P1D-m':'po4'}

	for pid in pids.keys():
		var=pids[pid]
		print('doing ' + var)
		command='python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id NWSHELF_ANALYSISFORECAST_BGC_004_002-TDS --product-id {:s} --date-min "{:s} 00:00:00" --date-max "{:s} 23:59:59" --depth-min 0 --depth-max 0 --variable {:s} --out-dir "./" --out-name "{:s}_{:s}.nc" --user "{:s}" --pwd "{:s}"'.format(pid,date,date,var,date.replace('-',''),var,user,passwd)
		subprocess.run(command, shell=True, check=True)


# conversion factors from mol/l   to mg/l
mmo2=31.99880 #g/mol
mmo2=31.99880/1000 #g/mmol

mmno3=62.0049 #g/mol
mmno3=62.0049/1000 #g/mmol

mmpo4=94.9714 #g/mol
mmpo4=94.9714/1000 #g/mmol

mmC=12.01070 #g/mol
mmC=12.01070/1000 #g/mmol

mmN=14.0067 #g/mol
mmN=14.0067/1000 #g/mmol


mmSio2=60.08 #g/mol
mmSio2=60.08/1000 #g/mol


#mmol m-3

# is alread in gram per liter
#mmCHLa=893.509 #g·mol 
#mmCHLa=893.509/1000 #g·mmol 

# https://www.sciencedirect.com/topics/earth-and-planetary-sciences/dissolved-organic-nitrogen
#DON is composed of both labile and recalcitrant, high molecular weight and low molecular weight molecules


def read_cmemsbio(file,varname,factor):
	""" read bio gep paramters from CMEMs for ICM model converting mol/l to g/l via factor"""
	ds=xr.open_dataset(file)
	
	inputvals=ds[varname][0,0,:].values.flatten()
	ivalid=~np.isnan(inputvals)
	inputvals=inputvals[ivalid]
	#lat,lon=np.meshgrid(ds.latitude.values,ds.longitude.values)
	lon,lat=np.meshgrid(ds.longitude.values,ds.latitude.values)
	#coords=list(zip(ds.longitude.values.flatten()[ivalid],ds.latitude.values.flatten()[ivalid]))
	coords=list(zip(lon.flatten()[ivalid],lat.flatten()[ivalid]))
	nnTree=cKDTree(coords)
	nns=nnTree.query(list(zip(s.lon,s.lat)))[1]
	#inputvals[]=0
	return inputvals[nns]


	
# ospar	
# https://oap.ospar.org/en/ospar-assessments/intermediate-assessment-2017/pressures-human-activities/eutrophication/nutrients-concentrations/
#1 Mol/Liter [mol/l] = 1000000 Mikromolar [µM]	


#https://bg.copernicus.org/articles/16/1073/2019/  (Fig5)
# Sommer (Sollte winter sein) DON  5 µM 	 5/1000  mmol/Liter
CN=mmN*5/1000
	
config=[]
#
# using iZB= 0 no Zooplankton, using decay rate for phyto
# Plankton
# algal biomass in mg[C]/l = g[C]/m^3
#! 1  ZB1   :  1st zooplankton                            g/m^3
ZB1={'type':{'const':{'value':0}}}
#! 2  ZB2   :  2nd zooplankton                            g/m^3
ZB2={'type':{'const':{'value':0}}}
#! 3  PB1   :  Diatom                                     g/m^3
#PB1={'type':{'const':{'value':0}}}  # Diatoms gro at colder temperatures. higher concentrations earlier in year than phyto. I naivlly guess 3x the ammount to have something 
PB1={'type':{'file':{'fname':'./ICM/20220101_phyc.nc','varname':'phyc','factor':mmC*3}}} #  need to apply carbon to chlorphy ratuib
#! 4  PB2   :  Green Algae                                g/m^3
PB2={'type':{'file':{'fname':'./ICM/20220101_phyc.nc','varname':'phyc','factor':mmC}}} # Green Algae # need to apply carbon to chlorphy ratuib
# phyc as prox increas value using redfield ratio
# carbon equivalent cmems
#! 5  PB3   :  Cyanobacteria                              g/m^3
#PB3={'type':{'const':{'value':0}}}
PB3={'type':{'file':{'fname':'./ICM/20220101_phyc.nc','varname':'phyc','factor':mmC/10}}}

#refractory:= not accessibale to rapid microbial degradartion
#Carbon
#! 6  RPOC  :  Refractory Particulate Organic Carbon      g/m^3
#RPOC={'type':{'const':{'value':0}}}
RPOC={'type':{'file':{'fname':'./ICM/20220101_chl.nc','varname':'chl','factor':1/1000*0.5}}} 
#! 7  LPOC  :  Labile Particulate Organic Carbon          g/m^3
#LPOC={'type':{'const':{'value':0}}}
LPOC={'type':{'file':{'fname':'./ICM/20220101_chl.nc','varname':'chl','factor':1/1000*0.5}}} 
#! 8  DOC   :  Dissolved Orgnaic Carbon                   g/m^3
DOC={'type':{'file':{'fname':'./ICM/20220101_chl.nc','varname':'chl','factor':1/1000}}} #not that   CMEMS is mg/m3  ICM g/m3

# Nitrate
#DON is composed of both labile and recalcitrant, high molecular weight and low molecular weight molecules
#!9  RPON  :  Refractory Particulate Organic Nitrogen    g/m^3
RPON={'type':{'const':{'value':CN*0.5}}}   # one can relate RPON=0.5 DON
#! 10 LPON  :  Labile Particulate Organic Nitrogen        g/m^3
LPON={'type':{'const':{'value':CN*0.5}}}
#! 11 DON   :  Dissolved Orgnaic Nitrogen                 g/m^3
DON={'type':{'const':{'value':CN}}}
#! 12 NH4   :  Ammonium Nitrogen                          g/m^3
#NHR={'type':{'const':{'value':0}}}
NHR={'type':{'file':{'fname':'./ICM/20220101_po4.nc','varname':'po4','factor':mmpo4*5}}}
#https://www.researchgate.net/publication/257538675_Microbial_biogeography_of_the_North_Sea_during_summer/figures?lo=1
# I take ammonium 5 times phosphate for which I have data
#! 13 NO3   :  Nitrate Nitrogen                           g/m^3
NO3={'type':{'file':{'fname':'./ICM/20220101_no3.nc','varname':'no3','factor':mmno3}}}


#https://link.springer.com/article/10.1007/bf02764035




# Phosphorus
#! 14 RPOP  :  Refractory Particulate Organic Phosphorus  g/m^3
PROP={'type':{'const':{'value':0}}}
#! 15 LPOP  :  Labile Particulate Organic Phosphorus      g/m^3
LPOP={'type':{'const':{'value':0}}}
#! 16 DOP   :  Dissolved Orgnaic Phosphorus               g/m^3
DOP={'type':{'const':{'value':0}}}
#! 17 PO4t  :  Total Phosphate                            g/m^3
PO4t={'type':{'file':{'fname':'./ICM/20220101_po4.nc','varname':'po4','factor':mmpo4}}}





#Silica
#! 18 SU    :  Particulate Biogenic Silica                g/m^3
#1 µ mol / l # https://www.researchgate.net/publication/236628526_Diatom_succession_silicification_and_silicic_acid_availability_in_Belgian_coastal_waters_Southern_North_Sea fig2
SU={'type':{'const':{'value':mmSio2/1000 * 1}}}
#! 19 SAt   :  Available Silica                           g/m^3
#SAt={'type':{'const':{'value':0}}}
SAt={'type':{'const':{'value':mmSio2/1000 * 1}}}

# Oxygen
#! 20 COD   :  Chemical Oxygen Demand                     g/m^3
#Der Chemische Sauerstoffbedarf ist als Summenparameter ein Maß für die Summe aller im Wasser vorhandenen, unter bestimmten Bedingungen oxidierbaren Stoffe. Er gibt die Menge an Sauerstoff an, die zu ihrer Oxidation benötigt würde, wenn Sauerstoff das Oxidationsmittel wäre
COD={'type':{'const':{'value':0}}}
#! 21 DOX   :  Dissolved Oxygen                           g/m^3
#CMEMS Mole concentration of dissolved molecular oxygen in sea water (O2)
DOX={'type':{'file':{'fname':'./ICM/20220101_o2.nc','varname':'o2','factor':mmo2}}}

s=schism_setup()
ones=np.ones(s.nnodes)

# loop over varibale configs and create respective ic file
config=[ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,LPON,DON,NHR,NO3,PROP,LPOP,DOP,PO4t,SU,SAt,COD,DOX]
for nr,varib in enumerate(config,1):
	key=list(varib['type'].keys())[0]
	if key =='const':
		values=varib['type'][key]['value']*ones
	#s.dump_gr3_spat_var('ICM_hvar_{:d}.ic'.format(nr),values,comment='created by genICMinuputs.py')
	elif key =='file':
		values=read_cmemsbio(file=varib['type'][key]['fname'],varname=varib['type'][key]['varname'],factor=varib['type'][key]['factor'])
	s.dump_gr3_spat_var('ICM_hvar_{:d}.ic'.format(nr),values,comment='created by genICMinuputs.py')
# plot maps
# load schism setup
#plt.ion()





