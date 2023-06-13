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

# conversion factors from mol/l   to mg/l
mmo2=31.99880 #g/mol
mmo2=31.99880/1000 #g/mmol

mmno3=62.0049 #g/mol
mmno3=62.0049/1000 #g/mmol

mmpo4=94.9714 #g/mol
mmpo4=94.9714/1000 #g/mmol

mmC=12.01070 #g/mol
mmC=12.01070/1000 #g/mmol


# is alread in gram per liter
#mmCHLa=893.509 #g·mol 
#mmCHLa=893.509/1000 #g·mmol 

def read_cmemsbio(file,varname,factor):
	""" read bio gep paramters from CMEMs for ICM model converting mol/l to g/l via factor"""
	ds=xr.open_dataset(file)
	
	inputvals=ds[varname][0,0,:].values.flatten()
	ivalid=!np.isnan(inputvals)
	inputvals=inputvals[ivalid]
	coords=list(zip(ds.longitude.values.flatten()[ivalid],ds.latitude.values.flatten()[ivalid]))
	nnTree=cKDTree(coords)
	nns=nnTree.query(list(zip(s.lon,s.lat)))[1]
	
	inputvals[]=0
	return inputvals[nns]


config=[]
# Plankton
#! 1  ZB1   :  1st zooplankton                            g/m^3
ZB1={'type':{'const':{'value':0}}}
#! 2  ZB2   :  2nd zooplankton                            g/m^3
ZB2={'type':{'const':{'value':0}}}
#! 3  PB1   :  Diatom                                     g/m^3
PB1={'type':{'const':{'value':0}}}
#! 4  PB2   :  Green Algae                                g/m^3
PB2={'type':{'file':{'fname':'./ICM/20220101_phyc.nc','varname':'phyc','factor':1}}}
# phyc as prox increas value using redfield ratio
# carbon equivalent cmems
#! 5  PB3   :  Cyanobacteria                              g/m^3
PB3={'type':{'const':{'value':0}}}

#Carbon
#! 6  RPOC  :  Refractory Particulate Organic Carbon      g/m^3
RPOC={'type':{'const':{'value':0}}}
#! 7  LPOC  :  Labile Particulate Organic Carbon          g/m^3
LPOC={'type':{'const':{'value':0}}}
#! 8  DOC   :  Dissolved Orgnaic Carbon                   g/m^3
DOC={'type':{'file':{'fname':'./ICM/20220101_chl.nc','varname':'chl','factor':1}}} #not that

# Nitrate
#! 10 LPON  :  Labile Particulate Organic Nitrogen        g/m^3
LPON={'type':{'const':{'value':0}}}
#! 11 DON   :  Dissolved Orgnaic Nitrogen                 g/m^3
DON={'type':{'const':{'value':0}}}
#! 12 NH4   :  Ammonium Nitrogen                          g/m^3
NHR={'type':{'const':{'value':0}}}
#! 13 NO3   :  Nitrate Nitrogen                           g/m^3
NO3={'type':{'file':{'fname':'./ICM/20220101_no3.nc','varname':'no3','factor':mmno3}}}

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
SU={'type':{'const':{'value':0}}}
#! 19 SAt   :  Available Silica                           g/m^3
SAt={'type':{'const':{'value':0}}}

# Oxygen
#! 20 COD   :  Chemical Oxygen Demand                     g/m^3
#Der Chemische Sauerstoffbedarf ist als Summenparameter ein Maß für die Summe aller im Wasser vorhandenen, unter bestimmten Bedingungen oxidierbaren Stoffe. Er gibt die Menge an Sauerstoff an, die zu ihrer Oxidation benötigt würde, wenn Sauerstoff das Oxidationsmittel wäre
COD={'type':{'const':{'value':0}}}
#! 21 DOX   :  Dissolved Oxygen                           g/m^3
#CMEMS Mole concentration of dissolved molecular oxygen in sea water (O2)
DOX={'type':{'file':{'fname':'./ICM/20220101_po2.nc','varname':'o2','factor':mmo2}}}

s=schism_setup()
ones=np.ones(s.nnodes)

# loop over varibale configs and create respective ic file
config=[ZB1,ZB2,PB1,PB2,PB3,PROC,LPOC,DOC,LPON,DON,NHR,NO3,PROP,LPOP,DOP,PO4t,SU,SAt,COD,DOX]
for nr,varib in enumerate(config,1):
	key=list(varib['type'].keys())[0]
	if key =='const':
		values=varib['type'][key]['value']*ones
		s.dump_gr3_spat_var('ICM_hvar_{:d}.ic'.format(nr),values,comment='created by genICMinuputs.py')
	elif key =='file':
		values=read_cmemsbio(file=varib['type'][key]['fname'],varname=varib['type'][key]['varname'])

# plot maps
# load schism setup
#plt.ion()





######################### motu download of CMEMS files #######################################################
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
