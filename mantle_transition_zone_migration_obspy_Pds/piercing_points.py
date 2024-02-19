#!/usr/bin/python -u
"""
This script reads receiver function data from a set of stations, and
calculates the ray path between all pairs station-event. The data MUST be 
located in some root folder (RF_DIR). The station-event informations 
(coordinates) MUST be in the header of the files. In the current version 
of the program, files should be organized inside their directory as you prefer, 
but they need to finish with some extension name (RF_EXT), e.g.:
*.mseed, *.sac, *.itr, *.eqr
The implemented algorithm follows the lines of Gao and Liu (2014),
"Imaging mantle discontinuities using multiply-reflected P-to-S
conversions", Earth and Planetary Science Letters 402 (2014) 99–106.
Note that all the parameters are defined in the configuration file.
"""
import warnings
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool
import time
import glob
from tqdm import tqdm


# =====================================
# Importing piercing points scritp_py 
# =====================================

from piercing_points_py.piercing_points_phase import arrivals_calculation

# ==================================================
# Importing some parameters from configuration file 
# ==================================================

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,MP_PROCESSES,OUTPUT_DIR,DEPTH_TARGET
				   )

# ============================================
# Importing station dictionary from JSON file 
# ============================================

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

print('\n')
print('Looking for receiver functions data in JSON file in '+STA_DIR)
print('\n')

filename_STA = STA_DIR+'sta_dic.json'

sta_dic = json.load(open(filename_STA))

event_depth = sta_dic['event_depth']
event_lat = sta_dic['event_lat']
event_long = sta_dic['event_long']
event_dist = sta_dic['event_dist']
event_gcarc = sta_dic['event_gcarc']
event_sta = sta_dic['event_sta']
event_ray = sta_dic['event_ray']
sta_lat = sta_dic['sta_lat']
sta_long = sta_dic['sta_long']
sta_data = sta_dic['sta_data']
sta_time = sta_dic['sta_time']

# ==================================================
# Importing earth model from obspy.taup.TauPyModel 
# ==================================================

print('Importing earth model from obspy.taup.TauPyModel')
print('\n')
model_THICKNESS_km = TauPyModel(model=MODEL_FILE_NPZ)


# =========================
# Creating the earth layers
# =========================

print('Creating the earth layers')
print('\n')

camadas_terra_10_km = np.arange(MIN_DEPTH,MAX_DEPTH+INTER_DEPTH,INTER_DEPTH)

# ====================
# Creating Pds list
# ====================

print('Creating Pds list')
print('\n')

PHASES = 'P410s','P'+str(DEPTH_TARGET)+'s','P660s'

# ========================
# Creating output Folder
# ========================

print('Creating output folders')

PP_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Piercing_Points'+'/'

Phase_P410s_folder = PP_DIR+'P410s/'
os.makedirs(Phase_P410s_folder,exist_ok=True)
print('Phase_P410s_folder = '+Phase_P410s_folder)

Phase_TARGET_Ps_folder = PP_DIR+'P'+str(DEPTH_TARGET)+'s/'
print('Phase_'+'P'+str(DEPTH_TARGET)+'s'+'_folder = '+Phase_TARGET_Ps_folder)
os.makedirs(Phase_TARGET_Ps_folder,exist_ok=True)

Phase_P660s_folder = PP_DIR+'P660s/'
print('Phase_P660s_folder = '+Phase_P660s_folder)
os.makedirs(Phase_P660s_folder,exist_ok=True)


# =========================
# Creating Pds Input lists
# ========================

print('Creating Pds Input lists')
print('\n')

input_list_410 = [[i+1,PHASES[0],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],Phase_P410s_folder] for i,j in enumerate(event_depth)]

input_list_TARGET = [[i+1,PHASES[1],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],Phase_TARGET_Ps_folder] for i,j in enumerate(event_depth)]

input_list_660 = [[i+1,PHASES[2],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],Phase_P660s_folder] for i,j in enumerate(event_depth)]



# ============================================================================================================================
# PDS PHASES
# ============================================================================================================================


# ===================================
# Calculating P410s Piercing Points
# ===================================

print('Calculating Piercing Points to '+PHASES[0])
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES) as p:
	max_ = len(input_list_410)
	with tqdm(total=max_,desc='Calculating Piercing Points') as pbar:
		for i, _ in enumerate(p.imap_unordered(arrivals_calculation,input_list_410)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(PHASES[0]+' Piercing Points estimated!')
print('\n')


# =============================
# Saving Piercing Points P410s
# =============================

filename_pds_json_P410s = sorted(glob.glob(Phase_P410s_folder+'*'))

PP_dic_P410s_files = [json.load(open(i)) for i in filename_pds_json_P410s]

PP_dic_P410s = {'depth':[],'time':[],'lat':[],'lon':[]}

for i,j in enumerate(PP_dic_P410s_files):
	PP_dic_P410s['time'].append(j['time'])
	PP_dic_P410s['depth'].append(j['depth'])
	PP_dic_P410s['lat'].append(j['lat']) 
	PP_dic_P410s['lon'].append(j['lon']) 


print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[0]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_P410s, fp)

# ===================================
# Calculating TARGET Piercing Points
# ===================================
 
print('Calculating Piercing Points to '+PHASES[1])
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES) as p:
	max_ = len(input_list_TARGET)
	with tqdm(total=max_,desc='Calculating Piercing Points') as pbar:
		for i, _ in enumerate(p.imap_unordered(arrivals_calculation,input_list_TARGET)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(PHASES[1]+' Piercing Points estimated!')
print('\n')

# =============================
# Saving Piercing Points TARGET
# =============================

filename_pds_json_TARGET_Ps = sorted(glob.glob(Phase_TARGET_Ps_folder+'*'))

PP_dic_TARGET_Ps_files = [json.load(open(i)) for i in filename_pds_json_TARGET_Ps]

PP_dic_TARGET_Ps = {'depth':[],'time':[],'lat':[],'lon':[]}

for i,j in enumerate(PP_dic_TARGET_Ps_files):
	PP_dic_TARGET_Ps['time'].append(j['time'])
	PP_dic_TARGET_Ps['depth'].append(j['depth'])
	PP_dic_TARGET_Ps['lat'].append(j['lat']) 
	PP_dic_TARGET_Ps['lon'].append(j['lon']) 


print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[1]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_TARGET_Ps, fp)

# ===================================
# Calculating P660s Piercing Points
# ===================================

print('Calculating Piercing Points to '+PHASES[2])
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES) as p:
	max_ = len(input_list_660)
	with tqdm(total=max_,desc='Calculating Piercing Points') as pbar:
		for i, _ in enumerate(p.imap_unordered(arrivals_calculation,input_list_660)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(PHASES[2]+' Piercing Points estimated!')
print('\n')

# =============================
# Saving Piercing Points P660s
# =============================

filename_pds_json_P660s = sorted(glob.glob(Phase_P660s_folder+'*'))

PP_dic_P660s_files = [json.load(open(i)) for i in filename_pds_json_P660s]

PP_dic_P660s = {'depth':[],'time':[],'lat':[],'lon':[]}

for i,j in enumerate(PP_dic_P660s_files):
	PP_dic_P660s['time'].append(j['time'])
	PP_dic_P660s['depth'].append(j['depth'])
	PP_dic_P660s['lat'].append(j['lat']) 
	PP_dic_P660s['lon'].append(j['lon']) 


print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[2]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_P660s, fp)
