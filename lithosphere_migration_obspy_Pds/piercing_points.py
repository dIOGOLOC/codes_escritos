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
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,MP_PROCESSES,OUTPUT_DIR,DEPTH_MOHO,DEPTH_LAB
				   )
# =====================================
# Function to estimate piercing points  
# =====================================

def parallel_piercing_points(number,PHASE,ev_depth,ev_lat,ev_long,st_lat,st_long,phase_folder):
	arrivals_calculation(number=number,fase=PHASE,ev_depth=ev_depth,ev_lat=ev_lat,ev_long=ev_long,st_lat=st_lat,st_long=st_long,phase_folder=phase_folder)
	print('Event '+str(number)+' of '+str(len(event_depth)))
	print('\n')

# ============================================
# Importing station dictionary from JSON file 
# ============================================

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Stations'+'/'

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

PHASES = 'P'+str(DEPTH_MOHO)+'s','P'+str(DEPTH_LAB)+'s'

# ========================
# Creating output Folder
# ========================

print('Creating output folders')

PP_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Piercing_Points'+'/'

Phase_MOHO_folder = PP_DIR+'P'+str(DEPTH_MOHO)+'s/'
print('Phase_MOHO_folder = '+Phase_MOHO_folder)
os.makedirs(Phase_MOHO_folder,exist_ok=True)

Phase_LAB_folder = PP_DIR+'P'+str(DEPTH_LAB)+'s/'
print('Phase_LAB_folder = '+Phase_LAB_folder)
os.makedirs(Phase_LAB_folder,exist_ok=True)

# =========================
# Creating Pds Input lists
# ========================

print('Creating Pds Input lists')
print('\n')

input_list_MOHO = [(i+1,PHASES[0],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],Phase_MOHO_folder) for i,j in enumerate(event_depth)]

input_list_LAB = [(i+1,PHASES[1],event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],Phase_LAB_folder) for i,j in enumerate(event_depth)]

# ============================================================================================================================
# PDS PHASES
# ============================================================================================================================


# ===================================
# Calculating Moho Piercing Points
# ===================================

print('Calculating Piercing Points to '+PHASES[0])
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES,maxtasksperchild=2) as pool_MOHO:
	pool_MOHO.starmap(parallel_piercing_points, input_list_MOHO,chunksize=1)
pool_MOHO.close()
pool_MOHO.join()
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(PHASES[0]+' Piercing Points estimated!')
print('\n')


# =============================
# Saving Piercing Points Moho
# =============================

filename_pds_json_MOHO = sorted(glob.glob(Phase_MOHO_folder+'*'))

PP_dic_MOHO_files = [json.load(open(i)) for i in filename_pds_json_MOHO]

PP_dic_MOHO = {'depth':[],'time':[],'lat':[],'lon':[]}

for i,j in enumerate(PP_dic_MOHO_files):
	PP_dic_MOHO['time'].append(j['time'])
	PP_dic_MOHO['depth'].append(j['depth'])
	PP_dic_MOHO['lat'].append(j['lat']) 
	PP_dic_MOHO['lon'].append(j['lon']) 


print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[0]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_MOHO, fp)

# =================================
# Calculating LAB Piercing Points
# =================================
 
print('Calculating Piercing Points to '+PHASES[1])
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES,maxtasksperchild=2) as pool_LAB:
	pool_LAB.starmap(parallel_piercing_points, input_list_LAB,chunksize=1)
pool_LAB.close()
pool_LAB.join()
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(PHASES[1]+' Piercing Points estimated!')
print('\n')

# ===========================
# Saving LAB Piercing Points 
# ===========================

filename_pds_json_LAB = sorted(glob.glob(Phase_LAB_folder+'*'))

PP_dic_LAB_files = [json.load(open(i)) for i in filename_pds_json_LAB]

PP_dic_LAB = {'depth':[],'time':[],'lat':[],'lon':[]}

for i,j in enumerate(PP_dic_LAB_files):
	PP_dic_LAB['time'].append(j['time'])
	PP_dic_LAB['depth'].append(j['depth'])
	PP_dic_LAB['lat'].append(j['lat']) 
	PP_dic_LAB['lon'].append(j['lon']) 


print('Saving Piercing Points in JSON file')

with open(PP_DIR+'PP_'+PHASES[1]+'_dic.json', 'w') as fp:
	json.dump(PP_dic_LAB, fp)