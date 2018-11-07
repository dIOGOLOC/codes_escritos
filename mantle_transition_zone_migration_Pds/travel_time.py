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
Here we utilize both the P-to-S converted phase (Pds) and the multiply reflected 
and converted phase (Ppds) at the discontinuities to simultaneously determine the 
depth of mantle discontinuities and velocity anomalies in the overlying layer. 
Note that all the parameters are defined in the configuration file.
"""
from parameters_py import mgconfig,get_header_data_RF
import warnings
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool

# =====================================
# Importing travel times scritp_py 
# =====================================

from time_py.time_P_Pds import travel_time_calculation_Pds

# ==================================================
# Importing some parameters from configuration file 
# ==================================================

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,PdS_DIR,
					PP_DIR,PP_SELEC_DIR,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,STA_DIR,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_1_SHP_NAME,BOUNDARY_2_SHP,BOUNDARY_2_SHP_NAME,					
					RAY_PATH_FIGURE,PP_FIGURE,EXT_FIG,DPI_FIG,MP_PROCESSES
				   )


# ======================================
# Function to estimate Pds travel times  
# ======================================

def parallel_travel_times_Pds(number,ev_depth,ev_lat,ev_long,st_lat,st_long):
	TT_dic_Pds = {'dist':[],'depth':[],'time':[],'number':[]}
	travel_time_calculation_Pds_lst = travel_time_calculation_Pds(ev_depth=ev_depth,ev_lat=ev_lat,ev_long=ev_long,st_lat=st_lat,st_long=st_long)
	TT_dic_Pds['time'].append(travel_time_calculation_Pds_lst[0])
	TT_dic_Pds['depth'].append(travel_time_calculation_Pds_lst[1])
	TT_dic_Pds['dist'].append(travel_time_calculation_Pds_lst[2])
	print('Number = '+str(number)+' of '+str(len(event_depth)))
	TT_dic_Pds['number'].append(number)

	return TT_dic_Pds


# ============================================
# Importing station dictionary from JSON file 
# ============================================

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
# Creating Pds Input lists
# ========================

print('Creating Input list')
print('\n')


input_list = [[i,event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

# =============================
# Calculating Pds Travel Times
# =============================

print('Calculating Pds Travel Times')
print('\n')

pool_Pds = Pool(MP_PROCESSES)
TT_dic_Pds = pool_Pds.starmap(parallel_travel_times_Pds, input_list)
pool_Pds.close()
#Saving Piercing Points in JSON file
print('Saving in JSON file')
os.makedirs(PdS_DIR,exist_ok=True)
with open(PdS_DIR+'Pds_dic.json', 'w') as fp:
	json.dump(TT_dic_Pds, fp)

print('Pds Travel Times estimated!')
print('\n')