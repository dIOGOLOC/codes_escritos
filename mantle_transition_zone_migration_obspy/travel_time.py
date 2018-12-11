#!/usr/bin/python -u

from parameters_py import mgconfig,get_header_data_RF
import warnings
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool
import glob
import time


# =====================================
# Importing travel times scritp_py 
# =====================================

from time_py.time_P_Pds_Ppds import travel_time_calculation_Pds,travel_time_calculation_Ppds

# ==================================================
# Importing some parameters from configuration file 
# ==================================================

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,MP_PROCESSES,OUTPUT_DIR
				   )


# ======================================
# Function to estimate Pds travel times  
# ======================================

def parallel_travel_times_Pds(number,ev_depth,ev_lat,ev_long,st_lat,st_long):
	travel_time_calculation_Pds(number=number,ev_depth=ev_depth,ev_lat=ev_lat,ev_long=ev_long,st_lat=st_lat,st_long=st_long,JSON_FOLDER=travel_times_pds)
	print('Saving travel times Pds number = '+str(number)+' of '+str(len(event_depth)))
	print('\n')

def parallel_travel_times_Ppds(number,ev_depth,ev_lat,ev_long,st_lat,st_long):
	travel_time_calculation_Ppds(number=number,ev_depth=ev_depth,ev_lat=ev_lat,ev_long=ev_long,st_lat=st_lat,st_long=st_long,JSON_FOLDER=travel_times_ppds)
	print('Saving travel times Ppds number = '+str(number)+' of '+str(len(event_depth)))
	print('\n')

# ============================================
# Importing station dictionary from JSON file 
# ============================================

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'/'+'Stations'+'/'

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

# ========================
# Creating output Folder
# ========================

print('Creating travel times output folders')

Travel_times_folder = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'/'+'Phases'+'/'
os.makedirs(Travel_times_folder,exist_ok=True)
print('travel_times_dic = '+Travel_times_folder)

travel_times_pds = Travel_times_folder+'travel_times_pds/'
os.makedirs(travel_times_pds,exist_ok=True)
print('travel_times_pds = '+travel_times_pds)

travel_times_ppds =  Travel_times_folder+'travel_times_ppds/'
print('travel_times_ppds = '+travel_times_ppds)
os.makedirs(travel_times_ppds,exist_ok=True)

print('\n')

# =========================
# Creating Pds Input lists
# ========================

print('Creating Input list')
print('\n')

input_list = [[i+1,event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i]] for i,j in enumerate(event_depth)]

# =============================
# Calculating Pds Travel Times
# =============================

print('Calculating Pds Travel Times')
print('\n')

start_time = time.time()
pool_Pds = Pool(MP_PROCESSES)
pool_Pds.starmap(parallel_travel_times_Pds, input_list)
pool_Pds.close()
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('Pds Travel Times estimated!')
print('\n')

# =============================
# Saving Pds Travel Times
# =============================

# Importing STATION-EVENT Pds dictionary from JSON file 

filename_pds_json = sorted(glob.glob(travel_times_pds+'*'))
filename_pds_json_raw_float = [float(i.split('.')[0].split('Pds_dic_')[1]) for i in filename_pds_json]  
idx = np.argsort(filename_pds_json_raw_float) 	
arrivals = [json.load(open(filename_pds_json[i]))['arrivals'] for i in idx]

TT_dic_Pds = {'time':[],'depth':[],'rayparam':[],'ev_lat': [],'ev_long': [],'st_lat': [],'st_long': []}
for i,j in enumerate(arrivals):
	TT_dic_Pds['time'].append([arrivals[i][k]['time']-arrivals[i][0]['time'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1])
	TT_dic_Pds['depth'].append([arrivals[i][k]['phase'][1:-1] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1])
	TT_dic_Pds['rayparam'].append([arrivals[i][k]['rayparam'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['ev_lat'].append([arrivals[i][k]['ev_lat'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['ev_long'].append([arrivals[i][k]['ev_long'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['st_lat'].append([arrivals[i][k]['st_lat'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['st_long'].append([arrivals[i][k]['st_long'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 


#Saving Piercing Points in JSON file
print('Saving Pds Travel Times in JSON file')
print('\n')

with open(Travel_times_folder+'Pds_dic.json', 'w') as fp:
	json.dump(TT_dic_Pds, fp)

# =============================
# Calculating Ppds Travel Times
# =============================

print('Calculating Ppds Travel Times')
print('\n')

start_time = time.time()
pool_Ppds = Pool(MP_PROCESSES)
pool_Ppds.starmap(parallel_travel_times_Ppds, input_list)
pool_Ppds.close()
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('Ppds Travel Times estimated!')
print('\n')

# =============================
# Saving Ppds Travel Times
# =============================

# Importing STATION-EVENT Pds dictionary from JSON file 

filename_ppds_json = sorted(glob.glob(travel_times_ppds+'*'))
filename_ppds_json_raw_float = [float(i.split('.')[0].split('PPvs_dic_')[1]) for i in filename_ppds_json]  
idx_ppds = np.argsort(filename_ppds_json_raw_float) 	
arrivals = [json.load(open(filename_ppds_json[i]))['arrivals'] for i in idx_ppds]

TT_dic_Ppds = {'time':[],'depth':[],'rayparam':[],'ev_lat': [],'ev_long': [],'st_lat': [],'st_long': []}

for i,j in enumerate(arrivals):
	TT_dic_Ppds['time'].append([arrivals[i][k]['time']-arrivals[i][0]['time'] for k,l in enumerate(j)  if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1])
	TT_dic_Ppds['depth'].append([arrivals[i][k]['phase'][3:-1] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1])
	TT_dic_Ppds['rayparam'].append([arrivals[i][k]['rayparam'] for k,l in enumerate(j)  if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Ppds['ev_lat'].append([arrivals[i][k]['ev_lat'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Ppds['ev_long'].append([arrivals[i][k]['ev_long'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Ppds['st_lat'].append([arrivals[i][k]['st_lat'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Ppds['st_long'].append([arrivals[i][k]['st_long'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 

#Saving Piercing Points in JSON file
print('Saving Ppds Travel Times in JSON file')
print('\n')

with open(Travel_times_folder+'PPvs_dic.json', 'w') as fp:
	json.dump(TT_dic_Ppds, fp)