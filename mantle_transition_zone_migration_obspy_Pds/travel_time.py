#!/usr/bin/python -u

from parameters_py import mgconfig,get_header_data_RF
import warnings
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
from multiprocessing import Pool
import glob
import time
from tqdm import tqdm
import pandas as pd
import pyarrow.feather as feather


# =====================================
# Importing travel times scritp_py 
# =====================================

from time_py.time_P_Pds import travel_time_calculation_Pds

# ==================================================
# Importing some parameters from configuration file 
# ==================================================

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,MP_PROCESSES,OUTPUT_DIR,DEPTH_TARGET
				   )

# ==============================================
# Importing station dictionary from FEATHER file 
# ==============================================

STA_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Stations'+'/'

print('\n')
print('Looking for receiver functions data in FEATHER file in '+STA_DIR)
print('\n')


filename_STA = STA_DIR+'sta_dic.feather'

sta_dic = pd.read_feather(filename_STA)  

event_depth = sta_dic['event_depth'].tolist()
event_lat = sta_dic['event_lat'].tolist()
event_long = sta_dic['event_long'].tolist()
event_dist = sta_dic['event_dist'].tolist()
event_gcarc = sta_dic['event_gcarc'].tolist()
event_sta = sta_dic['event_sta'].tolist()
event_ray = sta_dic['event_ray'].tolist()
sta_lat = sta_dic['sta_lat'].tolist()
sta_long = sta_dic['sta_long'].tolist()
sta_data = sta_dic['sta_data'].tolist()
sta_time = sta_dic['sta_time'].tolist()

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

Travel_times_folder = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_TARGET_'+str(DEPTH_TARGET)+'/'+'Phases'+'/'
os.makedirs(Travel_times_folder,exist_ok=True)
print('travel_times_dic = '+Travel_times_folder)

travel_times_pds = Travel_times_folder+'travel_times_pds/'
os.makedirs(travel_times_pds,exist_ok=True)
print('travel_times_pds = '+travel_times_pds)

print('\n')

# =========================
# Creating Pds Input lists
# ========================

print('Creating Input list')
print('\n')

input_list = [[i+1,event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],travel_times_pds] for i,j in enumerate(event_depth)]

# =============================
# Calculating Pds Travel Times
# =============================

print('Calculating Pds Travel Times')
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES) as p:
	max_ = len(input_list)
	with tqdm(total=max_,desc='Calculating Pds travel times') as pbar:
		for i, _ in enumerate(p.imap_unordered(travel_time_calculation_Pds,input_list)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('Pds Travel Times estimated!')
print('\n')

# =============================
# Saving Pds Travel Times
# =============================

# Importing STATION-EVENT Pds dictionary from FEATHER file 

filename_pds_feather = sorted(glob.glob(travel_times_pds+'*'))
filename_pds_feather_raw_float = [float(i.split('.')[0].split('Pds_dic_')[1]) for i in filename_pds_feather]  
idx = np.argsort(filename_pds_feather_raw_float)

arrivals = [pd.read_feather(filename_pds_feather[i])['arrivals'] for i in idx]


TT_dic_Pds = {'time':[],'depth':[],'rayparam':[],'ev_lat': [],'ev_long': [],'st_lat': [],'st_long': []}
for i,j in enumerate(arrivals):
	TT_dic_Pds['time'].append([arrivals[i][k]['time']-arrivals[i][0]['time'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1])
	TT_dic_Pds['depth'].append([arrivals[i][k]['phase'][1:-1] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1])
	TT_dic_Pds['rayparam'].append([arrivals[i][k]['rayparam'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['ev_lat'].append([arrivals[i][k]['ev_lat'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['ev_long'].append([arrivals[i][k]['ev_long'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['st_lat'].append([arrivals[i][k]['st_lat'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 
	TT_dic_Pds['st_long'].append([arrivals[i][k]['st_long'] for k,l in enumerate(j) if arrivals[i][k]['time']-arrivals[i][0]['time'] > 1]) 

#Saving Piercing Points in FEATHER file
print('Saving Pds Travel Times in FEATHER file')
print('\n')

TT_df_Pds = pd.DataFrame.from_dict(TT_dic_Pds)

file_feather_name = Travel_times_folder+'Pds_dic.feather'
feather.write_feather(TT_df_Pds, file_feather_name)