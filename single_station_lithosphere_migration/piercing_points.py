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

from parameters_py import mgconfig,get_header_data_RF
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
from multiprocessing import Pool
import time
import glob
from tqdm import tqdm
import pyarrow.feather as feather
import pandas as pd
import matplotlib.pyplot as plt

# =====================================
# Importing piercing points scritp_py 
# =====================================

from piercing_points_py.piercing_points_phase import get_pierce_points_geo_calculation

# ==================================================
# Importing some parameters from configuration file 
# ==================================================

from parameters_py.mgconfig import (
					RF_DIR,RF_EXT,MODEL_FILE_NPZ,MIN_DEPTH,MAX_DEPTH,INTER_DEPTH,NUMBER_PP_PER_BIN,
					LLCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLON_LARGE,URCRNRLAT_LARGE,LLCRNRLON_SMALL,
					URCRNRLON_SMALL,LLCRNRLAT_SMALL,URCRNRLAT_SMALL,PROJECT_LAT,PROJECT_LON,Ps_OR_Sp_PHASE,
					BOUNDARY_1_SHP,BOUNDARY_2_SHP,EXT_FIG,DPI_FIG,MP_PROCESSES,OUTPUT_DIR,DEPTH_MOHO,DEPTH_LAB
				   )

# ============================================
# Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for receiver functions data in FEATHER file in '+OUTPUT_DIR)
print('\n')

filename_STA = glob.glob(OUTPUT_DIR+'*/*/sta_dic.feather')[0]
sta_dic = pd.read_feather(filename_STA)  

event_depth = sta_dic['event_depth'].tolist()
event_lat = sta_dic['event_lat'].tolist()
event_long = sta_dic['event_long'].tolist()
event_dist = sta_dic['event_dist'].tolist()
event_gcarc = sta_dic['event_gcarc'].tolist()
event_ray = sta_dic['event_ray'].tolist()
sta_lat = sta_dic['sta_lat'].tolist()
sta_long = sta_dic['sta_long'].tolist()
sta_data = sta_dic['sta_data'].tolist()
sta_time = sta_dic['sta_time'].tolist()
sta_baz = sta_dic['sta_baz'].tolist()
sta_network = sta_dic['sta_knetwk'][0]
sta_station = sta_dic['sta_kstnm'][0]

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

# =======================
# Creating output Folder
# =======================

print('Creating piercing points output folders')

PP_DIR = OUTPUT_DIR+'MODEL_INTER_DEPTH_'+str(INTER_DEPTH)+'_DEPTH_MOHO_'+str(DEPTH_MOHO)+'_DEPTH_LAB_'+str(DEPTH_LAB)+'/'+'Piercing_Points'+'/'

os.makedirs(PP_DIR,exist_ok=True)
print('piercing_points_dic = '+PP_DIR)

piercing_points_Pds_Sdp = PP_DIR+'piercing_points_Pds_Sdp/'
os.makedirs(piercing_points_Pds_Sdp,exist_ok=True)
print('piercing_points_Pds_Sdp = '+piercing_points_Pds_Sdp)

print('\n')

# ====================
# Creating Input lists
# ====================

print('Creating '+Ps_OR_Sp_PHASE+' input lists')
print('\n')


input_list = [[i+1,event_depth[i],event_lat[i],event_long[i],sta_lat[i],sta_long[i],piercing_points_Pds_Sdp] for i,j in enumerate(event_depth)]

# ============================================================================================================================
# PdS_Sdp PHASE
# ============================================================================================================================


# ===================================
# Calculating PdS_Sdp Piercing Points
# ===================================

print('Calculating Piercing Points')
print('\n')

start_time = time.time()
with Pool(processes=MP_PROCESSES, maxtasksperchild=1) as p:
	max_ = len(input_list)
	with tqdm(total=max_,desc='Calculating Piercing Points') as pbar:
		for i, _ in enumerate(p.imap_unordered(get_pierce_points_geo_calculation,input_list)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print(' Piercing Points estimated!')
print('\n')

# ==============================
# Saving Piercing Points PdS_Sdp
# ==============================

filename_pds_json_PdS_Sdp = sorted(glob.glob(piercing_points_Pds_Sdp+'*'))

PP_dic_PdS_Sdp_files = [pd.read_feather(i)for i in filename_pds_json_PdS_Sdp]

PP_dic_PdS_Sdp = {'depth':[],'time':[],'lat':[],'lon':[]}


for i,j in enumerate(PP_dic_PdS_Sdp_files):
	time_reverse = np.abs(j['time'].values - j['time'].values[-1])

	PP_dic_PdS_Sdp['time'].append(time_reverse[::-1][:200].tolist())
	PP_dic_PdS_Sdp['depth'].append(j['depth'][::-1][:200].tolist())
	PP_dic_PdS_Sdp['lat'].append(j['lat'][::-1][:200].tolist()) 
	PP_dic_PdS_Sdp['lon'].append(j['lon'][::-1][:200].tolist()) 

print('Saving Piercing Points in FEATHER file')

PP_dic_PdS_Sdp_df = pd.DataFrame.from_dict(PP_dic_PdS_Sdp)
file_feather_name_PdS_Sdp = PP_DIR+'PP_dic.feather'
feather.write_feather(PP_dic_PdS_Sdp_df, file_feather_name_PdS_Sdp)
