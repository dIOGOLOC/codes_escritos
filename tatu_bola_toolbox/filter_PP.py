#!/usr/bin/python -u
"""
This script reads sac data from a set of stations and trim the data.
"""
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool


# =====================================
# Importing trim data script_py 
# =====================================

from pre_processing_py.filter_PP_seismogram import cut_PP_data_by_event

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
					DIR_SAC,DIR_EVENT,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES
				   )



# ==================================
#  Function to call cut data script
# ==================================

def parallel_trim_data(data_folder,kstnm):
			
	station_data_result = cut_PP_data_by_event(data_folder=data_folder,kstnm=kstnm)

	return station_data_result


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['KSTNM']
stla = sta_dic['STLA']
stlo = sta_dic['STLO']
stel = sta_dic['STEL']
sensor_keys = sta_dic['SENSOR_KEYS']
datalogger_keys = sta_dic['DATALOGGER_KEYS']

# ==============================
#  Creating stations Input lists
# ==============================


print('========================= Searching .SAC files: ========================= ')
print('Looking for EVENT files in '+DIR_EVENT)

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	print('Station = '+kstnm[i])
	datafile_lst = [] 
	for root, dirs, files in os.walk(DIR_EVENT):
		for datadirs in dirs:
			datafile_name = os.path.join(root, datadirs)
			if '/'+kstnm[i]+'/' in datafile_name and len(datafile_name.split('/')) >= 12:
				datafile_lst.append(datafile_name)
	datafile_lstS = sorted(datafile_lst)

	print(' Number of files = '+ str(len(datafile_lstS)))

print('Creating stations input lists')
print('\n')

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	print('Creating input list: '+j)
	print('\n')
	input_list[i] = [
			[l,kstnm[i]]
			 for k,l in enumerate(datafile_lstS)
			]
print('\n')

# ==============
#  Cutting data
# ==============

print('Cutting data and saving events folder for each station')
print('\n')

pool_trim = Pool(MP_PROCESSES)
error_station_data_dic = {'kstnm':[],'error_dir':[]}
for i,j in enumerate(input_list):
	print('Station: '+kstnm[i])
	trim_data_result = pool_trim.starmap(parallel_trim_data, j)
	
	error_station_data_dic['error_dir'].append(pool_trim.starmap(parallel_trim_data, j))
	error_station_data_dic['kstnm'].append(kstnm[i])

pool_trim.close()


#Saving stations data error in JSON file
print('Saving stations data error in JSON file')
with open(OUTPUT_JSON_FILE_DIR+'STA_data_error.json', 'w') as fp:
	json.dump(error_station_data_dic, fp)

print('Cutting finished!')