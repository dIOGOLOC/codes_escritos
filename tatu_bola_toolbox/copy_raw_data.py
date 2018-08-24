#!/usr/bin/python -u
"""
This script copy raw data.
"""
import numpy as np
import obspy
import os
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
import json
from multiprocessing import Pool
import time

start_time = time.time()

# =====================================
# Importing copy data script_py 
# =====================================

from pre_processing_py.copy_convert_data import copy_convert_raw_data

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import ( 
					DIR_RAW_DATA,DIR_SAC,KCMPNM_N,KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,
					NAME_SUFFIX_Z,MP_PROCESSES,OUTPUT_JSON_FILE_DIR
				   )



# ==================================
#  Function to call copy data script
# ==================================

def parallel_copy_data(folder_to_send,name_raw_data,name_station,name_sensor):

	result = copy_convert_raw_data(FOLDER_OUTPUT=folder_to_send,DATA_RAW=name_raw_data,STA_NAME=name_station,SENSOR_TYPE=name_sensor)

	return print(result)


# ============================================
#  Importing station dictionary from JSON file 
# ============================================

print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['kstnm']
stla = sta_dic['stla']
stlo = sta_dic['stlo']
sensor_type = sta_dic['sensor_type']

print('\n')
print('Looking for Data files in '+DIR_RAW_DATA)
print('\n')

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	print('Station = '+j)
	if sensor_type[i] == 'GCF':
		datafile_lst = [] 
		for root, dirs, files in os.walk(DIR_RAW_DATA):
			for datafile in files:
				datafile_name = os.path.join(root, datafile)
				if datafile.endswith('.gcf') and '/'+j+'/' in datafile_name:
					datafile_lst.append(datafile_name)

		datafile_lstS = sorted(datafile_lst)

		print(' Number of files = '+ str(len(datafile_lstS)))
		# ==============================
		#  Creating stations Input lists
		# ==============================

		print('Creating stations input lists')
		print('\n')

		input_list[i] = [
					[DIR_SAC,l,j,sensor_type[i]] for l in datafile_lstS
					]

	if sensor_type[i] == 'MSEED':
		datafile_lst = [] 
		for root, dirs, files in os.walk(DIR_RAW_DATA):
			for datafile in files:
				datafile_name = os.path.join(root, datafile)
				if datafile.endswith('.miniseed') and '/'+j+'/' in datafile_name and 'SOH' not in datafile_name:
					datafile_lst.append(datafile_name)

		datafile_lstS = sorted(datafile_lst)

		print(' Number of files = '+ str(len(datafile_lstS)))

		# ==============================
		#  Creating stations Input lists
		# ==============================

		print('Creating stations input lists')
		print('\n')

		input_list[i] = [
					[DIR_SAC,l,j,sensor_type[i]] for l in datafile_lstS
					]

	if sensor_type[i] == 'REFTEK130':
		datafile_lst = [] 
		for root, dirs, files in os.walk(DIR_RAW_DATA):
			for datafile in files:
				datafile_name = os.path.join(root, datafile)
				if '/'+j+'/' in datafile_name and '/1/' in datafile_name:
					datafile_lst.append(datafile_name)

		datafile_lstS = sorted(datafile_lst)

		print(' Number of files = '+ str(len(datafile_lstS)))

		# ==============================
		#  Creating stations Input lists
		# ==============================

		print('Creating stations input lists')
		print('\n')

		input_list[i] = [
					[DIR_SAC,l,j,sensor_type[i]] for l in datafile_lstS
					]
		print('==============================')

# ==============
#  Copying Data 
# ==============

print('Copying, converting and saving data for each station')
print('\n')

pool = Pool(MP_PROCESSES)
for x,y in enumerate(input_list):
	pool.starmap(parallel_copy_data, y)
pool.close()

print('Time for processing:')
print(str((time.time() - start_time))+" seconds")
print('Copy finished!')
print('Have a nice day (or night)!')