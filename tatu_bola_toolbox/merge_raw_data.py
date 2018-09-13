#!/usr/bin/python -u
"""
This script reads sac data from a set of stations and merge the data per day.
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

from pre_processing_py.merge_data import merge_data_ZNE

# ==================================================
#  Importing some parameters from configuration file 
# ==================================================

from parameters_py.config import (
					DIR_SAC,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES,OUTPUT_JSON_FILE_DIR,FILE_SUFFIX,FILE_TYPE
				   )


# ==================================
#  Function to call cut data script
# ==================================

def parallel_merge_data(folder_name,knetwk,kstnm,NAME_SUFFIX_E,NAME_SUFFIX_N,NAME_SUFFIX_Z,FILE_FORMAT):
			
	merge_data_result = merge_data_ZNE(folder_name=folder_name,knetwk=knetwk,kstnm=kstnm,
						NAME_SUFFIX_E=NAME_SUFFIX_E,NAME_SUFFIX_N=NAME_SUFFIX_N,NAME_SUFFIX_Z=NAME_SUFFIX_Z,FILE_FORMAT=FILE_TYPE)

	return print(merge_data_result)


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
print('Looking for raw files in '+DIR_SAC)

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	print('Station = '+kstnm[i])
	datafile_lst = [] 
	for root, dirs, files in os.walk(DIR_SAC):
		for datadirs in dirs:
			datafile_name = os.path.join(root, datadirs)
			if '/'+kstnm[i]+'/' in datafile_name and len(datafile_name.split('/')) >= 10:
				datafile_lst.append(datafile_name)
	datafile_lstS = sorted(datafile_lst)

	print(' Number of files = '+ str(len(datafile_lstS)))


	# ==============================
	#  Creating stations Input lists
	# ==============================

	print('Creating stations input lists')
	print('\n')

	input_list[i] = [
					[l,knetwk,kstnm[i],NAME_SUFFIX_E,NAME_SUFFIX_N,NAME_SUFFIX_Z,FILE_TYPE] for l in datafile_lstS
					]

# ==============
#  Merging data
# ==============
print('Merging data')
print('\n')

pool = Pool(MP_PROCESSES)
for x,y in enumerate(input_list):
	pool.starmap(parallel_merge_data, y)
pool.close()

print('Merging finished!')

