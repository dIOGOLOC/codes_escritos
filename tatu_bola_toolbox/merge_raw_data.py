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
					knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES,DETREND_TYPE,
					TAPER_TYPE,TAPER_MAX_PERCENTAGE,LOWPASS_FREQ,LOWPASS_CORNER,DIR_SAC,
					LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,HIGHPASS_ZEROPHASE,RMEAN_TYPE,
					SAMPLING_RATE,OUTPUT_JSON_FILE_DIR,FILTERS
				   )


# ==================================
#  Function to call cut data script
# ==================================

def parallel_merge_data(folder_name,knetwk,kstnm,year,month,day,NAME_SUFFIX_E,NAME_SUFFIX_N,NAME_SUFFIX_Z,FILTERS,RMEAN_TYPE,DETREND_TYPE,
			TAPER_TYPE,TAPER_MAX_PERCENTAGE,LOWPASS_FREQ,LOWPASS_CORNER,LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,
			HIGHPASS_ZEROPHASE,SAMPLING_RATE):
			
	merge_data_result = merge_data_ZNE(folder_name=folder_name,knetwk=knetwk,kstnm=kstnm,year=year,month=month,day=day,
						NAME_SUFFIX_E=NAME_SUFFIX_E,NAME_SUFFIX_N=NAME_SUFFIX_N,NAME_SUFFIX_Z=NAME_SUFFIX_Z,
						FILTERS=FILTERS,RMEAN_TYPE=RMEAN_TYPE,DETREND_TYPE=DETREND_TYPE,TAPER_TYPE=TAPER_TYPE,
						TAPER_MAX_PERCENTAGE=TAPER_MAX_PERCENTAGE,LOWPASS_FREQ=LOWPASS_FREQ,
						LOWPASS_CORNER=LOWPASS_CORNER,LOWPASS_ZEROPHASE=LOWPASS_ZEROPHASE,
						HIGHPASS_FREQ=HIGHPASS_FREQ,HIGHPASS_CORNER=HIGHPASS_CORNER,
						HIGHPASS_ZEROPHASE=HIGHPASS_ZEROPHASE,SAMPLING_RATE=SAMPLING_RATE)

	return print(merge_data_result)


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

for i in kstnm:
	print('Station = '+i)
print('\n')

# ==============================
#  Creating stations Input lists
# ==============================


print('========================= Searching .SAC files: ========================= ')

datalist = []
datalistS = []
folderlist = []
dirname = []
for root, dirs, files in os.walk(DIR_SAC+'BP02/'):
    for datafile in files:
        if datafile.endswith('.sac'):
            datalist.append(os.path.join(root, datafile))

datalistS = sorted(datalist)

dir_name = [i.split(knetwk+'.')[0] for i in datalistS]

folder_real = sorted(dir_name)

folder_name = sorted(list(set(dir_name)))

STA_name = [i.split('/')[-5] for i in folder_name]
year_name = [i.split('/')[-4] for i in folder_name]
month_name = [i.split('/')[-3] for i in folder_name]
day_name = [i.split('/')[-2] for i in folder_name]


print('Creating stations input lists')
print('\n')


print('Creating input list:')
print('\n')
input_list = 	[
		[folder_name[k],knetwk,STA_name[k],year_name[k],month_name[k],day_name[k],NAME_SUFFIX_E,NAME_SUFFIX_N,NAME_SUFFIX_Z,FILTERS,RMEAN_TYPE,DETREND_TYPE,
			TAPER_TYPE,TAPER_MAX_PERCENTAGE,LOWPASS_FREQ,LOWPASS_CORNER,LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,HIGHPASS_ZEROPHASE,SAMPLING_RATE]
		 for k,l in enumerate(year_name)
		]
print('\n')

# ==============
#  Merging data
# ==============

print('Merging data')
print('\n')

pool_trim = Pool(MP_PROCESSES)
pool_trim.starmap(parallel_merge_data, input_list)
pool_trim.close()

print('Merging finished!')

