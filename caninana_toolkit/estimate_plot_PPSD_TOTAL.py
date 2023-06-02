#!/usr/bin/python -u
"""
Script to estimate and save PPSD data of the raw data
"""
import os
import json
import obspy
import datetime
import random
from multiprocessing import Pool
import tqdm
# ==============================
# Generating DATA availability
# ==============================

from PSD_py.PSD_save import calc_PSD

from visual_py.plot_PSD_DATA import plot_PPSD_TOTAL_data



from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,NETWORK_CODE,OUTPUT_PSD_DIR,DAY_PERCENTAGE
				   )

# Importing stations list

# ============================
# Importing stations datetime
# ============================
print('\n')
print('Looking for STATIONS date in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'TIME_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['kstnm']
data = sta_dic['data']

input_list_day = [[]]*len(data)

for i,j in enumerate(data):
	input_list_day[i] = [l['input_list'] for l in j if isinstance(l,dict) == True]

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	input_list[i] = [
					[l,kstnm[i]] for l in input_list_day[i]
					]

# ==============
#  Get PPSD Data 
# ==============

print('Get Data for calculating PPSD for each station')
print('\n')
pool_trim = Pool(MP_PROCESSES)
for i,j in enumerate(input_list):
	calc_PSD_result = pool_trim.starmap(calc_PSD,tqdm.tqdm(j, total=len(j)))

# ===========================
# Finding stations PPSD data
# ===========================

print('\n')
print('Looking for PPSD STATIONS data in'+OUTPUT_PSD_DIR)
print('\n')

datafile_lst = [] 
for root, dirs, files in os.walk(OUTPUT_PSD_DIR):
	for directories in dirs:
		datafile_name = os.path.join(root, directories)
		if '.PPSD' in datafile_name:
			datafile_lst.append(datafile_name)
datafile_lstS = sorted(datafile_lst)

# ================
#  plot PPSD Data 
# ================

for i,j in enumerate(datafile_lstS):
	plot_PPSD_TOTAL_data(j)