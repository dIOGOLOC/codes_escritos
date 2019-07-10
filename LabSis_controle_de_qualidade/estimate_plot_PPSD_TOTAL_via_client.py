#!/usr/bin/python -u
"""
Script to estimate and save PPSD data of the raw data
"""
import os
import json
import obspy
import datetime
import random



# ==============================
# Generating DATA availability
# ==============================

from PSD_py.PSD_save import calc_PSD_client

from visual_py.plot_PSD_DATA import plot_PPSD_TOTAL_data



from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,NETWORK_CODE,OUTPUT_PSD_DIR,DAY_PERCENTAGE
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
input_list_day = sta_dic['data']["input_list"]

# ==============
#  Get PPSD Data 
# ==============

print('Get Data for calculating PPSD via Client Arclink')
print('\n')

calc_PSD_result = calc_PSD_client(kstnm)




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