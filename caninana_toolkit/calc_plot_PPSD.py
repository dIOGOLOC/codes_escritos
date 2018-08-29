#!/usr/bin/python -u
"""
Scritp to generate DATA availability plots
"""
import os
import json
import obspy
import datetime
from multiprocessing import Pool




# ==============================
# Generating DATA availability
# ==============================

from PSD_py.PSD import calc_PSD


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,NETWORK_CODE
				   )

# ===================================
#  Function to call get data function
# ===================================

def parallel_copy_data(lst):

	result = get_date_file(input_list=lst)

	if result != 1:
		return result


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

data_time = [[]]*len(data)
data_time_day = [[]]*len(data)
input_list_day = [[]]*len(data)

for i,j in enumerate(data):
	data_time[i] = [[l['date_time'],l['input_list']] for l in j if isinstance(l,dict) == True]
	data_time_day[i] = [l['date_time'] for l in j if isinstance(l,dict) == True]
	input_list_day[i] = [l['input_list'] for l in j if isinstance(l,dict) == True]


# ==============
#  Get PPSD Data 
# ==============

print('Get Data for calculating PPSD for each station')
print('\n')

datetime_day = [sorted(list(set(i))) for i in data_time_day]

for i,j in enumerate(input_list_day):
	calc_PSD(j,kstnm[i])

