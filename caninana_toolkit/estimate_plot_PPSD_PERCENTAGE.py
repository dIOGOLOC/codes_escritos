#!/usr/bin/python -u
'''
Script to estimate and save PPSD plot of the raw data
'''

import os
import json
import obspy
import datetime
import random
from multiprocessing import Pool




# ==============================
# Generating DATA availability
# ==============================

from PSD_py.PSD import calc_PSD


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,NETWORK_CODE,OUTPUT_PSD_DIR,DAY_PERCENTAGE
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

input_list_day = [[]]*len(data)

for i,j in enumerate(data):
	input_list_day[i] = [l['input_list'] for l in j if isinstance(l,dict) == True]


# ==============
#  Get PPSD Data 
# ==============

chosen_files =  [[]]*len(input_list_day)
for i,j in enumerate(input_list_day):
	chosen_files[i] = random.choices(j,k=round((len(j)*DAY_PERCENTAGE)/100))

input_list = [[]]*len(kstnm)
for i,j in enumerate(kstnm):
	input_list[i] = [
					[l,kstnm[i]] for l in chosen_files[i]
					]

print('Get Data for calculating PPSD for each station')
print('\n')
pool_trim = Pool(MP_PROCESSES)
for i,j in enumerate(input_list):
	calc_PSD_result = pool_trim.starmap(calc_PSD, j)
