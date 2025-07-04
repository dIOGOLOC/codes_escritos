#!/usr/bin/python -u
"""
Script to get raw DATA availability
"""


import os
import json
import obspy
from obspy import UTCDateTime
import datetime
from multiprocessing import Pool
import glob
import tqdm


# ==============================
# Generating DATA availability
# ==============================

from visual_py.data_availability import get_date_file,plot_data_mosaic


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,OUTPUT_FIGURE_DIR
				   )
				   

# ===================================
#  Function to call get data function
# ===================================

def parallel_get_data(lst):
	result = get_date_file(input_list=lst)
	
	if isinstance(result,dict):
		return result

# Importing stations list
print('\n')
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

print('Creating input list with endtime of each raw file')
print('\n')

input_list = [[]]*len(kstnm)

for i,j in enumerate(kstnm):
	datafileS = glob.glob('*/*/*/*/*', root_dir=DIR_DATA)
	datafile_lst = []
	for datafile_name in datafileS:
		if kstnm[i] in datafile_name:
			datafile_lst.append(DIR_DATA+datafile_name)
	
	
	input_list[i] = [
					[l] for l in sorted(datafile_lst)
					]	
	print(' Number of days = '+ str(len(datafile_lst)))		
	print('==============================')

# ==============
#  Get time Data 
# ==============

print('Get time Data  for each station')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
pool = Pool(MP_PROCESSES)
result_dic = {'kstnm':[],'data':[]}
for x,y in enumerate(input_list):
	result_dic['data'].append(pool.starmap(parallel_get_data, tqdm.tqdm(y,total=len(y))))
	result_dic['kstnm'].append(kstnm[x])

with open(OUTPUT_JSON_FILE_DIR+'TIME_dic.json', 'w') as fp:
		json.dump(result_dic, fp)


# ======================
#  Importing time Data 
# ======================


filename_STA = OUTPUT_JSON_FILE_DIR+'TIME_dic.json'

sta_dic = json.load(open(filename_STA))

data = sta_dic['data']

data_time = [[]]*len(data)
for i,j in enumerate(data):
	data_time[i] = [datetime.date(UTCDateTime(l['endtime']).year,UTCDateTime(l['endtime']).month,UTCDateTime(l['endtime']).day) for l in j if isinstance(l,dict) == True]

print('\n')
print('Plotting Data Availability')
print('\n')

plot_data_mosaic(data_time,kstnm)