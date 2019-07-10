#!/usr/bin/python -u
"""
Script to get raw DATA availability
"""


import os
import json
import obspy
import datetime
from multiprocessing import Pool



# ==============================
# Generating DATA availability
# ==============================

from visual_py.data_availability import get_date_file_via_client,plot_data_mosaic


from parameters_py.config import (
					CHANNEL_CODE,OUTPUT_JSON_FILE_DIR,DIR_DATA,OUTPUT_FIGURE_DIR,OUTPUT_PSD_DIR
				   )
				   

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

		
# ==============
#  Get time Data 
# ==============

print('Get time Data  for each station')
print('\n')

for k,l in enumerate(kstnm):
	if l in DIR_DATA:
		sta_name = l 

print('STATION NAME = '+sta_name)

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

result_dic = {'kstnm':[],'data':[]}
result_dic['data'] = get_date_file_via_client(directory_data=OUTPUT_PSD_DIR+'2019'+'/'+sta_name+'/'+CHANNEL_CODE+'.PPSD/')
result_dic['kstnm'] = sta_name

with open(OUTPUT_JSON_FILE_DIR+'TIME_dic.json', 'w') as fp:
		json.dump(result_dic, fp)


# ======================
#  Importing time Data 
# ======================


filename_STA = OUTPUT_JSON_FILE_DIR+'TIME_dic.json'

sta_dic = json.load(open(filename_STA))

data_starttime = sta_dic['data']['starttime']

data_time = []
for i,j in enumerate(data_starttime):
	data_time += [datetime.date(obspy.UTCDateTime(j).year,obspy.UTCDateTime(j).month,obspy.UTCDateTime(j).day)]

print('\n')
print('Plotting Data Availability')
print('\n')

plot_data_mosaic(data_time,sta_name,CHANNEL_CODE)