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

from visual_py.data_availability import get_date_file_via_client


from parameters_py.config import (
					OUTPUT_FIGURE_DIR,INITIAL_DATE,FINAL_DATE,CHANNEL_CODE,DIR_DATA,
					USER,HOST,PORT,INSTITUTION,NETWORK_CODE,OUTPUT_JSON_FILE_DIR
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

get_date_file_via_client(FIG_FOLDER_OUTPUT=OUTPUT_FIGURE_DIR,STATION_NAME=sta_name)