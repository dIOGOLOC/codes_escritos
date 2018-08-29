#!/usr/bin/python -u
"""
Scritp to generate DATA availability plots
"""
import os
import json
import obspy
import datetime
from multiprocessing import Pool
from obspy import UTCDateTime
import numpy as np


# ==============================
# Generating DATA availability
# ==============================

from visual_py.data_availability import plot_data_mosaic


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,OUTPUT_FIGURE_DIR
				   )

# Importing stations list
print('\n')
print('Looking for STATIONS date in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'TIME_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['kstnm']
data = sta_dic['data']

data_time = [[]]*len(data)
for i,j in enumerate(data):
	data_time[i] = [datetime.date(UTCDateTime(l['endtime']).year,UTCDateTime(l['endtime']).month,UTCDateTime(l['endtime']).day) for l in j if isinstance(l,dict) == True]

print('\n')
print('Plotting Data Availability')
print('\n')

plot_data_mosaic(data_time,kstnm)