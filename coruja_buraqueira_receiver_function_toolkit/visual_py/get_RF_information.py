'''
Script to collect information from NEIC (National Earthquake Information Center) 
csv file downloaded in https://earthquake.usgs.gov/earthquakes/search/
'''
#Importando módulos

import numpy as np
from obspy import UTCDateTime
import os
import json

from visual_py.config import (
					OUTPUT_JSON_FILE_DIR
				   )



print('Get RF Parameters')
print('\n')

datalist = []
datalistS = []
for root, dirs, files in os.walk(DIR_SAC):
    for datafile in files:
        if datafile.endswith(GAUSSIAN_FILTER):
            datalist.append(os.path.join(root, datafile))

datalistS = sorted(datalist)

dic_RF = {
		'ev_timeUTC':[],
		'ev_year':[],
		'ev_month':[],
		'ev_day':[],
		'ev_julday':[],
		'ev_hour':[],
		'ev_minute':[],
		'ev_second':[],
		'ev_microsecond':[],
		'evla':[],
		'evlo':[],
		'evdp':[],
		'mag':[]}

for i,j in enumerate(ev_time):
	temp = UTCDateTime(j)
	dic_RF['ev_year'].append('{:04}'.format(temp.year))
	dic_RF['ev_month'].append('{:02}'.format(temp.month))
	dic_RF['ev_julday'].append('{:03}'.format(temp.julday))
	dic_RF['ev_day'].append('{:02}'.format(temp.day))
	dic_RF['ev_hour'].append('{:02}'.format(temp.hour))


print('Saving Event Parameters in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

with open(OUTPUT_JSON_FILE_DIR+'RF_dic.json', 'w') as fp:
	json.dump(dic_RF, fp)

