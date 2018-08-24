'''
Script to collect information from NEIC (National Earthquake Information Center) 
csv file downloaded in https://earthquake.usgs.gov/earthquakes/search/
'''
#Importando módulos

import numpy as np
import os
import json

from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,STA_CSV_FILE
				   )



print('Get Station Information')
print('\n')

sta_info = np.genfromtxt(STA_CSV_FILE,skip_header=1)

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[1])

sensor_type =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[0])

sta_event = {

		'kstnm':[],
		'stla':[],
		'stlo':[],
		'sensor_type':[]

	    }

for i,j in enumerate(sta_info):
	sta_event['sensor_type'].append(sensor_type[i])
	sta_event['kstnm'].append(sta_name[i])
	sta_event['stla'].append(j[2])
	sta_event['stlo'].append(j[3])

print('Number of Stations: '+str(len(sta_event['kstnm'])))
for i,j in enumerate(sta_event['kstnm']):
	print('Station: '+j+' Sensor: '+sta_event['sensor_type'][i])
print('\n')

print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)

