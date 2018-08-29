'''
Script to collect information from stations and create the XML file
#Sensor keys names (see http://ds.iris.edu/NRL/)
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

sta_lat_lon = np.genfromtxt(STA_CSV_FILE,skip_header=1,usecols=[1,2,3],delimiter=';')

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[0,4,5],delimiter=';')

sta_event = {

		'KSTNM':[],
		'STLA':[],
		'STLO':[],
		'STEL':[],
		'SENSOR_KEYS':[],
		'DATALOGGER_KEYS':[]
	    }

for i,j in enumerate(sta_name):
	sta_event['STLA'].append(sta_lat_lon[i][0])
	sta_event['STLO'].append(sta_lat_lon[i][1])
	sta_event['STEL'].append(sta_lat_lon[i][2])
	sta_event['KSTNM'].append(j[0])
	sta_event['SENSOR_KEYS'].append(j[1])
	sta_event['DATALOGGER_KEYS'].append(j[2])

print('Number of Stations: '+str(len(sta_event['KSTNM'])))
for i,j in enumerate(sta_event['KSTNM']):
	print('Station: '+j)
	print('SENSOR_KEYS: '+sta_event['SENSOR_KEYS'][i])
	print('DATALOGGER_KEYS: '+sta_event['DATALOGGER_KEYS'][i])
	print('\n')


print('\n')

print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)

