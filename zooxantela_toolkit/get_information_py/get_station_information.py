'''
Script to collect information from stations

An example of STA_CSV_FILE is shown bellow:

LOC;SENSOR;NAME;LAT;LON;ELEV;FDAY;EDAY
RIO;1456;OBS19;----;----;----;2019-07-27;2020-06-14
SAOPAULO;1456;OBS23;----;----;----;2019-07-28;2020-06-15;
'''
#Importing modules

import numpy as np
import os
import json

from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,STA_CSV_FILE
				   )



print('Get Station Information')
print('\n')

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,delimiter=';')

sta_event = {

		'LOC':[],
		'SENSOR':[],
		'NAME':[],
		'LAT':[],
		'LON':[],
		'ELEV':[],
		'FDAY':[],
		'EDAY':[],
	    }

for i,j in enumerate(sta_name):
	sta_event['LOC'].append(j[0])
	sta_event['SENSOR'].append(j[1])
	sta_event['NAME'].append(j[2])
	sta_event['LAT'].append(float(j[3]))
	sta_event['LON'].append(float(j[4]))
	sta_event['ELEV'].append(float(j[5]))
	sta_event['FDAY'].append(j[6])
	sta_event['EDAY'].append(j[7])



print('Number of Stations: '+str(len(sta_event['NAME'])))
print('\n')

for i,j in enumerate(sta_event['NAME']):
	print('Station: '+j)
	print('\n')

print('\n')
print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)