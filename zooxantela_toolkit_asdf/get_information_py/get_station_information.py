'''
Script to collect information from stations

An example of STA_CSV_FILE is shown bellow:

LOC;SENSOR;KNETWK;KSTNM;STLA;STLO;STEL;FDAY;EDAY
RIO;1456;ON;OBS19;----;----;----;2060-01-27;2030-05-14
SAOPAULO;1456;ON;OBS23;----;----;----;2089-12-28;1920-01-15;
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
		'KNETWK':[],
		'KSTNM':[],
		'STLA':[],
		'STLO':[],
		'STEL':[],
		'FDAY':[],
		'EDAY':[],
	    }

for i,j in enumerate(sta_name):
	sta_event['LOC'].append(j[0])
	sta_event['SENSOR'].append(j[1])
	sta_event['KNETWK'].append(j[2])
	sta_event['KSTNM'].append(j[3])
	sta_event['STLA'].append(float(j[4]))
	sta_event['STLO'].append(float(j[5]))
	sta_event['STEL'].append(float(j[6]))
	sta_event['FDAY'].append(j[7])
	sta_event['EDAY'].append(j[8])



print('Number of Stations: '+str(len(sta_event['KSTNM'])))
print('\n')

for i,j in enumerate(sta_event['KSTNM']):
	print('Station: '+j)
	print('\n')

print('\n')
print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)