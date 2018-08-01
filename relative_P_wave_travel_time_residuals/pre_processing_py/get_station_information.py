'''
Script to collect information from NEIC (National Earthquake Information Center) 
csv file downloaded in https://earthquake.usgs.gov/earthquakes/search/
'''
#Importando módulos

import numpy as np
import os
import json

from parameters_py.rrconfig import (
					OUTPUT_JSON_FILE_DIR,STA_CSV_FILE
				   )



print('Get Station Information')
print('\n')

sta_info = np.genfromtxt(STA_CSV_FILE,delimiter='\t',skip_header=1)

sta_event = {

		'kstnm':[],
		'stla':[],
		'stlo':[],

	    }


sta_event['kstnm'] = event_info[:,1]
sta_event['stla'] = event_info[:,2]
sta_event['stlo'] = event_info[:,3]

print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'EVENT_dic.json', 'w') as fp:
	json.dump(sta_event, fp)

