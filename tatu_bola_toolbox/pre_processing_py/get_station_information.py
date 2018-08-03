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

sta_event = {

		'kstnm':[],
		'stla':[],
		'stlo':[],

	    }

for i,j in enumerate(sta_info):
	sta_event['kstnm'].append(sta_name[i])
	sta_event['stla'].append(j[2])
	sta_event['stlo'].append(j[3])

print('Saving Station Information in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)

