'''
Script to collect information from NEIC (National Earthquake Information Center) 
csv file downloaded in https://earthquake.usgs.gov/earthquakes/search/
'''
#Importando módulos

import numpy as np
from obspy import UTCDateTime
import os
import json

from parameters_py.rrconfig import (
					OUTPUT_JSON_FILE_DIR,NEIC_CSV_FILE
				   )



print('Get Event Parameters')
print('\n')

event_info = np.genfromtxt(NEIC_CSV_FILE,delimiter=',',skip_header=1,usecols=[1,2,3,4])

#Guardando as variáves cada evento:

ev_time =  np.genfromtxt(NEIC_CSV_FILE,delimiter=',',dtype='str',skip_header=1,usecols=[0])

dic_event = {
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
	dic_event['ev_year'].append('{:02}'.format(temp.year))
	dic_event['ev_month'].append('{:02}'.format(temp.month))
	dic_event['ev_julday'].append('{:02}'.format(temp.julday))
	dic_event['ev_day'].append('{:02}'.format(temp.day))
	dic_event['ev_hour'].append('{:02}'.format(temp.hour))
	dic_event['ev_minute'].append('{:02}'.format(temp.minute))
	dic_event['ev_second'].append('{:02}'.format(temp.second))
	dic_event['ev_microsecond'].append('{:02}'.format(temp.microsecond))
	dic_event['ev_timeUTC'].append(str(temp))

	dic_event['evla'] = event_info[:,0].toalist()
	dic_event['evlo'] = event_info[:,1].toalist()
	dic_event['evdp'] = event_info[:,2].toalist()
	dic_event['mag'] = event_info[:,3].toalist()

print('Saving Event Parameters in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

with open(OUTPUT_JSON_FILE_DIR+'EVENT_dic.json', 'w') as fp:
	json.dump(dic_event, fp)

