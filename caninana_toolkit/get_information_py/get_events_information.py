'''
Script to get events information from CATALOG obspy
(https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_events.html)
'''

import numpy as np
from obspy import UTCDateTime
import os
import json
from obspy import read_events
from obspy.clients.fdsn import Client

from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,INITIAL_DATE_EVENT,FINAL_DATE_EVENT,EV_MAGNITUDE_MB
				   )



print('Get Event Parameters')
print('\n')

irisclient=Client("IRIS")

starttime = UTCDateTime(INITIAL_DATE_EVENT)
endtime = UTCDateTime(FINAL_DATE_EVENT)

events = irisclient.get_events(starttime=starttime, endtime=endtime,minmagnitude=EV_MAGNITUDE_MB)

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

for i,j in enumerate(events):
	temp = j['origins'][0]['time']
	dic_event['ev_year'].append('{:04}'.format(temp.year))
	dic_event['ev_month'].append('{:02}'.format(temp.month))
	dic_event['ev_julday'].append('{:03}'.format(temp.julday))
	dic_event['ev_day'].append('{:02}'.format(temp.day))
	dic_event['ev_hour'].append('{:02}'.format(temp.hour))
	dic_event['ev_minute'].append('{:02}'.format(temp.minute))
	dic_event['ev_second'].append('{:02}'.format(temp.second))
	dic_event['ev_microsecond'].append('{:04}'.format(temp.microsecond))
	dic_event['ev_timeUTC'].append(str(temp))
	dic_event['evla'].append(j['origins'][0]['latitude'])
	dic_event['evlo'].append(j['origins'][0]['longitude'])
	dic_event['evdp'].append(j['origins'][0]['depth']/1000)
	dic_event['mag'].append(j['magnitudes'][0]['mag'])

print('Number of Events: '+str(len(dic_event['mag'])))
print('\n')

print('Saving Event Parameters in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

with open(OUTPUT_JSON_FILE_DIR+'EVENT_dic.json', 'w') as fp:
	json.dump(dic_event, fp)