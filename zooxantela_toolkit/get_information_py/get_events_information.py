'''
--------------------------------------------------------------------------------
            Function to collect information of a regional events
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a starttime and endtime, this code will return a JSON file with a list of
events downloaded from Data Centers using OBSPY

More information in:
https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_events.html
and
https://docs.obspy.org/tutorial/code_snippets/retrieving_data_from_datacenters.html

Keep in mind that data centers and web services are constantly changing so this recommendation
might not be valid anymore at the time you read this.

Inputs:
INITIAL_DATE_EVENT: Initial date for looking for events
FINAL_DATE_EVENT: Final date for looking for events
EV_MAGNITUDE_MB: Event magnitude threshold


Outputs:
JSON file with event description:
	ev_timeUTC: event time in UTC (str)
	ev_year: year of the event
	ev_month: month of the event
	ev_day: day of the event
	ev_julday: julian day of the event
	ev_hour: hour of the event
	ev_minute: minute of the event
	ev_second: second of the event
	ev_microsecond: microsecond of the event
	evla: latitude of the event
	evlo: longitude of the event
	evdp: depth of the event
	mag: magnitude of the event
'''

import numpy as np
from obspy import UTCDateTime
import os
import json
from obspy import read_events
from obspy.clients.fdsn import Client

from parameters_py.config import (
                    OUTPUT_JSON_FILE_DIR,INITIAL_DATE_EVENT,FINAL_DATE_EVENT,EV_MAGNITUDE_MB,LABEL_LANG
                    )


if LABEL_LANG == 'br':
    print('Obtendo Parâmetros dos eventos')
    print('\n')

else:
    print('Getting Parameters of the events')
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
            'mag':[]
            }

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

if LABEL_LANG == 'br':
    print('Número de eventos: '+str(len(dic_event['mag'])))
    print('\n')

    print('Salvando os arquivos')
    print('\n')
else:
    print('Number of Events: '+str(len(dic_event['mag'])))
    print('\n')

    print('Saving files.')
    print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

with open(OUTPUT_JSON_FILE_DIR+'EVENT_dic.json', 'w') as fp:
	json.dump(dic_event, fp)
