'''
--------------------------------------------------------------------------------
    Function for collecting information of a selected group of stations
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
Given a CSV file in a specific format, this code will return a JSON file that
will be used as input in other programs.


Inputs:
An example of STA_CSV_FILE:
LOC;SENSOR;KNETWK;KSTNM;STLA;STLO;STEL;FDAY;EDAY
RIO;1456;ON;OBS55;----;----;----;2060-01-27;2030-05-14
SAOPAULO;1456;ON;OBS97;----;----;----;2089-12-28;1920-01-15;

Header explanation:
		LOC: Location of the station (str)
		SENSOR: Serial number of the sensor (int)
		KNETWK: Network name (str)
		KSTNM: Network name (str)
		STLA: Latitude of the station (float)
		STLO: Longitude of the station (float)
		STEL: Elevation/Depth of the station (float)
		FDAY: Deployment day - First day (year-month-day)
		EDAY: Recovery day - End day (year-month-day)


Outputs:
JSON file with same structure of the input file


Examples of Usage (in command line):
   >> python get_STATTION_INFORMATION.py

--------------------------------------------------------------------------------
'''
# -----------------
# Importing modules
# -----------------

import numpy as np
import os
import json

from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,STA_CSV_FILE
				   )


if LABEL_LANG == 'br':
    print('\n')
    print('Obtendo Parâmetros das estações')
    print('\n')

else:
    print('\n')
    print('Getting Parameters of the stations')
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


if LABEL_LANG == 'br':
    print('Número de estações: '+str(len(sta_event['KSTNM'])))
    print('\n')

    for i,j in enumerate(sta_event['KSTNM']):
        print('Estação: '+j)
        print('\n')

    print('\n')
    print('Salvando os arquivo das estações (JSON):')
    print('\n')
else:
    print('Number of Stations: '+str(len(sta_event['KSTNM'])))
    print('\n')

    for i,j in enumerate(sta_event['KSTNM']):
        print('Station: '+j)
        print('\n')

    print('\n')
    print('Saving the station files (JSON):')
    print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
with open(OUTPUT_JSON_FILE_DIR+'STA_dic.json', 'w') as fp:
	json.dump(sta_event, fp)
