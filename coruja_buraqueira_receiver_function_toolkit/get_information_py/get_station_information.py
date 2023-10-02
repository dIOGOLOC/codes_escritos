'''
Script to collect information from stations csv file
'''
#Importando módulos

import numpy as np
import os
import pandas as pd

from parameters_py.config import (
					OUTPUT_FEATHER_FILE_DIR,STA_CSV_FILE
				   )



print('Get Station Information')
print('\n')

sta_lat_lon = np.genfromtxt(STA_CSV_FILE,skip_header=1,usecols=[2,3,4],delimiter=';')

sta_name =  np.genfromtxt(STA_CSV_FILE,dtype='str',skip_header=1,usecols=[0,1],delimiter=';')

sta_event = {
		'KNETWK':[],
		'KSTNM':[],
		'STLA':[],
		'STLO':[],
		'STEL':[]
	    }

for i,j in enumerate(sta_name):
	sta_event['STLA'].append(sta_lat_lon[i][0])
	sta_event['STLO'].append(sta_lat_lon[i][1])
	sta_event['STEL'].append(sta_lat_lon[i][2])
	sta_event['KNETWK'].append(j[0])
	sta_event['KSTNM'].append(j[1])


print('Number of Stations: '+str(len(sta_event['KSTNM'])))
for i,j in enumerate(sta_event['KSTNM']):
	print('Station: '+j)
	print('\n')


print('\n')

print('Saving Station Information in FEATHER file')
print('\n')

os.makedirs(OUTPUT_FEATHER_FILE_DIR,exist_ok=True)
df = pd.DataFrame.from_dict(sta_event)
df.to_feather(OUTPUT_FEATHER_FILE_DIR+'STA_dic.feather')
