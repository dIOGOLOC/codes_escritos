'''
Script to collect information from Receiver Functions
'''
#Importing modules

import numpy as np
from obspy import UTCDateTime
import os
import json
import obspy

from parameters_py.config import (
					DIR_SAC,OUTPUT_JSON_FILE_DIR,GAUSSIAN_FILTER,RADIAL_EXT,TRANSVERSAL_EXT,RF_PERCENT,CODA_TRACE_CHECK_AMP,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,SAMPLING_RATE,CODA_TRACE_CHECK,CODA_RATIO_AMP
				   )



print('Get RF Parameters')
print('\n')


#### Radial ####

datalist = []
datalistS = []
for root, dirs, files in os.walk(DIR_SAC):
    for datafile in files:
        if datafile.endswith(RADIAL_EXT):
            datalist.append(os.path.join(root, datafile))

datalistS = sorted(datalist)

#### Transversal ####

datalistT = []
datalistST = []
for root, dirs, files in os.walk(DIR_SAC):
    for datafile in files:
        if datafile.endswith(TRANSVERSAL_EXT):
            datalistT.append(os.path.join(root, datafile))

datalistST = sorted(datalistT)



dic_RF = {
		'dataR':[],
		'dataR_time':[],
		'dataT':[],
		'npts':[],
		'kstnm':[],
		'nzyear':[],
		'nzjday':[],
		'nzhour':[],
		'nzmin':[],
		'nzmsec':[],
		'evla':[],
		'evlo':[],
		'evdp':[],
		'mag':[],
		'stla':[],
		'stlo':[],
		'user0':[],
		'user5':[],
		'user8':[],
		'dist':[],
		'az':[],
		'baz':[],
		'gcarc':[],
		'b':[],
		'e':[]
}

for i,j in enumerate(datalistS):


	data_RF = obspy.read(j)
	data_RF_T = obspy.read(datalistST[i])	

	if len(data_RF[0].data) > CUT_AFTER_P*SAMPLING_RATE:

		print('Station: '+data_RF[0].stats.sac.kstnm)

		#RF P-arrival amplitudes check
		P_arrival_start = round((CUT_BEFORE_P*SAMPLING_RATE)-SAMPLING_RATE)
		P_arrival_mid = round(CUT_BEFORE_P*SAMPLING_RATE)
		P_arrival_end = round((CUT_BEFORE_P*SAMPLING_RATE)+SAMPLING_RATE)

		#RF Coda amplitudes check

		amp_Coda = data_RF[0].data[int((CODA_TRACE_CHECK+CUT_BEFORE_P)*SAMPLING_RATE):int(CUT_AFTER_P*SAMPLING_RATE)]
		std_Coda = CODA_RATIO_AMP*np.std(amp_Coda)



		if data_RF[0].stats.sac.user0 == GAUSSIAN_FILTER and data_RF[0].stats.sac.user5 > RF_PERCENT and min(data_RF[0].data) > CODA_TRACE_CHECK_AMP and data_RF[0].data[P_arrival_start] > 0 and data_RF[0].data[P_arrival_end] > 0 and data_RF[0].data[P_arrival_mid] > 0 and std_Coda > amp_Coda.max() and -std_Coda < amp_Coda.min():
			dic_RF['dataR'].append(data_RF[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)].tolist())
			dic_RF['dataR_time'].append(np.linspace(-CUT_BEFORE_P, CUT_AFTER_P,len(data_RF[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)])).tolist())
			dic_RF['dataT'].append(data_RF_T[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)].tolist())
			dic_RF['npts'].append(float(data_RF[0].stats.sac.npts))
			dic_RF['kstnm'].append(data_RF[0].stats.sac.kstnm)
			dic_RF['nzyear'].append(float(data_RF[0].stats.sac.nzyear))
			dic_RF['nzjday'].append(float(data_RF[0].stats.sac.nzjday))
			dic_RF['nzhour'].append(float(data_RF[0].stats.sac.nzhour))
			dic_RF['nzmin'].append(float(data_RF[0].stats.sac.nzmin))
			dic_RF['nzmsec'].append(float(data_RF[0].stats.sac.nzmsec))
			dic_RF['evla'].append(float(data_RF[0].stats.sac.evla))
			dic_RF['evlo'].append(float(data_RF[0].stats.sac.evlo))
			dic_RF['evdp'].append(float(data_RF[0].stats.sac.evdp))
			dic_RF['mag'].append(float(data_RF[0].stats.sac.mag))
			dic_RF['stla'].append(float(data_RF[0].stats.sac.stla))
			dic_RF['stlo'].append(float(data_RF[0].stats.sac.stlo))
			dic_RF['user0'].append(float(data_RF[0].stats.sac.user0))
			dic_RF['user5'].append(float(data_RF[0].stats.sac.user5))
			dic_RF['user8'].append(float(data_RF[0].stats.sac.user8))
			dic_RF['dist'].append(float(data_RF[0].stats.sac.dist))
			dic_RF['az'].append(float(data_RF[0].stats.sac.az))
			dic_RF['baz'].append(float(data_RF[0].stats.sac.baz))
			dic_RF['gcarc'].append(float(data_RF[0].stats.sac.gcarc))
			dic_RF['b'].append(float(data_RF[0].stats.sac.b))
			dic_RF['e'].append(float(data_RF[0].stats.sac.e))


print('Saving Event Parameters in JSON file')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)

JSON_FILE = 'RF_dic_'+str(GAUSSIAN_FILTER)

JSON_FILE_NAME = JSON_FILE.replace(".", "_")


with open(OUTPUT_JSON_FILE_DIR+JSON_FILE_NAME+'.json', 'w') as fp:
	json.dump(dic_RF, fp)