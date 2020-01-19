'''
Script to collect information from Receiver Functions
'''
#Importing modules

import numpy as np
from obspy import UTCDateTime
import os
import json
import obspy
import shutil

from parameters_py.config import (
					DIR_SAC,OUTPUT_JSON_FILE_DIR,GAUSSIAN_FILTER,RADIAL_EXT,TRANSVERSAL_EXT,RF_PERCENT,CODA_TRACE_CHECK_AMP_MIN,CODA_TRACE_CHECK_AMP_MAX,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,SAMPLING_RATE,CODA_TRACE_CHECK,DIR_SEL_SAC,ZERO_AMP_MIN,
					CODA_TRACE_MAX_AMP,CODA_TRACE_MIN_AMP
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

	data_name = j.split('/')[-1]


	if len(data_RF[0].data) > CUT_AFTER_P*SAMPLING_RATE and data_RF[0].stats.sac.mag > EV_MAGNITUDE_MB:

		print('Station: '+data_RF[0].stats.sac.kstnm)

		#RF P-arrival amplitudes check
		P_arrival_start = round((CUT_BEFORE_P*SAMPLING_RATE)-SAMPLING_RATE)
		P_arrival_mid = round(CUT_BEFORE_P*SAMPLING_RATE)
		P_arrival_end = round((CUT_BEFORE_P*SAMPLING_RATE)+SAMPLING_RATE)
		
		amp_mid = data_RF[0].data[P_arrival_mid]
		#RF Coda amplitudes check
		amp_Coda = data_RF[0].data[int((CODA_TRACE_CHECK+CUT_BEFORE_P)*SAMPLING_RATE):int(CUT_AFTER_P*SAMPLING_RATE)]

		mean_amp_Coda = np.mean(amp_Coda)
		std_amp_Coda = np.std(amp_Coda)
		
		#Gaussian Filter
		if data_RF[0].stats.sac.user0 == GAUSSIAN_FILTER and data_RF[0].stats.sac.user5 > RF_PERCENT:
			#Reconstruction value 
			#data_RF[0].stats.sac.user5 > RF_PERCENT and
			
			#Minimum amplitude threshold
			#min(data_RF[0].data) > CODA_TRACE_CHECK_AMP_MIN and

			#Maximum amplitude threshold 
			#max(data_RF[0].data) < CODA_TRACE_CHECK_AMP_MAX and

			#Maximum coda amplitude threshold 
			#and amp_Coda.max() < CODA_TRACE_MAX_AMP
	 		#amp_Coda.max() < mean_amp_Coda + 2*std_amp_Coda and

			#Minimum coda amplitude threshold 
			#amp_Coda.min() > mean_amp_Coda - 2*std_amp_Coda and

			#Origin amplitude larger than zero
			#amp_mid > ZERO_AMP_MIN and all(elem > 0 for elem in data_RF[0].data[P_arrival_start:P_arrival_end]):
				
				RF_directory = DIR_SEL_SAC+str(GAUSSIAN_FILTER)+'/'+data_RF[0].stats.sac.kstnm+'/'
				os.makedirs(RF_directory,exist_ok=True)
				shutil.copy2(j,RF_directory+data_name)

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
