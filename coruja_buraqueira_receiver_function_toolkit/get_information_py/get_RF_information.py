'''
Script to collect information from Receiver Functions
'''
#Importing modules

import numpy as np
from obspy import UTCDateTime
import os
import obspy
import shutil
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool

from parameters_py.config import (
					DIR_SAC,OUTPUT_DIR,GAUSSIAN_FILTER,RADIAL_EXT,TRANSVERSAL_EXT,CODA_TRACE_CHECK_AMP_MIN,CODA_TRACE_CHECK_AMP_MAX,
					EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,SAMPLING_RATE,CODA_TRACE_CHECK,ZERO_AMP_MIN,ZERO_AMP_MAX,
					CODA_TRACE_CHECK_MULT
				   )


###########################################
################ FUNCTIONS ################
###########################################

def selection_rf(in_lst):
      
    '''
    Function to select RF waveforms
    '''

    # Reading RF files (Radial and Transversal)

    data_RF = obspy.read(in_lst[0])
    data_RF_T = obspy.read(in_lst[1])

    dir_name_RFR, data_name_RF = os.path.split(in_lst[0])
    dir_name_RFT, data_name_RF_T = os.path.split(in_lst[1])

    if len(data_RF[0].data) > CUT_AFTER_P*SAMPLING_RATE and data_RF[0].stats.sac.mag > EV_MAGNITUDE_MB:

        #RF P-arrival amplitudes check
        P_arrival_start = round((CUT_BEFORE_P*SAMPLING_RATE)-SAMPLING_RATE)
        P_arrival_mid = round(CUT_BEFORE_P*SAMPLING_RATE)
        P_arrival_end = round((CUT_BEFORE_P*SAMPLING_RATE)+SAMPLING_RATE)
		
        amp_mid = data_RF[0].data[P_arrival_mid]
		
        #RF Coda amplitudes check
        amp_Coda = data_RF[0].data[int((CODA_TRACE_CHECK+CUT_BEFORE_P)*SAMPLING_RATE):int(CUT_AFTER_P*SAMPLING_RATE)]

        amp_Coda_mask = np.ma.masked_outside(amp_Coda,CODA_TRACE_CHECK_AMP_MIN,CODA_TRACE_CHECK_AMP_MAX)

        mean_amp_Coda = np.mean(amp_Coda_mask)
        std_amp_Coda = abs(np.std(amp_Coda_mask))

        if (
            #Gaussian Filter
            #data_RF[0].stats.sac.internal0 == GAUSSIAN_FILTER and
                  
            #Minimum data amplitude threshold 
            #data_RF[0].data.min() >= -6*ZERO_AMP_MIN and

            #Maximum data amplitude threshold 
            #data_RF[0].data.max() <= ZERO_AMP_MAX and

            #Origin amplitude larger than zero
            #ZERO_AMP_MIN <= amp_mid <= ZERO_AMP_MAX and all(elem > 0 for elem in data_RF[0].data[P_arrival_start:P_arrival_end]) and
            ZERO_AMP_MIN <= amp_mid <= ZERO_AMP_MAX and all(elem > 0 for elem in data_RF[0].data[P_arrival_start:P_arrival_end])

            #Maximum coda amplitude threshold 
            #amp_Coda.max() <= mean_amp_Coda + CODA_TRACE_CHECK_MULT*std_amp_Coda and

            #Minimum coda amplitude threshold 
            #amp_Coda.min() >= mean_amp_Coda - CODA_TRACE_CHECK_MULT*std_amp_Coda
            ):
                    
                RF_directory = OUTPUT_DIR+'RF_SELECT_STATION_EVENT/'+str(GAUSSIAN_FILTER)+'/'+data_RF[0].stats.sac.knetwk+'.'+data_RF[0].stats.sac.kstnm+'/'
                os.makedirs(RF_directory,exist_ok=True)

                RF_name_RFR = data_name_RF.split(data_RF[0].stats.sac.knetwk+'.'+data_RF[0].stats.sac.kstnm+'.')[1].split('.RFR')[0][:-4]+'_P_R.sac'
                RF_name_RFT = data_name_RF_T.split(data_RF[0].stats.sac.knetwk+'.'+data_RF[0].stats.sac.kstnm+'.')[1].split('.RFT')[0][:-4]+'_P_T.sac'

                shutil.copy2(in_lst[0],RF_directory+RF_name_RFR)
                shutil.copy2(in_lst[1],RF_directory+RF_name_RFT)
                                    
                df = pd.DataFrame([[data_RF[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)].tolist()],
                                        [np.linspace(-CUT_BEFORE_P, CUT_AFTER_P,len(data_RF[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)])).tolist()],
                                        [data_RF_T[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)].tolist()],
                                        [float(data_RF[0].stats.sac.npts)],
                                        [data_RF[0].stats.sac.kstnm],
                                        [float(data_RF[0].stats.sac.nzyear)],
                                        [float(data_RF[0].stats.sac.nzjday)],
                                        [float(data_RF[0].stats.sac.nzhour)],
                                        [float(data_RF[0].stats.sac.nzmin)],
                                        [float(data_RF[0].stats.sac.nzmsec)],
                                        [float(data_RF[0].stats.sac.evla)],
                                        [float(data_RF[0].stats.sac.evlo)],
                                        [float(data_RF[0].stats.sac.evdp)],
                                        [float(data_RF[0].stats.sac.mag)],
                                        [float(data_RF[0].stats.sac.stla)],
                                        [float(data_RF[0].stats.sac.stlo)],
                                        [GAUSSIAN_FILTER],
                                        [float(data_RF[0].stats.sac.dist)],
                                        [float(data_RF[0].stats.sac.az)],
                                        [float(data_RF[0].stats.sac.baz)],
                                        [float(data_RF[0].stats.sac.gcarc)],
                                        [float(data_RF[0].stats.sac.b)],
                                        [float(data_RF[0].stats.sac.e)],
                                        [True]], 
                                        index=['dataR','dataR_time','dataT','npts',
                                              'kstnm','nzyear','nzjday','nzhour',
                                              'nzmin','nzmsec','evla','evlo','evdp',
                                              'mag','stla','stlo','internal0','dist',
                                              'az','baz','gcarc','b','e','selection']).T
        else:     
            df = pd.DataFrame([[data_RF[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)].tolist()],
                                        [np.linspace(-CUT_BEFORE_P, CUT_AFTER_P,len(data_RF[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)])).tolist()],
                                        [data_RF_T[0].data[:int(CUT_AFTER_P*SAMPLING_RATE)].tolist()],
                                        [float(data_RF[0].stats.sac.npts)],
                                        [data_RF[0].stats.sac.kstnm],
                                        [float(data_RF[0].stats.sac.nzyear)],
                                        [float(data_RF[0].stats.sac.nzjday)],
                                        [float(data_RF[0].stats.sac.nzhour)],
                                        [float(data_RF[0].stats.sac.nzmin)],
                                        [float(data_RF[0].stats.sac.nzmsec)],
                                        [float(data_RF[0].stats.sac.evla)],
                                        [float(data_RF[0].stats.sac.evlo)],
                                        [float(data_RF[0].stats.sac.evdp)],
                                        [float(data_RF[0].stats.sac.mag)],
                                        [float(data_RF[0].stats.sac.stla)],
                                        [float(data_RF[0].stats.sac.stlo)],
                                        [GAUSSIAN_FILTER],
                                        [float(data_RF[0].stats.sac.dist)],
                                        [float(data_RF[0].stats.sac.az)],
                                        [float(data_RF[0].stats.sac.baz)],
                                        [float(data_RF[0].stats.sac.gcarc)],
                                        [float(data_RF[0].stats.sac.b)],
                                        [float(data_RF[0].stats.sac.e)],
                                        [False]], 
                                        index=['dataR','dataR_time','dataT','npts',
                                              'kstnm','nzyear','nzjday','nzhour',
                                              'nzmin','nzmsec','evla','evlo','evdp',
                                              'mag','stla','stlo','internal0','dist',
                                              'az','baz','gcarc','b','e','selection']).T
                     
        return df


print('Get RF Parameters')
print('\n')

#### Radial ####

datalist = []
datalistS = []
for root, dirs, files in os.walk(DIR_SAC):
    for datafile in files:
        if RADIAL_EXT+str(GAUSSIAN_FILTER)+'.' in datafile:
            datalist.append(os.path.join(root, datafile))

datalistS = sorted(datalist)

#### Transversal ####

datalistT = []
datalistST = []
for root, dirs, files in os.walk(DIR_SAC):
    for datafile in files:
        if TRANSVERSAL_EXT+str(GAUSSIAN_FILTER)+'.' in datafile:
            datalistT.append(os.path.join(root, datafile))

datalistST = sorted(datalistT)

#### Creatint the input list for multiprocessing ####

input_lst = []
for i,j in enumerate(datalistS):
	input_lst.append([j,datalistST[i]])
      
pandas_RF_lst = []
for i in tqdm(input_lst):
    a = selection_rf(i)
    pandas_RF_lst.append(a)


dataframe_RF_final = pd.concat(pandas_RF_lst, ignore_index=True)

print('Saving RF selecting result:')
print('\n')

os.makedirs(OUTPUT_DIR+'FEATHER/',exist_ok=True)

FEATHER_FILE = 'RF_dic_'+str(GAUSSIAN_FILTER)

FEATHER_FILE_NAME = FEATHER_FILE.replace(".", "_")

dataframe_RF_final.to_feather(OUTPUT_DIR+'FEATHER/'+FEATHER_FILE_NAME+'.feather')

'''

### NEED TO IMPLEMENT #SOME CRAZY ERROR 
print('Multiprocessing the dataset:')
print('\n')

pandas_RF_lst = []

with Pool(processes=1) as p:
    max_ = len(input_lst)
    with tqdm(total=max_) as pbar:
        for result in p.imap(selection_rf,input_lst):
            pbar.update()
            print(result)
            pandas_RF_lst.append(result)

print('Saving RF selecting result:')
print('\n')

dataframe_RF_final = pd.concat(pandas_RF_lst, ignore_index=True)

os.makedirs(OUTPUT_FEATHER_FILE_DIR,exist_ok=True)

FEATHER_FILE = 'RF_dic_'+str(GAUSSIAN_FILTER)

FEATHER_FILE_NAME = FEATHER_FILE.replace(".", "_")

dataframe_RF_final.to_feather(OUTPUT_FEATHER_FILE_DIR+FEATHER_FILE_NAME+'.feather')
'''