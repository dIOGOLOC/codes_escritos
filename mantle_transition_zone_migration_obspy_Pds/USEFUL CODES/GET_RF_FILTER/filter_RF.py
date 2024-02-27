'''
--------------------------------------------------------------------------------
 Function to filter P-wave Receiver Function (PRF)
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 08/2023


Description:
This code will filter and plot the receiver functions estimated.

More information in:
https://seispy.xumijian.me/

Inputs:
Receiver functions files
'''


import os
import glob
import obspy as op
import json
from tqdm import tqdm
from multiprocessing import Pool

import seispy
from seispy.decon import RFTrace

from obspy import read,read_inventory, UTCDateTime, Stream
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace
from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta,plot_trigger

import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

# =====
# paths
# =====

#Directory to get RF data
INPUT_RF_DIR = '/home/sysop/dados_posdoc/MTZ_2024/DATA_MTZ/DATA_2024_YES_PP_FILTER/'

#Directory to save RF data
OUTPUT_RF_DIR = '/home/sysop/dados_posdoc/MTZ_2024/PRF_SEISPY_DATA_YES_PP_FILTER/'

# ======
# inputs
# ======

GAUSSIAN_FILTER = 0.5

CUT_BEFORE = -10

# =======
# quality
# =======

#Trace check (coda amplitude starts in (seconds))
CODA_TRACE_CHECK = 30

#Trace check (coda minimium amplitude of trace)
CODA_TRACE_CHECK_AMP_MIN = -0.3

#Trace check (coda maximum amplitude of trace)
CODA_TRACE_CHECK_AMP_MAX = 0.3

#Trace check (Standard Deviation Multiplier)
CODA_TRACE_CHECK_MULT = 5

#Minimum data amplitude threshold
ZERO_AMP_MIN = 0.1

#Maximum data amplitude threshold 
ZERO_AMP_MAX = 1


# =====================
# Function to save PRF:
# =====================

def RF_calc(input):
    if len(glob.glob(input+'/*')) > 0:

        #name of the RFs
        #%Y.%j.%H.%M.%S_P_[RT].sac
        
        # Decimate to 10 Hz and trim the PRF_R
        rf_R = read(input+'/*_P_R.sac')

        rf_R[0].stats.sac.user1 = GAUSSIAN_FILTER
        rf_R[0].trim(rf_R[0].stats.starttime,rf_R[0].stats.starttime+170)
        rf_R[0].decimate(factor=int(rf_R[0].stats.sampling_rate / 10), strict_length=False, no_filter=True)

        # Decimate to 10 Hz and trim the PRF_T
        rf_T = read(input+'/*_P_T.sac')

        rf_T[0].stats.sac.user1 = GAUSSIAN_FILTER
        rf_T[0].trim(rf_T[0].stats.starttime,rf_T[0].stats.starttime+170)
        rf_T[0].decimate(factor=int(rf_T[0].stats.sampling_rate / 10), strict_length=False, no_filter=True)

        #Check lenght
        if len(rf_R[0].data) > 170*10:
            #RF P-arrival amplitudes check
            P_arrival_start = round((10*10)-10)
            P_arrival_mid = round(10*10)
            P_arrival_end = round((10*10)+10)
            
            amp_mid = rf_R[0].data[P_arrival_mid]
            
            #RF Coda amplitudes check
            amp_Coda = rf_R[0].data[int((CODA_TRACE_CHECK+10)*10):int(170*10)]

            amp_Coda_mask = np.ma.masked_outside(amp_Coda,CODA_TRACE_CHECK_AMP_MIN,CODA_TRACE_CHECK_AMP_MAX)

            mean_amp_Coda = np.mean(amp_Coda_mask)
            std_amp_Coda = abs(np.std(amp_Coda_mask))
            
            if (
                #Minimum data amplitude threshold 
                rf_R[0].data.min() >= -6*ZERO_AMP_MIN and

                #Maximum data amplitude threshold 
                rf_R[0].data.max() <= ZERO_AMP_MAX and

                #Origin amplitude larger than zero
                ZERO_AMP_MIN <= amp_mid <= ZERO_AMP_MAX and all(elem > 0 for elem in rf_R[0].data[P_arrival_start:P_arrival_end])  and

                #Maximum coda amplitude threshold 
                amp_Coda.max() <= mean_amp_Coda + CODA_TRACE_CHECK_MULT*std_amp_Coda and

                #Minimum coda amplitude threshold 
                amp_Coda.min() >= mean_amp_Coda - CODA_TRACE_CHECK_MULT*std_amp_Coda
                ):
                
                os.makedirs(OUTPUT_RF_DIR+'BP.'+rf_R[0].stats.station,exist_ok=True)

                RF_name = rf_R[0].stats.starttime.strftime('%Y.%j.%H.%M.%S')
                    
                rf_R.write(OUTPUT_RF_DIR+'BP.'+rf_R[0].stats.station+'/'+RF_name+'_P_R.sac')

                rf_T.write(OUTPUT_RF_DIR+'BP.'+rf_R[0].stats.station+'/'+RF_name+'_P_T.sac')

    else:
	    pass
    #----------------------------------

# =========================
#  Looking for event folder 
# =========================

PRF_folders = sorted(glob.glob(INPUT_RF_DIR+'*/*/*/*/PRF_'+str(GAUSSIAN_FILTER)+'/'))

# =======
# RF CALC
# =======

with Pool(processes=20) as p:
	max_ = len(PRF_folders)
	with tqdm(total=max_,desc='RF processing') as pbar:
		for i, _ in enumerate(p.imap_unordered(RF_calc,PRF_folders)):
			pbar.update()

print('Finished!')