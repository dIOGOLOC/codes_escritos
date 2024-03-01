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
INPUT_RF_DIR = '/home/sysop/dados_posdoc/MTZ_2024/DATA_MTZ/DATA_2024_SYNTHETIC/'

#Directory to save RF data
OUTPUT_RF_DIR = '/home/sysop/dados_posdoc/MTZ_2024/PRF_SEISPY_DATA_YES_PP_FILTER_SYNTHETIC/'

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