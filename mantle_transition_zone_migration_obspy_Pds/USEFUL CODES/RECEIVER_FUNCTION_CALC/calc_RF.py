'''
--------------------------------------------------------------------------------
 Function to Calculate a P-wave Receiver Function (PRF)
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 08/2023


Description:
This code will trim and plot the local datase according to a given an event time
and a list of stations.

More information in:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html


Inputs:
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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter,LatitudeLocator,LongitudeLocator


# =====
# paths
# =====

#Directory to save EVENT data
OUTPUT_EV_DIR = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/MTZ_2024/DATA_MTZ/DATA_2024_YES_PP_FILTER/'

# =====
# event
# =====

# ==================
# Calculate PRF
# ==================

# time shift
tshift=10

# Gauss factor
f0 = 0.5

# The maximum number of iterations
itmax = 400

# Minimum error
minderr = 0.001

# ===============================
# Function to cut and plot event:
# ===============================

def RF_calc(input):
    if len(glob.glob(input+'/*')) > 0:
        try:
            st = read(input+'/*')

            st.detrend()

            st.filter("bandpass", freqmin=0.05, freqmax=2.0, zerophase=True)
          
            st_TRZ = st.copy().rotate('NE->RT', back_azimuth=st[0].stats.sac.baz)

            rf_R = RFTrace.deconvolute(st_TRZ[1], st_TRZ[2], method='iter',tshift=tshift, f0=f0, itmax=itmax, minderr=minderr)
            
            rf_T = RFTrace.deconvolute(st_TRZ[0], st_TRZ[2], method='iter',tshift=tshift, f0=f0, itmax=itmax, minderr=minderr)
            
            #name of the RFs
            #%Y.%j.%H.%M.%S_P_[RT].sac
			
            RF_name = st_TRZ[2].stats.starttime.strftime('%Y.%j.%H.%M.%S')
			
            os.makedirs(input+'/PRF_'+str(f0),exist_ok=True)
			
            rf_R.write(input+'/PRF_'+str(f0)+'/'+RF_name+'_P_R.sac')

            rf_T.write(input+'/PRF_'+str(f0)+'/'+RF_name+'_P_T.sac')

        except:
            pass
    else:
	    pass
    #----------------------------------

# =========================
#  Looking for event folder 
# =========================

EVENT_folders = sorted(glob.glob(OUTPUT_EV_DIR+'*/*/*/*'))

# =======
# RF CALC
# =======

with Pool(processes=8) as p:
	max_ = len(EVENT_folders)
	with tqdm(total=max_,desc='RF processing') as pbar:
		for i, _ in enumerate(p.imap_unordered(RF_calc,EVENT_folders)):
			pbar.update()

print('Finished!')