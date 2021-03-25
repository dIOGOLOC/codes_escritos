#!/usr/bin/env python
# coding: utf-8

import time
from tqdm import tqdm
from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
from matplotlib.patches import Ellipse
import matplotlib.cbook as cbook
from matplotlib.patches import Rectangle


import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream, Trace
from obspy.io.xseed import Parser
from obspy.signal.cross_correlation import correlate
from obspy.signal.filter import bandpass,lowpass
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.util import prev_pow_2
from obspy.signal.tf_misfit import cwt
from scipy.stats import moment,kurtosis,skew

import json
import glob
import os
import numpy as np
from numpy.fft import rfft, irfft, fft, ifft, fftfreq
from itertools import combinations
from numpy.lib.stride_tricks import as_strided
import pandas as pd
from scipy.signal import spectrogram, detrend, resample,savgol_filter
from scipy.linalg import norm
from sklearn.preprocessing import normalize as normalize_matrix

import random
import collections
from copy import copy
import datetime
import matplotlib.dates as mdates
from itertools import compress


import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature

from pyasdf import ASDFDataSet

from mtspec import mtspec, wigner_ville_spectrum

from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta


#Configuration file

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_OUTPUT/ASDF_FILES/'

FILTER_DATA = [1,40]

NETWORK = 'ON'

STATION = 'OBS17'

CHANNEL = 'HHZ'

# ========================
# Constants and parameters
# ========================

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)


# ================
# MULTIPROCESSING
# ================
num_processes = 12

# =========
# Functions
# =========


# =======
# Classes
# =======

#Creating ASDF preprocessed files folder
#output_PREPROCESS_DATA_DAY = ASDF_FILES+'PREPROCESS_DATA_DAY_FILES/'+NETWORK+'.'+STATION+'.'+CHANNEL+'/'


#    ds = ASDFDataSet(output_PREPROCESS_DATA_DAY+'PREPROCESS_DATA_DAY_'+NETWORK+'_'+STATION+'_'+channel+'_'+year_day+'_'+julday_day+".h5", compression="gzip-3")
#    ds.add_waveforms(st, tag="preprocessed_recording")
#    ds.add_stationxml(STATIONXML_DIR+'.'.join([NETWORK,STATION,'xml']))

# The type always should be camel case.
#data_type_f = "Frequencies"

# Name to identify the particular piece of data.
#path_f = sta_channel_id+'.f'

# Any additional parameters as a Python dictionary which will end up as
# attributes of the array.
#parameters_f = {'f':'Array of sample frequencies.','sampling_rate_in_hz': st[0].stats.sampling_rate,'station_id': sta_channel_id}

#ds.add_auxiliary_data(data=f, data_type=data_type_f, path=path_f, parameters=parameters_f)

# The type always should be camel case.
#data_type_t = "Times"

# Name to identify the particular piece of data.
#path_t = sta_channel_id+'.t'

# Any additional parameters as a Python dictionary which will end up as
# attributes of the array.
#parameters_t = {'t':'Array of segment times.','sampling_rate_in_hz': st[0].stats.sampling_rate,'station_id': sta_channel_id}


#ds.add_auxiliary_data(data=t, data_type=data_type_t, path=path_t, parameters=parameters_t)

# The type always should be camel case.
#data_type_s = "Spectrogram"

# Name to identify the particular piece of data.
#path_s = sta_channel_id+'.S'

# Any additional parameters as a Python dictionary which will end up as
# attributes of the array.
#parameters_s = {'S':'Spectrogram of x. By default, the last axis of S corresponds to the segment times.','sampling_rate_in_hz': st[0].stats.sampling_rate,'station_id': sta_channel_id}
#ds.add_auxiliary_data(data=psd, data_type=data_type_s, path=path_s, parameters=parameters_s)

# ============
# Main program
# ============

print('========================')
print('Opening ASDF Event files')
print('========================')
print('\n')

EVENT_FILES_H5 = sorted(glob.glob(ASDF_FILES+'EVENT_DATA_FILES/'+NETWORK+'.'+STATION+'.'+CHANNEL+'/*.h5'))

for data_h5 in tqdm(EVENT_FILES_H5,desc='Events loop'):
    event_data = ASDFDataSet(data_h5)

    for i in event_data.waveforms[NETWORK+'.'+STATION]['event_recording']:
        tr = i
        npts = tr.stats.npts
        dt = tr.stats.delta
        df = tr.stats.sampling_rate
        datetime_window = tr.stats.starttime
        #tr.filter(type="bandpass",freqmin=FILTER_DATA[0],freqmax=FILTER_DATA[1],zerophase=True)

        #----------------------------------------------------------------------------

        tr_data = i
        trim_data = tr_data.copy()
        trim_data.filter("bandpass", freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])

        ## Get frequencies corresponding to signal PSD
        freq_cwt_lst = fftfreq(npts,d=dt)

        ## Get positive half of frequencies
        positive_frequencies = freq_cwt_lst>0
        freq_cwt_lst = freq_cwt_lst[positive_frequencies]

        t = trim_data.times('matplotlib')
        f_min = freq_cwt_lst[0]
        f_max = freq_cwt_lst[-1]

        scalogram = cwt(trim_data.data, dt, 8, f_min, f_max)
        x, y = np.meshgrid(
        t,
        np.linspace(f_min, f_max, scalogram.shape[0]))

        #---------------------------------------------------------------------------------------------------------------------------

        # Plotting the results
        axis_major = SecondLocator(interval=5)   # every 5-second
        axis_minor = SecondLocator(interval=1)   # every 1-second
        axis_Fmt = DateFormatter('%H:%M:%S')

        plt.rcParams.update({'font.size': 20})
        fig, axes = plt.subplots(ncols=1, nrows=2,figsize=(20,20),sharex=True)

        ax1 = axes[0]
        ax4 = axes[1]

        tr_raw_data = i
        ax1.set_title('Raw Data')
        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5)
        ax1.tick_params(axis='both',which='minor',width=2,length=3)
        ax1.plot(tr_raw_data.times('matplotlib'),tr_raw_data.data, color='k', linewidth=2)

        #----------------------------------------------------------------------------
        ax4.set_title('Continuous Wavelet Transform')
        ax4.xaxis.set_major_formatter(axis_Fmt)
        ax4.xaxis.set_major_locator(axis_major)
        ax4.xaxis.set_minor_locator(axis_minor)
        ax4.tick_params(axis='both',which='major',width=2,length=5)
        ax4.tick_params(axis='both',which='minor',width=2,length=3)

        im = ax4.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
        ax4.set_ylabel("Frequency [Hz]")
        ax4.set_xlabel('Date: '+datetime_window.strftime('%d/%m/%Y'))

        axins = inset_axes(ax4,
                           width="15%",
                           height="5%",
                           loc='upper left',
                           bbox_to_anchor=(0.8, 0.1, 1, 1),
                           bbox_transform=ax4.transAxes,
                           borderpad=0,
                           )

        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')
        plt.show()
        #daily_event_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'.'+STATION+'/Selected_event_data/'+CHANNEL+'/'
        #os.makedirs(daily_event_output,exist_ok=True)
        #fig.savefig(daily_event_output+NETWORK+'_'+STATION+'_'+CHANNEL+'_'+datetime_window.strftime('%d_%m_%Y_%H_%M_%S_%f')+'_trim.png', dpi=300, facecolor='w', edgecolor='w')
        #plt.close()
