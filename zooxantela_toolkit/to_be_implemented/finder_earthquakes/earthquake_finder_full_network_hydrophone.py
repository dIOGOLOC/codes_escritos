#!/usr/bin/env python
# coding: utf-8

from pprint import pprint

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
from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta

# ====================================================================================================
# Configuration file
# ====================================================================================================

MSEED_DIR_OBS = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_HYDROPHONE_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_HYDROPHONE_OUTPUT/ASDF_FILES/'

FIRSTDAY = '2019-08-01'
LASTDAY = '2019-12-31'

FILTER_DATA = [2,16]

NETWORK = 'ON'

CHANNEL = 'HHX'

WINDOW_LENGTH = 600

VERBOSE_MODE = True

STA = 1
LTA = 60

THRON = 3
THROFF = 2

EVENT_LENGTH_MIN = 1.0

PEM = 5
PET = 5

# ========================
# Constants and parameters
# ========================

DTINY = np.finfo(0.0).tiny

WINDOWSEC = datetime.timedelta(seconds=WINDOW_LENGTH)
ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)


# ================
# MULTIPROCESSING
# ================

num_processes = 12

# =================
# Filtering by date
# =================

fday = UTCDateTime(FIRSTDAY)
lday = UTCDateTime(LASTDAY)

INTERVAL_PERIOD = [UTCDateTime(x.astype(str)) for x in np.arange(fday.datetime,lday.datetime+ONEDAY,ONEDAY)]
INTERVAL_PERIOD_DATE = [str(x.year)+'.'+"%03d" % x.julday for x in INTERVAL_PERIOD]

# =========
# Functions
# =========


def filelist(basedir,interval_period_date):
    """
    Returns the list of files in *basedir* whose are in the specified period
    """
    files = []
    files_list = glob.glob(basedir+'/*')
    for s in files_list:
    	if any(day_s in s for day_s in interval_period_date):
    		files.append(s)

    files = [i for i in files if CHANNEL in i]

    return sorted(files)

#-------------------------------------------------------------------------------
def classic_sta_lta_py(a, nsta, nlta):
    """
    Computes the standard STA/LTA from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples. Written in Python by Obspy.
    """
    m = len(a)
    #
    # compute the short time average (STA) and long time average (LTA)
    sta = np.zeros(m, dtype=np.float64)
    lta = np.zeros(m, dtype=np.float64)

    for i in range(m):
        sta[i] = np.mean(a[i:int(i+nsta)]**2)
        lta[i] = np.mean(a[i:int(i+nlta)]**2)

    # Pad zeros
    sta[:nlta - 1] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(0.0).tiny
    idx = lta < dtiny
    lta[idx] = dtiny

    return sta / lta

# =======
# Classes
# =======

def find_events(input_list_FIND_EVENT):
        st = input_list_FIND_EVENT[0]
        st1 = st.copy()

        #===============
        #Preprocess Data
        #===============

        st1 = st.slice(starttime=input_list_FIND_EVENT[1], endtime=input_list_FIND_EVENT[2], keep_empty_traces=False, nearest_sample=True)
        st_good = Stream()
        for tr in st1:
            try:
                network = tr.stats.network
                name = tr.stats.station
                channel = tr.stats.channel
                npts = tr.stats.npts
                dt = tr.stats.delta
                df = tr.stats.sampling_rate

                inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))
                pre_filt = [0.001, 0.005, 45., 50.]
                tr.detrend('demean')
                tr.detrend('linear')
                st_good.append(tr)
            except:
                pass

        st_good.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering
        st2 = st_good.copy()
        trig = coincidence_trigger("recstalta", THRON, THROFF, st2, 3, sta=STA, lta=LTA,details=True)
        for i in trig:
            if any('OBS' in string for string in i['stations']):
                st_selected = Stream()
                for sta_sel in i['stations']:
                    st3 = st2.copy()
                    st4 = st3.select(station=sta_sel)
                    st_selected.append(st4[0])

                    st_selected.trim(i['time']-PEM,i['time']+i['duration']+PET)

                st_selected_time = st_selected[0].stats.starttime

                #----------------------------------------------------------------------------

                ampli_max = []
                for trac in st_selected:
                    ampli_max.append(np.mean(abs(trac.data)))

                if sum(np.array(ampli_max) > 10**-9) > 1:

                    # Plotting the results
                    axis_major = SecondLocator(interval=10)   # every 5-second
                    axis_minor = SecondLocator(interval=1)   # every 1-second
                    axis_Fmt = DateFormatter('%H:%M:%S')

                    plt.rcParams.update({'font.size': 20})
                    fig, axes = plt.subplots(ncols=2, nrows=len(st_selected),figsize=(20,20),sharex='col')
                    fig.suptitle('Date: '+st_selected_time.strftime('%d/%m/%Y'))

                    #----------------------------------------------------------------------------

                    for ind,traces in enumerate(st_selected):
                        try:

                            #Creating ASDF preprocessed files folder
                            output_EVENT_DATA = ASDF_FILES+'EVENT_DATA_FILES/'+NETWORK+'/'+NETWORK+'_'+st_selected_time.strftime('%Y_%m_%d_%H_%M_%S_%f')+'/'
                            os.makedirs(output_EVENT_DATA,exist_ok=True)

                            event_asdf = ASDFDataSet(output_EVENT_DATA+'EVENT_DATA_'+NETWORK+'_'+traces.stats.station+'_'+traces.stats.channel+'_'+traces.stats.starttime.strftime('%Y_%m_%d_%H_%M_%S_%f')+"_event.h5", compression="gzip-3")
                            tr_2_save = traces
                            event_asdf.add_waveforms(tr_2_save, tag="event_recording")
                            #----------------------------------------------------------------------------

                            ax1 = axes[ind,0]

                            trace_amp_max = np.max(np.abs(traces.data))

                            ax1.xaxis.set_major_locator(axis_major)
                            ax1.xaxis.set_major_formatter(axis_Fmt)
                            ax1.xaxis.set_minor_locator(axis_minor)
                            ax1.tick_params(axis='both',which='major',width=2,length=5)
                            ax1.tick_params(axis='both',which='minor',width=2,length=3)
                            ax1.plot(traces.times('matplotlib'),traces.data, 'k')
                            ax1.set_ylim(-trace_amp_max,trace_amp_max)
                            ax1.set_ylabel(traces.stats.station)
                            ax1.ticklabel_format(axis='y', style='sci')

                            #----------------------------------------------------------------------------

                            dt = traces.stats.delta
                            df = traces.stats.sampling_rate

                            t = traces.times('matplotlib')
                            f_min = FILTER_DATA[0]
                            f_max = round(df/2)

                            scalogram = cwt(traces.data, dt, 8, f_min, f_max)
                            x, y = np.meshgrid(t,np.linspace(f_min, f_max, scalogram.shape[0]))

                            ax2 = axes[ind,1]
                            #ax2.set_title('Continuous Wavelet Transform')
                            ax2.xaxis.set_major_formatter(axis_Fmt)
                            ax2.xaxis.set_major_locator(axis_major)
                            ax2.xaxis.set_minor_locator(axis_minor)
                            ax2.tick_params(axis='both',which='major',width=2,length=5)
                            ax2.tick_params(axis='both',which='minor',width=2,length=3)

                            im = ax2.pcolormesh(x, y, np.abs(scalogram), shading='auto', cmap='viridis')
                            ax2.set_ylabel("Frequency [Hz]")
                            ax2.set_ylim(f_min, f_max)

                            axins = inset_axes(ax2,
                                               width="40%",
                                               height="5%",
                                               loc='upper left',
                                               bbox_to_anchor=(0.6, 0.1, 1, 1),
                                               bbox_transform=ax2.transAxes,
                                               borderpad=0,
                                              )

                            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

                        #----------------------------------------------------------------------------
                        except:
                            pass

                    daily_event_output = EARTHQUAKE_FINDER_OUTPUT+'/Daily_event_data_windows/'+'/'+st_selected_time.strftime('%Y-%m-%d')+'/'
                    os.makedirs(daily_event_output,exist_ok=True)
                    fig.savefig(daily_event_output+NETWORK+'_'+st_selected_time.strftime('%Y_%m_%d_%H_%M_%S_%f')+'.png', dpi=300, facecolor='w', edgecolor='w')
                    plt.close()

# ============
# Main program
# ============

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

# initializing list of stations by scanning name of miniseed files

lst_OBS = glob.glob(MSEED_DIR_OBS+NETWORK+'/*')
print('Total of OBS stations = '+str(len(lst_OBS)))

files_OBS = []
for i,j in enumerate(lst_OBS):
    files_OBS.append(filelist(basedir=j+'/'+CHANNEL+'.D/',interval_period_date=INTERVAL_PERIOD_DATE))

files_INTERVAL_PERIOD_DATE = []
for i,j in enumerate(INTERVAL_PERIOD_DATE):
    temp2 = []
    for k in files_OBS:
        temp1 = [w for w in k if j in w]
        temp2.append(temp1)
    temp = filter(lambda x: len(x) > 0, temp2)
    flat_list = [item for sublist in temp for item in sublist]
    files_INTERVAL_PERIOD_DATE.append(list(flat_list))


print('\n')
print('==============')
print('Finding events')
print('==============')
print('\n')

start_time = time.time()
for i,j in enumerate(files_INTERVAL_PERIOD_DATE):
    print('Dia: '+INTERVAL_PERIOD[i].strftime('%d/%m/%Y'))
    st = Stream()
    for filename in j:
        st += read(filename)

    input_list_FIND_EVENT = []

    INTERVAL_PERIOD_WINDOW_START = [UTCDateTime(x.astype(str)) for x in np.arange(INTERVAL_PERIOD[i].datetime,INTERVAL_PERIOD[i].datetime+ONEDAY-WINDOWSEC/2,WINDOWSEC/2)]
    INTERVAL_PERIOD_WINDOW_END = [UTCDateTime(x.astype(str)) for x in np.arange(INTERVAL_PERIOD[i].datetime+WINDOWSEC,INTERVAL_PERIOD[i].datetime+ONEDAY+WINDOWSEC/2,WINDOWSEC/2)]

    for k,l in enumerate(INTERVAL_PERIOD_WINDOW_START):
        input_list_FIND_EVENT.append([st,INTERVAL_PERIOD_WINDOW_START[k],INTERVAL_PERIOD_WINDOW_END[k]])

    with Pool(processes=num_processes) as p:
    	max_ = len(input_list_FIND_EVENT)
    	with tqdm(total=max_, desc='Slice loop') as pbar:
    		for i, _ in enumerate(p.imap_unordered(find_events, input_list_FIND_EVENT)):
    			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
