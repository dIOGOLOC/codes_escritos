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
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, MicrosecondLocator, DateFormatter
from matplotlib.patches import Ellipse
import matplotlib.cbook as cbook
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

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
from obspy.signal.trigger import plot_trigger

# ====================================================================================================
# Configuration file
# ====================================================================================================

MSEED_DIR_OBS = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

MSEED_DIR_STA = '/home/diogoloc/dados_posdoc/ON_MAR/data/'

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/Figuras/EVENTS_SELECTED/'

CHANNEL_LST = ['HHX.D']

DATE_DAY = '2019-12-07T12:18:40.00'

FDAY = UTCDateTime(DATE_DAY)
INTERVAL_PERIOD_DATE = str(FDAY.year)+'.'+"%03d" % FDAY.julday

NETWORK = 'ON'

STATIONS_LST = []

OBS_LST = ['OBS17','OBS18','OBS20','OBS22']

PEM = 10
PET = 60

#Low-cut frequency (Hz) and High-cut frequency (Hz) for bandpass filter
FILTER_DATA = [2,20]

#Length of the short-term average window (seconds)
stalen = 0.5

#Length of the long-term average window (seconds)
ltalen = 10

#sta/lta ratio to trigger a detection/pick
trig_on = 2

#sta/lta ratio to turn the trigger off - no further picks\
trig_off = 2

#Show picks on waveform.
VERBOSE_MODE = True

#-------------------------------------------------------------------------------

def filelist(basedir,interval_period_date,channel_list):
    """
    Returns the list of files in *basedir* whose are in the specified period
    """
    files = []
    files_list = glob.glob(basedir+'/*')
    files_list_ch = []
    for s in files_list:
        if any(day_s in s for day_s in channel_list):
            files_list_ch.append(s)
    day_files = []
    for ch_folder in files_list_ch:
        files = glob.glob(ch_folder+'/*')
        date_file = [file for file in files if interval_period_date in file]
        if date_file != []:
            day_files.append(date_file[0])
    return sorted(day_files)
#-------------------------------------------------------------------------------
def smooth(data, nd, axis=0):
    """
    Function to smooth power spectral density functions from the convolution
    of a boxcar function with the PSD

    Parameters
    ----------
    data : :class:`~numpy.ndarray`
        Real-valued array to smooth (PSD)
    nd : int
        Number of samples over which to smooth
    axis : int
        axis over which to perform the smoothing
    Returns
    -------
    filt : :class:`~numpy.ndarray`, optional
        Filtered data
    """
    if np.any(data):
        if data.ndim > 1:
            filt = np.zeros(data.shape)
            for i in range(data.shape[::-1][axis]):
                if axis == 0:
                    filt[:, i] = np.convolve(
                        data[:, i], np.ones((nd,))/nd, mode='same')
                elif axis == 1:
                    filt[i, :] = np.convolve(
                        data[i, :], np.ones((nd,))/nd, mode='same')
        else:
            filt = np.convolve(data, np.ones((nd,))/nd, mode='same')
        return filt
    else:
        return None
#-------------------------------------------------------------------------------
# ============
# Main program
# ============

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

# initializing list of stations by scanning name of miniseed files

lst_OBS1 = glob.glob(MSEED_DIR_OBS+NETWORK+'/*')
lst_OBS = []
for sta in OBS_LST:
    b = filelist(basedir=[i for i in lst_OBS1 if sta in i][0],interval_period_date=INTERVAL_PERIOD_DATE,channel_list=CHANNEL_LST)
    if b != []:
        lst_OBS.append(b)

files_INTERVAL_PERIOD_DATE = lst_OBS

print('\n')
print('==============')
print('Finding events')
print('==============')
print('\n')

st = Stream()
for k in files_INTERVAL_PERIOD_DATE:
    for filename in k:
        st += read(filename)

#===============
#Preprocess Data
#===============
st.trim(starttime=FDAY-PEM, endtime=FDAY+PET)

stX = Stream()
for tr in st:

    network = tr.stats.network
    name = tr.stats.station
    channel = tr.stats.channel
    npts = tr.stats.npts
    dt = tr.stats.delta
    df = tr.stats.sampling_rate

    if channel == 'HHX':
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(max_percentage=0.05)
        tr.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering
        stX += tr
stX2 = stX.copy()

st_selected_time = st[0].stats.starttime

print('\n')
print('===============')
print('Plotting events')
print('===============')
print('\n')

# Plotting the results
axis_major = SecondLocator(interval=30)   # every 5-second
axis_minor = SecondLocator(interval=5) # every 1-second
axis_Fmt = DateFormatter('%H:%M:%S')

#----------------------------------------------------------------------------

plt.rcParams.update({'font.size': 15})
fig = plt.figure(figsize=(30,15))
gs = gridspec.GridSpec(nrows=3, ncols=len(stX2))
fig.suptitle('Event Date: '+FDAY.strftime('%d/%m/%Y'),y=0.95,fontsize=20)
#----------------------------------------------------------------------------
for ind,traces in enumerate(stX2):
        ax1 = fig.add_subplot(gs[0,ind])
        trace_amp_max = np.max(np.abs(traces.data))

        ax1.set_title(traces.stats.station)
        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5)
        ax1.tick_params(axis='both',which='minor',width=2,length=3)
        ax1.plot(traces.times('matplotlib'),traces.data, 'k')
        ax1.set_ylim(-trace_amp_max,trace_amp_max)
        ax1.text(0.1, 0.85, traces.stats.channel, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
        ax1.ticklabel_format(axis='y', style='sci')
        if ind == 0:
            ax1.set_ylabel("Amplitude")
        plt.setp(ax1.get_xticklabels(), visible=False)

        #----------------------------------------------------------------------------

        t = traces.times('matplotlib')
        f_min = FILTER_DATA[0]
        f_max = FILTER_DATA[1]

        scalogram = cwt(traces.data, dt, 8, f_min, f_max)
        x, y = np.meshgrid(t,np.linspace(f_min, f_max, scalogram.shape[0]))

        ax3 = fig.add_subplot(gs[1,ind],sharex=ax1)
        ax3.xaxis.set_major_formatter(axis_Fmt)
        ax3.xaxis.set_major_locator(axis_major)
        ax3.xaxis.set_minor_locator(axis_minor)
        ax3.tick_params(axis='both',which='major',width=2,length=5)
        ax3.tick_params(axis='both',which='minor',width=2,length=3)


        im = ax3.pcolormesh(x, y, np.abs(scalogram), shading='auto', cmap='viridis')
        if ind == 0:
            ax3.set_ylabel("Frequency [Hz]")
        ax3.set_ylim(f_min, f_max)
        plt.setp(ax3.get_xticklabels(), visible=False)

        axins = inset_axes(ax3,
                           width="40%",
                           height="5%",
                           loc='upper left',
                           bbox_to_anchor=(0.60, 0.1, 1, 1),
                           bbox_transform=ax3.transAxes,
                           borderpad=0,
                           )

        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

        #-------------------------------------------------------------------------------------------------------------------------------
        ax5 = fig.add_subplot(gs[2,ind],sharex=ax1)

        cft = classic_sta_lta(traces.data, int(stalen * df), int(ltalen * df))
        triggers = trigger_onset(cft, trig_on, trig_off)
        t = np.arange(npts, dtype=np.float32)/df

        on_off = np.array(trigger_onset(cft, trig_on, trig_off))
        time_on = traces.stats.starttime+float((on_off[:,0]/df)[0])
        time_off = traces.stats.starttime+float((on_off[:,1]/df)[0])

        ax1.axvline(time_on.matplotlib_date, color='blue', lw=1, ls='--')
        ax1.axvline(time_off.matplotlib_date, color='red', lw=1, ls='--')

        ax3.axvline(time_on.matplotlib_date, color='blue', lw=1, ls='--')
        ax3.axvline(time_off.matplotlib_date, color='red', lw=1, ls='--')

        ax5.plot(traces.times('matplotlib'), cft, 'k')
        ax5.axhline(trig_on, color='red', lw=1, ls='--')
        ax5.axhline(trig_off, color='blue', lw=1, ls='--')
plt.show()
#-------------------------------------------------------------------------------
daily_event_output = EARTHQUAKE_FINDER_OUTPUT+'/EVENTS_SELECTED_PLOT_STA/'
os.makedirs(daily_event_output,exist_ok=True)
fig.savefig(daily_event_output+NETWORK+'_'+traces.stats.starttime.strftime('%Y_%m_%d_%H_%M_%S_%f')+'_hydrophone_only.png', dpi=300, facecolor='w', edgecolor='w')
plt.close()
#-------------------------------------------------------------------------------
