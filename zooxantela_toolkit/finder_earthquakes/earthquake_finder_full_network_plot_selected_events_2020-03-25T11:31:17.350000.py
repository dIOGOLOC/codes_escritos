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

MSEED_DIR_STA = '/home/diogoloc/dados_posdoc/ON_MAR/data/'

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_OUTPUT/EVENTS_SELECTED/'

CHANNEL_LST = ['HHZ.D','HHN.D','HHE.D','HH1.D','HH2.D']

DATE_DAY = '2020-03-25T11:31:17.350000'

FDAY = UTCDateTime(DATE_DAY)
INTERVAL_PERIOD_DATE = str(FDAY.year)+'.'+"%03d" % FDAY.julday

NETWORK = 'ON'

STATIONS_LST = ['DUB01','CAM01','RIB01','GUA01','PET01','TIJ01']

OBS_LST = ['OBS18']

PEM = 30
PET = 120

FILTER_DATA = [1,20]

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
    return day_files

# ============
# Main program
# ============

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

# initializing list of stations by scanning name of miniseed files

lst_INLAND1 = glob.glob(MSEED_DIR_STA+str(FDAY.year)+'/'+NETWORK+'/*')
lst_INLAND = []
for sta in STATIONS_LST:
    a = filelist(basedir=[i for i in lst_INLAND1 if sta in i][0],interval_period_date=INTERVAL_PERIOD_DATE,channel_list=CHANNEL_LST)
    if a != []:
        lst_INLAND.append(a)

lst_OBS1 = glob.glob(MSEED_DIR_OBS+NETWORK+'/*')
lst_OBS = []
for sta in OBS_LST:
    b = filelist(basedir=[i for i in lst_OBS1 if sta in i][0],interval_period_date=INTERVAL_PERIOD_DATE,channel_list=CHANNEL_LST)
    if b != []:
        lst_OBS.append(b)

files_INTERVAL_PERIOD_DATE = lst_OBS+lst_INLAND
print('\n')
print('==============')
print('Finding events')
print('==============')
print('\n')

stE = Stream()
stN = Stream()
stZ = Stream()
for k in files_INTERVAL_PERIOD_DATE:
    for filename in k:
        if 'HHZ.D' in filename:
            stZ += read(filename)
        elif 'HHN.D' in filename or 'HH1.D' in filename:
            stN += read(filename)
        else:
            stE += read(filename)

#===============
#Preprocess Data
#===============
stZ.trim(starttime=FDAY-PEM, endtime=FDAY+PET)

for tr in stZ:
    network = tr.stats.network
    name = tr.stats.station
    channel = tr.stats.channel
    npts = tr.stats.npts
    dt = tr.stats.delta
    df = tr.stats.sampling_rate

    inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))
    pre_filt = [0.001, 0.005, 45., 50.]
    tr.remove_response(inventory=inv,pre_filt=pre_filt,output="VEL",water_level=60)
    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.05)
    tr.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering
stZ2 = stZ.copy()
st_selected_time = stZ2[0].stats.starttime
# Plotting the results
axis_major = SecondLocator(interval=30)   # every 5-second
axis_minor = SecondLocator(interval=5)   # every 1-second
axis_Fmt = DateFormatter('%H:%M:%S')

plt.rcParams.update({'font.size': 15})
fig, axes = plt.subplots(ncols=2, nrows=len(stZ),figsize=(20,20),sharex='col')
fig.suptitle('Event Date: '+FDAY.strftime('%d/%m/%Y'),y=0.95,fontsize=20)

#----------------------------------------------------------------------------

for ind,traces in enumerate(stZ2):
    #----------------------------------------------------------------------------

    ax1 = axes[ind,0]

    trace_amp_max = np.max(np.abs(traces.data))

    if ind == 0:
        ax1.set_title('Filter: '+str(FILTER_DATA[0])+'-'+str(FILTER_DATA[1])+' Hz')

    ax1.xaxis.set_major_locator(axis_major)
    ax1.xaxis.set_major_formatter(axis_Fmt)
    ax1.xaxis.set_minor_locator(axis_minor)
    ax1.tick_params(axis='both',which='major',width=2,length=5)
    ax1.tick_params(axis='both',which='minor',width=2,length=3)
    ax1.plot(traces.times('matplotlib'),traces.data, 'k')
    ax1.set_ylim(-trace_amp_max,trace_amp_max)
    ax1.set_ylabel(traces.stats.station)
    ax1.ticklabel_format(axis='y', style='sci')

    t = traces.times('matplotlib')
    f_min = FILTER_DATA[0]
    f_max = round(df/2)

    scalogram = cwt(traces.data, dt, 8, f_min, f_max)
    x, y = np.meshgrid(t,np.linspace(f_min, f_max, scalogram.shape[0]))

    ax4 = axes[ind,1]
    #ax4.set_title('Continuous Wavelet Transform')
    ax4.xaxis.set_major_formatter(axis_Fmt)
    ax4.xaxis.set_major_locator(axis_major)
    ax4.xaxis.set_minor_locator(axis_minor)
    ax4.tick_params(axis='both',which='major',width=2,length=5)
    ax4.tick_params(axis='both',which='minor',width=2,length=3)

    im = ax4.pcolormesh(x, y, np.abs(scalogram), shading='auto', cmap='viridis')
    ax4.set_ylabel("Frequency [Hz]",fontsize=12)
    ax4.set_ylim(f_min, f_max)

    axins = inset_axes(ax4,
                       width="25%",
                       height="5%",
                       loc='upper left',
                       bbox_to_anchor=(0.75, 0.1, 1, 1),
                       bbox_transform=ax4.transAxes,
                       borderpad=0,
                       )

    plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #----------------------------------------------------------------------------

    daily_event_output = EARTHQUAKE_FINDER_OUTPUT+'/EVENTS_SELECTED_PLOT/'
    os.makedirs(daily_event_output,exist_ok=True)
    fig.savefig(daily_event_output+NETWORK+'_'+st_selected_time.strftime('%Y_%m_%d_%H_%M_%S_%f')+'.png', dpi=300, facecolor='w', edgecolor='w')
    plt.close()
