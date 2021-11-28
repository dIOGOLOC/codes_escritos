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
from matplotlib.widgets import Button, RadioButtons, CheckButtons

import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream, Trace
from obspy.io.xseed import Parser
from obspy.signal.cross_correlation import correlate
from obspy.signal.filter import bandpass,lowpass
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.util import prev_pow_2
from obspy.signal.tf_misfit import cwt
import pywt
import scipy.stats as stats

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
from sklearn.preprocessing import normalize

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

from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta,plot_trigger

# ==================
# Configuration file
# ==================

EARTHQUAKE_FINDER_OUTPUT = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'
#EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'

#ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'
ASDF_FILES = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'

#STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'
STATIONXML_DIR = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

FILTER_DATA = [5,15]

NETWORK = 'ON'

OBS_NAME = 'OBS17'

CHANNEL_LST = ['HHX','HHZ','HHN','HHE']

EVENT_PATTERN_DATE = '2019,08,04,20,46,58'

# ========================
# Constants and parameters
# ========================

WINDOW_LENGTH = 30

STA = 0.5
LTA = 10

THRON = 3
THROFF = 0.5

PEM = 1
PET = 10


DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# ======
# Print?
# ======

VERBOSE = False

# ============
# Main program
# ============

print('=========================')
print('Loading HHZ Event Pattern')
print('=========================')
print('\n')

EVENT_STANDARD_TIME = UTCDateTime(EVENT_PATTERN_DATE)
EVENT_STANDARD_TIME_STR = str(EVENT_STANDARD_TIME.year)+'.'+"%03d" % EVENT_STANDARD_TIME.julday

#obs_HHZ_standard_pattern = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/ON/'+OBS_NAME+'/HHZ.D/ON.'+OBS_NAME+'..HHZ.D.'+EVENT_STANDARD_TIME_STR
obs_HHZ_standard_pattern = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/obs_data_MSEED/ON/'+OBS_NAME+'/HHZ.D/ON.'+OBS_NAME+'..HHZ.D.'+EVENT_STANDARD_TIME_STR

obs_HHZ_standard_pattern_waveform = read(obs_HHZ_standard_pattern)
tr = obs_HHZ_standard_pattern_waveform[0]

tr.trim(EVENT_STANDARD_TIME-PEM,EVENT_STANDARD_TIME+PET)
inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,OBS_NAME,'xml']))
pre_filt = [0.001, 0.005, 45., 50.]
tr.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
tr.detrend('demean')
tr.detrend('linear')
tr.taper(max_percentage=0.1, type="hann")
tr.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])
#----------------------------------------------------------------------------
npts = tr.stats.npts
dt = tr.stats.delta
df = tr.stats.sampling_rate
datetime_window = tr.stats.starttime
EVENT_LENGTH_MIN = tr.stats.endtime-tr.stats.starttime

# ----------------------------------------------------------------------------
# Time Series --> Spectrogram
# ----------------------------------------------------------------------------

tr_data = obs_HHZ_standard_pattern_waveform[0]
trim_data = tr_data.copy()
## Get frequencies corresponding to signal PSD
freq_cwt_lst = fftfreq(npts,d=dt)

## Get positive half of frequencies
positive_frequencies = freq_cwt_lst>0
freq_cwt_lst = freq_cwt_lst[positive_frequencies]

t = trim_data.times('matplotlib')
f_min = FILTER_DATA[0]
f_max = FILTER_DATA[-1]

scalogram_pattern = cwt(trim_data.data, dt, 8, f_min, f_max)

x, y = np.meshgrid(t,np.linspace(f_min, f_max, scalogram_pattern.shape[0]))

# Normalizing data
#event_pattern = normalize(np.abs(scalogram_pattern), axis=1, norm='l2')
event_pattern = np.abs(scalogram_pattern)
# ----------------------------------------------------------------------------
# Spectral Image --> Wavelet Transform
# ----------------------------------------------------------------------------
#
# 2D Discrete Wavelet Transform.
# compress image
#                            -------------------
#                            |        |        |
#                            | cA(LL) | cH(LH) |
#                            |        |        |
# (cA, (cH, cV, cD))  <--->   -------------------
#                            |        |        |
#                            | cV(HL) | cD(HH) |
#                            |        |        |
#                            -------------------
# ------------------------------------------------------------------------------
coeffs = pywt.dwt2(event_pattern, 'db1')
cA, lst_c = coeffs
# normalize each coefficient array independently for better visibility
cA /= np.abs(cA).max()

lst_c_ = []
for coef in lst_c:
    lst_c_.append([d/np.abs(d).max() for d in coef])

cH, cV, cD = lst_c_
haar_image = np.concatenate((np.concatenate((cA, cV),axis= 1),np.concatenate((cH, cD),axis = 1)),axis=0)

# ----------------------------------------------------------------------------
# Wavelet Transform --> Top Coefficients
# ----------------------------------------------------------------------------
# Key discriminative features are concentrated in a few # wavelet coefficients
# with highest deviation:
# – Deviation defined by median/MAD over entire data set
# – Keep only sign (+ or -) of these coefficients, set rest to 0
#  Data compression, robust to noise

shape = haar_image.shape

medians = []
for i in range(shape[1]):
    medians.append(np.median(haar_image[:, i]))
haar_medians = np.array(medians)

mad = []
for i in range(shape[1]):
    tmp = abs(haar_image[:, i] - medians[i])
    mad.append(np.median(tmp))
haar_absdevs = np.array(mad)

haar_image_top = (haar_image - haar_medians)/haar_absdevs
#haar_image_top = np.sign(haar_image_top)

# ----------------------------------------------------------------------------
# Top Coefficients è Binary Fingerprint
# ----------------------------------------------------------------------------
# Fingerprint must be compact and sparse to store in database
# – Convert top coefficients to a binary sequence of 0’s, 1’s
# • Negative: 01, Zero: 00, Positive: 10
# Plotting the results

binaryFingerprints_bool = haar_image_top > 0
binaryFingerprint = 1*binaryFingerprints_bool

# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

axis_major = SecondLocator(interval=5)   # every 5-second
axis_minor = SecondLocator(interval=1)   # every 1-second
axis_Fmt = DateFormatter('%H:%M:%S')

plt.rcParams.update({'font.size': 12})
fig = plt.figure(figsize=(30,30))
gs = fig.add_gridspec(3, 3)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
#----------------------------------------------------------------------------
ax1 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, :],sharex=ax1)
ax3 = fig.add_subplot(gs[2, 0])
ax4 = fig.add_subplot(gs[2, 1])
ax5 = fig.add_subplot(gs[2, 2])
#-----------------------
ax1.set_title('Channel: HHZ')
ax1.xaxis.set_major_formatter(axis_Fmt)
ax1.xaxis.set_major_locator(axis_major)
ax1.xaxis.set_minor_locator(axis_minor)
ax1.tick_params(axis='both',which='major',width=2,length=5)
ax1.tick_params(axis='both',which='minor',width=2,length=3)
ax1.plot(trim_data.times('matplotlib'),trim_data.data, color='k', linewidth=2)
ax1.text(0.05, 0.9, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 0.2, 'boxstyle': 'round'})

#----------------------------------------------------------------------------
ax2.set_title('Continuous Wavelet Transform')
fig.suptitle('Event Date: '+datetime_window.strftime('%d/%m/%Y'),y=0.95,fontsize=25)

ax2.xaxis.set_major_formatter(axis_Fmt)
ax2.xaxis.set_major_locator(axis_major)
ax2.xaxis.set_minor_locator(axis_minor)
ax2.tick_params(axis='both',which='major',width=2,length=5)
ax2.tick_params(axis='both',which='minor',width=2,length=3)

im = ax2.pcolormesh(x, y, event_pattern,shading='auto', cmap='plasma')
ax2.set_ylabel("Frequency [Hz]")
ax2.text(0.05, 0.9, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 0.2, 'boxstyle': 'round'})

axins = inset_axes(ax2,
                    width="30%",
                    height="5%",
                    loc='upper left',
                    bbox_to_anchor=(0.650, 0.1, 1, 1),
                    bbox_transform=ax2.transAxes,
                    borderpad=0,
                    )

plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

#----------------------------------------------------------------------------
ax3.set_title('|Haar transform|')

im3 = ax3.pcolormesh(haar_image, cmap='plasma')
ax3.set_xlabel('wavelet transform x index')
ax3.set_ylabel('wavelet transform y index')
axins3 = inset_axes(ax3,
                    width="20%",
                    height="5%",
                    loc='upper left',
                    bbox_to_anchor=(0.75, 0.1, 1, 1),
                    bbox_transform=ax3.transAxes,
                    borderpad=0,
                    )

plt.colorbar(im3, cax=axins3, orientation="horizontal", ticklocation='top')
#----------------------------------------------------------------------------

ax4.set_title('Top Coefficients')

im4 = ax4.pcolormesh(haar_image_top, cmap='gray')
ax4.set_xlabel('wavelet transform x index')
ax4.set_ylabel('wavelet transform y index')
axins4 = inset_axes(ax4,
                    width="20%",
                    height="5%",
                    loc='upper left',
                    bbox_to_anchor=(0.75, 0.1, 1, 1),
                    bbox_transform=ax4.transAxes,
                    borderpad=0,
                    )

plt.colorbar(im4, cax=axins4, orientation="horizontal", ticklocation='top')

#----------------------------------------------------------------------------

ax5.set_title('Binary Fingerprint')

im5 = ax5.pcolormesh(binaryFingerprint, cmap='gray')
ax5.set_xlabel('wavelet transform x index')
ax5.set_ylabel('wavelet transform y index')
axins5 = inset_axes(ax5,
                    width="20%",
                    height="5%",
                    loc='upper left',
                    bbox_to_anchor=(0.75, 0.1, 1, 1),
                    bbox_transform=ax5.transAxes,
                    borderpad=0,
                    )

plt.colorbar(im5, cax=axins5, orientation="horizontal", ticklocation='top')
#----------------------------------------------------------------------------
plt.show()
'''

EVENT_FOLDER_H5 = [i for i in sorted(glob.glob(ASDF_FILES+'EVENT_DATA_FILES/'+NETWORK+'/**/**')) if OBS_NAME in i]
EVENT_FOLDER_H5 = EVENT_FOLDER_H5[:5]

for data_h5 in tqdm(EVENT_FOLDER_H5,desc='Events loop'):
    _date_str = data_h5.split('_')[-8:-4]
    file_hour =  UTCDateTime(','.join(_date_str))
    file__TIME_STR = str(file_hour.year)+'.'+"%03d" % file_hour.julday
    obs_HHZ_day_files = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/ON/'+OBS_NAME+'/HHZ.D/ON.'+OBS_NAME+'..HHZ.D.'+file__TIME_STR
    obs_HHZ_standard_pattern_waveform = read(obs_HHZ_day_files)

    tr_hour = obs_HHZ_standard_pattern_waveform[0]
    tr_hour.trim(file_hour,file_hour+3600)

    slide_data = [k for k in tr_hour.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH/2,include_partial_windows=False, nearest_sample=True)]

    for tr in slide_data:
        npts = tr.stats.npts
        dt = tr.stats.delta
        df = tr.stats.sampling_rate
        datetime_window = tr.stats.starttime

        tr_data = tr
        trim_data = tr_data.copy()
        inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,OBS_NAME,'xml']))
        pre_filt = [0.001, 0.005, 45., 50.]
        trim_data.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
        trim_data.detrend('demean')
        trim_data.detrend('linear')
        trim_data.taper(max_percentage=0.05, type="hann")
        trim_data.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering        ## Get frequencies corresponding to signal PSD

        #----------------------------------------------------------------------------
        ## Get frequencies corresponding to signal PSD
        freq_cwt_lst = fftfreq(npts,d=dt)

        ## Get positive half of frequencies
        positive_frequencies = freq_cwt_lst>0
        freq_cwt_lst = freq_cwt_lst[positive_frequencies]

        t = trim_data.times('matplotlib')
        f_min = FILTER_DATA[0]
        f_max = FILTER_DATA[-1]

        scalogram = cwt(trim_data.data, dt, 8, f_min, f_max)
        x, y = np.meshgrid(
        t,
        np.linspace(f_min, f_max, scalogram.shape[0]))

        scalogram = cwt(trim_data.data, dt, 8, f_min, f_max)
        scalogram_norm = normalize_matrix(np.abs(scalogram), axis=1, norm='max')

        #---------------------------------------------------------------------------
        # Characteristic function and trigger onsets
        cft = classic_sta_lta(trim_data.data, int(STA*df), int(LTA*df))
        on_off = np.array(trigger_onset(charfct=cft, thres1=THRON, thres2=THROFF))

        #cft = classic_sta_lta(trim_data.data, int(STA * df), int(LTA * df))
        #plot_trigger(trim_data, cft,THRON, THROFF)


        #try:
        for trg in on_off:
                trigger_on = datetime_window+int(trg[0])/df
                trigger_off = datetime_window+int(trg[1])/df

                #----------------------------------------------------------------------------
                # Plotting the results
                axis_major = SecondLocator(interval=5)   # every 5-second
                axis_minor = SecondLocator(interval=1)   # every 1-second
                axis_Fmt = DateFormatter('%H:%M:%S')

                plt.rcParams.update({'font.size': 12})
                fig, axes = plt.subplots(ncols=1, nrows=2,figsize=(20,20),sharex='col')
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
                #----------------------------------------------------------------------------

                ax1 = axes[0]
                ax4 = axes[1]

                ax1.set_title('Channel: HHZ')

                ax1.xaxis.set_major_formatter(axis_Fmt)
                ax1.xaxis.set_major_locator(axis_major)
                ax1.xaxis.set_minor_locator(axis_minor)
                ax1.tick_params(axis='both',which='major',width=2,length=5)
                ax1.tick_params(axis='both',which='minor',width=2,length=3)
                ax1.plot(trim_data.times('matplotlib'),trim_data.data, color='k', linewidth=2)
                ax1.text(0.05, 0.9, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 0.2, 'boxstyle': 'round'})
                ymin, ymax = ax1.get_ylim()
                ax1.vlines(trigger_on.matplotlib_date, ymin, ymax, color='r', linewidth=1)
                ax1.vlines(trigger_off.matplotlib_date, ymin, ymax, color='b', linewidth=1)
                #----------------------------------------------------------------------------
                ax4.set_title('Continuous Wavelet Transform')
                fig.suptitle('Event Date: '+datetime_window.strftime('%d/%m/%Y'),y=0.95,fontsize=25)

                ax4.xaxis.set_major_formatter(axis_Fmt)
                ax4.xaxis.set_major_locator(axis_major)
                ax4.xaxis.set_minor_locator(axis_minor)
                ax4.tick_params(axis='both',which='major',width=2,length=5)
                ax4.tick_params(axis='both',which='minor',width=2,length=3)

                im = ax4.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
                ax4.set_ylabel("Frequency [Hz]")
                ax4.text(0.05, 0.9, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 0.2, 'boxstyle': 'round'})


                axins = inset_axes(ax4,
                                    width="30%",
                                    height="5%",
                                    loc='upper left',
                                    bbox_to_anchor=(0.650, 0.1, 1, 1),
                                    bbox_transform=ax4.transAxes,
                                    borderpad=0,
                                    )

                plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

                plt.show()

                #if trigger_off-trigger_on >= EVENT_LENGTH_MIN:
                    #a
        #except:
            #pass
    #---------------------------------------------------------------------------




    Ashape = event_pattern.shape[1]
    Bshape = scalogram_norm.shape[1]

    if Ashape > Bshape:
        A = event_pattern[:,:Bshape]
        B = scalogram_norm[:,:Bshape]
    else:
        A = event_pattern[:,:Ashape]
        B = scalogram_norm[:,:Ashape]

    diff_A_B = np.subtract(B,A)

    close_matrix = np.mean(np.mean(diff_A_B,axis=1)) < 0.05

    if close_matrix:



processed_EVENT_FOLDER_H5 = []
for data_h5 in tqdm(EVENT_FOLDER_H5,desc='Events loop'):
    EVENT_FILES = sorted(glob.glob(data_h5+'/*.h5'))

for ich in CHANNEL_LST:
    obs_standard_folder_prefix = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/ON/'+OBS_NAME+'/'+ich+'.D/ON.'+OBS_NAME+'..'+ich+'.D.'

    EVENT_STANDARD_TIME = UTCDateTime(EVENT_PATTERN_DATE)
    EVENT_STANDARD_TIME_STR = str(EVENT_STANDARD_TIME.year)+'.'+"%03d" % EVENT_STANDARD_TIME.julday

    #----------------------------------------------------------------------------
    if VERBOSE:
        a = read(obs_standard_folder_prefix+EVENT_STANDARD_TIME_STR)
        a.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering
        a.plot(type="dayplot", interval=60, right_vertical_labels=False,
                vertical_scaling_range=1e4, one_tick_per_line=True,
                color=['k'], show_y_UTC_label=False,
                events={'min_magnitude': 6})
    #----------------------------------------------------------------------------

    obs_waveform = read(obs_standard_folder_prefix+EVENT_STANDARD_TIME_STR)
    tr = obs_waveform[0]

    tr.trim(EVENT_STANDARD_TIME-PEM,EVENT_STANDARD_TIME+PET)
    inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,OBS_NAME,'xml']))
    pre_filt = [0.001, 0.005, 45., 50.]
    tr.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.1, type="hann")
    tr.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering

    #----------------------------------------------------------------------------
    # Plotting the results
    axis_major = SecondLocator(interval=5)   # every 5-second
    axis_minor = SecondLocator(interval=1)   # every 1-second
    axis_Fmt = DateFormatter('%H:%M:%S')

    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(ncols=1, nrows=2,figsize=(20,20),sharex='col')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    #----------------------------------------------------------------------------

    npts = tr.stats.npts
    dt = tr.stats.delta
    df = tr.stats.sampling_rate
    datetime_window = tr.stats.starttime

    #----------------------------------------------------------------------------

    tr_data = obs_waveform[0]
    trim_data = tr_data.copy()
    ## Get frequencies corresponding to signal PSD
    freq_cwt_lst = fftfreq(npts,d=dt)

    ## Get positive half of frequencies
    positive_frequencies = freq_cwt_lst>0
    freq_cwt_lst = freq_cwt_lst[positive_frequencies]

    t = trim_data.times('matplotlib')
    f_min = FILTER_DATA[0]
    f_max = FILTER_DATA[-1]

    scalogram = cwt(trim_data.data, dt, 8, f_min, f_max)
    x, y = np.meshgrid(
    t,
    np.linspace(f_min, f_max, scalogram.shape[0]))

    #---------------------------------------------------------------------------------------------------------------------------

    ax1 = axes[0]
    ax4 = axes[1]

    tr_raw_data = tr
    max_amp = max(tr_raw_data.data)

    ax1.set_title('Channel: '+ich)

    ax1.xaxis.set_major_formatter(axis_Fmt)
    ax1.xaxis.set_major_locator(axis_major)
    ax1.xaxis.set_minor_locator(axis_minor)
    ax1.tick_params(axis='both',which='major',width=2,length=5)
    ax1.tick_params(axis='both',which='minor',width=2,length=3)
    ax1.plot(tr_raw_data.times('matplotlib'),tr_raw_data.data, color='k', linewidth=2)
    ax1.text(0.05, 0.9, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,bbox={'facecolor': 'white', 'alpha': 0.2, 'pad': 1, 'boxstyle': 'round'})

    #----------------------------------------------------------------------------
    ax4.set_title('Continuous Wavelet Transform')
    fig.suptitle('Event Date: '+datetime_window.strftime('%d/%m/%Y'),y=0.95,fontsize=25)

    ax4.xaxis.set_major_formatter(axis_Fmt)
    ax4.xaxis.set_major_locator(axis_major)
    ax4.xaxis.set_minor_locator(axis_minor)
    ax4.tick_params(axis='both',which='major',width=2,length=5)
    ax4.tick_params(axis='both',which='minor',width=2,length=3)

    im = ax4.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
    ax4.set_ylabel("Frequency [Hz]")


    axins = inset_axes(ax4,
                    width="30%",
                    height="5%",
                    loc='upper left',
                    bbox_to_anchor=(0.650, 0.1, 1, 1),
                    bbox_transform=ax4.transAxes,
                    borderpad=0,
                    )

    plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

    plt.show()







for ich in CHANNEL_LST:
    obs_standard_folder_prefix = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/ON/'+OBS_NAME+'/'+ich+'.D/ON.'+OBS_NAME+'..'+ich+'.D.'

    EVENT_STANDARD_TIME = UTCDateTime(EVENT_PATTERN_DATE)
    EVENT_STANDARD_TIME_STR = str(EVENT_STANDARD_TIME.year)+'.'+"%03d" % EVENT_STANDARD_TIME.julday

    #----------------------------------------------------------------------------
    if VERBOSE:
        a = read(obs_standard_folder_prefix+EVENT_STANDARD_TIME_STR)
        a.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering
        a.plot(type="dayplot", interval=60, right_vertical_labels=False,
                vertical_scaling_range=1e4, one_tick_per_line=True,
                color=['k'], show_y_UTC_label=False,
                events={'min_magnitude': 6})
    #----------------------------------------------------------------------------

    obs_waveform = read(obs_standard_folder_prefix+EVENT_STANDARD_TIME_STR)
    tr = obs_waveform[0]

    tr.trim(EVENT_STANDARD_TIME-PEM,EVENT_STANDARD_TIME+PET)
    inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,OBS_NAME,'xml']))
    pre_filt = [0.001, 0.005, 45., 50.]
    tr.remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)
    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.1, type="hann")
    tr.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering

    #----------------------------------------------------------------------------
    # Plotting the results
    axis_major = SecondLocator(interval=5)   # every 5-second
    axis_minor = SecondLocator(interval=1)   # every 1-second
    axis_Fmt = DateFormatter('%H:%M:%S')

    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(ncols=1, nrows=2,figsize=(20,20),sharex='col')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    #----------------------------------------------------------------------------

    npts = tr.stats.npts
    dt = tr.stats.delta
    df = tr.stats.sampling_rate
    datetime_window = tr.stats.starttime

    #----------------------------------------------------------------------------

    tr_data = obs_waveform[0]
    trim_data = tr_data.copy()
    ## Get frequencies corresponding to signal PSD
    freq_cwt_lst = fftfreq(npts,d=dt)

    ## Get positive half of frequencies
    positive_frequencies = freq_cwt_lst>0
    freq_cwt_lst = freq_cwt_lst[positive_frequencies]

    t = trim_data.times('matplotlib')
    f_min = FILTER_DATA[0]
    f_max = FILTER_DATA[-1]

    scalogram = cwt(trim_data.data, dt, 8, f_min, f_max)
    x, y = np.meshgrid(
    t,
    np.linspace(f_min, f_max, scalogram.shape[0]))

    #---------------------------------------------------------------------------------------------------------------------------

    ax1 = axes[0]
    ax4 = axes[1]

    tr_raw_data = tr
    max_amp = max(tr_raw_data.data)

    ax1.set_title('Channel: '+ich)

    ax1.xaxis.set_major_formatter(axis_Fmt)
    ax1.xaxis.set_major_locator(axis_major)
    ax1.xaxis.set_minor_locator(axis_minor)
    ax1.tick_params(axis='both',which='major',width=2,length=5)
    ax1.tick_params(axis='both',which='minor',width=2,length=3)
    ax1.plot(tr_raw_data.times('matplotlib'),tr_raw_data.data, color='k', linewidth=2)
    ax1.text(0.05, 0.9, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,bbox={'facecolor': 'white', 'alpha': 0.2, 'pad': 1, 'boxstyle': 'round'})

    #----------------------------------------------------------------------------
    ax4.set_title('Continuous Wavelet Transform')
    fig.suptitle('Event Date: '+datetime_window.strftime('%d/%m/%Y'),y=0.95,fontsize=25)

    ax4.xaxis.set_major_formatter(axis_Fmt)
    ax4.xaxis.set_major_locator(axis_major)
    ax4.xaxis.set_minor_locator(axis_minor)
    ax4.tick_params(axis='both',which='major',width=2,length=5)
    ax4.tick_params(axis='both',which='minor',width=2,length=3)

    im = ax4.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
    ax4.set_ylabel("Frequency [Hz]")


    axins = inset_axes(ax4,
                    width="30%",
                    height="5%",
                    loc='upper left',
                    bbox_to_anchor=(0.650, 0.1, 1, 1),
                    bbox_transform=ax4.transAxes,
                    borderpad=0,
                    )

    plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

    plt.show()

print('========================')
print('Opening ASDF Event files')
print('========================')
print('\n')

if os.path.isfile(ASDF_FILES+'EVENT_DATA_FILES/processed_folders.npy'):
    EVENT_FOLDER_H5 = np.load(ASDF_FILES+'EVENT_DATA_FILES/processed_folders.npy').tolist()
else:
    EVENT_FOLDER_H5 = sorted(glob.glob(ASDF_FILES+'EVENT_DATA_FILES/'+NETWORK+'/*'))

processed_EVENT_FOLDER_H5 = []
for data_h5 in tqdm(EVENT_FOLDER_H5,desc='Events loop'):
    EVENT_FILES = sorted(glob.glob(data_h5+'/*.h5'))

    # Plotting the results
    axis_major = SecondLocator(interval=5)   # every 5-second
    axis_minor = SecondLocator(interval=1)   # every 1-second
    axis_Fmt = DateFormatter('%H:%M:%S')

    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(ncols=2, nrows=len(EVENT_FILES),figsize=(20,20),sharex='col')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    for i,file_h5 in enumerate(EVENT_FILES):
        event_data = ASDFDataSet(file_h5)
        id_event = event_data.waveforms.list()[0]
        tr = event_data.waveforms[id_event]['event_recording'][0]

        npts = tr.stats.npts
        dt = tr.stats.delta
        df = tr.stats.sampling_rate
        datetime_window = tr.stats.starttime

        #----------------------------------------------------------------------------

        tr_data = tr
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

        ax1 = axes[i,0]
        ax4 = axes[i,1]

        tr_raw_data = tr
        max_amp = max(tr_raw_data.data)

        if i == 0:
            ax1.set_title('Raw Data')

        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5)
        ax1.tick_params(axis='both',which='minor',width=2,length=3)
        ax1.plot(tr_raw_data.times('matplotlib'),tr_raw_data.data, color='k', linewidth=2)
        ax1.text(0.1, 0.85, tr.stats.station, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)

        #----------------------------------------------------------------------------
        if i == 0:
            ax4.set_title('Continuous Wavelet Transform')
            fig.suptitle('Event Date: '+datetime_window.strftime('%d/%m/%Y'),y=0.95,fontsize=25)

        ax4.xaxis.set_major_formatter(axis_Fmt)
        ax4.xaxis.set_major_locator(axis_major)
        ax4.xaxis.set_minor_locator(axis_minor)
        ax4.tick_params(axis='both',which='major',width=2,length=5)
        ax4.tick_params(axis='both',which='minor',width=2,length=3)

        im = ax4.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
        ax4.set_ylabel("Frequency [Hz]")


        axins = inset_axes(ax4,
                       width="40%",
                       height="5%",
                       loc='upper left',
                       bbox_to_anchor=(0.60, 0.1, 1, 1),
                       bbox_transform=ax4.transAxes,
                       borderpad=0,
                       )

        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

    #---BUTTON----
    #Buttons
    NO_ax_button = plt.axes([0.15, 0.925, 0.05,0.05]) #xposition, yposition, width and height
    YES_ax_button = plt.axes([0.25, 0.925, 0.05,0.05]) #xposition, yposition, width and height
    MAYBE_ax_button = plt.axes([0.35, 0.925, 0.05,0.05]) #xposition, yposition, width and height

    #Properties of the button
    NO_button = Button(NO_ax_button, 'NOISE', color = 'red', hovercolor = 'grey')
    YES_button = Button(YES_ax_button, 'EVENT', color = 'green', hovercolor = 'grey')
    MAYBE_button = Button(MAYBE_ax_button, 'MAYBE', color = 'white', hovercolor = 'grey')

    #Saving informations if bottom if clicked:
    def NO_checked(val):
        plt.close(fig)

    def YES_checked(val):
        daily_event_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'/Selected_event_data_GOOD/'
        os.makedirs(daily_event_output,exist_ok=True)
        fig.savefig(daily_event_output+NETWORK+'_'+datetime_window.strftime('%Y_%d_%m_%H_%M_%S_%f')+'_trim.png', dpi=300, facecolor='w', edgecolor='w')
        plt.close(fig)

    def MAYBE_checked(val):
        daily_event_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'/Selected_event_data_MAYBE/'
        os.makedirs(daily_event_output,exist_ok=True)
        fig.savefig(daily_event_output+NETWORK+'_'+datetime_window.strftime('%Y_%d_%m_%H_%M_%S_%f')+'_trim.png', dpi=300, facecolor='w', edgecolor='w')
        plt.close(fig)

    #calling the function when the button gets clicked
    NO_button.on_clicked(NO_checked)
    YES_button.on_clicked(YES_checked)
    MAYBE_button.on_clicked(MAYBE_checked)
    plt.show()

    processed_EVENT_FOLDER_H5.append(data_h5)
    lst3 = [value for value in EVENT_FOLDER_H5 if value not in processed_EVENT_FOLDER_H5]
    np.save(ASDF_FILES+'EVENT_DATA_FILES/processed_folders',lst3)
'''
