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
import pyarrow.feather as feather
from scipy.signal import spectrogram, detrend, resample,savgol_filter
from skimage.transform import rescale, resize, downscale_local_mean
from scipy.linalg import norm
from sklearn.preprocessing import normalize
from sklearn.metrics import jaccard_score

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

# ========================
# Input and output folders
# ========================

#MSEED_FOLDER = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/obs_data_MSEED/'
MSEED_FOLDER = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

#EARTHQUAKE_FINDER_OUTPUT = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'
EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'
#ASDF_FILES = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'
#STATIONXML_DIR = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

BINARY_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/'
#BINARY_FILES = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/'

# ==========
# Parameters
# ==========

#Type of event:
EVENT_TYPE = 'airgun'

#Bandpass frequency (Hz) - minimum and maximum:
<<<<<<< HEAD
FILTER_DATA = [5,15]

NETWORK = 'ON'

OBS_NAME = 'OBS17'
=======
FILTER_DATA = [10,40]

NETWORK = 'ON'

OBS_NAME = 'OBS22'
>>>>>>> master

CHANNEL = 'HHZ'

#Datetime of the event: string with year,month,day,hour,minute,second:
<<<<<<< HEAD
EVENT_PATTERN_DATE = '2019,08,04,20,46,58' #airgun
#EVENT_PATTERN_DATE = '2020,02,28,19,53,20'
=======
#EVENT_PATTERN_DATE = '2019,08,04,20,46,58' #airgun
EVENT_PATTERN_DATE = '2019,12,07,12,19,40'
>>>>>>> master

#Time before the event starttime
PEM = 2

#Time after the event starttime
<<<<<<< HEAD
PET = 15
=======
PET = 3
>>>>>>> master

#Spectral image length (samples)
spectral_length = 32

#Number of wavelet coefficients to keep
k_coef = 200

# =========
# Constants
# =========

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# ======
# Print?
# ======

VERBOSE = True

# =========
# Functions
# =========

def fingerprint_creator(obspy_trace,freqmin,freqmax,fp_length):
    tr = obspy_trace
    #----------------------------------------------------------------------------

    inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,OBS_NAME,'xml']))
    pre_filt = [0.001, 0.005, 45., 50.]
    tr.remove_response(inventory=inv,output="DISP",pre_filt=pre_filt,water_level=60)
    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.2, type="hann")
    tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)

    npts = tr.stats.npts
    dt = tr.stats.delta
    df = tr.stats.sampling_rate
    datetime_window = tr.stats.starttime

    # --------------------------------------------------------------------------
    # Time Series --> Spectrogram
    # --------------------------------------------------------------------------

    tr_data = tr
    trim_data = tr_data.copy()

    t = trim_data.times('matplotlib')
    f_min = freqmin
    f_max = freqmax

    scalogram_pattern = cwt(st=trim_data.data, dt=dt, w0=5, fmin=f_min, fmax=f_max,wl='morlet')
    # Normalizing data
    event_pattern = normalize(np.abs(scalogram_pattern), axis=1, norm='max')

    # Resizing data
    new_event_pattern = resize(event_pattern, (fp_length, fp_length),anti_aliasing=False)

    # --------------------------------------------------------------------------
    # Spectral Image --> Wavelet Transform
    # --------------------------------------------------------------------------
    # 2D Discrete Wavelet Transform.
    # compress image
    # --------------------------------------------------------------------------

    coeffs = pywt.dwt2(new_event_pattern, 'haar')
    # reshapes output from PyWavelets 2d wavelet transform

    cA = coeffs[0]
    (cH, cV, cD) = coeffs[1]

    cA = normalize(np.abs(cA), axis=1, norm='max')
    cH = normalize(np.abs(cH), axis=1, norm='max')
    cV = normalize(np.abs(cV), axis=1, norm='max')
    cD = normalize(np.abs(cD), axis=1, norm='max')

    haar_image = np.concatenate((np.concatenate((cA, cH),axis= 1),np.concatenate((cV, cD),axis = 1)),axis=0)

    # ----------------------------------------------------------------------------
    # Wavelet Transform --> Top Coefficients --> Binary Fingerprint
    # ----------------------------------------------------------------------------
    # Top Coefficients are concentrated in features with value 95% of the maximum:
    # – Keep only sign (+ or -) of these coefficients, set rest to 0
    #  Data compression, robust to noise
    # Fingerprint must be compact and sparse to store in database
    # – Convert top coefficients to a binary sequence of 0’s, 1’s

    # ----------------------------------------------------------------------------
    # Wavelet Transform --> Top Coefficients
    # ----------------------------------------------------------------------------
    # Key discriminative features are concentrated in a few wavelet coefficients with
    # highest deviation
    # – Deviation defined by median/MAD over entire data set
    # – Keep only sign (+ or -) of these coefficients, set rest to 0
    # • Data compression, robust to noise

    #Top Coefficients --> Mean Absolute Deviation of a matrix (using PANDAS)
    haar_image_top_pd = pd.DataFrame(haar_image).mad().to_numpy()
    haar_image_top = haar_image - haar_image_top_pd

    #binarize vectors:
    K = k_coef
    N,M = np.shape(haar_image_top)
    binary_vectors = np.zeros((N,2*M), dtype=bool)
    for i in range(N):
        idx = np.argsort(abs(haar_image_top[i,:]))[-K:]
        binary_vectors[i,idx]   = haar_image_top[i,idx] > 0
        binary_vectors[i,idx+M] = haar_image_top[i,idx] < 0

    # ----------------------------------------------------------------------------
    # Converting the image into a array
    # ----------------------------------------------------------------------------

    binary_vector = binary_vectors.flatten()

    #-------------------------------------------------------------------------------
    if VERBOSE:
        # ----------------------------------------------------------------------------
        # Plotting
        # ----------------------------------------------------------------------------

        axis_major = SecondLocator(interval=2)   # every 5-second
        axis_minor = SecondLocator(interval=1)   # every 1-second
        axis_Fmt = DateFormatter('%H:%M:%S')

        plt.rcParams.update({'font.size': 12})
        fig = plt.figure(figsize=(15,15))
        gs = fig.add_gridspec(3, 3)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
        #----------------------------------------------------------------------------
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, :],sharex=ax1)
        ax3 = fig.add_subplot(gs[2, 0])
        ax4 = fig.add_subplot(gs[2, 1])
        ax5 = fig.add_subplot(gs[2, 2])
        #-----------------------
        ax1.set_title(trim_data.stats.network+'.'+trim_data.stats.station+'.'+trim_data.stats.channel)
        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5)
        ax1.tick_params(axis='both',which='minor',width=2,length=3)
        ax1.plot(trim_data.times('matplotlib'),trim_data.data, color='k', linewidth=2)

        #----------------------------------------------------------------------------
        ax2.set_title('Continuous Wavelet Transform')
        fig.suptitle('Event Date: '+datetime_window.strftime('%d/%m/%Y'),y=0.95,fontsize=25)

        ax2.xaxis.set_major_formatter(axis_Fmt)
        ax2.xaxis.set_major_locator(axis_major)
        ax2.xaxis.set_minor_locator(axis_minor)
        ax2.tick_params(axis='both',which='major',width=2,length=5)
        ax2.tick_params(axis='both',which='minor',width=2,length=3)

        im = ax2.imshow(event_pattern,extent= [t[0],t[-1],f_min,f_max],aspect='auto',cmap='plasma',origin='lower',vmin=-np.abs(event_pattern.max()),vmax=np.abs(event_pattern.max()))
        ax2.set_ylabel("Frequency [Hz]")
        axins = inset_axes(ax2,
                            width="20%",
                            height="5%",
                            loc='upper left',
                            bbox_to_anchor=(0.80, 0.1, 1, 1),
                            bbox_transform=ax2.transAxes,
                            borderpad=0,
                            )

        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

        #----------------------------------------------------------------------------

        ax3.set_title('|Haar transform|')

        im3 = ax3.imshow(haar_image,extent=[0,fp_length,0,fp_length],aspect='auto',cmap='plasma',origin='lower',vmin=-np.abs(haar_image.max()),vmax=np.abs(haar_image.max()))
        ax3.set_xlabel('wavelet transform x index')
        ax3.set_ylabel('wavelet transform y index')
        axins3 = inset_axes(ax3,
                            width="20%",
                            height="5%",
                            loc='upper left',
                            bbox_to_anchor=(0.8, 0.1, 1, 1),
                            bbox_transform=ax3.transAxes,
                            borderpad=0,
                            )

        plt.colorbar(im3, cax=axins3, orientation="horizontal", ticklocation='top')

        #----------------------------------------------------------------------------

        ax4.set_title('Top Coefficients')

        im4 = ax4.imshow(haar_image_top,extent=[0,fp_length,0,fp_length],aspect='auto',cmap='gray',origin='lower',vmin=-np.abs(haar_image_top.max()),vmax=np.abs(haar_image_top.max()))
        ax4.set_xlabel('wavelet transform x index')
        ax4.set_ylabel('wavelet transform y index')
        axins4 = inset_axes(ax4,
                            width="20%",
                            height="5%",
                            loc='upper left',
                            bbox_to_anchor=(0.8, 0.1, 1, 1),
                            bbox_transform=ax4.transAxes,
                            borderpad=0,
                            )

        plt.colorbar(im4, cax=axins4, orientation="horizontal", ticklocation='top')
        #----------------------------------------------------------------------------

        ax5.set_title('Binary Fingerprint')

        im5 = ax5.imshow(binary_vectors,aspect='auto',cmap='gray',origin='lower',vmin=0,vmax=1)
        ax5.set_xlabel('wavelet transform x index')
        ax5.set_ylabel('wavelet transform y index')
        axins5 = inset_axes(ax5,
                            width="20%",
                            height="5%",
                            loc='upper left',
                            bbox_to_anchor=(0.8, 0.1, 1, 1),
                            bbox_transform=ax5.transAxes,
                            borderpad=0,
                            )

        plt.colorbar(im5, cax=axins5, orientation="horizontal", ticklocation='top')
        #----------------------------------------------------------------------------
        plt.show()

    return binary_vector

# ============
# Main program
# ============

print('=====================')
print('Loading Event Pattern')
print('=====================')
print('\n')

EVENT_STANDARD_TIME = UTCDateTime(EVENT_PATTERN_DATE)
EVENT_STANDARD_TIME_STR = str(EVENT_STANDARD_TIME.year)+'.'+"%03d" % EVENT_STANDARD_TIME.julday

obs_HHZ_standard_pattern = MSEED_FOLDER+'/'+NETWORK+'/'+OBS_NAME+'/'+CHANNEL+'.D/'+NETWORK+'.'+OBS_NAME+'..'+CHANNEL+'.D.'+EVENT_STANDARD_TIME_STR

obs_HHZ_standard_pattern_waveform = read(obs_HHZ_standard_pattern)
tr = obs_HHZ_standard_pattern_waveform[0]

tr.trim(EVENT_STANDARD_TIME-PEM,EVENT_STANDARD_TIME+PET)

fingerprint_vector = fingerprint_creator(obspy_trace=tr,freqmin=FILTER_DATA[0],freqmax=FILTER_DATA[1],fp_length=spectral_length)

# Convert from pandas to Arrow and saving in feather formart file
os.makedirs(BINARY_FILES+EVENT_TYPE,exist_ok=True)
file_BINARY_name = BINARY_FILES+EVENT_TYPE+'/'+EVENT_TYPE+'_standard_pattern_date_'+'_'.join(EVENT_PATTERN_DATE.split(','))
np.save(file=file_BINARY_name, arr=fingerprint_vector)
