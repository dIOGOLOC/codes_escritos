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
from PIL import Image

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature

from pyasdf import ASDFDataSet

from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta,plot_trigger


# ========================
# Input and output folders
# ========================

MSEED_FOLDER = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/obs_data_MSEED/'
#MSEED_FOLDER = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

EARTHQUAKE_FINDER_OUTPUT = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'
#EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'

#ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'
ASDF_FILES = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'

#STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'
STATIONXML_DIR = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

#BINARY_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/'
BINARY_FILES = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/'

# ===================================
# Loading standard_pattern event file
# ===================================

#standard_pattern_binary = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/AIRGUN/AIRGUN_standard_pattern_date_2019_08_04_20_46_58.npy'
standard_pattern_binary = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/AIRGUN/AIRGUN_standard_pattern_date_2019_08_04_20_46_58.npy'

# ==========
# Parameters
# ==========

#Type of event
EVENT_TYPE = 'EVENT_HARMONIC'

#Bandpass frequency (Hz) - minimum and maximum
FILTER_DATA = [2,20]

NETWORK = 'ON'

OBS_NAME = 'OBS18'

CHANNEL = 'HHZ'

#Spectral image length (samples)
spectral_length = 32

#Spectral image lag (samples)
TIME_WINDOW_LAG = 2

#Time window length (same of the standard_pattern_binary)
TIME_WINDOW = 7

#Number of wavelet coefficients to keep
k_coef = 200

#The Jaccard index or Jaccard similarity coefficient threshold
#defined as the size of the intersection divided by the size of the union of two sets,
#is used to compare sets
JACC_COEFF_THR = 0.5

# =========
# Constants
# =========

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# =================
# Filtering by date
# =================

FIRSTDAY = '2020,02,20'
LASTDAY = '2020,03,08'

fday = UTCDateTime(FIRSTDAY)
lday = UTCDateTime(LASTDAY)

INTERVAL_PERIOD = [UTCDateTime(x.astype(str)) for x in np.arange(fday.datetime,lday.datetime+ONEDAY,ONEDAY)]
INTERVAL_PERIOD_DATE = [str(x.year)+'.'+"%03d" % x.julday for x in INTERVAL_PERIOD]


# ================
# MULTIPROCESSING
# ================

#Number of threads
num_processes = 4

# =========
# Functions
# =========

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

#----------------------------------------------------------------------------

def fingerprint_creator(obspy_trace,freqmin,freqmax,fp_length):
    tr = obspy_trace
    #----------------------------------------------------------------------------

    inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,OBS_NAME,'xml']))
    pre_filt = [0.001, 0.005, 45., 50.]
    tr.remove_response(inventory=inv,output="DISP",pre_filt=pre_filt,water_level=60)
    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.1, type="hann")
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
    bin_fingerprint = np.zeros((N,2*M), dtype=bool)
    for i in range(N):
        idx = np.argsort(abs(haar_image_top[i,:]))[-K:]
        bin_fingerprint[i,idx]   = haar_image_top[i,idx] > 0
        bin_fingerprint[i,idx+M] = haar_image_top[i,idx] < 0

    # ----------------------------------------------------------------------------
    # Converting the image into a array
    # ----------------------------------------------------------------------------

    binary_vector = bin_fingerprint.flatten()

    return [binary_vector,trim_data,event_pattern,[t[0],t[-1],f_min,f_max],haar_image,[0,fp_length,0,fp_length],haar_image_top,[0,fp_length,0,fp_length],bin_fingerprint]

# ----------------------------------------------------------------------------

def plot_binary_image(input):

        trim_data = input[0]
        event_pattern = input[1]
        event_pattern_extent = input[2]
        haar_image = input[3]
        haar_image_extent = input[4]
        haar_image_top = input[5]
        haar_image_top_extent = input[6]
        bin_fingerprint = input[7]
        jaccard_similarity = input[8]

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

        #----------------------------------------------------------------------------

        ax1.set_title(trim_data.stats.network+'.'+trim_data.stats.station+'.'+trim_data.stats.channel+' (Similarity='+str(jaccard_similarity)+')')
        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5)
        ax1.tick_params(axis='both',which='minor',width=2,length=3)
        ax1.plot(trim_data.times('matplotlib'),trim_data.data, color='k', linewidth=2)

        #----------------------------------------------------------------------------

        ax2.set_title('Continuous Wavelet Transform')
        fig.suptitle('Event Date: '+trim_data.stats.starttime.strftime('%d/%m/%Y'),y=0.95,fontsize=25)
        ax2.xaxis.set_major_formatter(axis_Fmt)
        ax2.xaxis.set_major_locator(axis_major)
        ax2.xaxis.set_minor_locator(axis_minor)
        ax2.tick_params(axis='both',which='major',width=2,length=5)
        ax2.tick_params(axis='both',which='minor',width=2,length=3)

        im = ax2.imshow(event_pattern,extent=event_pattern_extent,aspect='auto',cmap='plasma',origin='lower',vmin=-np.abs(event_pattern.max()),vmax=np.abs(event_pattern.max()))
        ax2.set_ylabel("Frequency [Hz]")
        axins = inset_axes(ax2,
                                width="20%",
                                height="5%",
                                loc='upper left',
                                bbox_to_anchor=(0.8, 0.1, 1, 1),
                                bbox_transform=ax2.transAxes,
                                borderpad=0,
                                )

        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

        #----------------------------------------------------------------------------

        ax3.set_title('|Haar transform|')
        im3 = ax3.imshow(haar_image,extent=haar_image_extent,aspect='auto',cmap='plasma',origin='lower',vmin=-np.abs(haar_image.max()),vmax=np.abs(haar_image.max()))
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
        im4 = ax4.imshow(haar_image_top,extent=haar_image_top_extent,aspect='auto',cmap='gray',origin='lower',vmin=-np.abs(haar_image_top.max()),vmax=np.abs(haar_image_top.max()))
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
        im5 = ax5.imshow(bin_fingerprint,aspect='auto',cmap='gray',origin='lower',vmin=0,vmax=1)
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

        folder_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'.'+OBS_NAME+'.'+CHANNEL+'/'+EVENT_TYPE+'_PATTERN_SELECTED/'
        os.makedirs(folder_output,exist_ok=True)
        fig.savefig(folder_output+NETWORK+'_'+OBS_NAME+'_'+CHANNEL+'_'+EVENT_TYPE+'_'+trim_data.stats.starttime.strftime('%Y_%d_%m_%H_%M_%S_%f')+'.png', dpi=300, facecolor='w', edgecolor='w')
        plt.close(fig)

#-------------------------------------------------------------------------------

def calculate_jaccard_similarity(input):

    data = input[0]
    freqmin = input[1]
    freqmax = input[2]
    fp_length = input[3]
    event_fp = input[4]
    threshold = input[5]

    output_lst = fingerprint_creator(obspy_trace=data,freqmin=freqmin,freqmax=freqmax,fp_length=fp_length)

    #-------------------------------------------------------------------------------

    binary_vector = output_lst[0]
    data_plot = output_lst[1]
    event_pattern = output_lst[2]
    event_pattern_extent = output_lst[3]
    haar_image = output_lst[4]
    haar_image_extent = output_lst[5]
    haar_image_top = output_lst[6]
    haar_image_top_extent = output_lst[7]
    bin_fingerprint = output_lst[8]

    jaccard_similarity = jaccard_score(event_fp, binary_vector)

    if jaccard_similarity > threshold:
        plot_binary_image([data_plot,event_pattern,event_pattern_extent,haar_image,haar_image_extent,haar_image_top,haar_image_top_extent,bin_fingerprint,jaccard_similarity])

#-------------------------------------------------------------------------------

# ============
# Main program
# ============

print('=========================')
print('Loading the Event pattern')
print('=========================')
print('\n')

binary_vector_float = np.load(standard_pattern_binary)
EVENT_FINGERPRINT = binary_vector_float

#----------------------------------------------------------------------------

print('===================================================')
print('Loading events files retrieved by LTA/STA procedure')
print('===================================================')
print('\n')

#----------------------------------------------------------------------------

for data_h5 in tqdm(INTERVAL_PERIOD_DATE,desc='File loop'):

    file__TIME_STR = data_h5
    obs_day_files = MSEED_FOLDER+'/'+NETWORK+'/'+OBS_NAME+'/'+CHANNEL+'.D/'+NETWORK+'.'+OBS_NAME+'..'+CHANNEL+'.D.'+file__TIME_STR
    obs_day_waveform = read(obs_day_files)[0]

    slide_data = [k for k in obs_day_waveform.slide(window_length=TIME_WINDOW, step=TIME_WINDOW_LAG,include_partial_windows=False, nearest_sample=True)]

    input_data_lst = []
    for slided_data in slide_data:
        input_data_lst.append([slided_data,FILTER_DATA[0],FILTER_DATA[1],spectral_length,EVENT_FINGERPRINT,JACC_COEFF_THR])

    with Pool(processes=num_processes) as p:
        max_ = len(input_data_lst)
        with tqdm(total=max_, desc='Slice loop') as pbar:
            for i, _ in enumerate(p.imap_unordered(calculate_jaccard_similarity, input_data_lst)):
                pbar.update()
