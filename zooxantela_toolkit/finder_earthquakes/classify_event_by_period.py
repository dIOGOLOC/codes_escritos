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
from pyrocko import obspy_compat
obspy_compat.plant()

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

# ===================================
# Loading standard_pattern event file
# ===================================

standard_pattern_binary = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/AIRGUN/AIRGUN_standard_pattern_date_2019_08_04_20_46_58.npy'
#standard_pattern_binary = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/BINARY_FILES/AIRGUN/AIRGUN_standard_pattern_date_2019_08_04_20_46_58.npy'

# ==========
# Parameters
# ==========

#Bandpass frequency (Hz) - minimum and maximum
FILTER_DATA = [10,40]

NETWORK = 'ON'

OBS_NAME = ['OBS17','OBS18','OBS20','OBS22']

CHANNEL = 'HHZ'

# =========
# Constants
# =========

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
MINUTES30 = datetime.timedelta(minutes=30)
ONEDAY = datetime.timedelta(days=1)

# =================
# Filtering by date
# =================

FIRSTDAY = '2019,12,07,11'
LASTDAY = '2019,12,08,11'

fday = UTCDateTime(FIRSTDAY)
lday = UTCDateTime(LASTDAY)

INTERVAL_PERIOD = [UTCDateTime(x.astype(str)) for x in np.arange(fday.datetime,lday.datetime+MINUTES30,MINUTES30)]
INTERVAL_PERIOD_DATE = [str(x.year)+'.'+"%03d" % x.julday for x in INTERVAL_PERIOD]

# =========
# Functions
# =========


# ============
# Main program
# ============

print('===================================================')
print('Loading events files retrieved by LTA/STA procedure')
print('===================================================')
print('\n')

xml_files = sorted(glob.glob(STATIONXML_DIR+'ON.OBS*'))
inv = read_inventory(xml_files[0])
for xml_file in xml_files[1:]:
    inv.extend(read_inventory(xml_file))

#----------------------------------------------------------------------------

for iperid,period_date in enumerate(tqdm(INTERVAL_PERIOD_DATE,desc='File loop')):
    obs_day_files = glob.glob(MSEED_FOLDER+'**/**/**/*'+str(period_date))

    st = Stream()
    for file in obs_day_files:
        if 'HHX' not in file and 'OBS19' not in file:
            st.append(read(file)[0])

    st.trim(starttime=INTERVAL_PERIOD[iperid], endtime=INTERVAL_PERIOD[iperid]+MINUTES30)

    # Start the Snuffler
    st.snuffle(inventory=inv)
