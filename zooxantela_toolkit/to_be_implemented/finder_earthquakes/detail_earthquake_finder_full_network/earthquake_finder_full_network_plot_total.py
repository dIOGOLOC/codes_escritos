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

# ====================================================================================================
# Configuration file
# ====================================================================================================

EARTHQUAKE_FINDER_FILES = '/run/user/1000/gvfs/smb-share:server=hatabackup.local,share=on_mar/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/ON.OBS17.HHZ/AIRGUN_PATTERN_SELECTED/'

MSEED_DIR_STA = '/home/diogoloc/dados_posdoc/ON_MAR/data/'

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/Figuras/EVENTS_FILTERED/'

#start date deployment
INITIAL_DATE = '2019-06-01'

#final date deployment
FINAL_DATE = '2020-08-01'

#---------------------------------------------------------------------------------------------------------

def filelist(basedir):
    """
    Returns the list of files in *basedir* whose are in the specified period
    """
    day_files = []
    for root, dirs, files in os.walk(basedir):
        for file in files:
            if file.endswith(".png"):
                day_files.append(os.path.join(file))
    dates_files = []

    for i in day_files:
        year = i.split('_')[1]
        day = i.split('_')[2]
        mounth = i.split('_')[3]
        hour = i.split('_')[4]
        dates_files.append(UTCDateTime(year+'-'+mounth+'-'+day+'T'+hour)-3)
    return sorted(dates_files)

# ============
# Main program
# ============

print('==============')
print('Scanning files')
print('==============')
print('\n')

# initializing list of stations by scanning name of miniseed files

dates_earthquakes = filelist(basedir=EARTHQUAKE_FINDER_FILES)
print(len(dates_earthquakes))

print('\n')
print('====================')
print('Plotting events date')
print('====================')
print('\n')

# ==========================================================
# Calculating datetime between INITIAL_DATE and  FINAL_DATE
# ==========================================================

datatime_initial = datetime.datetime(UTCDateTime(INITIAL_DATE).year,UTCDateTime(INITIAL_DATE).month,UTCDateTime(INITIAL_DATE).day)

datatime_final = datetime.datetime(UTCDateTime(FINAL_DATE).year,UTCDateTime(FINAL_DATE).month,UTCDateTime(FINAL_DATE).day)

datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)

xlim_initial = mdates.date2num(datatime_initial)
xlim_final = mdates.date2num(datatime_final)

#----------------------------
#Function to check if the dates in data set are inside the period chosen (INITIAL_DATE to FINAL_DATE)
#data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['NUMBER_HOUR'])


def check_datetime_in_period(datetime_lst,dates_earthquakes_lst):
    array_to_plot_by_xlim = []
    for x,c in enumerate(datetime_lst):
        lista_temp = []
        for t,y in enumerate(dates_earthquakes_lst):
            if datetime.datetime(y.year,y.month,y.day) == c:
                lista_temp.append([1 if y.hour == i  else 0 for i in np.arange(24)])
        array_to_plot_by_xlim.append(lista_temp)

    data_x_axis = []
    for x,c in enumerate(array_to_plot_by_xlim):
        if c != []:
            data_x_axis.append(c[0][::-1])
        else:
            data_x_axis.append(np.zeros_like(np.arange(24)))

    data_x_axis = np.array(data_x_axis).T

    return data_x_axis

data_x_axis = check_datetime_in_period(datetime_lista,dates_earthquakes)

# ====================================
# Function to plot DATA availability
# ====================================

#x axis parameters

days1 = DayLocator(interval=1)   # every day
days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every 5 day
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

days1.MAXTICKS = 10000

#Matplotlib parameters
fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,10))

data_x_axis = check_datetime_in_period(datetime_lista,dates_earthquakes)
ax.set_title('SDE Selecionados',fontsize=25)
im = ax.imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.binary,interpolation=None, vmin=0, vmax=1)
ax.set_xlim(datetime.datetime(UTCDateTime(INITIAL_DATE).year,UTCDateTime(INITIAL_DATE).month,UTCDateTime(INITIAL_DATE).day),datetime.datetime(UTCDateTime(FINAL_DATE).year,UTCDateTime(FINAL_DATE).month,UTCDateTime(FINAL_DATE).day))
ax.yaxis.set_major_locator(MultipleLocator(3))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_major_locator(days5)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days1)
ax.tick_params(which='minor', length=4)
ax.tick_params(which='major', length=10)
ax.set_aspect(5)
ax.set_ylim(0,24)
ax.set_ylabel('Hora Local',fontsize=15)
ax.grid(b=True, which='major', color='k', linestyle='-')
ax.grid(b=True, which='minor', color='k', linestyle='-')

plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)

os.makedirs(EARTHQUAKE_FINDER_OUTPUT,exist_ok=True)
fig.savefig(EARTHQUAKE_FINDER_OUTPUT+'Eventos_selecionados.png',dpi=300,bbox_inches='tight',pad_inches=0.1)
