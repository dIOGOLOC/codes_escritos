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


import obspy
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++ Configuration file ++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pasta com os arquivos em MSEED (MSEED_DIR+NETWORK+'/'+STATION+'/'+CHANNEL+'.D/')
MSEED_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

# Escolha uma rede:
NETWORK = 'ON'

# Escolha uma estação:
STATION = 'OBS19'

# Escolha um canal:
CHANNEL = 'HHZ'

# Número de processadores:
num_processes = 8

# Pasta para salvar as figuras
FIGURES_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/Figuras/'

# Data inicial e final da aquisição
INITIAL_DATE = '2019-08-02'
FINAL_DATE = '2020-06-17'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ========================
# Constants and parameters
# ========================

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# =================
# Filtering by date
# =================

fday = UTCDateTime(INITIAL_DATE)
lday = UTCDateTime(FINAL_DATE)
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

# =============================
# Function to store date in dic
# =============================

def get_date_file(FILE_CHANNEL):
    dataframe_lista = []
    #create a empty dataframe with pandas

    #Reading header from data
    st = obspy.read(FILE_CHANNEL)

    #----------------------------
    #Dataframe starting

    network = st[0].stats.network
    station = st[0].stats.station
    channel = st[0].stats.channel

    time_lst = []
    for t,trace in enumerate(st):
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime
        time_lst.append(np.arange(starttime,endtime,60))

    flat_time_lst = [item for sublist in time_lst for item in sublist]

    DATETIME = datetime.datetime(st[0].stats.starttime.year,st[0].stats.starttime.month,st[0].stats.starttime.day)

    #Contador da lista de horas
    NUMBER_HOUR = []
    for x,c in enumerate(flat_time_lst):
        NUMBER_HOUR.append(c.hour)

    dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[len(list(set(NUMBER_HOUR)))]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','NUMBER_HOUR']).T)
    #Dataframe ending
    #----------------------------
    df = pd.concat(dataframe_lista, ignore_index=True)
    # using dictionary to convert specific columns
    convert_dict = {'NETWORK':str, 'STATION':str, 'CHANNEL':str, 'DATETIME':str, 'NUMBER_HOUR':int}
    df = df.astype(convert_dict)

    return df

# ================================================================
# Function to check if the dates in data set are inside the period
# ================================================================

def check_datetime_in_period(datetime_lst,df_DATETIME,df_NUMBER_HOUR):
    array_to_plot_by_xlim = []
    for x,c in enumerate(datetime_lst):
        array_to_plot_by_xlim.append([df_NUMBER_HOUR.array[t] for t,y in enumerate(df_DATETIME.array) if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c])

    data_x_axis = []
    for x,c in enumerate(array_to_plot_by_xlim):
        if c != []:
            data_x_axis.append([c[0]])
        else:
            data_x_axis.append([0])
    data_x_axis = np.array(data_x_axis).T
    return data_x_axis

def plot_date_file(DATAFRAME):
    # ==========================================================
    # Calculating datetime between INITIAL_DATE and  FINAL_DATE
    # ==========================================================

    datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)
    datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)
    datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)

    xlim_initial = mdates.date2num(datatime_initial)
    xlim_final = mdates.date2num(datatime_final)

    # ====================================
    # Function to plot DATA availability
    # ====================================
    data_x_axis = check_datetime_in_period(datetime_lista,DATAFRAME['DATETIME'],DATAFRAME['NUMBER_HOUR'])

    #x axis parameters

    days1 = DayLocator(interval=1)   # every day
    days5 = DayLocator(interval=30)   # every day
    months = MonthLocator()  # every month
    yearsFmt = DateFormatter('%Y-%m-%d')

    days1.MAXTICKS = 10000

    #Percentage of data
    total_length = sum(len(row) for row in data_x_axis)
    total_no_data  = np.count_nonzero(data_x_axis)
    data_percentage = int(round((total_no_data*100)/total_length))


    #Matplotlib parameters
    fig, ax = plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(20,4))
    plt.rcParams.update({'font.size': 20})

    ax.set_title(STATION+' - '+str(data_percentage)+'% de dados',fontsize=20)

    im = ax.imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,1],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=24)
    ax.set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
    ax.set_aspect(10)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(days5)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(days1)
    ax.tick_params(which='minor', length=4)
    ax.tick_params(which='major', length=10)
    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='k', linestyle='-')
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=20, rotation=30)

    #criando a localização da barra de cores:
    axins = inset_axes(ax,
                        width="10%",  # width = 10% of parent_bbox width
                        height="10%",  # height : 50%
                        loc='upper left',
                        bbox_to_anchor=(0.85, 0.2, 1, 1),
                        bbox_transform=ax.transAxes,
                        borderpad=0,
                        )
    cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,12,24],label='Arquivos/dia')
    cbar.ax.set_xticklabels(['0%','50%','100%'])

    COMPLETUDE_FIG_FOLDER = FIGURES_OUTPUT+'COMPLETUDE/'
    os.makedirs(COMPLETUDE_FIG_FOLDER,exist_ok=True)
    fig.savefig(COMPLETUDE_FIG_FOLDER+STATION+'_'+'COMPLETUDE_BIG'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.jpg',dpi=300,bbox_inches='tight',pad_inches=0.3)

print('============')
print('Main program')
print('============')
print('\n')

#Procurando os arquivos em MSEED
data_lista = []
print('Data from = '+MSEED_DIR+NETWORK+'/'+STATION+'/'+CHANNEL+'.D/')
print('\n')

for root, dirs, files in os.walk(MSEED_DIR+NETWORK+'/'+STATION+'/'+CHANNEL+'.D/'):
    for name in files:
        data_lista.append(os.path.join(root, name))

data_lista = sorted(data_lista)

#Processando cada arquivos em MSEED

start_time = time.time()
df_lst = []
with Pool(processes=num_processes) as p:
    max_ = len(data_lista)
    with tqdm(total=max_, desc='Extracting data') as pbar:
        for i, df_data in enumerate(p.imap_unordered(get_date_file, data_lista)):
            pbar.update()
            df_lst.append(df_data)
df = pd.concat(df_lst)
df = df.sort_values(by="DATETIME")
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

#Plotando o resultado final
plot_date_file(df)

#FIM
