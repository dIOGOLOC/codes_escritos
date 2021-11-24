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


#Configuration file

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_NETWORK_OUTPUT/ASDF_FILES/'

FILTER_DATA = [2,40]

NETWORK = 'ON'

# ========================
# Constants and parameters
# ========================

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# ============
# Main program
# ============

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
