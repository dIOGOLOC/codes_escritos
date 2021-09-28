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
from obspy.signal import PPSD
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.core.util import AttribDict

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

#Directory to save PSD
OUTPUT_PSD_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_PSD_MSD/'

#Directory witg noise models
NOISE_MODEL_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/noise_models.npz'

#Channel list
CHANNEL_LST = ['HHZ.','HHN.','HHE.','HH1.','HH2.','HHX.']

#DATE_DAY = '2019-12-24T19:04:00.00'
#DATE_DAY = '2019-12-24T19:04:00.00'
DATE_DAY = '2019-12-24T16:47:00.00'

FDAY = UTCDateTime(DATE_DAY)
INTERVAL_PERIOD_DATE = str(FDAY.year)+'.'+"%03d" % FDAY.julday

NETWORK = 'ON'

#STATIONs
#STATIONS_LST = ['DUB01','CAM01','ABR01','GUA01']
STATIONS_LST = []
#OBSs
OBS_LST = ['OBS17','OBS18','OBS20','OBS22']

PEM = 20
PET = 210

#Low-cut frequency (Hz) and High-cut frequency (Hz) for bandpass filter
FILTER_DATA = [1,10]

#Length of the short-term average window (seconds)
stalen = 1

#Length of the long-term average window (seconds)
ltalen = 30

#sta/lta ratio to trigger a detection/pick
trig_on = 10

#sta/lta ratio to turn the trigger off - no further picks\
trig_off = 2

#Show picks on waveform.
VERBOSE_MODE = True

#Maximum and minimum period of the PSD (sensor).
PERIOD_LIM = (0.02, 200)

#Maximum and minimum amplitude of the PSD (sensor).
AMP_PSD_MIN = -200
AMP_PSD_MAX = -65

#Restricts the data that is included in the stack by time of day and weekday.
#Monday is 1, Sunday is 7, -1 for any day of week.
#For example, using time_of_weekday=[(-1, 22, 24)]
#only individual spectra that have a starttime in between 10pm and 12am are used in the stack for all days of week
#time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
TIME_OF_WEEKDAY_DAY = -1

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

files_INTERVAL_PERIOD_DATE = lst_INLAND+lst_OBS

# =====================
# Retrieving .NPZ files
# =====================

PPSD_files = []
lst_OBS1 = glob.glob(OUTPUT_PSD_DIR+'/*')
PPSD_dir = []
for sta in OBS_LST:
    b = filelist(basedir=[i for i in lst_OBS1 if sta in i][0],interval_period_date=INTERVAL_PERIOD_DATE,channel_list=CHANNEL_LST)
    if b != []:
        PPSD_files.append(b)

print('\n')
print('==================')
print('Loading MSEED data')
print('==================')
print('\n')

st = Stream()
for k in files_INTERVAL_PERIOD_DATE:
    for filename in k:
        st += read(filename)

print('\n')
print('=================')
print('Loading PPSD data')
print('=================')
print('\n')

TIME_OF_WEEKDAY_START_HOUR = FDAY.hour
TIME_OF_WEEKDAY_FINAL_HOUR = FDAY.hour+FDAY.minute

ppsd_lst_HHZ = []
ppsd_lst_HHE = []
ppsd_lst_HHN = []
for k in PPSD_files:
    for filename in k:
        ppsd = PPSD.load_npz(filename)
        #ppsd.calculate_histogram(starttime=UTCDateTime(FDAY),endtime=UTCDateTime(FDAY+PEM),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
        ppsd.calculate_histogram(time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
        ppsd_periods, ppsd_mean = ppsd.get_mean()

        if ppsd.channel == 'HHZ':
            ppsd_lst_HHZ.append([ppsd_periods, ppsd_mean])
        elif ppsd.channel == 'HHN' or ppsd.channel == 'HH1':
            ppsd_lst_HHN.append([ppsd_periods, ppsd_mean])
        elif ppsd.channel == 'HHE' or ppsd.channel == 'HH2':
            ppsd_lst_HHE.append([ppsd_periods, ppsd_mean])


#===============
#Preprocess Data
#===============

st.trim(starttime=FDAY-PEM, endtime=FDAY+PET)

stE = Stream()
stN = Stream()
stZ = Stream()
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
    else:
        inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))
        pre_filt = [0.001, 0.005, 45., 50.]
        tr.remove_response(inventory=inv,pre_filt=pre_filt,output="VEL",water_level=60)
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(max_percentage=0.05)
        tr.filter('bandpass', freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])  # optional prefiltering

        coordinates_lst = inv[0][0]
        tr.stats.coordinates = AttribDict({'latitude': coordinates_lst.latitude,'elevation': coordinates_lst.elevation,'longitude': coordinates_lst.longitude})

        if channel == 'HHZ':
            stZ += tr
        elif channel == 'HHN' or channel == 'HH1':
            stN += tr
        elif channel == 'HHE' or channel == 'HH2':
            stE += tr

stZ2 = stZ.copy()
stN2 = stN.copy()
stE2 = stE.copy()
stX2 = stX.copy()

st_selected_time = st[0].stats.starttime

print('\n')
print('==============')
print('Plotting event')
print('==============')
print('\n')

# Plotting the results
axis_major = SecondLocator(interval=60)   # every 5-second
axis_minor = SecondLocator(interval=5) # every 1-second
axis_Fmt = DateFormatter('%H:%M:%S')

#----------------------------------------------------------------------------

plt.rcParams.update({'font.size': 15})
fig = plt.figure(figsize=(30,20))
gs = gridspec.GridSpec(nrows=7, ncols=len(stZ2))


client = Client("IRIS")
ev = client.get_events(starttime=FDAY-20*60, endtime=FDAY+210, minmagnitude=6)
model = TauPyModel(model="iasp91")
O = ev[0].preferred_origin()
mag = ev[0].magnitudes[0].mag
mag_type = ev[0].magnitudes[0].magnitude_type

fig.suptitle('Event Date: '+FDAY.strftime('%d/%m/%Y')+'\n'+ev[0].event_descriptions[0].text.title()+'(M='+str(round(mag,1))+' '+mag_type+')',y=0.95,fontsize=25,color='white')
#----------------------------------------------------------------------------
for ind,traces in enumerate(stZ2):

        #-------------------------------------------------------------------------------------------
        d = op.geodetics.kilometer2degrees(op.geodetics.calc_vincenty_inverse(traces.stats.coordinates.latitude, traces.stats.coordinates.longitude, O.latitude, O.longitude)[0]/1000.)
        T = model.get_travel_times(O.depth / 1000., d)
        arrP = T[0]
        dist_km = op.geodetics.calc_vincenty_inverse(traces.stats.coordinates.latitude, traces.stats.coordinates.longitude, O.latitude, O.longitude)[0]/1000.

        #-------------------------------------------------------------------------------------------

        t = traces.times('matplotlib')
        f_min = FILTER_DATA[0]
        f_max = FILTER_DATA[1]

        #----------------------------------------------------------------------------
        data_model_noise = np.load(NOISE_MODEL_FILE)
        periods_model_noise_ln = data_model_noise['model_periods']
        nlnm_ln = data_model_noise['low_noise']

        periods_model_noise_hn = data_model_noise['model_periods']
        nlnm_hn = data_model_noise['high_noise']

        ax5 = fig.add_subplot(gs[0,ind])

        ax5.set_title(traces.stats.station+'('+str(round(dist_km))+' km)')
        ax5.plot(ppsd_lst_HHZ[ind][0],ppsd_lst_HHZ[ind][1], ls='solid',c='k', linewidth=3 ,zorder=10, label='HHZ')
        ax5.plot(ppsd_lst_HHN[ind][0],ppsd_lst_HHN[ind][1], ls='dotted',c='k', linewidth=3 ,zorder=10, label='HHN')
        ax5.plot(ppsd_lst_HHE[ind][0],ppsd_lst_HHE[ind][1], ls='--',c='k', linewidth=3 ,zorder=10, label='HHE')

        ax5.plot(periods_model_noise_ln,nlnm_ln, color='gray',ls='-', linewidth=2, zorder=10,alpha=0.7)
        ax5.plot(periods_model_noise_hn,nlnm_hn, color='gray',ls='-', linewidth=2, zorder=10,alpha=0.7)

        ax5.axvline(1/f_max, color='r',ls='--', linewidth=1, zorder=10)
        ax5.axvline(1/f_min, color='r',ls='--', linewidth=1, zorder=10)

        ax5.semilogx()
        ax5.set_xlim(PERIOD_LIM)
        ax5.set_ylim(AMP_PSD_MIN, AMP_PSD_MAX)
        ax5.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax5.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax5.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax5.xaxis.label.set_color('white')
        ax5.yaxis.label.set_color('white')
        ax5.title.set_color('white')

        if ind == 0:
            ax5.set_ylabel("Amp.(DB)")

        #----------------------------------------------------------------------------

        scalogram = cwt(traces.data*1000, dt, 8, f_min, f_max)
        x, y = np.meshgrid(t,np.linspace(f_min, f_max, scalogram.shape[0]))

        ax0 = fig.add_subplot(gs[1,ind])
        ax0.xaxis.set_major_formatter(axis_Fmt)
        ax0.xaxis.set_major_locator(axis_major)
        ax0.xaxis.set_minor_locator(axis_minor)
        ax0.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax0.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax0.xaxis.label.set_color('white')
        ax0.yaxis.label.set_color('white')
        ax0.title.set_color('white')

        im = ax0.pcolormesh(x, y, np.abs(scalogram), shading='auto', cmap='viridis')
        if ind == 0:
            ax0.set_ylabel("Freq.(Hz)")
        ax0.set_ylim(f_min, f_max)
        plt.setp(ax0.get_xticklabels(), visible=False)

        #axins = inset_axes(ax0,
                   #width="25%",
                   #height="5%",
                   #loc='upper left',
                   #bbox_to_anchor=(0.70, 0.1, 1, 1),
                   #bbox_transform=ax0.transAxes,
                   #borderpad=0,
                   #)

        #plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top', label='Vel[mm/s]')
        #-------------------------------------------------------------------------------------------

        ax1 = fig.add_subplot(gs[2,ind],sharex=ax0)
        trace_amp_max = np.max(np.abs(traces.data*1000))

        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax1.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax1.xaxis.label.set_color('white')
        ax1.yaxis.label.set_color('white')
        ax1.title.set_color('white')

        ax1.plot(traces.times('matplotlib'),traces.data*1000, 'k')
        ax1.text(x=((O.time+arrP.time)-10).matplotlib_date,y=trace_amp_max-(trace_amp_max*0.3),s='P',bbox=dict(facecolor='w', alpha=0.2,edgecolor='w'))

        ax1.set_ylim(-trace_amp_max,trace_amp_max)
        ax1.text(0.05, 0.1, traces.stats.channel, horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)
        ax1.ticklabel_format(axis='y', style='sci')

        if ind == 0:
            ax1.set_ylabel("Vel.(mm/s)")
        plt.setp(ax1.get_xticklabels(), visible=False)

        #----------------------------------------------------------------------------
        ax2 = fig.add_subplot(gs[3,ind],sharex=ax1)
        trace_HHN = stN2[ind]
        trace_HHN_amp_max = np.max(np.abs(trace_HHN.data*1000))

        ax2.xaxis.set_major_locator(axis_major)
        ax2.xaxis.set_major_formatter(axis_Fmt)
        ax2.xaxis.set_minor_locator(axis_minor)
        ax2.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax2.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax2.xaxis.label.set_color('white')
        ax2.yaxis.label.set_color('white')
        ax2.title.set_color('white')

        ax2.plot(trace_HHN.times('matplotlib'),trace_HHN.data*1000, 'k')
        ax2.text(x=((O.time+arrP.time)-10).matplotlib_date,y=trace_HHN_amp_max-(trace_HHN_amp_max*0.3),s='P',bbox=dict(facecolor='w', alpha=0.2,edgecolor='w'))

        ax2.set_ylim(-trace_HHN_amp_max,trace_HHN_amp_max)
        ax2.text(0.05, 0.1, trace_HHN.stats.channel, horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes)
        ax2.ticklabel_format(axis='y', style='sci')
        if ind == 0:
            ax2.set_ylabel("Vel.(mm/s)")
        plt.setp(ax2.get_xticklabels(), visible=False)

        #----------------------------------------------------------------------------

        ax3 = fig.add_subplot(gs[4,ind],sharex=ax2)
        trace_HHE = stE2[ind]
        trace_HHE_amp_max = np.max(np.abs(trace_HHE.data*1000))

        ax3.xaxis.set_major_locator(axis_major)
        ax3.xaxis.set_major_formatter(axis_Fmt)
        ax3.xaxis.set_minor_locator(axis_minor)
        ax3.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax3.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax3.xaxis.label.set_color('white')
        ax3.yaxis.label.set_color('white')
        ax3.title.set_color('white')

        ax3.plot(trace_HHN.times('matplotlib'),trace_HHE.data*1000, 'k')
        ax3.text(x=((O.time+arrP.time)-10).matplotlib_date,y=trace_HHE_amp_max-(trace_HHE_amp_max*0.3),s='P',bbox=dict(facecolor='w', alpha=0.2,edgecolor='w'))

        ax3.set_ylim(-trace_HHE_amp_max,trace_HHE_amp_max)
        ax3.text(0.05, 0.1, trace_HHE.stats.channel, horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
        ax3.ticklabel_format(axis='y', style='sci')
        if ind == 0:
            ax3.set_ylabel("Vel.(mm/s)")
        plt.setp(ax3.get_xticklabels(), visible=False)

        #----------------------------------------------------------------------------
        ax4 = fig.add_subplot(gs[5,ind],sharex=ax3)
        trace_HHX = stX2[ind]
        trace_HHX_amp_max = np.max(np.abs(trace_HHX.data/1000))

        ax4.xaxis.set_major_locator(axis_major)
        ax4.xaxis.set_major_formatter(axis_Fmt)
        ax4.xaxis.set_minor_locator(axis_minor)
        ax4.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax4.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax4.xaxis.label.set_color('white')
        ax4.yaxis.label.set_color('white')
        ax4.title.set_color('white')

        ax4.plot(trace_HHX.times('matplotlib'),trace_HHX.data/1000, 'k')
        ax4.text(x=((O.time+arrP.time)-10).matplotlib_date,y=trace_HHX_amp_max-(trace_HHX_amp_max*0.3),s='P',bbox=dict(facecolor='w', alpha=0.1,edgecolor='w'))

        ax4.set_ylim(-trace_HHX_amp_max,trace_HHX_amp_max)
        ax4.text(0.05, 0.1, trace_HHX.stats.channel, horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes)
        ax4.ticklabel_format(axis='y', style='sci')
        if ind == 0:
            ax4.set_ylabel("Digital Units (10³)")
        plt.setp(ax4.get_xticklabels(), visible=False)

        #------------------------------------------------------------------------------------------
        ax6 = fig.add_subplot(gs[6,ind],sharex=ax1)

        cft = classic_sta_lta(traces.data, int(stalen * df), int(ltalen * df))
        triggers = trigger_onset(cft, trig_on, trig_off)
        t = np.arange(npts, dtype=np.float32)/df

        on_off = np.array(trigger_onset(cft, trig_on, trig_off))
        time_on = traces.stats.starttime+float((on_off[:,0]/df)[0])
        time_off = traces.stats.starttime+float((on_off[:,1]/df)[0])

        #ax1.axvline(time_on.matplotlib_date, color='blue', lw=2, ls='--')
        ax1.axvline((O.time+arrP.time).matplotlib_date, color='r', lw=2, ls='--')

        #ax2.axvline(time_on.matplotlib_date, color='blue', lw=2, ls='--')
        ax2.axvline((O.time+arrP.time).matplotlib_date, color='r', lw=2, ls='--')

        #ax3.axvline(time_on.matplotlib_date, color='blue', lw=2, ls='--')
        ax3.axvline((O.time+arrP.time).matplotlib_date, color='r', lw=2, ls='--')

        #ax4.axvline(time_on.matplotlib_date, color='blue', lw=2, ls='--')
        ax4.axvline((O.time+arrP.time).matplotlib_date, color='r', lw=2, ls='--')


        ax6.plot(traces.times('matplotlib'), cft, 'k')
        ax6.axhline(trig_on, color='red', lw=1, ls='--')
        ax6.axhline(trig_off, color='blue', lw=1, ls='--')

        ax6.tick_params(axis='both',which='major',width=2,length=5,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax6.tick_params(axis='both',which='minor',width=2,length=3,bottom=True, top=True, left=True, right=True,color='white',labelcolor='white')
        ax6.xaxis.label.set_color('white')
        ax6.yaxis.label.set_color('white')
        ax6.title.set_color('white')

        if ind == 0:
            ax6.set_ylabel("STA/LTA")
        #-------------------------------------------------------------------------------------------------------------------------------
        fig.autofmt_xdate()

#-------------------------------------------------------------------------------
daily_event_output = EARTHQUAKE_FINDER_OUTPUT+'/EVENTS_SELECTED_PLOT_STA_PPSD/'
os.makedirs(daily_event_output,exist_ok=True)
fig.savefig(daily_event_output+NETWORK+'_'+traces.stats.starttime.strftime('%Y_%m_%d_%H_%M_%S_%f')+'.png', facecolor='None',dpi=300,bbox_inches='tight',pad_inches=0.3)
plt.close()

#-------------------------------------------------------------------------------

plt.rcParams.update({'font.size': 15})
# map
#--------------------------------------------------------
client = Client("IRIS")
ev = client.get_events(starttime=FDAY-20*60, endtime=FDAY+210, minmagnitude=6)
model = TauPyModel(model="iasp91")
O = ev[0].preferred_origin()
mag = ev[0].magnitudes[0].mag
mag_type = ev[0].magnitudes[0].magnitude_type
crs = ccrs.Orthographic(central_longitude=O.longitude, central_latitude=O.latitude)
fig = plt.figure(figsize=[10, 5])
ax1 = fig.add_subplot(1, 1, 1, projection=crs)
#-------------------------------------------------------------------------------------------

for ind,traces in enumerate(stZ2):
    d = op.geodetics.kilometer2degrees(op.geodetics.calc_vincenty_inverse(traces.stats.coordinates.latitude, traces.stats.coordinates.longitude, O.latitude, O.longitude)[0]/1000.)
    T = model.get_travel_times(O.depth / 1000., d)
    arrP = T[0]
    dist_km = op.geodetics.calc_vincenty_inverse(traces.stats.coordinates.latitude, traces.stats.coordinates.longitude, O.latitude, O.longitude)[0]/1000.
    #-------------------------------------------------------------------------------------------
    ax1.set_global()
    ax1.scatter(traces.stats.coordinates.longitude,traces.stats.coordinates.latitude,marker='^', color='k', s=50,transform=ccrs.PlateCarree(),zorder=10)
    ax1.scatter(O.longitude,O.latitude ,marker='*', color='y', s=50, transform=ccrs.PlateCarree(),zorder=10)
    ax1.plot([traces.stats.coordinates.longitude, O.longitude], [traces.stats.coordinates.latitude, O.latitude],c='grey', transform=ccrs.Geodetic(),zorder=1)
    ax1.coastlines()
    ax1.stock_img()
    #--------------------------------------------------------
#-------------------------------------------------------------------------------
daily_event_output = EARTHQUAKE_FINDER_OUTPUT+'/EVENTS_SELECTED_PLOT_STA_PPSD/'
os.makedirs(daily_event_output,exist_ok=True)
fig.savefig(daily_event_output+NETWORK+'map_event_'+O.time.strftime('%Y_%m_%d_%H_%M_%S_%f')+'.png', facecolor='None',dpi=300,bbox_inches='tight',pad_inches=0.3)
plt.close()
