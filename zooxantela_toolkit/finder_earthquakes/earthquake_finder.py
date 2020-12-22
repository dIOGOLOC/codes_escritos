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

from mtspec import mtspec, wigner_ville_spectrum

from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta


#Configuration file

MSEED_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_OBS/'

EARTHQUAKE_FINDER_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/EARTHQUAKE_FINDER_OUTPUT/ASDF_FILES/'

NOISE_MODEL_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/TRANSFER_FUNC/NOISE_MODEL_FILE/noise_models.npz'

FIRSTDAY = '2020-03-22'
LASTDAY = '2020-03-27'

FILTER_DATA = [2,10]

NETWORK = 'ON'

STATION = 'OBS18'

CHANNEL = 'HHZ'

MIN_WINDOWS = 24

WINDOW_LENGTH = 600

ALPHA = 0.05

TOL = 3.0

VERBOSE_MODE = True

STA = 0.5
LTA = 10
THRON = 15
THROFF = 0.5

# ========================
# Constants and parameters
# ========================

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)


# ================
# MULTIPROCESSING
# ================
num_processes = 12

# =================
# Filtering by date
# =================

fday = UTCDateTime(FIRSTDAY)
lday = UTCDateTime(LASTDAY)
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

#-------------------------------------------------------------------------------

def ftest(res1, pars1, res2, pars2):

    from scipy.stats import f as f_dist

    N1 = len(res1)
    N2 = len(res2)

    dof1 = N1 - pars1
    dof2 = N2 - pars2

    Ea_1 = np.sum(res1**2)
    Ea_2 = np.sum(res2**2)

    Fobs = (Ea_1/dof1)/(Ea_2/dof2)

    P = 1. - (f_dist.cdf(Fobs, dof1, dof2) - f_dist.cdf(1./Fobs, dof1, dof2))

    return P

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

def rotate_dir(tr1, tr2, direc):

    d = -direc*np.pi/180.+np.pi/2.
    rot_mat = np.array([[np.cos(d), -np.sin(d)],
                        [np.sin(d), np.cos(d)]])

    v12 = np.array([tr2, tr1])
    vxy = np.tensordot(rot_mat, v12, axes=1)
    tr_2 = vxy[0, :]
    tr_1 = vxy[1, :]

    return tr_1

#-------------------------------------------------------------------------------
def classic_sta_lta_py(a, nsta, nlta):
    """
    Computes the standard STA/LTA from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples. Written in Python by Obspy.
    """
    m = len(a)
    #
    # compute the short time average (STA) and long time average (LTA)
    sta = np.zeros(m, dtype=np.float64)
    lta = np.zeros(m, dtype=np.float64)

    for i in range(m):
        sta[i] = np.mean(a[i:int(i+nsta)]**2)
        lta[i] = np.mean(a[i:int(i+nlta)]**2)

    # Pad zeros
    sta[:nlta - 1] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(0.0).tiny
    idx = lta < dtiny
    lta[idx] = dtiny

    return sta / lta

#-------------------------------------------------------------------------------
def variance_sta_lta_py(a, nsta, nlta):
    """
    Computes the variance STA/LTA from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples. Written in Python by Kumomeperkoch et al.(2010).
    """
    m = len(a)
    #
    # compute the short time average (STA) and long time average (LTA)
    sta = np.zeros(m, dtype=np.float64)
    lta = np.zeros(m, dtype=np.float64)

    for i in range(m):
        sta[i] = np.var(a[i:int(i+nsta)])
        lta[i] = np.var(a[i:int(i+nlta)])

    # Pad zeros
    sta[:nlta - 1] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(0.0).tiny
    idx = lta < dtiny
    lta[idx] = dtiny

    return sta / lta

#-------------------------------------------------------------------------------

def skew_sta_lta_py(a, nlta):
    """
    Computes the skew from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples. Written in Python by Kumomeperkoch et al.(2010).
    """
    m = len(a)
    #
    # compute the short time average (STA) and long time average (LTA)
    sta = np.zeros(m, dtype=np.float64)
    lta = np.zeros(m, dtype=np.float64)

    for i in range(m):
        lta[i] = skew(a[i:int(i+nlta)])

    return lta

#-------------------------------------------------------------------------------

def kurtosis_sta_lta_py(a, nlta):
    """
    Computes the kurtosis from a given input array a. The length of
    the STA is given by nsta in samples, respectively is the length of the
    LTA given by nlta in samples. Written in Python by Kumomeperkoch et al.(2010).
    """
    m = len(a)
    #
    # compute the short time average (STA) and long time average (LTA)
    sta = np.zeros(m, dtype=np.float64)
    lta = np.zeros(m, dtype=np.float64)

    for i in range(m):
        lta[i] = kurtosis(a[i:int(i+nlta)])
    return lta

# =======
# Classes
# =======

def get_stations_data(f):
    """
    Gets stations daily data from miniseed file and convert in ASDF

    @type f: paht of the minissed file (str)
    @rtype: list of L{StationDayData}
    """

    # splitting subdir/basename
    subdir, filename = os.path.split(f)

    # network, station name and station channel in basename,
    # e.g., ON.TIJ01..HHZ.D.2020.002

    network, name = filename.split('.')[0:2]
    sta_channel_id = filename.split('.D.')[0]
    channel = sta_channel_id.split('..')[-1]
    time_day = filename.split('.D.')[-1]
    year_day = time_day.split('.')[0]
    julday_day = time_day.split('.')[1]

    st = read(f)
    inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))

    pre_filt = [0.001, 0.005, 45., 50.]
    st[0].remove_response(inventory=inv,pre_filt=pre_filt,output="DISP",water_level=60)

    st.detrend('demean')
    st.detrend('linear')
    st.taper(max_percentage=0.05, type="hann")

    #------------------------------------------------------------------
	# Points in window
    ws = int(WINDOW_LENGTH/st[0].stats.delta)

    # Number of points to overlap
    ss = int(WINDOW_LENGTH*0.5/st[0].stats.delta)

    # hanning window
    wind = np.ones(ws)

    # Calculating the spectrogram
    spec = spectrogram(x=st[0].data,fs=st[0].stats.sampling_rate, window=wind, nperseg=ws, noverlap=ss)
    f, t, psd = spec

    #-----------------------------------------------------------------
    #Creating ASDF preprocessed files folder
    output_PREPROCESS_DATA_DAY = ASDF_FILES+'PREPROCESS_DATA_DAY_FILES/'+NETWORK+'.'+STATION+'.'+CHANNEL+'/'
    os.makedirs(output_PREPROCESS_DATA_DAY,exist_ok=True)

    ds = ASDFDataSet(output_PREPROCESS_DATA_DAY+'PREPROCESS_DATA_DAY_'+NETWORK+'_'+STATION+'_'+channel+'_'+year_day+'_'+julday_day+".h5", compression="gzip-3")
    ds.add_waveforms(st, tag="preprocessed_recording")
    ds.add_stationxml(STATIONXML_DIR+'.'.join([NETWORK,STATION,'xml']))

    # The type always should be camel case.
    data_type_f = "Frequencies"

    # Name to identify the particular piece of data.
    path_f = sta_channel_id+'.f'

    # Any additional parameters as a Python dictionary which will end up as
    # attributes of the array.
    parameters_f = {'f':'Array of sample frequencies.',
        		  'sampling_rate_in_hz': st[0].stats.sampling_rate,
        		  'station_id': sta_channel_id}


    ds.add_auxiliary_data(data=f, data_type=data_type_f, path=path_f, parameters=parameters_f)

    # The type always should be camel case.
    data_type_t = "Times"

    # Name to identify the particular piece of data.
    path_t = sta_channel_id+'.t'

    # Any additional parameters as a Python dictionary which will end up as
    # attributes of the array.
    parameters_t = {'t':'Array of segment times.',
        		  'sampling_rate_in_hz': st[0].stats.sampling_rate,
        		  'station_id': sta_channel_id}


    ds.add_auxiliary_data(data=t, data_type=data_type_t, path=path_t, parameters=parameters_t)

    # The type always should be camel case.
    data_type_s = "Spectrogram"

    # Name to identify the particular piece of data.
    path_s = sta_channel_id+'.S'

    # Any additional parameters as a Python dictionary which will end up as
    # attributes of the array.
    parameters_s = {'S':'Spectrogram of x. By default, the last axis of S corresponds to the segment times.',
        		  'sampling_rate_in_hz': st[0].stats.sampling_rate,
        		  'station_id': sta_channel_id}


    ds.add_auxiliary_data(data=psd, data_type=data_type_s, path=path_s, parameters=parameters_s)

# ============
# Main program
# ============

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

# initializing list of stations by scanning name of miniseed files

files = filelist(basedir=MSEED_DIR+NETWORK+'/'+STATION+'/'+CHANNEL+'.D/',interval_period_date=INTERVAL_PERIOD_DATE)

print('Total of miniseed files = '+str(len(files)))
print('\n')

print('============================================================')
print('Opening miniseed files, preprocessing and converting to ASDF')
print('============================================================')
print('\n')

start_time = time.time()

with Pool(processes=num_processes) as p:
	max_ = len(files)
	with tqdm(total=max_) as pbar:
		for i, _ in enumerate(p.imap_unordered(get_stations_data, files)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

print('====================================')
print('Opening ASDF files')
print('====================================')
print('\n')

daily_lst_data = [[]]*len(INTERVAL_PERIOD_DATE)

DATA_DAY_FILES = sorted(glob.glob(ASDF_FILES+'PREPROCESS_DATA_DAY_FILES/'+NETWORK+'.'+STATION+'.'+CHANNEL+'/*.h5'))
DATA_DAYs = [i.split('.h5')[0].split('_')[-2]+'.'+i.split('.h5')[0].split('_')[-1] for i in DATA_DAY_FILES]

for l,k in enumerate(INTERVAL_PERIOD_DATE):
	daily_lst_data[l] = [j for i,j in enumerate(DATA_DAY_FILES) if DATA_DAYs[i] == k]

#-------------------------------------------------------------------------------
# Select bandpass frequencies

print('=======================')
print('Filtering daily windows')
print('=======================')
print('\n')

first_filter_day = []
good_windows_lst = []

for j in tqdm(daily_lst_data, desc='Daily loop'):

    spec = [ASDFDataSet(k) for k in j][0]

    sta_id = spec.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].id

    year_spec = str(spec.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.year)
    julday_spec = str(spec.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.julday)

    f = spec.auxiliary_data.Frequencies[sta_id+'.f'].data[:]
    t = spec.auxiliary_data.Times[sta_id+'.t'].data[:]
    psd = spec.auxiliary_data.Spectrogram[sta_id+'.S'].data[:]

    th = np.array([i/(60*60) for i in t])

    ff = (f > FILTER_DATA[0]) & (f < FILTER_DATA[1])
    freq_lst = f[ff]

	# avoid calculating log of zero
    idx = psd < DTINY
    psd[idx] = DTINY

    # go to dB
    log_psd = np.log10(psd)
    log_psd *= 10
    log_psd = smooth(log_psd,40)
    dsl = log_psd[ff, :] - np.mean(log_psd[ff, :], axis=0)

	#--------------------------------------------------------------

    good_windows = np.repeat([True], len(th))
    indwin = np.argwhere(good_windows == True)
    moveon = False

    while moveon == False:

        normvar = np.zeros(np.sum(good_windows))
        for ii, tmp in enumerate(indwin):
            ind = np.copy(indwin)
            ind = np.delete(ind, ii)
            normvar[ii] = norm(np.std(dsl[:, ind], axis=1), ord=2)
        ubernorm = np.median(normvar) - normvar

        eliminated_windows = ubernorm > TOL*np.std(ubernorm)

        if np.sum(eliminated_windows) == 0:
            moveon = True

        trypenalty = eliminated_windows[np.argwhere(eliminated_windows == False)].T[0]

        if ftest(eliminated_windows, 1, trypenalty, 1) < ALPHA:
            good_windows[indwin[eliminated_windows == True]] = False
            indwin = np.argwhere(good_windows == True)
            moveon = False
        else:
            moveon = True

    bad_windows = np.array([False if i == True else True for i in good_windows])

    #------------------------------------------------------------------------------

    if VERBOSE_MODE:

	        cmap = 'Spectral_r'
	        heights = [1,0.1]
	        widths = [1]

	        gs_kw = dict(width_ratios=widths, height_ratios=heights)
	        fig, (ax1,ax4) = plt.subplots(ncols=1, nrows=2,figsize=(15,10),sharex=True,gridspec_kw=gs_kw)

	        im = ax1.pcolor(th,freq_lst, dsl, cmap=cmap,shading='auto')
	        ax1.text(0.96, 0.8, CHANNEL, ha='center',bbox=dict(facecolor='w'),transform=ax1.transAxes)
	        ax1.xaxis.set_major_locator(MultipleLocator(4))
	        ax1.yaxis.set_major_locator(MultipleLocator(4))
	        ax1.set_ylabel('Frequency (Hz)')

	        axins = inset_axes(ax1,
	                       width="25%",  # width = 10% of parent_bbox width
	                       height="5%",  # height : 50%
	                       loc='upper left',
	                       bbox_to_anchor=(0.75, 0.1, 1, 1),
	                       bbox_transform=ax1.transAxes,
	                       borderpad=0,
	                       )
	        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

	        #-----------------------------------------------------


	        colors = ['darkred' if eli == True else 'white' for eli in bad_windows]
	        ax4.bar(th,height=1,width=WINDOW_LENGTH/3600,align='center',color=colors,edgecolor='k')

	        ax4.set_xlabel('Time (h)')
	        ax4.xaxis.set_major_locator(MultipleLocator(4))
	        ax4.set_ylim(0,1)
	        ax4.set_xlim(0,24)
	        ax4.set_yticks([])
	        #-----------------------------------------------------
	        daily_wind_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'.'+STATION+'/Daily_data_windows/'
	        os.makedirs(daily_wind_output,exist_ok=True)
	        fig.suptitle('Station = '+STATION+'.'+CHANNEL+' - Day = '+UTCDateTime(year=int(year_spec),julday=int(julday_spec)).strftime('%d/%m/%Y'),fontsize=18)
	        fig.savefig(daily_wind_output+NETWORK+'_'+STATION+'_'+CHANNEL+'_'+year_spec+'_'+julday_spec+'.png', dpi=300, facecolor='w', edgecolor='w')
	        plt.close()

    #----------------------------------------------------------------------------------------------------------------

    #==================
    #Looking for EVENTs
    #==================

    slide_data = [k for k in spec.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'].slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH/2,include_partial_windows=False, nearest_sample=True)]

    window_data_to_find = list(compress(slide_data, bad_windows))

    for i,wind_data in enumerate(tqdm(window_data_to_find, desc='Windows loop')):
        tr = wind_data[0]
        npts = tr.stats.npts
        dt = tr.stats.delta
        df = tr.stats.sampling_rate
        datetime_window = tr.stats.starttime
        tr.filter("bandpass", freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])

        # Characteristic function and trigger onsets
        cft = classic_sta_lta(tr.data, int(STA*df), int(LTA*df))
        on_off = np.array(trigger_onset(charfct=cft, thres1=THRON, thres2=THROFF, max_len=9e+99, max_len_delete=False))

        try:
            for trg in on_off:
                trigger_on = datetime_window+int(trg[0])/df
                trigger_off = datetime_window+int(trg[1])/df

                #----------------------------------------------------------------------------

                tr_data = wind_data[0]
                tr_trim = tr_data.copy()
                trim_data = tr_trim.trim(trigger_on-10,trigger_off+10)
                trim_data.filter("bandpass", freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])

                t = trim_data.times('matplotlib')
                f_min = FILTER_DATA[0]
                f_max = FILTER_DATA[1]

                scalogram = cwt(trim_data.data, dt, 8, f_min, f_max)
                x, y = np.meshgrid(
                t,
                np.linspace(f_min, f_max, scalogram.shape[0]))

                #----------------------------------------------------------------------------
                #Creating ASDF preprocessed files folder
                output_EVENT_DATA = ASDF_FILES+'EVENT_DATA_FILES/'+NETWORK+'.'+STATION+'.'+CHANNEL+'/'
                os.makedirs(output_EVENT_DATA,exist_ok=True)

                event_asdf = ASDFDataSet(output_EVENT_DATA+'EVENT_DATA_'+NETWORK+'_'+STATION+'_'+CHANNEL+'_'+datetime_window.strftime('%d_%m_%Y_%H_%M_%S_%f')+"_event.h5", compression="gzip-3")
                tr_2_save = wind_data[0]
                tr_2_save_trim = tr_2_save.copy()
                trim_tr_2_save = tr_2_save_trim.trim(trigger_on-10,trigger_off+10)
                event_asdf.add_waveforms(trim_tr_2_save, tag="event_recording")

                #----------------------------------------------------------------------------
                if VERBOSE_MODE:

                    # Plotting the results
                    axis_major = MinuteLocator(interval=2)   # every 2-minutes
                    axis_minor = SecondLocator(interval=10)   # every 10-second
                    axis_Fmt = DateFormatter('%H:%M:%S')

                    plt.rcParams.update({'font.size': 20})
                    fig, (ax1,ax2) = plt.subplots(ncols=1, nrows=2,figsize=(15,20),sharex=True)

                    ax1.set_title('Date: '+datetime_window.strftime('%d/%m/%Y'))
                    ax1.xaxis.set_major_locator(axis_major)
                    ax1.xaxis.set_major_formatter(axis_Fmt)
                    ax1.xaxis.set_minor_locator(axis_minor)
                    ax1.tick_params(axis='both',which='major',width=2,length=5)
                    ax1.tick_params(axis='both',which='minor',width=2,length=3)
                    ax1.plot(tr.times('matplotlib')[500:],tr.data[500:], 'k')
                    ymin, ymax = ax1.get_ylim()
                    ax1.vlines(trigger_on.matplotlib_date, ymin, ymax, color='r', linewidth=1)
                    ax1.vlines(trigger_off.matplotlib_date, ymin, ymax, color='b', linewidth=1)

                    #----------------------------------------------------------------------------

                    ax2.set_title('Absolute Mean Value STA/LTA trigger')
                    ax2.xaxis.set_major_locator(axis_major)
                    ax2.xaxis.set_major_formatter(axis_Fmt)
                    ax2.xaxis.set_minor_locator(axis_minor)
                    ax2.tick_params(axis='both',which='major',width=2,length=5)
                    ax2.tick_params(axis='both',which='minor',width=2,length=3)
                    ax2.set_ylim(0,20)
                    ax2.plot(tr.times('matplotlib'),cft, 'k')
                    ax2.axhline(THROFF, 0, 1, color='r', linestyle='--')
                    ax2.axhline(THRON, 0, 1, color='b', linestyle='--')

                    #----------------------------------------------------------------------------

                    daily_event_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'.'+STATION+'/Daily_event_data_windows/'+CHANNEL+'/'
                    os.makedirs(daily_event_output,exist_ok=True)
                    fig.savefig(daily_event_output+NETWORK+'_'+STATION+'_'+CHANNEL+'_'+trigger_on.strftime('%d_%m_%Y_%H_%M_%S_%f')+'.png', dpi=300, facecolor='w', edgecolor='w')
                    plt.close()


                    #==============================================================================================================================================================

                    # Plotting the results
                    axis_major = SecondLocator(interval=5)   # every 5-second
                    axis_minor = SecondLocator(interval=1)   # every 1-second

                    plt.rcParams.update({'font.size': 20})
                    fig, (ax1,ax2,ax3) = plt.subplots(ncols=1, nrows=3,figsize=(15,20),sharex=True)

                    ax1.set_title('Date: '+datetime_window.strftime('%d/%m/%Y')+' - Filter: '+str(FILTER_DATA[0])+'-'+str(FILTER_DATA[1])+' Hz')
                    ax1.xaxis.set_major_formatter(axis_Fmt)
                    ax1.xaxis.set_major_locator(axis_major)
                    ax1.xaxis.set_minor_locator(axis_minor)
                    ax1.tick_params(axis='both',which='major',width=2,length=5)
                    ax1.tick_params(axis='both',which='minor',width=2,length=3)
                    ax1.plot(trim_data.times('matplotlib'),trim_data.data, color='k', linewidth=2)

                    #----------------------------------------------------------------------------

                    ax2.set_title('Absolute Mean Value STA/LTA trigger')
                    ax2.xaxis.set_major_locator(axis_major)
                    ax2.xaxis.set_major_formatter(axis_Fmt)
                    ax2.xaxis.set_minor_locator(axis_minor)
                    ax2.tick_params(axis='both',which='major',width=2,length=5)
                    ax2.tick_params(axis='both',which='minor',width=2,length=3)
                    ax2.set_ylim(0,20)
                    ax2.plot(tr.times('matplotlib'),cft, 'k')
                    ax2.axhline(THROFF, 0, 1, color='r', linestyle='--')
                    ax2.axhline(THRON, 0, 1, color='b', linestyle='--')
                    #----------------------------------------------------------------------------

                    ax3.set_title('Continuous Wavelet Transform')
                    ax3.xaxis.set_major_formatter(axis_Fmt)
                    ax3.xaxis.set_major_locator(axis_major)
                    ax3.xaxis.set_minor_locator(axis_minor)
                    ax3.tick_params(axis='both',which='major',width=2,length=5)
                    ax3.tick_params(axis='both',which='minor',width=2,length=3)

                    im = ax3.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
                    ax3.set_ylabel("Frequency [Hz]")
                    ax3.set_ylim(f_min, f_max)

                    axins = inset_axes(ax3,
                            		   width="15%",
                            		   height="5%",
                            		   loc='upper left',
                            		   bbox_to_anchor=(0.8, 0.1, 1, 1),
                            		   bbox_transform=ax3.transAxes,
                            		   borderpad=0,
                            		  )

                    plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

                    daily_event_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'.'+STATION+'/Daily_event_data_windows/'+CHANNEL+'/'
                    os.makedirs(daily_event_output,exist_ok=True)
                    fig.savefig(daily_event_output+NETWORK+'_'+STATION+'_'+CHANNEL+'_'+trigger_on.strftime('%d_%m_%Y_%H_%M_%S_%f')+'_trim.png', dpi=300, facecolor='w', edgecolor='w')
                    plt.close()

        except IndexError:
            pass

print('========================')
print('Opening ASDF Event files')
print('========================')
print('\n')

EVENT_FILES_H5 = sorted(glob.glob(ASDF_FILES+'EVENT_DATA_FILES/'+NETWORK+'.'+STATION+'.'+CHANNEL+'/*.h5'))

for data_h5 in tqdm(EVENT_FILES_H5,desc='Events loop'):
    event_data = ASDFDataSet(data_h5)

    for i in event_data.waveforms[NETWORK+'.'+STATION]['event_recording']:
        tr = i
        npts = tr.stats.npts
        dt = tr.stats.delta
        df = tr.stats.sampling_rate
        datetime_window = tr.stats.starttime
        tr.filter(type="bandpass",freqmin=FILTER_DATA[0],freqmax=FILTER_DATA[1],zerophase=True)

        # Characteristic function and trigger onsets
        cft = classic_sta_lta_py(tr.data, int((STA/2)*df), int((LTA/2)*df))
        cft_var = variance_sta_lta_py(tr.data, int((STA/2)*df), int((LTA/2)*df))
        cft_skew = skew_sta_lta_py(tr.data, int((LTA/2)*df))
        cft_kur = kurtosis_sta_lta_py(tr.data, int((LTA/2)*df))

        #----------------------------------------------------------------------------

        tr_data = i
        trim_data = tr_data.copy()
        #trim_data.filter("bandpass", freqmin=FILTER_DATA[0], freqmax=FILTER_DATA[1])

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

        # Plotting the results
        axis_major = SecondLocator(interval=5)   # every 5-second
        axis_minor = SecondLocator(interval=1)   # every 1-second
        axis_Fmt = DateFormatter('%H:%M:%S')

        plt.rcParams.update({'font.size': 20})
        fig, axes = plt.subplots(ncols=1, nrows=5,figsize=(20,20),sharex=True)

        ax1 = axes[0]
        ax2 = axes[1]
        ax3 = axes[2]
        ax4 = axes[3]
        ax5 = axes[4]
        #ax6 = axes[2,0]

        tr_raw_data = i
        ax1.set_title('Raw Data')
        ax1.xaxis.set_major_formatter(axis_Fmt)
        ax1.xaxis.set_major_locator(axis_major)
        ax1.xaxis.set_minor_locator(axis_minor)
        ax1.tick_params(axis='both',which='major',width=2,length=5)
        ax1.tick_params(axis='both',which='minor',width=2,length=3)
        ax1.plot(tr_raw_data.times('matplotlib'),tr_raw_data.data, color='k', linewidth=2)

        #----------------------------------------------------------------------------

        ax2.set_title('Absolute Mean Value STA/LTA trigger')
        ax2.xaxis.set_major_locator(axis_major)
        ax2.xaxis.set_major_formatter(axis_Fmt)
        ax2.xaxis.set_minor_locator(axis_minor)
        ax2.tick_params(axis='both',which='major',width=2,length=5)
        ax2.tick_params(axis='both',which='minor',width=2,length=3)
        ax2.plot(tr_raw_data.times('matplotlib'),cft, 'k')

        #----------------------------------------------------------------------------

        ax3.set_title('Variance STA/LTA trigger')
        ax3.xaxis.set_major_locator(axis_major)
        ax3.xaxis.set_major_formatter(axis_Fmt)
        ax3.xaxis.set_minor_locator(axis_minor)
        ax3.tick_params(axis='both',which='major',width=2,length=5)
        ax3.tick_params(axis='both',which='minor',width=2,length=3)
        ax3.plot(tr_raw_data.times('matplotlib'),cft_var, 'k')

        #----------------------------------------------------------------------------

        ax4.set_title('Kurtosis trigger')
        ax4.xaxis.set_major_locator(axis_major)
        ax4.xaxis.set_major_formatter(axis_Fmt)
        ax4.xaxis.set_minor_locator(axis_minor)
        ax4.tick_params(axis='both',which='major',width=2,length=5)
        ax4.tick_params(axis='both',which='minor',width=2,length=3)
        ax4.plot(tr_raw_data.times('matplotlib'),cft_kur, 'k')

        #----------------------------------------------------------------------------

        ax5.set_title('Skew trigger')
        ax5.xaxis.set_major_locator(axis_major)
        ax5.xaxis.set_major_formatter(axis_Fmt)
        ax5.xaxis.set_minor_locator(axis_minor)
        ax5.tick_params(axis='both',which='major',width=2,length=5)
        ax5.tick_params(axis='both',which='minor',width=2,length=3)
        ax5.plot(tr_raw_data.times('matplotlib'),cft_skew, 'k')
        ax5.set_xlabel('Date: '+datetime_window.strftime('%d/%m/%Y'))


        #----------------------------------------------------------------------------
        '''
        ax6.set_title('Continuous Wavelet Transform')
        ax6.xaxis.set_major_formatter(axis_Fmt)
        ax6.xaxis.set_major_locator(axis_major)
        ax6.xaxis.set_minor_locator(axis_minor)
        ax6.tick_params(axis='both',which='major',width=2,length=5)
        ax6.tick_params(axis='both',which='minor',width=2,length=3)

        im = ax6.pcolormesh(x, y, np.abs(scalogram),shading='auto', cmap='viridis')
        ax6.set_ylabel("Frequency [Hz]")
        ax6.set_xlabel('Date: '+datetime_window.strftime('%d/%m/%Y'))

        axins = inset_axes(ax6,
                           width="15%",
                           height="5%",
                           loc='upper left',
                           bbox_to_anchor=(0.8, 0.1, 1, 1),
                           bbox_transform=ax6.transAxes,
                           borderpad=0,
                           )

        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')
        '''
        daily_event_output = EARTHQUAKE_FINDER_OUTPUT+NETWORK+'.'+STATION+'/Selected_event_data/'+CHANNEL+'/'
        os.makedirs(daily_event_output,exist_ok=True)
        fig.savefig(daily_event_output+NETWORK+'_'+STATION+'_'+CHANNEL+'_'+datetime_window.strftime('%d_%m_%Y_%H_%M_%S_%f')+'_trim.png', dpi=300, facecolor='w', edgecolor='w')
        plt.close()
