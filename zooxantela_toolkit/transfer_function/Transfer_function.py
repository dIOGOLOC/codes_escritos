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

#Configuration file

MSEED_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_OBS/'

TRANSFER_FUNC_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/TRANSFER_FUNC/FIGURAS/'

CORRECT_DATA_TRANSFER_FUNC_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/TRANSFER_FUNC/DATA_CORRECTION/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/TRANSFER_FUNC/ASDF_FILES/'

JSON_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/TRANSFER_FUNC/JSON_FILES/'

NOISE_MODEL_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/TRANSFER_FUNC/NOISE_MODEL_FILE/noise_models.npz'

FIRSTDAY = '2019-07-27'
LASTDAY = '2020-01-27'

NETWORK = 'ON'

STATION = 'OBS17'

MIN_WINDOWS = 24

WINDOW_LENGTH = 600

FILTER_DATA = [0.004,0.2]

NEW_SAMPLING_RATE = 100

tiltfreq=[0.005, 0.035]

ALPHA = 0.05
TOL = 2.0

VERBOSE_MODE = True

# ========================
# Constants and parameters
# ========================

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)


# ================
# MULTIPROCESSING
# ================
num_processes = 6

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
    flist = os.listdir(basedir)
    for d in flist:
        if os.path.isdir(os.path.join(basedir, d)):
            files_dir = os.path.join(basedir, d)
            files_list = glob.glob(files_dir+'/*')
            for s in files_list:
                if any(day_s in s for day_s in interval_period_date):
                    files.append(s)

    files = [i for i in files if not 'HHX' in i]

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

def Admittance_xy(G_xy,G_xx,G_yy,wind_number):

    A_xy = abs(G_xy)/G_xx

    Error_A_xy = np.sqrt(1 - (np.abs(G_xy)**2)/(G_xx*G_yy))/np.sqrt(2*wind_number)*abs(np.sqrt((np.abs(G_xy)**2)/(G_xx*G_yy)))

    return A_xy,Error_A_xy

#-------------------------------------------------------------------------------

def Coherence_xy(G_xy,G_xx,G_yy,wind_number):

    Gama2_xy = (np.abs(G_xy)**2)/(G_xx*G_yy)

    Error_Gama2_xy = np.sqrt(2)*(1 - np.abs(G_xy)**2)/(G_xx*G_yy)/np.sqrt(wind_number)*abs(np.sqrt((np.abs(G_xy)**2)/(G_xx*G_yy)))

    return Gama2_xy,Error_Gama2_xy

#-------------------------------------------------------------------------------

def Phase_xy(Q_xy,C_xy,G_xy,G_xx,G_yy,wind_number):

    Phase_xy = np.arctan2(Q_xy,C_xy)

    Error_phase_xy = np.sqrt(1 - (np.abs(G_xy)**2)/(G_xx*G_yy))/np.sqrt(2*wind_number)*abs(np.sqrt((np.abs(G_xy)**2)/(G_xx*G_yy)))


    return Phase_xy,Error_phase_xy

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

def daily_process_function(ASDF_FILES):
    """
    Gets ASDF files for each channel and estimate good windows

    @type ASDF_FILES:  list of ASDF files
    """
    spec_HHZ = [ASDFDataSet(k) for k in ASDF_FILES if 'HHZ' in k][0]
    spec_HHE = [ASDFDataSet(k) for k in ASDF_FILES if 'HHE' in k][0]
    spec_HHN = [ASDFDataSet(k) for k in ASDF_FILES if 'HHN' in k][0]

    year_spec = spec_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.year
    julday_spec = spec_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.julday

    sta_id_HHZ = spec_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].id

    fz = spec_HHZ.auxiliary_data.Frequencies[sta_id_HHZ+'.f'].data[:]
    tz = spec_HHZ.auxiliary_data.Times[sta_id_HHZ+'.t'].data[:]
    psdz = spec_HHZ.auxiliary_data.Spectrogram[sta_id_HHZ+'.S'].data[:]
    tzh = np.array([i/(60*60) for i in tz])

    #Frequency target
    ff = (fz > FILTER_DATA[0]) & (fz < FILTER_DATA[1])
    freq_lst = fz[ff]

    # avoid calculating log of zero
    idx = psdz < DTINY
    psdz[idx] = DTINY

    # go to dB
    log_psdz = np.log10(psdz)
    log_psdz *= 10
    log_psdz = smooth(log_psdz,40)
    dsl_HHZ = log_psdz[ff, :] - np.mean(log_psdz[ff, :], axis=0)

    #-----------------------------------------------------
    sta_id_HHN = spec_HHN.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].id

    fn = spec_HHN.auxiliary_data.Frequencies[sta_id_HHN+'.f'].data[:]
    tn = spec_HHN.auxiliary_data.Times[sta_id_HHN+'.t'].data[:]
    psdn = spec_HHN.auxiliary_data.Spectrogram[sta_id_HHN+'.S'].data[:]
    tnh = np.array([i/(60*60) for i in tn])

    # avoid calculating log of zero
    idx = psdn < DTINY
    psdn[idx] = DTINY

    # go to dB
    log_psdn = np.log10(psdn)
    log_psdn *= 10
    log_psdn = smooth(log_psdn,40)
    dsl_HHN = log_psdn[ff, :] - np.mean(log_psdn[ff, :], axis=0)
    #-----------------------------------------------------
    sta_id_HHE = spec_HHE.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].id

    fe = spec_HHE.auxiliary_data.Frequencies[sta_id_HHE+'.f'].data[:]
    te = spec_HHE.auxiliary_data.Times[sta_id_HHE+'.t'].data[:]
    psde = spec_HHE.auxiliary_data.Spectrogram[sta_id_HHE+'.S'].data[:]
    teh = np.array([i/(60*60) for i in te])

    # avoid calculating log of zero
    idx = psde < DTINY
    psde[idx] = DTINY

    # go to dB
    log_psde = np.log10(psde)
    log_psde *= 10
    log_psde = smooth(log_psde,40)
    dsl_HHE = log_psde[ff, :] - np.mean(log_psde[ff, :], axis=0)
    #-----------------------------------------------------

    if VERBOSE_MODE:
	    cmap = 'viridis_r'
	    heights = [1,1,1]
	    widths = [1]

	    gs_kw = dict(width_ratios=widths, height_ratios=heights)
	    fig, (ax1,ax2,ax3) = plt.subplots(ncols=1, nrows=3,figsize=(15,10),sharex=True,gridspec_kw=gs_kw)

	    im = ax1.pcolormesh(tzh,fz, log_psdz, cmap=cmap,shading='auto',vmin=-180,vmax=-60)
	    ax1.text(0.96, 0.8, 'HHZ', ha='center',bbox=dict(facecolor='w'),transform=ax1.transAxes)
	    ax1.xaxis.set_major_locator(MultipleLocator(4))
	    ax1.yaxis.set_major_locator(MultipleLocator(10))
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
	    ax2.pcolormesh(tnh,fn, log_psdn, cmap=cmap,shading='auto',vmin=-180,vmax=-60)
	    ax2.text(0.96, 0.8, 'HHN', ha='center',bbox=dict(facecolor='w'),transform=ax2.transAxes)
	    ax2.set_ylabel('Frequency (Hz)')
	    ax2.xaxis.set_major_locator(MultipleLocator(4))
	    ax2.yaxis.set_major_locator(MultipleLocator(10))
	    #-----------------------------------------------------
	    ax3.pcolormesh(tnh,fe, log_psde, cmap=cmap,shading='auto',vmin=-180,vmax=-60)
	    ax3.text(0.96, 0.8, 'HHE', ha='center',bbox=dict(facecolor='w'),transform=ax3.transAxes)
	    ax3.xaxis.set_major_locator(MultipleLocator(4))
	    ax3.yaxis.set_major_locator(MultipleLocator(10))
	    ax3.set_xlabel('Time (h)')
	    #-----------------------------------------------------
	    daily_wind_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_data_windows/'
	    os.makedirs(daily_wind_output,exist_ok=True)
	    fig.suptitle('Station = '+STATION+' - Day = '+UTCDateTime(year=year_spec,julday=julday_spec).strftime('%d/%m/%Y'),fontsize=18)
	    fig.savefig(daily_wind_output+NETWORK+'_'+STATION+'_'+str(year_spec)+'_'+str(julday_spec)+'_full.png', dpi=300, facecolor='w', edgecolor='w')
	    plt.close()
        #-------------------------------------------------------------------------------------------------------------

    good_windows = np.repeat([True], len(teh))
    indwin = np.argwhere(good_windows == True)
    moveon = False

    while moveon == False:

        normvar_HHZ = np.zeros(np.sum(good_windows))
        for ii, tmp in enumerate(indwin):
            ind = np.copy(indwin)
            ind = np.delete(ind, ii)
            normvar_HHZ[ii] = norm(np.std(dsl_HHZ[:, ind], axis=1), ord=2)
        ubernorm_HHZ = np.median(normvar_HHZ) - normvar_HHZ

        normvar_HHN = np.zeros(np.sum(good_windows))
        for ii, tmp in enumerate(indwin):
            ind = np.copy(indwin)
            ind = np.delete(ind, ii)
            normvar_HHN[ii] = norm(np.std(dsl_HHN[:, ind], axis=1), ord=2)
        ubernorm_HHN = np.median(normvar_HHN) - normvar_HHN

        normvar_HHE = np.zeros(np.sum(good_windows))
        for ii, tmp in enumerate(indwin):
            ind = np.copy(indwin)
            ind = np.delete(ind, ii)
            normvar_HHE[ii] = norm(np.std(dsl_HHE[:, ind], axis=1), ord=2)
        ubernorm_HHE = np.median(normvar_HHE) - normvar_HHE

        norm_allchannels = ubernorm_HHZ+ubernorm_HHN+ubernorm_HHE

        eliminated_windows = norm_allchannels > TOL*np.std(norm_allchannels)

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

    len_good = int(sum(good_windows))
    len_bad = int(sum(bad_windows))

    #------------------------------------------------------------------------------

    if len_good > MIN_WINDOWS:

        first_filter_day = True

        dsl_HHZs = dsl_HHZ[:,good_windows]
        dsl_HHNs = dsl_HHN[:,good_windows]
        dsl_HHEs = dsl_HHE[:,good_windows]

        #------------------------------------------------------------------------------
        if VERBOSE_MODE:

	        figPSD_mean = plt.figure(figsize=(15,10))
	        ax1 = figPSD_mean.add_subplot(3, 2, 1)
	        ax1.semilogx(fz, log_psdz[:,good_windows], 'k', lw=0.5,alpha=0.5)
	        ax1.text(0.96, 0.8, 'HHZ', ha='center',bbox=dict(facecolor='w'),transform=ax1.transAxes)
	        ax1.set_xlim(1/200,max(fe))
	        ax1.set_ylim(-180,-20)

	        ax2 = figPSD_mean.add_subplot(3, 2, 2)
	        ax2.semilogx(fz, log_psdz[:,bad_windows], 'r', lw=0.5,alpha=0.5)
	        ax2.text(0.96, 0.8, 'HHZ', ha='center',bbox=dict(facecolor='w'),transform=ax2.transAxes)
	        ax2.set_xlim(1/200,max(fe))
	        ax2.set_ylim(-180,-20)

	        ax3 = figPSD_mean.add_subplot(3, 2, 3)
	        ax3.semilogx(fn, log_psdn[:,good_windows], 'k', lw=0.5,alpha=0.5)
	        ax3.text(0.96, 0.8, 'HHN', ha='center',bbox=dict(facecolor='w'),transform=ax3.transAxes)
	        ax3.set_xlim(1/200,max(fe))
	        ax3.set_ylim(-180,-20)

	        ax4 = figPSD_mean.add_subplot(3, 2, 4)
	        ax4.semilogx(fn, log_psdn[:,bad_windows], 'r', lw=0.5,alpha=0.5)
	        ax4.text(0.96, 0.8, 'HHN', ha='center',bbox=dict(facecolor='w'),transform=ax4.transAxes)
	        ax4.set_xlim(1/200,max(fe))
	        ax4.set_ylim(-180,-20)

	        ax5 = figPSD_mean.add_subplot(3, 2, 5)
	        ax5.semilogx(fe, log_psde[:,good_windows], 'k', lw=0.5,alpha=0.5)
	        ax5.set_xlim(1/200,max(fe))
	        ax5.text(0.96, 0.8, 'HHE', ha='center',bbox=dict(facecolor='w'),transform=ax5.transAxes)
	        ax5.set_xlabel('Frequency (Hz)',fontsize=15)
	        ax5.set_ylim(-180,-20)


	        ax6 = figPSD_mean.add_subplot(3, 2, 6)
	        ax6.semilogx(fe, log_psde[:,bad_windows], 'r', lw=0.5,alpha=0.5)
	        ax6.set_xlim(1/100,max(fe))
	        ax6.text(0.96, 0.8, 'HHE', ha='center',bbox=dict(facecolor='w'),transform=ax6.transAxes)
	        ax6.set_xlabel('Frequency (Hz)',fontsize=15)
	        ax6.set_ylim(-180,-20)

	        figPSD_mean.suptitle('Station = '+STATION+' - Day = '+UTCDateTime(year=year_spec,julday=julday_spec).strftime('%d/%m/%Y')+' - '+'good_wind='+str(len_good)+'/'+'bad_wind='+str(len_bad),fontsize=18)
	        daily_PSD_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_PSD_windows/'
	        os.makedirs(daily_PSD_output,exist_ok=True)
	        figPSD_mean.savefig(daily_PSD_output+'PSD_windows_good_bad_'+str(year_spec)+'_'+str(julday_spec)+'.png', dpi=300, facecolor='w', edgecolor='w')
	        plt.close()

	        #------------------------------------------------------------------------------

	        cmap = 'viridis_r'
	        heights = [1,1,1,0.1]
	        widths = [1]

	        gs_kw = dict(width_ratios=widths, height_ratios=heights)
	        fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=1, nrows=4,figsize=(15,10),sharex=True,gridspec_kw=gs_kw)

	        im = ax1.pcolormesh(tzh,freq_lst, dsl_HHZ, cmap=cmap,shading='auto')
	        ax1.text(0.96, 0.8, 'HHZ', ha='center',bbox=dict(facecolor='w'),transform=ax1.transAxes)
	        ax1.xaxis.set_major_locator(MultipleLocator(4))
	        ax1.yaxis.set_major_locator(MultipleLocator(0.05))

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
	        ax2.pcolormesh(tzh,freq_lst, dsl_HHN, cmap=cmap,shading='auto')
	        ax2.text(0.96, 0.8, 'HHN', ha='center',bbox=dict(facecolor='w'),transform=ax2.transAxes)
	        ax2.set_ylabel('Frequency (Hz)')
	        ax2.xaxis.set_major_locator(MultipleLocator(4))
	        ax2.yaxis.set_major_locator(MultipleLocator(0.05))
	        #-----------------------------------------------------
	        ax3.pcolormesh(tzh,freq_lst, dsl_HHE, cmap=cmap,shading='auto')
	        ax3.text(0.96, 0.8, 'HHE', ha='center',bbox=dict(facecolor='w'),transform=ax3.transAxes)
	        ax3.xaxis.set_major_locator(MultipleLocator(4))
	        ax3.yaxis.set_major_locator(MultipleLocator(0.05))
	        #-----------------------------------------------------
	        colors = ['yellowgreen' if eli == True else 'firebrick' for eli in good_windows]
	        ax4.bar(teh,height=1,width=WINDOW_LENGTH/3600,align='center',color=colors,edgecolor='k')

	        ax4.set_xlabel('Time (h)')
	        ax4.xaxis.set_major_locator(MultipleLocator(4))
	        ax4.set_ylim(0,1)
	        ax4.set_xlim(0,24)
	        ax4.set_yticks([])
	        #-----------------------------------------------------
	        daily_wind_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_data_windows/'
	        os.makedirs(daily_wind_output,exist_ok=True)
	        fig.suptitle('Station = '+STATION+' - Day = '+UTCDateTime(year=year_spec,julday=julday_spec).strftime('%d/%m/%Y'),fontsize=18)
	        fig.savefig(daily_wind_output+NETWORK+'_'+STATION+'_'+str(year_spec)+'_'+str(julday_spec)+'.png', dpi=300, facecolor='w', edgecolor='w')
	        plt.close()

        #----------------------------------------------------------------------------------------------------------------

        SELECTED_DAYS_dic = {
                            'PSDe': dsl_HHEs.tolist(),
                            'PSDn': dsl_HHNs.tolist(),
                            'PSDz': dsl_HHZs.tolist(),
                            'f': freq_lst.tolist(),
                            't': tzh.tolist()
                            }
        output_FOLDER_SPECTROGRAM = JSON_FILES+'SELECTED_SPECTROGRAM_WINDOWS_FILES/'+NETWORK+'.'+STATION+'/'
        os.makedirs(output_FOLDER_SPECTROGRAM,exist_ok=True)
        with open(output_FOLDER_SPECTROGRAM+'SELECTED_SPECTROGRAM_'+str(year_spec)+'_'+str(julday_spec)+'.json', 'w') as fp:
            json.dump(SELECTED_DAYS_dic, fp)

    else:
        first_filter_day = False

    return [good_windows.tolist(),first_filter_day]

#--------------------------------------------------------------------------------------------------------------------------------------------------------

def calculate_tilt_function(ASDF_FILES,GOOD_WINDOWS_LIST):

    data_HHZ = [ASDFDataSet(k) for k in ASDF_FILES if 'HHZ' in k][0]
    data_HHE = [ASDFDataSet(k) for k in ASDF_FILES if 'HHE' in k][0]
    data_HHN = [ASDFDataSet(k) for k in ASDF_FILES if 'HHN' in k][0]

    year_day = data_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.year
    julday_day = data_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.julday

    tr_HHZ = data_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0]
    tr_HHN = data_HHN.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0]
    tr_HHE = data_HHE.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0]

    slide_HHZ = np.array([k.data for k in tr_HHZ.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)],dtype=object)[GOOD_WINDOWS_LIST]
    slide_HHE = np.array([k.data for k in tr_HHN.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)],dtype=object)[GOOD_WINDOWS_LIST]
    slide_HHN = np.array([k.data for k in tr_HHE.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)],dtype=object)[GOOD_WINDOWS_LIST]

    #-------------------------------------------------------------------------------------------------------------
    n2 = prev_pow_2(int(WINDOW_LENGTH*NEW_SAMPLING_RATE))
    f_FFT = NEW_SAMPLING_RATE/2. * np.linspace(0., 1., int(n2/2) + 1)

    FFT_ZN = []
    FFT_ZE = []
    FFT_NE = []

    FFT_ZZ = []
    FFT_NN = []
    FFT_EE = []

    cFFT_ZN = []
    cFFT_ZE = []
    cFFT_NE = []

    qFFT_ZN = []
    qFFT_ZE = []
    qFFT_NE = []

    Z_lst = []
    N_lst = []
    E_lst = []
    for c, v in enumerate(slide_HHZ):
        ftZ = fft(slide_HHZ[c][0:len(f_FFT)], n=None, axis=-1, norm=None)
        ftN = fft(slide_HHE[c][0:len(f_FFT)], n=None, axis=-1, norm=None)
        ftE = fft(slide_HHN[c][0:len(f_FFT)], n=None, axis=-1, norm=None)

        Z = ftZ
        N = ftN
        E = ftE

        Z_lst.append(ftZ)
        N_lst.append(ftN)
        E_lst.append(ftE)

        FFT_ZN.append(np.conj(Z)*N)
        FFT_ZE.append(np.conj(Z)*E)
        FFT_NE.append(np.conj(N)*E)
        FFT_ZZ.append(np.conj(Z)*Z)
        FFT_NN.append(np.conj(N)*N)
        FFT_EE.append(np.conj(E)*E)

        cFFT_ZN.append(Z.real*N.real+Z.imag*N.imag)
        cFFT_ZE.append(Z.real*E.real+Z.imag*E.imag)
        cFFT_NE.append(N.real*E.real+N.imag*E.imag)

        qFFT_ZN.append(Z.real*N.imag+Z.imag*N.real)
        qFFT_ZE.append(Z.real*E.imag+Z.imag*E.real)
        qFFT_NE.append(N.real*E.imag+N.imag*E.real)

    G_ZN = smooth(np.abs(np.mean(FFT_ZN,axis=0)),40)
    G_ZE = smooth(np.abs(np.mean(FFT_ZE,axis=0)),40)
    G_NE = smooth(np.abs(np.mean(FFT_NE,axis=0)),40)

    G_ZZ = smooth(np.abs(np.mean(FFT_ZZ,axis=0)),40)
    G_NN = smooth(np.abs(np.mean(FFT_NN,axis=0)),40)
    G_EE = smooth(np.abs(np.mean(FFT_EE,axis=0)),40)

    C_ZN = smooth(np.abs(np.mean(cFFT_ZN,axis=0)),40)
    C_ZE = smooth(np.abs(np.mean(cFFT_ZE,axis=0)),40)
    C_NE = smooth(np.abs(np.mean(cFFT_NE,axis=0)),40)

    Q_ZN = smooth(np.abs(np.mean(qFFT_ZN,axis=0)),40)
    Q_ZE = smooth(np.abs(np.mean(qFFT_ZE,axis=0)),40)
    Q_NE = smooth(np.abs(np.mean(qFFT_NE,axis=0)),40)

    #---------------------------------------------

    wind_number = len(slide_HHZ)

    A_ZN,Error_A_ZN = Admittance_xy(G_ZN,G_ZZ,G_NN,wind_number)

    A_ZE,Error_A_ZE = Admittance_xy(G_ZE,G_ZZ,G_EE,wind_number)

    A_NE,Error_A_NE = Admittance_xy(G_NE, G_NN,G_EE,wind_number)

    #---------------------------------------------

    Gama2_ZN,Error_Gama2_ZN = Coherence_xy(G_ZN,G_ZZ,G_NN,wind_number)

    Gama2_ZE,Error_Gama2_ZE = Coherence_xy(G_ZE,G_ZZ,G_EE,wind_number)

    Gama2_NE,Error_Gama2_NE = Coherence_xy(G_NE,G_NN,G_EE,wind_number)

    #---------------------------------------------

    phase_ZN,Error_phase_ZN = Phase_xy(Q_ZN,C_ZN,G_ZN,G_ZZ,G_NN,wind_number)

    phase_ZE,Error_phase_ZE = Phase_xy(Q_ZE,C_ZE,G_ZE,G_ZZ,G_EE,wind_number)

    phase_NE,Error_phase_NE = Phase_xy(Q_NE,C_NE,G_NE,G_NN,G_EE,wind_number)

    #---------------------------------------------
    if VERBOSE_MODE:
        figA = plt.figure(figsize=(15,15))
        ax = figA.add_subplot(3, 3, 1)
        ax.semilogx(f_FFT,np.log(A_ZN),'ok',ms=3,alpha=0.5)
        ax.errorbar(f_FFT,np.log(A_ZN), yerr=Error_A_ZN, ecolor='k',fmt='none',elinewidth=0.5,capsize=2)
        ax.set_title('log(Admittance)',fontsize=15)
        ax.set_ylabel('HHZ-HHN',fontsize=14)
        ax.set_xlim(1/200,1/10)

        ax1 = figA.add_subplot(3, 3, 2)
        ax1.set_ylim(0,1)
        ax1.semilogx(f_FFT,Gama2_ZN, 'ok',ms=3,alpha=0.5)
        ax1.set_title('Coherence',fontsize=15)
        ax1.set_xlim(1/200,1/10)

        ax2 = figA.add_subplot(3, 3, 3)
        ax2.semilogx(f_FFT,phase_ZN, 'ok',ms=3,alpha=0.5)
        ax2.errorbar(f_FFT,phase_ZN, yerr=Error_phase_ZN, ecolor='k',fmt='none',elinewidth=0.5,capsize=2)
        ax2.set_ylim(-0.4,0.4)
        ax2.set_title('Phase',fontsize=15)
        ax2.set_xlim(1/200,1/10)

	    #----------------------------------------------------------------------------------------------------
        ax3 = figA.add_subplot(3, 3, 4)
        ax3.semilogx(f_FFT,np.log(A_ZE), 'ok',ms=3,alpha=0.5)
        ax3.errorbar(f_FFT,np.log(A_ZE), yerr=Error_A_ZE, ecolor='k',fmt='none',elinewidth=0.5,capsize=2)
        ax3.set_ylabel('HHZ-HHE',fontsize=14)
        ax3.set_xlim(1/200,1/10)

        ax4 = figA.add_subplot(3, 3, 5)
        ax4.set_ylim(0,1)
        ax4.semilogx(f_FFT,Gama2_ZE, 'ok',ms=3,alpha=0.5)
        ax4.set_xlim(1/200,1/10)

        ax5 = figA.add_subplot(3, 3, 6)
        ax5.semilogx(f_FFT,phase_ZE, 'ok',ms=3,alpha=0.5)
        ax5.errorbar(f_FFT,phase_ZE, yerr=Error_phase_ZE, ecolor='k',fmt='none',elinewidth=0.5,capsize=2)
        ax5.set_ylim(-0.4,0.4)
        ax5.set_xlim(1/200,1/10)

	    #-----------------------------------------------------------------------------------------------------
        ax6 = figA.add_subplot(3, 3, 7)
        ax6.semilogx(f_FFT,np.log(A_NE), 'ok',ms=3,alpha=0.5)
        ax6.errorbar(f_FFT,np.log(A_NE), yerr=Error_A_NE, ecolor='k',fmt='none',elinewidth=0.5,capsize=2)
        ax6.set_xlabel('Frequency (Hz)',fontsize=14)
        ax6.set_ylabel('HHN-HHE',fontsize=14)
        ax6.set_xlim(1/200,1/10)

        ax7 = figA.add_subplot(3, 3, 8)
        ax7.set_ylim(0,1)
        ax7.semilogx(f_FFT,Gama2_NE, 'ok',ms=3,alpha=0.5)
        ax7.set_xlabel('Frequency (Hz)',fontsize=14)
        ax7.set_xlim(1/200,1/10)

        ax8 = figA.add_subplot(3, 3, 9)
        ax8.semilogx(f_FFT,phase_NE, 'ok',ms=3,alpha=0.5)
        ax8.errorbar(f_FFT,phase_NE, yerr=Error_phase_NE, ecolor='k',fmt='none',elinewidth=0.5,capsize=2)
        ax8.set_ylim(-0.4,0.4)
        ax8.set_xlabel('Frequency (Hz)',fontsize=14)
        ax8.set_xlim(1/200,1/10)

	    #-----------------------------------------------------------------------------------------------------

        figA.suptitle('Station = '+STATION+' - Day = '+UTCDateTime(year=int(year_day),julday=int(julday_day)).strftime('%d/%m/%Y'),fontsize=18)
        daily_A_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_Admittance_Coherence_Phase/'
        os.makedirs(daily_A_output,exist_ok=True)
        figA.savefig(daily_A_output+'Admittance_Coherence_Phase_'+str(year_day)+'_'+str(julday_day)+'.png', dpi=300, facecolor='w', edgecolor='w')
        plt.close()

    #-----------------------------------------------------------------------------------------------------

    #Calculating TILT

    direc = np.arange(0., 360., 10.)
    coh = np.zeros(len(direc))
    ph = np.zeros(len(direc))
    cZZ = np.abs(np.mean(np.array(Z_lst)*np.conj(np.array(Z_lst)), axis=0))

    for i, d in enumerate(direc):

        # Rotate horizontals
        ftH = rotate_dir(np.array(N_lst), np.array(E_lst),d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH*np.conj(ftH), axis=0))
        cHZ = np.mean(np.conj(ftH)*np.array(Z_lst), axis=0)

        Co = (np.abs(cHZ)**2)/(cZZ*cHH)
        Ph = 180/np.pi*np.arctan2(cHZ.imag,cHZ.real)

        # Calculate coherence over frequency band
        coh[i] = np.mean(Co[(f_FFT > tiltfreq[0]) & (f_FFT < tiltfreq[1])])
        ph[i] = np.pi/2. - np.mean(Ph[(f_FFT > tiltfreq[0]) & (f_FFT < tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(coh == coh.max())

    # Phase and direction at maximum coherence
    phase_value = ph[ind[0]][0]
    coh_value = coh[ind[0]][0]
    tilt = direc[ind[0]][0]

    # Refine search
    rdirec = np.arange(direc[ind[0]][0]-10., direc[ind[0]][0]+10., 1.)
    rcoh = np.zeros(len(direc))
    rph = np.zeros(len(direc))

    for i, d in enumerate(rdirec):

        # Rotate horizontals
        ftH = rotate_dir(np.array(N_lst),np.array(E_lst), d)

        # Get transfer functions
        cHH = np.abs(np.mean(np.conj(ftH)*ftH, axis=0))
        cHZ = np.mean(np.conj(ftH)*np.array(Z_lst), axis=0)

        Co = (np.abs(cHZ)**2)/(cZZ*cHH)
        Ph = 180/np.pi*np.arctan2(cHZ.imag,cHZ.real)

        # Calculate coherence over frequency band
        rcoh[i] = np.mean(Co[(f_FFT > tiltfreq[0]) & (f_FFT < tiltfreq[1])])
        rph[i] = np.pi/2. - np.mean(Ph[(f_FFT > tiltfreq[0]) & (f_FFT < tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(rcoh == rcoh.max())

    # Phase and direction at maximum coherence
    phase_value = rph[ind[0]][0]
    coh_value = rcoh[ind[0]][0]
    tilt = rdirec[ind[0]][0]

    # Phase has to be close to zero - otherwise add pi
    if phase_value > 0.5*np.pi:
        tilt += 180.
    if tilt > 360.:
        tilt -= 360.

    # Now calculate spectra at tilt direction
    ftH1 = rotate_dir(np.array(N_lst), np.array(E_lst), tilt)

    # Get transfer functions
    cHH1 = np.abs(np.mean(np.conj(ftH1)*ftH1, axis=0))

    cHZ1 = np.mean(np.conj(ftH1)*np.array(Z_lst), axis=0)

    tf_ZH1 = np.conj(cHZ1)/cHH1

    #-----------------------------------------------------------------------------------------------------
    if VERBOSE_MODE:
	    colors = plt.cm.cividis(np.linspace(0, 1, coh.shape[0]))

	    figureTILT, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,5))
	    for i, (co, p) in enumerate(zip(coh, ph)):
	        ax1.plot(direc[i], co, 'ok')
	        ax2.plot(direc[i], p, 'ok')
	    ax1.set_ylabel('Coherence')
	    ax1.set_xlabel('Angle from HHE')
	    ax1.set_ylim((0, 1.))
	    ax2.set_ylabel('Phase')
	    ax2.set_xlabel('Angle from HHE')
	    ax1.set_title('Maximum coherence = '+str(coh_value))
	    ax2.set_title('Tilt = '+str(tilt))
	    figureTILT.suptitle('Station = '+STATION+' - Day = '+UTCDateTime(year=int(year_day),julday=int(julday_day)).strftime('%d/%m/%Y'),fontsize=18)
	    figureTILT.savefig(daily_A_output+'Tilt_Coherence_'+str(year_day)+'_'+str(julday_day)+'.png', dpi=300, facecolor='w', edgecolor='w')
	    plt.close()

    #-------------------------------------------------------------------------------

    return [tilt,[int(year_day),int(julday_day)],coh_value,phase_value,tf_ZH1]

    #-------------------------------------------------------------------------------

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
    ss = int(WINDOW_LENGTH/st[0].stats.delta)

    # hanning window
    wind = np.ones(ws)

    # Calculating the spectrogram
    spec = spectrogram(x=st[0].data,fs=st[0].stats.sampling_rate, window=wind, nperseg=ws, noverlap=0)
    f, t, psd = spec

    #-----------------------------------------------------------------
    #Creating ASDF preprocessed files folder
    output_DATA_DAY = ASDF_FILES+'DATA_DAY_FILES/'+NETWORK+'.'+STATION+'/'
    os.makedirs(output_DATA_DAY,exist_ok=True)

    ds = ASDFDataSet(output_DATA_DAY+'DATA_DAY_'+NETWORK+'_'+STATION+'_'+channel+'_'+year_day+'_'+julday_day+".h5", compression="gzip-3")
    ds.add_waveforms(st, tag="preprocessed_recording")
    ds.add_stationxml(STATIONXML_DIR+'.'.join([NETWORK,STATION,'xml']))

    #-----------------------------------------------------------------
    # The type always should be camel case.
    data_type_f = "Frequencies"
    # Name to identify the particular piece of data.
    path_f = sta_channel_id+'.f'
    # Any additional parameters as a Python dictionary which will end up as attributes of the array.
    parameters_f = {'f':'Array of sample frequencies.',
        		  'sampling_rate_in_hz': st[0].stats.sampling_rate,
        		  'station_id': sta_channel_id}

    ds.add_auxiliary_data(data=f, data_type=data_type_f, path=path_f, parameters=parameters_f)

    # The type always should be camel case.
    data_type_t = "Times"
    # Name to identify the particular piece of data.
    path_t = sta_channel_id+'.t'
    # Any additional parameters as a Python dictionary which will end up as attributes of the array.
    parameters_t = {'t':'Array of segment times.',
        		  'sampling_rate_in_hz': st[0].stats.sampling_rate,
        		  'station_id': sta_channel_id}


    ds.add_auxiliary_data(data=t, data_type=data_type_t, path=path_t, parameters=parameters_t)

    # The type always should be camel case.
    data_type_s = "Spectrogram"
    # Name to identify the particular piece of data.
    path_s = sta_channel_id+'.S'
    # Any additional parameters as a Python dictionary which will end up as attributes of the array.
    parameters_s = {'S':'Spectrogram of x. By default, the last axis of S corresponds to the segment times.',
        		  'sampling_rate_in_hz': st[0].stats.sampling_rate,
        		  'station_id': sta_channel_id}

    ds.add_auxiliary_data(data=psd, data_type=data_type_s, path=path_s, parameters=parameters_s)

    #-------------------------------------------------------------------------------

# ============
# Main program
# ============

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

# initializing list of stations by scanning name of miniseed files
files = filelist(basedir=MSEED_DIR+NETWORK+'/'+STATION+'/',interval_period_date=INTERVAL_PERIOD_DATE)

print('Total of miniseed files = '+str(len(files)))
print('\n')
'''
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
'''
print('=======================')
print('Filtering daily windows')
print('=======================')
print('\n')

daily_lst_data = [[]]*len(INTERVAL_PERIOD_DATE)

DATA_DAY_FILES = sorted(glob.glob(ASDF_FILES+'DATA_DAY_FILES/'+NETWORK+'.'+STATION+'/*.h5'))
DATA_DAYs = [i.split('.h5')[0].split('_')[-2]+'.'+i.split('.h5')[0].split('_')[-1] for i in DATA_DAY_FILES]

for l,k in enumerate(INTERVAL_PERIOD_DATE):
	daily_lst_data[l] = [j for i,j in enumerate(DATA_DAY_FILES) if DATA_DAYs[i] == k]
'''
#----------------------------------------------------------------------------------------------------------

first_filter_day = []
good_windows_lst = []

start_time = time.time()

with Pool(processes=int(num_processes/2)) as p:
    max_ = len(daily_lst_data)
    with tqdm(total=max_) as pbar:
        for i,j in enumerate(p.imap_unordered(daily_process_function, daily_lst_data)):
            good_windows_lst.append(j[0])
            first_filter_day.append(j[1])
            pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

#-------------------------------------------------------------------------------

FILTER_DAYS_dic = {
					'daily_lst_data':daily_lst_data,
					'first_filter_day': first_filter_day,
					'good_windows_lst': good_windows_lst
					}

output_FOLDER_WINDOWS_FILTER_FILES = JSON_FILES+'SELECTED_SPECTROGRAM_WINDOWS_FILTER_FILES/'+NETWORK+'.'+STATION+'/'
os.makedirs(output_FOLDER_WINDOWS_FILTER_FILES,exist_ok=True)
with open(output_FOLDER_WINDOWS_FILTER_FILES+'FIRST_FILTER_SPECTROGRAM_WINDOWS.json', 'w') as fp:
	json.dump(FILTER_DAYS_dic, fp)
#-------------------------------------------------------------------------------
'''
print('\n')
print('==============')
print('Filtering days')
print('==============')
print('\n')

output_FOLDER_WINDOWS_FILTER_FILES = JSON_FILES+'SELECTED_SPECTROGRAM_WINDOWS_FILTER_FILES/'+NETWORK+'.'+STATION+'/'
FIRST_FILTER_WINDOWS = json.load(open(output_FOLDER_WINDOWS_FILTER_FILES+'FIRST_FILTER_SPECTROGRAM_WINDOWS.json'))

second_daily_lst_data = list(compress(FIRST_FILTER_WINDOWS['daily_lst_data'], FIRST_FILTER_WINDOWS['first_filter_day']))
second_good_windows_lst = list(compress(FIRST_FILTER_WINDOWS['good_windows_lst'], FIRST_FILTER_WINDOWS['first_filter_day']))

output_FOLDER_SPECTROGRAM = JSON_FILES+'SELECTED_SPECTROGRAM_WINDOWS_FILES/'+NETWORK+'.'+STATION+'/'

json_daily_files = sorted(glob.glob(output_FOLDER_SPECTROGRAM+'*'))

PSDe = []
PSDn = []
PSDz = []
f = []
PSD_data_lst = []

for j in tqdm(json_daily_files):
    sta_dic = json.load(open(j))
    STA_NAME = j.split('/')[-2]
    PSD_julday = j.split('/')[-1].split('.')[0].split('_')[-2:][1]
    PSD_year = j.split('/')[-1].split('.')[0].split('_')[-2:][0]

    t = np.array(sta_dic['t'])
    f = np.array(sta_dic['f'])

    PSDe.append(np.mean(sta_dic['PSDe'],axis=1))
    PSDn.append(np.mean(sta_dic['PSDn'],axis=1))
    PSDz.append(np.mean(sta_dic['PSDz'],axis=1))

PSDe = np.array(PSDe).T
PSDn = np.array(PSDn).T
PSDz = np.array(PSDz).T


#-------------------------------------------------------------------------------

good_days = np.repeat([True], len(json_daily_files))
indwin = np.argwhere(good_days == True)
moveon = False

while moveon == False:

    normvar_HHZ = np.zeros(np.sum(good_days))
    for ii, tmp in enumerate(indwin):
        ind = np.copy(indwin)
        ind = np.delete(ind, ii)
        normvar_HHZ[ii] = norm(np.std(PSDz[:, ind], axis=1), ord=2)
    ubernorm_HHZ = np.median(normvar_HHZ) - normvar_HHZ

    normvar_HHN = np.zeros(np.sum(good_days))
    for ii, tmp in enumerate(indwin):
        ind = np.copy(indwin)
        ind = np.delete(ind, ii)
        normvar_HHN[ii] = norm(np.std(PSDn[:, ind], axis=1), ord=2)
    ubernorm_HHN = np.median(normvar_HHN) - normvar_HHN

    normvar_HHE = np.zeros(np.sum(good_days))
    for ii, tmp in enumerate(indwin):
        ind = np.copy(indwin)
        ind = np.delete(ind, ii)
        normvar_HHE[ii] = norm(np.std(PSDe[:, ind], axis=1), ord=2)
    ubernorm_HHE = np.median(normvar_HHE) - normvar_HHE

    norm_allchannels = ubernorm_HHZ+ubernorm_HHN+ubernorm_HHE

    eliminated_days = norm_allchannels > TOL*np.std(norm_allchannels)

    if np.sum(eliminated_days) == 0:
        moveon = True

    trypenalty = eliminated_days[np.argwhere(eliminated_days == False)].T[0]

    if ftest(eliminated_days, 1, trypenalty, 1) < ALPHA:
        good_days[indwin[eliminated_days == True]] = False
        indwin = np.argwhere(good_days == True)
        moveon = False
    else:
        moveon = True

bad_days = np.array([False if i == True else True for i in good_days])

len_good_days = int(sum(good_days))
len_bad_days = int(sum(bad_days))

#-------------------------------------------------------------------------------
if VERBOSE_MODE:
	figPSD_days = plt.figure(figsize=(15,10))
	ax1 = figPSD_days.add_subplot(3, 2, 1)
	if len_good_days > 0:
	    ax1.semilogx(f, PSDz[:,good_days], 'k', lw=0.5,alpha=0.8)
	ax1.text(0.96, 0.8, 'HHZ', ha='center',bbox=dict(facecolor='w'),transform=ax1.transAxes)
	ax1.set_xlim(1/200,max(f))

	ax2 = figPSD_days.add_subplot(3, 2, 2)
	if len_bad_days > 0:
	    ax2.semilogx(f, PSDz[:,bad_days], 'r', lw=0.5,alpha=0.5)
	ax2.text(0.96, 0.8, 'HHZ', ha='center',bbox=dict(facecolor='w'),transform=ax2.transAxes)
	ax2.set_xlim(1/200,max(f))

	ax3 = figPSD_days.add_subplot(3, 2, 3)
	if len_good_days > 0:
	    ax3.semilogx(f, PSDn[:,good_days], 'k', lw=0.5,alpha=0.5)
	ax3.text(0.96, 0.8, 'HHN', ha='center',bbox=dict(facecolor='w'),transform=ax3.transAxes)
	ax3.set_xlim(1/200,max(f))

	ax4 = figPSD_days.add_subplot(3, 2, 4)
	if len_bad_days > 0:
	    ax4.semilogx(f, PSDn[:,bad_days], 'r', lw=0.5,alpha=0.5)
	ax4.text(0.96, 0.8, 'HHN', ha='center',bbox=dict(facecolor='w'),transform=ax4.transAxes)
	ax4.set_xlim(1/200,max(f))

	ax5 = figPSD_days.add_subplot(3, 2, 5)
	if len_good_days > 0:
	    ax5.semilogx(f, PSDe[:,good_days], 'k', lw=0.5,alpha=0.5)
	ax5.set_xlim(1/200,max(f))
	ax5.text(0.96, 0.8, 'HHE', ha='center',bbox=dict(facecolor='w'),transform=ax5.transAxes)
	ax5.set_xlabel('Frequency (Hz)',fontsize=15)


	ax6 = figPSD_days.add_subplot(3, 2, 6)
	if len_bad_days > 0:
	    ax6.semilogx(f, PSDe[:,bad_days], 'r', lw=0.5,alpha=0.5)
	ax6.set_xlim(1/100,max(f))
	ax6.text(0.96, 0.8, 'HHE', ha='center',bbox=dict(facecolor='w'),transform=ax6.transAxes)
	ax6.set_xlabel('Frequency (Hz)',fontsize=15)

	figPSD_days.suptitle('Station = '+STA_NAME+' - good_days='+str(len_good_days)+'/'+'bad_days='+str(len_bad_days),fontsize=18)
	daily_PSD_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_PSD_windows/'
	os.makedirs(daily_PSD_output,exist_ok=True)
	figPSD_days.savefig(daily_PSD_output+'PSD_days_good_bad'+'.png', dpi=300, facecolor='w', edgecolor='w')
	plt.close()

#-------------------------------------------------------------------------------
print('\n')
print('====================')
print('Calculating the tilt')
print('====================')
print('\n')

last_daily_lst_data = list(compress(second_daily_lst_data, good_days.tolist()))
last_good_windows_lst = list(compress(second_good_windows_lst, good_days.tolist()))

tilt_lst = []
day_time_tilt = []
coh_lst = []
phase_lst = []
tf_ZH1_lst = []

#-------------------------
start_time = time.time()
with Pool(processes=num_processes) as p:
    max_ = len(last_daily_lst_data)
    with tqdm(total=max_) as pbar:
        for i,j in enumerate(p.starmap(calculate_tilt_function,  zip(last_daily_lst_data,last_good_windows_lst))):
            tilt_lst.append(j[0])
            day_time_tilt.append(j[1])
            coh_lst.append(j[2])
            phase_lst.append(j[3])
            tf_ZH1_lst.append(j[4])
            pbar.update()
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

#----------------------------------------------------------------------------------------------------------------------------------------------------
if VERBOSE_MODE:
    daily_data = [mdates.num2date(UTCDateTime(year=i[0],julday=i[1]).matplotlib_date) for i in day_time_tilt]

    figureTILTfinal, (ax1, ax2) = plt.subplots(2, 1,sharex=True,figsize=(10,5))
    ax1.plot(daily_data, coh_lst, 'ok')
    ax1.set_ylabel('Coherence')
    ax1.set_ylim(0, 1)
    ax1.set_xlim(mdates.num2date(INTERVAL_PERIOD[0].matplotlib_date),mdates.num2date(INTERVAL_PERIOD[1].matplotlib_date))

    ax2.plot(daily_data, tilt_lst, 'ok')
    ax2.set_ylabel('Tilt')
    ax2.set_xlabel('Day')
    ax2.set_ylim(0, 360)
    ax2.yaxis.set_major_locator(MultipleLocator(90))
    ax2.yaxis.set_minor_locator(MultipleLocator(10))
    ax2.set_xlim(mdates.num2date(INTERVAL_PERIOD[0].matplotlib_date),mdates.num2date(INTERVAL_PERIOD[-1].matplotlib_date))
    plt.gcf().autofmt_xdate()

    figureTILTfinal.suptitle('Station = '+STATION+' - Period = '+UTCDateTime(year=int(INTERVAL_PERIOD[0].year),julday=int(INTERVAL_PERIOD[0].julday)).strftime('%d/%m/%Y')+'--'+UTCDateTime(year=int(INTERVAL_PERIOD[-1].year),julday=int(INTERVAL_PERIOD[-1].julday)).strftime('%d/%m/%Y'),fontsize=18)
    daily_A_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_Admittance_Coherence_Phase/'
    os.makedirs(daily_A_output,exist_ok=True)
    figureTILTfinal.savefig(daily_A_output+'Tilt.Coherence.total.png', dpi=300, facecolor='w', edgecolor='w')
    plt.close()

#-------------------------------------------------------------------------------

print('\n')
print('===================')
print('Correcting the data')
print('===================')
print('\n')

# Method to apply transfer functions between multiple components to produce corrected/cleaned vertical components.

for i,j in enumerate(tqdm(daily_lst_data)):

    data_HHZ = [ASDFDataSet(k) for k in j if 'HHZ' in k][0]
    data_HHE = [ASDFDataSet(k) for k in j if 'HHE' in k][0]
    data_HHN = [ASDFDataSet(k) for k in j if 'HHN' in k][0]

    year_day = data_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.year
    julday_day = data_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0].stats.starttime.julday

    tr_HHZ = data_HHZ.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0]
    tr_HHN = data_HHN.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0]
    tr_HHE = data_HHE.waveforms[NETWORK+'.'+STATION]['preprocessed_recording'][0]

    slide_HHZ = np.array([k.data for k in tr_HHZ.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])
    slide_HHE = np.array([k.data for k in tr_HHN.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])
    slide_HHN = np.array([k.data for k in tr_HHE.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])

    #-------------------------------------------------------------------------------------------------------------
    n2 = prev_pow_2(int(WINDOW_LENGTH*NEW_SAMPLING_RATE))
    f_FFT = NEW_SAMPLING_RATE/2. * np.linspace(0., 1., int(n2/2) + 1)

    FFT_ZN = []
    FFT_ZE = []
    FFT_NE = []

    FFT_ZZ = []
    FFT_NN = []
    FFT_EE = []

    cFFT_ZN = []
    cFFT_ZE = []
    cFFT_NE = []

    qFFT_ZN = []
    qFFT_ZE = []
    qFFT_NE = []

    Z_lst = []
    N_lst = []
    E_lst = []
    for c, v in enumerate(slide_HHZ):
        ftZ = fft(slide_HHZ[c][0:len(f_FFT)], n=None, axis=-1, norm=None)
        ftN = fft(slide_HHE[c][0:len(f_FFT)], n=None, axis=-1, norm=None)
        ftE = fft(slide_HHN[c][0:len(f_FFT)], n=None, axis=-1, norm=None)

        Z = ftZ
        N = ftN
        E = ftE

        Z_lst.append(ftZ)
        N_lst.append(ftN)
        E_lst.append(ftE)

    tilt_mean = sum(tilt_lst)/len(tilt_lst)

	# Now calculate spectra at tilt direction
    ftH1 = rotate_dir(np.array(N_lst), np.array(E_lst), tilt_mean)

    # Get transfer functions and correcting data
    fTF_ZH1 = smooth(np.array(sum(tf_ZH1_lst)/len(tf_ZH1_lst)),50)
    corrspec1 = np.array([np.array(Z_lst[i]) - fTF_ZH1*ftH1[i] for i,j in enumerate(ftH1)])
    corrtime1 = np.real(np.fft.ifft(corrspec1))

    #-------------------------------------------------------------------------------------------------------------

    #Saving data
    slide_HHZ_CORRECTED_H1 = Stream([k for k in tr_HHZ.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])

    for i,j in enumerate(slide_HHZ_CORRECTED_H1):
        j.data = corrtime1[i]

    slide_HHZ_CORRECTED_H1.merge(method=0, fill_value=0)
    inv = read_inventory(STATIONXML_DIR+'.'.join([NETWORK,STATION,'xml']))
    slide_HHZ_CORRECTED_H1.attach_response(inv)

    CORRECT_DATA_OUTPUT = CORRECT_DATA_TRANSFER_FUNC_OUTPUT+'/TF_ZH1/'+NETWORK+'/'+STATION+'/HHZ.D/'
    os.makedirs(CORRECT_DATA_OUTPUT,exist_ok=True)

    CORRECT_DATA_OUTPUT_STR = NETWORK+'.'+STATION+'..HHZ.D.'+str(year_day)+'.'+str(julday_day)
    slide_HHZ_CORRECTED_H1.write(CORRECT_DATA_OUTPUT+CORRECT_DATA_OUTPUT_STR,format='MSEED')
    #----------------------------------------------------------------------------------------------

    #Importing noise models:
    NOISE_MODEL_DATA = np.load(NOISE_MODEL_FILE)
    nhnm = NOISE_MODEL_DATA['high_noise']
    nlnm = NOISE_MODEL_DATA['low_noise']
    periods = NOISE_MODEL_DATA['model_periods']
    xdata = 1/periods

    #------------------------------------------------------------------
    # Points in window
    ws = int(len(f_FFT))

    # hanning window
    wind = np.ones(ws)

    figPSD_mean, (ax0,ax1) = plt.subplots(1, 2,figsize=(15,5),sharey=True,sharex=False)

    for i in slide_HHZ:
        spec0 = spectrogram(x=i,fs=NEW_SAMPLING_RATE, window=wind, nperseg=ws, noverlap=None)
        f0, t0, psd0 = spec0

        # avoid calculating log of zero
        idx = psd0 < DTINY
        psd0[idx] = DTINY

        idf = f0 < DTINY
        f0[idf] = DTINY

        # go to dB
        log_psdz0 = np.log10(psd0)
        log_psdz0 *= 10

        periods_0 =  1.0 / f0
        #log_psdz0 = smooth(log_psdz0,8)
        ax0.semilogx(periods_0,log_psdz0, 'k', lw=0.5,alpha=0.5)

    ax0.plot(periods,nhnm, '0.4', linewidth=2, zorder=10)
    ax0.plot(periods,nlnm, '0.4', linewidth=2, zorder=10)
    ax0.text(0.9, 0.9, 'Raw', ha='center',bbox=dict(facecolor='w'),transform=ax0.transAxes)
    ax0.set_xlim(0.01, 50)
    ax0.set_ylim(-180,-30)
    ax0.set_ylabel('Amplitude [dB]')
    ax0.set_xlabel('Period [s]')
    ax0.set_title('Number of Segments: '+str(len(psd0)))


    for i in corrtime1:
        spec1 = spectrogram(x=i,fs=NEW_SAMPLING_RATE, window=wind, nperseg=ws, noverlap=None)
        f1, t1, psd1 = spec1

        # avoid calculating log of zero
        idx = psd1 < DTINY
        psd1[idx] = DTINY

        idf = f1 < DTINY
        f1[idf] = DTINY

        # go to dB
        log_psdz1 = np.log10(psd1)
        log_psdz1 *= 10

        periods_1 =  1.0 / f1
        #log_psdz1 = smooth(log_psdz1,8)
        ax1.semilogx(periods_1,log_psdz1, 'k', lw=0.5,alpha=0.5)

    ax1.plot(periods,nhnm, '0.4', linewidth=2, zorder=10)
    ax1.plot(periods,nlnm, '0.4', linewidth=2, zorder=10)
    ax1.text(0.9, 0.9, 'Corrected', ha='center',bbox=dict(facecolor='w'),transform=ax1.transAxes)
    ax1.set_xlim(0.01, 50)
    ax1.set_ylim(-180,-30)
    ax1.set_xlabel('Period [s]')
    ax1.set_title('Number of Segments: '+str(len(psd1)))



    figPSD_mean.suptitle('Station = '+STATION+' - Day = '+UTCDateTime(year=int(year_day),julday=int(julday_day)).strftime('%d/%m/%Y'),fontsize=18)
    daily_TRANSFER_FUNC_CORRECTION_output = TRANSFER_FUNC_OUTPUT+NETWORK+'.'+STATION+'/Daily_TRANSFER_FUNC_CORRECTION/'
    os.makedirs(daily_TRANSFER_FUNC_CORRECTION_output,exist_ok=True)
    figPSD_mean.savefig(daily_TRANSFER_FUNC_CORRECTION_output+'DAILY_TRANSFER_FUNC_CORRECTION_'+str(year_day)+'_'+str(julday_day)+'.png', dpi=300, facecolor='w', edgecolor='w')
    plt.close()
