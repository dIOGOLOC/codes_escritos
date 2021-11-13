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
import matplotlib.dates as mdates
from matplotlib.patches import Ellipse
import matplotlib.cbook as cbook
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream, Trace
from obspy.io.xseed import Parser
from obspy.signal.filter import bandpass,lowpass
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.util import prev_pow_2
from obspy.signal.cross_correlation import correlate as obscorr
from obspy.signal.cross_correlation import xcorr_max
from obspy.core.util import AttribDict

import glob
import os
import numpy as np
from numpy.fft import rfft, irfft, fft, ifft, fftfreq
from itertools import combinations,product,compress
from numpy.lib.stride_tricks import as_strided
import pandas as pd
import pyarrow.feather as feather
from scipy.signal import spectrogram, detrend, resample,savgol_filter,decimate
from scipy.linalg import norm
import pickle
import random
import collections
from copy import copy
import datetime

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

from pyasdf import ASDFDataSet

# ==================
# Configuration file
# ==================

# Folders input

MSEED_DIR_OBS = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/obs_data_MSEED/'

MSEED_DIR_STA = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/data/'

# -------------------------------

# Shapefile  boundary states input

BOUNDARY_STATES_SHP = '/media/diogoloc/Backup/dados_posdoc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/UFEBRASIL.shp'

# -------------------------------

# Stations and OBSs information

OBS_LST = ['OBS17','OBS18','OBS20','OBS22']

STATIONS_LST = ['ABR01','DUB01','MAN01','OBS20','OBS22','TER01','ALF01','GDU01','NAN01','TIJ01','CAJ01','GUA01','OBS17','PET01','TRI01','CAM01','JAC01','OBS18','RIB01','VAS01','CMC01','MAJ01','SLP01','PARB','CNLB','BSFB']
STATIONS_LST = sorted(STATIONS_LST)

STATIONXML_DIR = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

CHANNEL_LST = ['HHZ.D','HHN.D','HHE.D','HH1.D','HH2.D']

# -------------------------------

# Folders output

CLOCK_DRIFT_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/ASDF_FILES/'

FEATHER_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/FEATHER_FILES/'

# -------------------------------
#create figures?
VERBOSE = False

# Input parameters

FIRSTDAY = '2019-08-01'
LASTDAY = '2020-06-01'

#Each hour-long seismogram is amplitude clipped at twice its standard deviation of that hour-long time window.
CLIP_FACTOR = 2

MIN_WINDOWS = 30

WINDOW_LENGTH = 3600

#max time window (s) for cross-correlation
SHIFT_LEN = 1800

PERIOD_BANDS = [[2, 5], [5, 10], [10, 30]]
# (these bands focus on periods ~3, 7, 20 seconds)

# default parameters to define the signal and noise windows used to
# estimate the SNR:
# - the signal window is defined according to a min and a max velocity as:
#   dist/vmax < t < dist/vmin
# - the noise window has a fixed size and starts after a fixed trailing
#   time from the end of the signal window
SIGNAL_WINDOW_VMIN = 2.0
SIGNAL_WINDOW_VMAX = 4.0
SIGNAL2NOISE_TRAIL = 700.0
NOISE_WINDOW_SIZE = 700.0

#Returns pairs and spectral SNR array whose spectral SNRs are all >= minspectSNR
minspectSNR = 1

#RESAMPLING
NEW_SAMPLING_RATE = 2

# -------------------------------

# Constants and parameters

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# -------------------------------

# MULTIPROCESSING

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
    files_lst = []
    for root, dirs, files in os.walk(basedir):
        for file in files:
            files_path  = os.path.join(root, file)
            if any(day_s in files_path for day_s in interval_period_date):
                files_lst.append(files_path)

    file_lsts = sorted(files_lst)

    return file_lsts

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

# Calculating signal-to-noise ratio
def obscorr_window(data1,time_data1,data2,time_data2,dist,vmin,vmax):
    """
    Calculate the CrossCorrelation according to the distance between two stations.

    The signal window is defined by *vmin* and *vmax*:
        dist/*vmax* < t < dist/*vmin*

    @type data1: numpy array
    @type time_data1: numpy array
    @type data2: numpy array
    @type time_data2: numpy array
    @type dist: float
    @type vmin: float
    @type vmax: float
    """

    # signal window
    tmin_signal = dist/vmax
    tmax_signal = dist/vmin

    signal_window1 = (time_data1 >= tmin_signal) & (time_data1 <= tmax_signal)
    signal_window2 = (time_data2 >= tmin_signal) & (time_data2 <= tmax_signal)

    trace1 = data1[signal_window1]
    trace2 = data2[signal_window2]

    cc = obscorr(trace1,trace2,1000)
    shift, coefficient = xcorr_max(cc)

    return shift, coefficient

#-------------------------------------------------------------------------------

# Calculating Cross-Correlation between two traces
def SNR(data,dist,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX):
	"""
    The signal window is defined by *vmin* and *vmax*:
    	dist/*vmax* < t < dist/*vmin*

    @type data: numpy array
    @type dist: float
    @type vmin: float
    @type vmax: float
    """

    # signal window
	tmin_signal = dist/vmax
	tmax_signal = dist/vmin

    # noise window
	tmin_noise = tmax_signal + signal2noise_trail
	tmax_noise = tmin_noise + noise_window_size

	signal_window = (time_data >= tmin_signal) & (time_data <= tmax_signal)
	noise_window = (time_data >= tmin_noise) & (time_data <= tmax_noise)

	peak = np.abs(data[signal_window]).max()
	noise = data[noise_window].std()

    # appending SNR
	SNR = peak / noise

    # returning SNR
	return SNR

#-------------------------------------------------------------------------------

def Normalize(data):
	"""
	z(i)=2*(x(i)−min(x)/max(x)-min(x))−1

	where x=(x1,...,xn) and z(i) is now your ith normalized data between -1 and 1.

	@type data: list
	"""

	normalized_data = [2*(i-data.min()/(data.max()-data.min()))-1 for i in data]

	return normalized_data

#-------------------------------------------------------------------------------

def get_stations_data(f,Amp_clip=True,onebit_norm=True,white_spectral=True):
    """
    Gets stations daily data from miniseed file

    @type f: path of the minissed file (str)
    """

    # splitting subdir/basename
    subdir, filename = os.path.split(f)

    # network, station name and station channel in basename,
    # e.g., ON.TIJ01..HHZ.D.2020.002

    network, name = filename.split('.')[0:2]
    sta_channel_id = filename.split('.D.')[0]
    sta_channel = sta_channel_id.split('..')[1]
    time_day = filename.split('.D.')[-1]
    year_day = time_day.split('.')[0]
    julday_day = time_day.split('.')[1]

    if sta_channel == 'HHZ':
        sta_channel = 'HHZ'
    elif sta_channel == 'HHN' or sta_channel == 'HH1':
        sta_channel = 'HHN'
    elif sta_channel == 'HHE' or sta_channel == 'HH2':
        sta_channel = 'HHE'

    output_DATA_DAY = ASDF_FILES+'DATA_DAY_FILES/'+year_day+'.'+julday_day+'/'

    if os.path.isfile(output_DATA_DAY+'DATA_DAY_'+network+'_'+name+'_'+sta_channel+'_'+year_day+'_'+julday_day+'.h5'):
        pass

    else:

        try:
            st = read(f)
            st_starttime = st[0].stats.starttime
            st_endtime = st[0].stats.endtime

            if len(st[0].data) > WINDOW_LENGTH*100:
        	    st_traces = [k for k in st.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH/2)]

        	    st_traces_check = []
        	    st_hours = []
        	    for k in st_traces:
        	    	if len(k[0].data) >= WINDOW_LENGTH*100:
        	    		k[0].data = k[0].data[:WINDOW_LENGTH*100]
        	    		st_traces_check.append(k)
        	    		st_hours.append(str(k[0].stats.starttime.hour)+':'+str(k[0].stats.starttime.minute))
        	    	else:
        	    		pass

        	    if len(st_hours) > MIN_WINDOWS:
        		    inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))
        		    coordinates_lst = inv[0][0]

        		    traces_resp = [tr.remove_response(inventory=inv,output="DISP",water_level=60) for tr in st_traces_check]
        		    traces_demean = [tr.detrend('demean') for tr in traces_resp]
        		    traces_detrend = [tr.detrend('linear') for tr in traces_demean]
        		    traces_filter = [tr.filter('bandpass', freqmin=0.01,freqmax=10,corners=2, zerophase=True) for tr in traces_detrend]
        		    traces_resample = [tr.resample(NEW_SAMPLING_RATE) for tr in traces_filter]

        		    # ===================
        		    # Amplitude  clipping
        		    # ===================

        		    if Amp_clip:
        			    for i,tr in enumerate(traces_resample):
        			    	lim = CLIP_FACTOR * np.std(tr[0].data)
        			    	tr[0].data[tr[0].data > lim] = lim
        			    	tr[0].data[tr[0].data < -lim] = -lim

        		    # ======================
        		    # One-bit normalization
        		    # ======================

        		    if onebit_norm:
        		    	for i,tr in enumerate(traces_resample):
        		    		tr[0].data = np.sign(tr[0].data)

        		    # ==================
        		    # Spectral whitening
        		    # ==================

        		    if white_spectral:
        		    	freqmin=0.05
        		    	freqmax=0.5
        		    	for i,tr in enumerate(traces_resample):
        		    		n = len(tr[0].data)
        		    		nsamp = tr[0].stats.sampling_rate
        		    		frange = float(freqmax) - float(freqmin)
        		    		nsmo = int(np.fix(min(0.01, 0.5 * (frange)) * float(n) / nsamp))
        		    		f = np.arange(n) * nsamp / (n - 1.)
        		    		JJ = ((f > float(freqmin)) & (f<float(freqmax))).nonzero()[0]

        		    		# signal FFT
        		    		FFTs = fft(tr[0].data)
        		    		FFTsW = np.zeros(n) + 1j * np.zeros(n)

        		    		# Apodization to the left with cos^2 (to smooth the discontinuities)
        		    		smo1 = (np.cos(np.linspace(np.pi/2, np.pi, nsmo+1))**2)
        		    		FFTsW[JJ[0]:JJ[0]+nsmo+1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0]+nsmo+1]))

        		    		# boxcar
        		    		FFTsW[JJ[0]+nsmo+1:JJ[-1]-nsmo] = np.ones(len(JJ) - 2 * (nsmo+1))\
        		    		* np.exp(1j * np.angle(FFTs[JJ[0]+nsmo+1:JJ[-1]-nsmo]))

        		    		# Apodization to the right with cos^2 (to smooth the discontinuities)
        		    		smo2 = (np.cos(np.linspace(0, np.pi/2, nsmo+1))**2)
        		    		espo = np.exp(1j * np.angle(FFTs[JJ[-1]-nsmo:JJ[-1]+1]))
        		    		FFTsW[JJ[-1]-nsmo:JJ[-1]+1] = smo2 * espo

        		    		whitedata = 2. * ifft(FFTsW).real

        		    		tr[0].data = np.require(whitedata, dtype="float32")
        		    	traces_white_spectral = traces_resample

        		    traces_data_day = np.array([ton[0].data for ton in traces_white_spectral])

        		    os.makedirs(output_DATA_DAY,exist_ok=True)

        		    ds = ASDFDataSet(output_DATA_DAY+'DATA_DAY_'+network+'_'+name+'_'+sta_channel+'_'+year_day+'_'+julday_day+'.h5', compression="gzip-3")

                    #Adding Waveforms
        		    ds.add_waveforms(st,tag='raw_recording')

                    # Adding Auxiliary Data
        		    # Name to identify the particular piece of data.
        		    path_SD = sta_channel_id

        		    # Any additional parameters as a Python dictionary which will end up as
        		    # attributes of the array.
        		    parameters_SD = {'SD':'Array with 1 hour slices.', 'latitude': coordinates_lst.latitude, 'longitude': coordinates_lst.longitude, 'time_day':time_day, 'hours_day': st_hours}

        		    ds.add_auxiliary_data(data=traces_data_day,data_type='StationOneHourData',path=path_SD, parameters=parameters_SD)

        except:
            print('Problem with the file: '+f)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Creating Dictionaries to allocate results ###
def nested_dict():
    return collections.defaultdict(nested_dict)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

def crosscorr_func(stationtrace_pairs,name_suffix='CROSS_CORR_DAY_FILES'):
    """
    Gets Cross-correlation daily data

    @type stationtrace_pairs: path of the asdf file (str)
    """
    CrossCorrelation_dic = nested_dict()

    sta1_asdf = ASDFDataSet(stationtrace_pairs[0],mode='r')
    sta2_asdf = ASDFDataSet(stationtrace_pairs[1],mode='r')

    sta1 = sta1_asdf.auxiliary_data.StationOneHourData.list()[0]
    sta2 = sta2_asdf.auxiliary_data.StationOneHourData.list()[0]

    sta1_parameters = sta1_asdf.auxiliary_data.StationOneHourData[sta1].parameters
    sta2_parameters = sta2_asdf.auxiliary_data.StationOneHourData[sta2].parameters

    sta1_data_day = sta1_asdf.auxiliary_data.StationOneHourData[sta1].data[::]
    sta2_data_day = sta2_asdf.auxiliary_data.StationOneHourData[sta2].data[::]

    year_day = sta1_parameters['time_day'].split('.')[0]
    julday_day = sta1_parameters['time_day'].split('.')[1]

    # ----------------------------------------------------------------------------------------------------------------------------------------------
    #Check if file exists
    output_CrossCorrelation_DAY = ASDF_FILES+name_suffix+'/'+year_day+'.'+julday_day+'/'
    if os.path.isfile(output_CrossCorrelation_DAY+name_suffix+'_'+sta1+'_'+sta2+'_'+year_day+'_'+julday_day+'.h5'):
        pass

    else:

        day_crosscor_causal = CrossCorrelation(name1=sta1,name2=sta2,lat1=sta1_parameters['latitude'],lon1=sta1_parameters['longitude'],lat2=sta2_parameters['latitude'],lon2=sta2_parameters['longitude'],pair_time_day=sta1_parameters['time_day'])
        day_crosscor_acausal = CrossCorrelation(name1=sta2,name2=sta1,lat1=sta2_parameters['latitude'],lon1=sta2_parameters['longitude'],lat2=sta1_parameters['latitude'],lon2=sta1_parameters['longitude'],pair_time_day=sta1_parameters['time_day'])

        day_crosscor_causal.add(sta1_data_day,sta2_data_day,sta1_parameters['hours_day'].tolist(),sta2_parameters['hours_day'].tolist())
        day_crosscor_acausal.add(sta2_data_day,sta1_data_day,sta2_parameters['hours_day'].tolist(),sta1_parameters['hours_day'].tolist())

        if len(day_crosscor_acausal.dataarray) > 1:

            day_time_crosscor_all = day_crosscor_acausal.timearray+day_crosscor_causal.timearray
            day_data_crosscor_all = day_crosscor_acausal.dataarray+day_crosscor_causal.dataarray
            raw_SNR = SNR(day_data_crosscor_all,day_time_crosscor_all,day_crosscor_causal.dist(),vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)

            if raw_SNR > minspectSNR:

                output_CrossCorrelation_DAY = ASDF_FILES+name_suffix+'/'+year_day+'.'+julday_day+'/'
                os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)

                # ----------------------------------------------------------------------------------------------------------------------------------------------
                cc_asdf = ASDFDataSet(output_CrossCorrelation_DAY+name_suffix+'_'+sta1+'_'+sta2+'_'+year_day+'_'+julday_day+'.h5', compression="gzip-3")

                # Causal part of the CrossCorrelation
                path_CC_causal = sta1+'/'+sta2+'/'

                # Additional parameters of the causal part of the CrossCorrelation
                parameters_CC_causal = {
                            'CrossCorrelation':'Cross-correlation data between '+sta1+' and '+sta2+'.',
                            'dist': round(day_crosscor_causal.dist()),
                            'date': sta1_parameters['time_day'],
                            'sta1_loc': [sta1_parameters['latitude'],sta1_parameters['longitude']],
                            'sta1_name': sta1,
                            'sta2_loc': [sta2_parameters['latitude'],sta2_parameters['longitude']],
                            'sta2_name': sta2,
                            'crosscorr_daily_causal_time':day_crosscor_causal.timearray
                            }

                cc_asdf.add_auxiliary_data(data=day_crosscor_causal.dataarray,data_type='CrossCorrelation',path=path_CC_causal, parameters=parameters_CC_causal)

                # ----------------------------------------------------------------------------------------------------------------------------------------------

                # Acausal part of the CrossCorrelation
                path_CC_acausal = sta2+'/'+sta1+'/'

                # Additional parameters of the acausal part of the CrossCorrelation
                parameters_CC_acausal = {
                            'CrossCorrelation':'Cross-correlation data between '+sta2+' and '+sta1+'.',
                            'dist': round(day_crosscor_causal.dist()),
                            'date': sta1_parameters['time_day'],
                            'sta1_loc': [sta1_parameters['latitude'],sta1_parameters['longitude']],
                            'sta1_name': sta1,
                            'sta2_loc': [sta2_parameters['latitude'],sta2_parameters['longitude']],
                            'sta2_name': sta2,
                            'crosscorr_daily_acausal_time': day_crosscor_acausal.timearray
                            }

                cc_asdf.add_auxiliary_data(data=day_crosscor_acausal.dataarray,data_type='CrossCorrelation',path=path_CC_acausal, parameters=parameters_CC_acausal)

                # ============================
                # Plot: map and pair crosscorr
                # ============================

                fig = plt.figure(figsize=(15, 15))
                fig.suptitle(sta1+'-'+sta2+' - Day - '+UTCDateTime(year=int(year_day),julday=int(julday_day)).strftime('%d/%m/%Y'),fontsize=20)

                gs = gridspec.GridSpec(2, 1,wspace=0.2, hspace=0.5)

                #-------------------------------------------

                map_loc = fig.add_subplot(gs[0],projection=ccrs.PlateCarree())

                LLCRNRLON_LARGE = -52
                URCRNRLON_LARGE = -28
                LLCRNRLAT_LARGE = -30
                URCRNRLAT_LARGE = -12

                map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
                map_loc.yaxis.set_ticks_position('both')
                map_loc.xaxis.set_ticks_position('both')

                map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,3), crs=ccrs.PlateCarree())
                map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,3), crs=ccrs.PlateCarree())
                map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)

                map_loc.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

                reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
                shape_1_SHP = list(reader_1_SHP.geometries())
                plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
                map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)

                # Use the cartopy interface to create a matplotlib transform object
                # for the Geodetic coordinate system. We will use this along with
                # matplotlib's offset_copy function to define a coordinate system which
                # translates the text by 25 pixels to the left.

                geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
                text_transform = offset_copy(geodetic_transform, units='dots', y=0,x=80)
                text_transform_mag = offset_copy(geodetic_transform, units='dots', y=15,x=15)

                map_loc.scatter(sta1_parameters['longitude'],sta1_parameters['latitude'], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
                map_loc.scatter(sta2_parameters['longitude'],sta2_parameters['latitude'], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
                map_loc.plot([sta1_parameters['longitude'],sta2_parameters['longitude']],[sta1_parameters['latitude'],sta2_parameters['latitude']],c='k', transform=ccrs.PlateCarree())

                map_loc.text(sta1_parameters['longitude'],sta1_parameters['latitude'], sta1,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
                map_loc.text(sta2_parameters['longitude'],sta2_parameters['latitude'], sta2,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)

                #-------------------------------------------

                ax = fig.add_subplot(gs[1])
                data_to_plot = np.flip(day_crosscor_acausal.dataarray)+day_crosscor_causal.dataarray
                time_to_plot = np.flip(day_crosscor_acausal.timearray)*-1 + day_crosscor_causal.timearray
                ax.plot(time_to_plot,data_to_plot,color='k')
                ax.set_xlabel('time (s)',fontsize=14)
                ax.set_title('Dist = '+str(round(day_crosscor_causal.dist()))+' km',fontsize=14)

                output_figure_CrossCorrelation_DAY = CLOCK_DRIFT_OUTPUT+name_suffix+'_FIGURES/'+year_day+'.'+julday_day+'/'
                os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
                fig.savefig(output_figure_CrossCorrelation_DAY+name_suffix+'_FIG_'+sta1+'_'+sta2+'_'+year_day+'_'+julday_day+'.png')
                plt.close('all')

            return sta1_parameters['time_day']

        else:
            print("Problem: CrossCorrelation between "+sta1+" and "+sta2+" in "+sta1_parameters['time_day'])

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

def crosscorr_10_days_stack_func(input_lst):

        """
        Stacking 10 days of cross-correlation daily data

        @type stationtrace_pairs: path of the asdf file (str)
        @type stack_date: date of the day - julday.year (str)
        @type name_suffix: name of the output folder (str)
        """

        stationtrace_pairs = input_lst[0]
        stack_date = input_lst[1]
        name_suffix = input_lst[2]

        year_day = stack_date.split('.')[0]
        julday_day = stack_date.split('.')[1]

        #Reading data
        sta1_sta2_asdf_files = [ASDFDataSet(i, mode='r') for i in stationtrace_pairs]
        sta1 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation.list()[0]
        sta2 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation.list()[1]

        dist_pair = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[sta1][sta2].parameters['dist']

        loc_sta1 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[sta1][sta2].parameters['sta1_loc']
        loc_sta2 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[sta1][sta2].parameters['sta2_loc']

        #Stacking data
        causal_lst = np.mean(np.array([i.auxiliary_data.CrossCorrelation[sta1][sta2].data for i in sta1_sta2_asdf_files]),axis=0)
        acausal_lst = np.mean(np.array([i.auxiliary_data.CrossCorrelation[sta2][sta1].data for i in sta1_sta2_asdf_files]),axis=0)

        causal_time = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[sta1][sta2].parameters['crosscorr_daily_causal_time']
        acausal_time = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[sta2][sta1].parameters['crosscorr_daily_acausal_time']

        # ----------------------------------------------------------------------------------------------------------------------------------------------
        #Check if file exists
        output_CrossCorrelation_DAY = ASDF_FILES+name_suffix+'/'+year_day+'.'+julday_day+'/'

        if os.path.isfile(output_CrossCorrelation_DAY+name_suffix+'_'+sta1+'_'+sta2+'_'+year_day+'_'+julday_day+'.h5'):
            pass

        else:

            output_CrossCorrelation_DAY = ASDF_FILES+name_suffix+'/'+year_day+'.'+julday_day+'/'
            os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)
            cc_asdf = ASDFDataSet(output_CrossCorrelation_DAY+name_suffix+'_'+sta1+'_'+sta2+'_'+year_day+'_'+julday_day+'.h5', compression="gzip-3")

            # -----------------------------------------------------------------------------------------------------------------------------------------------
            # Causal part of the CrossCorrelation
            path_CC_causal = sta1+'/'+sta2+'/'

            # Additional parameters of the causal part of the CrossCorrelation
            parameters_CC_causal = {
                                    'CrossCorrelation':'Cross-correlation data between '+sta1+' and '+sta2+'.',
                                    'dist': dist_pair,
                                    'date': stack_date,
                                    'sta1_loc': loc_sta1,
                                    'sta1_name': sta1,
                                    'sta2_loc': loc_sta2,
                                    'sta2_name': sta2,
                                    'crosscorr_daily_causal_time': causal_time
                                    }

            cc_asdf.add_auxiliary_data(data=causal_lst,data_type='CrossCorrelation',path=path_CC_causal, parameters=parameters_CC_causal)

            # -----------------------------------------------------------------------------------------------------------------------------------------------
            # Acausal part of the CrossCorrelation
            path_CC_acausal = sta2+'/'+sta1+'/'

            # Additional parameters of the acausal part of the CrossCorrelation
            parameters_CC_acausal = {
                                    'CrossCorrelation':'Cross-correlation data between '+sta2+' and '+sta1+'.',
                                    'dist': dist_pair,
                                    'date': stack_date,
                                    'sta1_loc': loc_sta1,
                                    'sta1_name': sta1,
                                    'sta2_loc': loc_sta2,
                                    'sta2_name': sta2,
                                    'crosscorr_daily_acausal_time': acausal_time
                                    }

            cc_asdf.add_auxiliary_data(data=acausal_lst,data_type='CrossCorrelation',path=path_CC_acausal, parameters=parameters_CC_acausal)

            return stack_date

# ----------------------------------------------------------------------------------------------------------------------------------------------

def crosscorr_stack_asdf(input):
    """
    Stacking crosscorrelation data
    @type crosscorr_pairs_data: list of ASDF files
    @type cross_name_suffix: name of the folder (str)
    """
    crosscorr_pairs_data = input[0]
    cross_name_suffix = input[1]

    #Reading data
    sta1_sta2_asdf_files = []
    for i in crosscorr_pairs_data:
        try:
            sta1_sta2_asdf_files.append(ASDFDataSet(i, mode='r'))
        except:
            print('Problem in file: '+i)

    name_sta1 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation.list()[0]
    name_sta2 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation.list()[1]

    dist_pair = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[name_sta1][name_sta2].parameters['dist']

    loc_sta1 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[name_sta1][name_sta2].parameters['sta1_loc']
    loc_sta2 = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[name_sta1][name_sta2].parameters['sta2_loc']

    #Stacking data
    causal_lst = np.mean(np.array([i.auxiliary_data.CrossCorrelation[name_sta1][name_sta2].data for i in sta1_sta2_asdf_files]),axis=0)
    acausal_lst = np.mean(np.array([i.auxiliary_data.CrossCorrelation[name_sta2][name_sta1].data for i in sta1_sta2_asdf_files]),axis=0)

    causal_time = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[name_sta1][name_sta2].parameters['crosscorr_daily_causal_time']
    acausal_time = sta1_sta2_asdf_files[0].auxiliary_data.CrossCorrelation[name_sta2][name_sta1].parameters['crosscorr_daily_acausal_time']

    data_causal = causal_lst
    causal_time = causal_time

    data_acausal = acausal_lst[::-1]
    acausal_time = acausal_time[::-1]*-1

    # ----------------------------------------------------------------------------------------------------------------------------------------------

    output_CrossCorrelation_DAY = ASDF_FILES+cross_name_suffix+'/'+name_sta1+'.'+name_sta2+'/'
    os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)

    cc_asdf = ASDFDataSet(output_CrossCorrelation_DAY+cross_name_suffix+'_'+name_sta1+'_'+name_sta2+'.h5', compression="gzip-3")

    # Causal part of the CrossCorrelation
    path_CC_stacked_causal = name_sta1+'/'+name_sta2+'/'

    # Additional parameters of the causal part of the CrossCorrelation
    parameters_CC_stacked_causal = {
                            'CrossCorrelation':'Causal part of the cross-correlation data stacked between '+name_sta1+' and '+name_sta2+'.',
                            'dist': dist_pair,
                            'number_days': len(crosscorr_pairs_data),
                            'sta1_name': name_sta1,
                            'sta2_name': name_sta2,
                            'sta1_loc': loc_sta1,
                            'sta2_loc': loc_sta2,
                            'crosscorr_stack_time': causal_time
                            }

    cc_asdf.add_auxiliary_data(data=data_causal,data_type='CrossCorrelationStacked',path=path_CC_stacked_causal, parameters=parameters_CC_stacked_causal)

    # ----------------------------------------------------------------------------------------------------------------------------------------------------

    # Acausal part of the CrossCorrelation
    path_CC_stacked_acausal = name_sta2+'/'+name_sta1+'/'

    # Additional parameters of the acausal part of the CrossCorrelation
    parameters_CC_stacked_acausal = {
                            'CrossCorrelation':'Acausal part of the cross-correlation data stacked between '+name_sta2+' and '+name_sta1+'.',
                            'dist': dist_pair,
                            'number_days': len(crosscorr_pairs_data),
                            'sta1_name': name_sta1,
                            'sta2_name': name_sta2,
                            'sta1_loc': loc_sta1,
                            'sta2_loc': loc_sta2,
                            'crosscorr_stack_time': acausal_time
                            }

    cc_asdf.add_auxiliary_data(data=data_acausal,data_type='CrossCorrelationStacked',path=path_CC_stacked_acausal, parameters=parameters_CC_stacked_acausal)

    return [name_sta1,name_sta2]

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_stacked_cc_interstation_distance(folder_name):
    '''
    Plotting Stacked Cross-correlations according to interstation distance
    @type folder_name: folder of the cross-correlations files (str)
    '''
    chan_lst = ['HHE-HHE','HHN-HHN','HHZ-HHZ','HHE-HHN','HHE-HHZ','HHN-HHZ']

    #Collecting daily list of cross-correlations
    crosscorr_days_lst = sorted(glob.glob(ASDF_FILES+folder_name+'/*'))

    crosscorr_pairs_lst = []
    for i,j in enumerate(crosscorr_days_lst):
    	crosscorr_file = sorted(glob.glob(j+'/*'))
    	crosscorr_pairs_lst.append(crosscorr_file)

    #Make a list of list flat
    crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

    HHE_HHE_lst = []
    HHN_HHN_lst = []
    HHZ_HHZ_lst = []

    HHE_HHN_lst = []
    HHE_HHZ_lst = []
    HHN_HHZ_lst = []

    for i in crosscorr_pairs:
        # splitting subdir/basename
        subdir, filename = os.path.split(i)
        name_pair1 = filename.split('.h5')[0].split('_')[-1]
        name_pair2 = filename.split('.h5')[0].split('_')[-2]

        if 'HHE' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HH2' in name_pair2 or 'HH2' in name_pair1 and 'HHE' in name_pair2  or 'HHE' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHE_lst.append(i)

        if 'HHN' in name_pair1 and 'HHN' in name_pair2 or 'HH1' in name_pair1 and 'HH1' in name_pair2 or 'HHN' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHN' in name_pair2:
            HHN_HHN_lst.append(i)

        if 'HHZ' in name_pair1 and 'HHZ' in name_pair2:
            HHZ_HHZ_lst.append(i)

        if 'HHE' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HH2' in name_pair2 or 'HHE' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHN_lst.append(i)

        if 'HHE' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHZ_lst.append(i)

        if 'HHZ' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHZ' in name_pair2:
            HHN_HHZ_lst.append(i)

    CHANNEL_fig_lst = [HHE_HHE_lst,HHN_HHN_lst,HHZ_HHZ_lst,HHE_HHN_lst,HHE_HHZ_lst,HHN_HHZ_lst]

    for ipairs,crosscorr_pairs in enumerate(CHANNEL_fig_lst):
        #Creating the figure
        fig = plt.figure(figsize=(10, 20))
        fig.suptitle('Cross-correlations according to interstation distance: '+chan_lst[ipairs],fontsize=20)

        gs = gridspec.GridSpec(5, 2,wspace=0.2, hspace=0.5)

        #-------------------------------------------

        map_loc = fig.add_subplot(gs[:2,:],projection=ccrs.PlateCarree())

        LLCRNRLON_LARGE = -52
        URCRNRLON_LARGE = -38
        LLCRNRLAT_LARGE = -30
        URCRNRLAT_LARGE = -12

        map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
        map_loc.yaxis.set_ticks_position('both')
        map_loc.xaxis.set_ticks_position('both')

        map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,3), crs=ccrs.PlateCarree())
        map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,3), crs=ccrs.PlateCarree())
        map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)

        map_loc.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

        reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
        shape_1_SHP = list(reader_1_SHP.geometries())
        plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
        map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)

        crosscorr_stack_style_lst = []
        crosscorr_stack_name_lst = []
        crosscorr_stack_data_normalized_lst = []
        crosscorr_stack_data_normalized_dist_lst = []
        crosscorr_stack_data_normalized_vmin_lst = []
        crosscorr_stack_data_normalized_vmax_lst = []
        time_to_plot = []

        for i in tqdm(crosscorr_pairs, desc='Reading '+chan_lst[ipairs]+' data'):

            #Reading data
            sta1_sta2_asdf_file = ASDFDataSet(i, mode='r')

            name_sta1 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked.list()[0]
            name_sta2 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked.list()[1]

            if 'OBS' in name_sta1 or 'OBS' in name_sta2:
                crosscorr_stack_style_lst.append('k')
                crosscorr_stack_name_lst.append(name_sta1+'-'+name_sta2)
                dist_pair = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['dist']

                loc_sta1 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['sta1_loc']
                loc_sta2 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['sta2_loc']

        	    #Stacked data
                data_causal = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].data[::]
                causal_time = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['crosscorr_stack_time']

                data_acausal = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta2][name_sta1].data[::]
                acausal_time = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta2][name_sta1].parameters['crosscorr_stack_time']

                crosscorr_stack_data = data_acausal + data_causal
                crosscorr_stack_time = acausal_time + causal_time

                try:
                    OBS_SNR = SNR(crosscorr_stack_data,crosscorr_stack_time,dist_pair,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)
                    if OBS_SNR >= 1:

                        crosscorr_stack_data_normalized_lst.append(crosscorr_stack_data)
                        time_to_plot.append(crosscorr_stack_time)

                        crosscorr_stack_data_normalized_dist_lst.append(dist_pair)
                        crosscorr_stack_data_normalized_vmin_lst.append(dist_pair/SIGNAL_WINDOW_VMIN)
                        crosscorr_stack_data_normalized_vmax_lst.append(dist_pair/SIGNAL_WINDOW_VMAX)

                        # Use the cartopy interface to create a matplotlib transform object
                        # for the Geodetic coordinate system. We will use this along with
                        # matplotlib's offset_copy function to define a coordinate system which
                        # translates the text by 25 pixels to the left.
                        geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
                        text_transform = offset_copy(geodetic_transform, units='dots', y=0,x=80)
                        text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-15,x=15)

                        map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k',alpha=0.5, transform=ccrs.PlateCarree())
                        map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
                        map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

                except:
                    pass


        #-------------------------------------------
        orglst = np.argsort(crosscorr_stack_data_normalized_dist_lst)
        crosscorr_stack_name_lst
        crosscorr_stack_name_org_lst = [crosscorr_stack_name_lst[i] for i in orglst]
        crosscorr_stack_style_org_lst = [crosscorr_stack_style_lst[i] for i in orglst]
        crosscorr_stack_data_normalized_dist_org_lst = [crosscorr_stack_data_normalized_dist_lst[i] for i in orglst]
        crosscorr_stack_data_normalized_vmax_org_lst = [crosscorr_stack_data_normalized_vmax_lst[i] for i in orglst]
        crosscorr_stack_data_normalized_vmin_org_lst = [crosscorr_stack_data_normalized_vmin_lst[i] for i in orglst]
        crosscorr_stack_data_normalized_org_lst = [crosscorr_stack_data_normalized_lst[i] for i in orglst]


        #--------------------------------------------------------------------------------------------------------------------
        ax2 = fig.add_subplot(gs[2,0])
        crosscorr_stack_data_normalized_org_lsts = [bandpass(data_2_plot,1.0/PERIOD_BANDS[0][1], 1.0/PERIOD_BANDS[0][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
        crosscorr_stack_data_normalized_org_lst = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lsts]

        y_factor = 1
        sav_f = 51

        for i,j in enumerate(crosscorr_stack_data_normalized_org_lst):
        	ax2.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

        ax2.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))][::50])
        ax2.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::50])

        ax2.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')
        ax2.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')

        ax2.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')
        ax2.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')

        ax2.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=0.5)
        ax2.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

        # adding labels
        ax2.set_xlabel('Lapse time (s)',fontsize=14)
        ax2.set_ylabel('Distance (km)',fontsize=14)
        ax2.set_title('Filter: '+str(PERIOD_BANDS[0][0])+'-'+str(PERIOD_BANDS[0][1])+'s')

        #---------------------------------------------------------
        ax3 = fig.add_subplot(gs[2,1])

        vector_plot = np.array(crosscorr_stack_data_normalized_org_lst)

        extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst)]
        im = ax3.imshow(vector_plot,extent=extent,origin='lower', interpolation='sinc',aspect='auto',cmap='RdBu')
        ax3.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='--')
        ax3.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

        ax3.set_yticks([i for i in range(len(crosscorr_stack_data_normalized_org_lst))][::50])
        ax3.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::50])

        ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')
        ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')

        ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')
        ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')

        # adding labels
        ax3.set_xlabel('Lapse time (s)',fontsize=14)

        axins = inset_axes(ax3,
                           width="30%",  # width = 10% of parent_bbox width
                           height="2%",  # height : 5%
                           loc='upper left',
                           bbox_to_anchor=(0.65, 0.03, 1, 1),
                           bbox_transform=ax3.transAxes,
                           borderpad=0,
                           )
        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

        #-----------------------------------------------------------------------------------------------
        crosscorr_stack_data_normalized_org_lst_20_50s = [bandpass(data_2_plot, 1.0/PERIOD_BANDS[1][1], 1.0/PERIOD_BANDS[1][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
        crosscorr_stack_data_normalized_org_lst_20_50 = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lst_20_50s]

        ax4 = fig.add_subplot(gs[3,0])

        for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_20_50):
        	ax4.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst_20_50[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

        ax4.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))][::50])
        ax4.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::50])

        ax4.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')
        ax4.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')

        ax4.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')
        ax4.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')

        ax4.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
        ax4.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

        # adding labels
        ax4.set_xlabel('Lapse time (s)',fontsize=14)
        ax4.set_ylabel('Distance (km)',fontsize=14)
        ax4.set_title('Filter: '+str(PERIOD_BANDS[1][0])+'-'+str(PERIOD_BANDS[1][1])+'s')

        #---------------------------------------------------------
        ax5 = fig.add_subplot(gs[3,1])

        vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_20_50)

        extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst_20_50)]
        im = ax5.imshow(vector_plot,extent=extent,origin='lower', interpolation='sinc',aspect='auto',cmap='RdBu')
        ax5.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='--')
        ax5.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

        ax5.set_yticks([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))][::50])
        ax5.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::50])

        ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')
        ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')

        ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')
        ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')

        # adding labels
        ax5.set_xlabel('Lapse time (s)',fontsize=14)

        axins = inset_axes(ax5,
                           width="30%",  # width = 10% of parent_bbox width
                           height="2%",  # height : 5%
                           loc='upper left',
                           bbox_to_anchor=(0.65, 0.03, 1, 1),
                           bbox_transform=ax5.transAxes,
                           borderpad=0,
                           )
        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

        #-----------------------------------------------------------------------------------------------
        crosscorr_stack_data_normalized_org_lst_50_100s = [bandpass(data_2_plot, 1.0/PERIOD_BANDS[2][1], 1.0/PERIOD_BANDS[2][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
        crosscorr_stack_data_normalized_org_lst_50_100 = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lst_50_100s]

        ax6 = fig.add_subplot(gs[4,0])
        for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_50_100):
        	ax6.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst_50_100[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

        ax6.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))][::50])
        ax6.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::50])

        ax6.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')
        ax6.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')

        ax6.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')
        ax6.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')

        ax6.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
        ax6.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

        # adding labels
        ax6.set_xlabel('Lapse time (s)',fontsize=14)
        ax6.set_ylabel('Distance (km)',fontsize=14)
        ax6.set_title('Filter: '+str(PERIOD_BANDS[2][0])+'-'+str(PERIOD_BANDS[2][1])+'s')

        #---------------------------------------------------------
        ax7 = fig.add_subplot(gs[4,1])

        vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_50_100)

        extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst_50_100)]
        im = ax7.imshow(vector_plot,extent=extent,origin='lower', interpolation='sinc',aspect='auto',cmap='RdBu')
        ax7.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='--')
        ax7.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

        ax7.set_yticks([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))][::50])
        ax7.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::50])

        ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')
        ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')

        ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')
        ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')


        # adding labels
        ax7.set_xlabel('Lapse time (s)',fontsize=14)

        axins = inset_axes(ax7,
                           width="30%",  # width = 10% of parent_bbox width
                           height="2%",  # height : 5%
                           loc='upper left',
                           bbox_to_anchor=(0.65, 0.03, 1, 1),
                           bbox_transform=ax7.transAxes,
                           borderpad=0,
                           )
        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

        output_figure_CrossCorrelation_DAY = CLOCK_DRIFT_OUTPUT+'CROSS_CORR_STACK_INTERSTATION_DISTANCE_FIGURES/'
        os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
        fig.savefig(output_figure_CrossCorrelation_DAY+folder_name+'_INTERSTATION_DISTANCE_FIG_'+chan_lst[ipairs]+'.png',dpi=300)
        plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_stacked_cc_interstation_distance_per_obs(folder_name):
    '''
    Plotting Stacked Cross-correlations according to interstation distance
    @type folder_name: folder of the cross-correlations files (str)
    '''

    #Collecting daily list of cross-correlations
    crosscorr_days_lst = sorted(glob.glob(ASDF_FILES+folder_name+'/*'))

    crosscorr_pairs_lst = []
    for i,j in enumerate(crosscorr_days_lst):
    	crosscorr_file = sorted(glob.glob(j+'/*'))
    	crosscorr_pairs_lst.append(crosscorr_file)

    #Make a list of list flat
    crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

    HHE_HHE_lst = []
    HHN_HHN_lst = []
    HHZ_HHZ_lst = []

    HHE_HHN_lst = []
    HHE_HHZ_lst = []
    HHN_HHZ_lst = []

    for i in crosscorr_pairs:
        # splitting subdir/basename
        subdir, filename = os.path.split(i)
        name_pair1 = filename.split('.h5')[0].split('_')[-1]
        name_pair2 = filename.split('.h5')[0].split('_')[-2]

        if 'HHE' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HH2' in name_pair2 or 'HH2' in name_pair1 and 'HHE' in name_pair2  or 'HHE' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHE_lst.append(i)

        if 'HHN' in name_pair1 and 'HHN' in name_pair2 or 'HH1' in name_pair1 and 'HH1' in name_pair2 or 'HHN' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHN' in name_pair2:
            HHN_HHN_lst.append(i)

        if 'HHZ' in name_pair1 and 'HHZ' in name_pair2:
            HHZ_HHZ_lst.append(i)

        if 'HHE' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HH2' in name_pair2 or 'HHE' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHN_lst.append(i)

        if 'HHE' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHZ_lst.append(i)

        if 'HHZ' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHZ' in name_pair2:
            HHN_HHZ_lst.append(i)

    CHANNEL_fig_lst = [HHE_HHE_lst,HHN_HHN_lst,HHZ_HHZ_lst,HHE_HHN_lst,HHE_HHZ_lst,HHN_HHZ_lst]
    chan_lst = ['HHE-HHE','HHN-HHN','HHZ-HHZ','HHE-HHN','HHE-HHZ','HHN-HHZ']

    for iOBS in OBS_LST:

        for ipairs,crosscorr_pairs in enumerate(CHANNEL_fig_lst):
            #Creating the figure
            fig = plt.figure(figsize=(10, 20))
            fig.suptitle(iOBS+' cross-correlations according to interstation distance: '+chan_lst[ipairs],fontsize=20)

            gs = gridspec.GridSpec(5, 2,wspace=0.2, hspace=0.5)

            #-------------------------------------------

            map_loc = fig.add_subplot(gs[:2,:],projection=ccrs.PlateCarree())

            LLCRNRLON_LARGE = -52
            URCRNRLON_LARGE = -38
            LLCRNRLAT_LARGE = -30
            URCRNRLAT_LARGE = -12

            map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
            map_loc.yaxis.set_ticks_position('both')
            map_loc.xaxis.set_ticks_position('both')

            map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,3), crs=ccrs.PlateCarree())
            map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,3), crs=ccrs.PlateCarree())
            map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)

            map_loc.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

            reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
            shape_1_SHP = list(reader_1_SHP.geometries())
            plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
            map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)

            crosscorr_stack_style_lst = []
            crosscorr_stack_name_lst = []
            crosscorr_stack_data_normalized_lst = []
            crosscorr_stack_data_normalized_dist_lst = []
            crosscorr_stack_data_normalized_vmin_lst = []
            crosscorr_stack_data_normalized_vmax_lst = []
            time_to_plot = []

            for i in tqdm(crosscorr_pairs, desc='Reading '+chan_lst[ipairs]+' data'):

                #Reading data
                sta1_sta2_asdf_file = ASDFDataSet(i, mode='r')

                name_sta1 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked.list()[0]
                name_sta2 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked.list()[1]

                if iOBS in name_sta1 or iOBS in name_sta2:
                    crosscorr_stack_style_lst.append('k')
                    crosscorr_stack_name_lst.append(name_sta1+'-'+name_sta2)
                    dist_pair = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['dist']

                    loc_sta1 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['sta1_loc']
                    loc_sta2 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['sta2_loc']

            	    #Stacked data
                    data_causal = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].data[::]
                    causal_time = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['crosscorr_stack_time']

                    data_acausal = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta2][name_sta1].data[::]
                    acausal_time = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta2][name_sta1].parameters['crosscorr_stack_time']

                    crosscorr_stack_data = data_acausal + data_causal
                    crosscorr_stack_time = acausal_time + causal_time

                    try:
                        OBS_SNR = SNR(crosscorr_stack_data,crosscorr_stack_time,dist_pair,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)
                        if OBS_SNR >= 1:

                            crosscorr_stack_data_normalized_lst.append(crosscorr_stack_data)
                            time_to_plot.append(crosscorr_stack_time)

                            crosscorr_stack_data_normalized_dist_lst.append(dist_pair)
                            crosscorr_stack_data_normalized_vmin_lst.append(dist_pair/SIGNAL_WINDOW_VMIN)
                            crosscorr_stack_data_normalized_vmax_lst.append(dist_pair/SIGNAL_WINDOW_VMAX)

                            # Use the cartopy interface to create a matplotlib transform object
                            # for the Geodetic coordinate system. We will use this along with
                            # matplotlib's offset_copy function to define a coordinate system which
                            # translates the text by 25 pixels to the left.
                            geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
                            text_transform = offset_copy(geodetic_transform, units='dots', y=0,x=80)
                            text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-15,x=15)

                            map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k',alpha=0.5, transform=ccrs.PlateCarree())
                            map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
                            map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

                    except:
                        pass


            #-------------------------------------------
            orglst = np.argsort(crosscorr_stack_data_normalized_dist_lst)
            crosscorr_stack_name_org_lst = [crosscorr_stack_name_lst[i] for i in orglst]
            crosscorr_stack_style_org_lst = [crosscorr_stack_style_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_dist_org_lst = [crosscorr_stack_data_normalized_dist_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_vmax_org_lst = [crosscorr_stack_data_normalized_vmax_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_vmin_org_lst = [crosscorr_stack_data_normalized_vmin_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_org_lst = [crosscorr_stack_data_normalized_lst[i] for i in orglst]
            #--------------------------------------------------------------------------------------------------------------------
            #plot parameters
            y_factor = 1
            sav_f = 11
            cmap_interpolation = 'RdBu'
            type_interpolation = 'sinc'
            #--------------------------------------------------------------------------------------------------------------------

            ax2 = fig.add_subplot(gs[2,0])
            crosscorr_stack_data_normalized_org_lsts = [bandpass(data_2_plot,1.0/PERIOD_BANDS[0][1], 1.0/PERIOD_BANDS[0][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
            crosscorr_stack_data_normalized_org_lst = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lsts]

            for i,j in enumerate(crosscorr_stack_data_normalized_org_lst):
            	ax2.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

            ax2.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')
            ax2.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')

            ax2.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')
            ax2.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax2.set_ylabel('Distance (km)',fontsize=14)
            ax2.set_title('Filter: '+str(PERIOD_BANDS[0][0])+'-'+str(PERIOD_BANDS[0][1])+'s')
            ax2.set_xlabel('Lapse time (s)',fontsize=14)
            ax2.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))])
            ax2.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst])
            ax2.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax2.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            #---------------------------------------------------------
            ax3 = fig.add_subplot(gs[2,1])

            vector_plot = np.array(crosscorr_stack_data_normalized_org_lst)

            extent = [-SHIFT_LEN,SHIFT_LEN,min(crosscorr_stack_data_normalized_dist_org_lst),max(crosscorr_stack_data_normalized_dist_org_lst)]
            im = ax3.imshow(vector_plot,extent=extent,origin='lower', interpolation=type_interpolation,aspect='auto',cmap=cmap_interpolation)

            ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax3.set_xlabel('Lapse time (s)',fontsize=14)
            ax3.set_yticks(crosscorr_stack_data_normalized_dist_org_lst)
            ax3.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst])
            ax3.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax3.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            axins = inset_axes(ax3,
                               width="30%",  # width = 10% of parent_bbox width
                               height="2%",  # height : 5%
                               loc='upper left',
                               bbox_to_anchor=(0.65, 0.03, 1, 1),
                               bbox_transform=ax3.transAxes,
                               borderpad=0,
                               )
            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            #-----------------------------------------------------------------------------------------------
            crosscorr_stack_data_normalized_org_lst_20_50s = [bandpass(data_2_plot, 1.0/PERIOD_BANDS[1][1], 1.0/PERIOD_BANDS[1][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
            crosscorr_stack_data_normalized_org_lst_20_50 = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lst_20_50s]

            ax4 = fig.add_subplot(gs[3,0])

            for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_20_50):
            	ax4.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst_20_50[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

            ax4.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')
            ax4.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')

            ax4.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')
            ax4.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax4.set_xlabel('Lapse time (s)',fontsize=14)
            ax4.set_ylabel('Distance (km)',fontsize=14)
            ax4.set_title('Filter: '+str(PERIOD_BANDS[1][0])+'-'+str(PERIOD_BANDS[1][1])+'s')
            ax4.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))])
            ax4.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst])
            ax4.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax4.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            #---------------------------------------------------------
            ax5 = fig.add_subplot(gs[3,1])

            vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_20_50)

            extent = [-SHIFT_LEN,SHIFT_LEN,min(crosscorr_stack_data_normalized_dist_org_lst),max(crosscorr_stack_data_normalized_dist_org_lst)]
            im = ax5.imshow(vector_plot,extent=extent,origin='lower', interpolation=type_interpolation,aspect='auto',cmap=cmap_interpolation)

            ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax5.set_xlabel('Lapse time (s)',fontsize=14)
            ax5.set_yticks(crosscorr_stack_data_normalized_dist_org_lst)
            ax5.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst])
            ax5.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax5.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            axins = inset_axes(ax5,
                               width="30%",  # width = 10% of parent_bbox width
                               height="2%",  # height : 5%
                               loc='upper left',
                               bbox_to_anchor=(0.65, 0.03, 1, 1),
                               bbox_transform=ax5.transAxes,
                               borderpad=0,
                               )
            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            #-----------------------------------------------------------------------------------------------
            crosscorr_stack_data_normalized_org_lst_50_100s = [bandpass(data_2_plot, 1.0/PERIOD_BANDS[2][1], 1.0/PERIOD_BANDS[2][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
            crosscorr_stack_data_normalized_org_lst_50_100 = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lst_50_100s]

            ax6 = fig.add_subplot(gs[4,0])
            for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_50_100):
            	ax6.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst_50_100[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

            ax6.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')
            ax6.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')

            ax6.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')
            ax6.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax6.set_xlabel('Lapse time (s)',fontsize=14)
            ax6.set_ylabel('Distance (km)',fontsize=14)
            ax6.set_title('Filter: '+str(PERIOD_BANDS[1][0])+'-'+str(PERIOD_BANDS[1][1])+'s')
            ax6.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))])
            ax6.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst])
            ax6.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax6.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            #---------------------------------------------------------
            ax7 = fig.add_subplot(gs[4,1])

            vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_50_100)

            extent = [-SHIFT_LEN,SHIFT_LEN,min(crosscorr_stack_data_normalized_dist_org_lst),max(crosscorr_stack_data_normalized_dist_org_lst)]
            im = ax7.imshow(vector_plot,extent=extent,origin='lower', interpolation=type_interpolation,aspect='auto',cmap=cmap_interpolation)

            ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax7.set_xlabel('Lapse time (s)',fontsize=14)
            ax7.set_yticks(crosscorr_stack_data_normalized_dist_org_lst)
            ax7.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst])
            ax7.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax7.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            axins = inset_axes(ax7,
                               width="30%",  # width = 10% of parent_bbox width
                               height="2%",  # height : 5%
                               loc='upper left',
                               bbox_to_anchor=(0.65, 0.03, 1, 1),
                               bbox_transform=ax7.transAxes,
                               borderpad=0,
                               )
            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            output_figure_CrossCorrelation_DAY = CLOCK_DRIFT_OUTPUT+'CROSS_CORR_STACK_INTERSTATION_DISTANCE_FIGURES/'
            os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
            fig.savefig(output_figure_CrossCorrelation_DAY+folder_name+'_INTERSTATION_DISTANCE_FIG_'+chan_lst[ipairs]+'_'+iOBS+'.png',dpi=300)
            plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def Calculating_clock_drift_func(iOBS):
    '''
    Calculating clock drift from cross-correlation data
    @type input: name of the OBS (str)
    '''
    # -------------------------------------------
    # Collecting daily list of cross-correlations
    # -------------------------------------------

    crosscorr_pairs = sorted(glob.glob(ASDF_FILES+'CROSS_CORR_10_DAYS_FILES/**/*.h5', recursive=True))

    # --------------------------------
    # Separating according to OBS name
    # --------------------------------

    crosscorr_pairs_obs = [i for i in crosscorr_pairs if iOBS in i]

    # ------------------
    # Separating by pair
    # ------------------

    crosscorr_pairs_name_lst = []
    for i in crosscorr_pairs_obs:

        # splitting subdir/basename
        subdir, filename = os.path.split(i)
        nameslst = filename.split("_20")[0]

        name_sta1 = nameslst.split('_')[-2].split('..')[0]
        name_sta2 = nameslst.split('_')[-1].split('..')[0]

        if name_sta1 != name_sta2:
            crosscorr_pairs_name_lst.append(name_sta1+'_'+name_sta2)

    crosscorr_pairs_names = sorted(list(set(crosscorr_pairs_name_lst)))

    for ipair in crosscorr_pairs_names:

        pair_sta_1 = ipair.split('_')[0].split('..')[0]
        pair_sta_2 = ipair.split('_')[1].split('..')[0]

        # ---------------------
        # Separating by channel
        # ---------------------

        HHE_HHE_lst = []
        HHN_HHN_lst = []
        HHZ_HHZ_lst = []

        HHE_HHN_lst = []
        HHN_HHE_lst = []

        HHE_HHZ_lst = []
        HHZ_HHE_lst = []

        HHN_HHZ_lst = []
        HHZ_HHN_lst = []

        for i in tqdm(crosscorr_pairs_obs,desc='pair: '+pair_sta_1+'-'+pair_sta_2):

            if pair_sta_1 and pair_sta_2 in i:

                # splitting subdir/basename
                subdir, filename = os.path.split(i)
                nameslst = filename.split("_20")[0]

                name_pair1 = nameslst.split('_')[-2]
                name_pair2 = nameslst.split('_')[-1]

                name1 = nameslst.split('_')[-2].split('..')[0]
                name2 = nameslst.split('_')[-1].split('..')[0]

                channel_sta1 = nameslst.split('_')[-2].split('..')[1]
                channel_sta2 = nameslst.split('_')[-1].split('..')[1]

                if pair_sta_1 == name1 and pair_sta_2 == name2:

                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHE' and channel_sta2 == 'HHE' or channel_sta1 == 'HH2' and channel_sta2 == 'HH2' or channel_sta1 == 'HH2' and channel_sta2 == 'HHE' or  channel_sta1 == 'HHE' and channel_sta2 == 'HH2':
                        HHE_HHE_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHN' and channel_sta2 == 'HHN' or channel_sta1 == 'HH1' and channel_sta2 == 'HH1' or channel_sta1 == 'HHN' and channel_sta2 == 'HH1' or channel_sta1 == 'HH1' and channel_sta2 == 'HHN':
                        HHN_HHN_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHZ' and channel_sta2 == 'HHZ':
                        HHZ_HHZ_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHE' and channel_sta2 == 'HHN' or channel_sta1 == 'HH2' and channel_sta2 == 'HHN' or channel_sta1 == 'HHE' and channel_sta2 == 'HH1' or channel_sta1 == 'HH2' and channel_sta2 == 'HH1':
                        HHE_HHN_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHN' and channel_sta2 == 'HHE' or channel_sta1 == 'HH1' and channel_sta2 == 'HHE' or channel_sta1 == 'HHN' and channel_sta2 == 'HH2' or channel_sta1 == 'HH1' and channel_sta2 == 'HH2':
                        HHN_HHE_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHE' and channel_sta2 == 'HHZ' or channel_sta1 == 'HH2' and channel_sta2 == 'HHZ':
                        HHE_HHZ_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHZ' and channel_sta2 == 'HHE' or channel_sta1 == 'HHZ' and channel_sta2 == 'HH2':
                        HHZ_HHE_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHN' and channel_sta2 == 'HHZ' or channel_sta1 == 'HH1' and channel_sta2 == 'HHZ':
                        HHN_HHZ_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------
                    if channel_sta1 == 'HHZ' and channel_sta2 == 'HHN' or channel_sta1 == 'HHZ' and channel_sta2 == 'HH1':
                        HHZ_HHN_lst.append(i)
                    # ------------------------------------------------------------------------------------------------------

        CHANNEL_fig_lst = [HHE_HHE_lst,HHE_HHN_lst,HHE_HHZ_lst,HHN_HHN_lst,HHN_HHE_lst,HHN_HHZ_lst,HHZ_HHE_lst,HHZ_HHN_lst,HHZ_HHZ_lst]
        chan_lst = ['HHE-HHE','HHE-HHN','HHE-HHZ','HHN-HHN','HHN-HHE','HHN-HHZ','HHZ-HHE','HHZ-HHN','HHZ-HHZ']

        # ------------------------------------------------------------------------------------------------------
        # Starting the value list and the columns_headers:
        columns_headers = ['sta_1','sta_2']

        data_drift_lst = [pair_sta_1,pair_sta_2]

        # ---------------------
        # Calculating the drift
        # ---------------------

        for idch, i in enumerate(CHANNEL_fig_lst):

            column_date_n = chan_lst[idch]+' date'
            column_coefficient_n = chan_lst[idch]+' coefficient'
            column_shift_n = chan_lst[idch]+' shift'

            columns_headers.extend([column_date_n,column_coefficient_n,column_shift_n])

            if len(i) > 1:

                # ------------
                # Reading data
                # ------------

                crosscorr_pair_date_filename = [filename.split('/')[-1] for filename in i]
                crosscorr_pair_date = [datetime.datetime.strptime(filename.split('/')[-1].split('_')[-2]+'.'+filename.split('/')[-1].split('_')[-1].split('.')[0], '%Y.%j') for filename in crosscorr_pair_date_filename]

                sta1_sta2_asdf_file = [ASDFDataSet(k, mode='r') for k in i]

                name_sta1 = [file.auxiliary_data.CrossCorrelation.list()[0] for file in sta1_sta2_asdf_file]
                name_sta2 = [file.auxiliary_data.CrossCorrelation.list()[1] for file in sta1_sta2_asdf_file]

                dist_pair = [sta1_sta2_asdf_file[id].auxiliary_data.CrossCorrelation[name_sta1[id]][name_sta2[id]].parameters['dist'] for id,jd in enumerate(name_sta1)]
                loc_sta1 = [sta1_sta2_asdf_file[id].auxiliary_data.CrossCorrelation[name_sta1[id]][name_sta2[id]].parameters['sta1_loc'] for id,jd in enumerate(name_sta1)]
                loc_sta2 = [sta1_sta2_asdf_file[id].auxiliary_data.CrossCorrelation[name_sta1[id]][name_sta2[id]].parameters['sta2_loc'] for id,jd in enumerate(name_sta1)]

                causal_time = [sta1_sta2_asdf_file[id].auxiliary_data.CrossCorrelation[name_sta1[id]][name_sta2[id]].parameters['crosscorr_daily_causal_time']  for id,jd in enumerate(name_sta1)]
                acausal_time = [sta1_sta2_asdf_file[id].auxiliary_data.CrossCorrelation[name_sta2[id]][name_sta1[id]].parameters['crosscorr_daily_acausal_time'] for id,jd in enumerate(name_sta1)]

                # ------------
                # Stacked data
                # ------------

                causal_lst = [jd.auxiliary_data.CrossCorrelation[name_sta1[id]][name_sta2[id]].data[::] for id,jd in enumerate(sta1_sta2_asdf_file)]
                acausal_lst = [jd.auxiliary_data.CrossCorrelation[name_sta2[id]][name_sta1[id]].data[::] for id,jd in enumerate(sta1_sta2_asdf_file)]

                # ----------------------------------------------------------------------------------------------------

                date_to_plot_clock = []
                data_to_plot_coefficient_clock_drift = []
                data_to_plot_shift_clock_drift = []

                for k in tqdm(range(len(causal_lst)), desc=chan_lst[idch]+' drift'):

                    data_acausal_causal = np.array(acausal_lst[k] + causal_lst[k])
                    data_normalized = np.array(Normalize(data_acausal_causal))
                    time_acausal_causal = np.array(acausal_time[k] + causal_time[k])

                    dist_pair_norm = dist_pair[k]
                    # --------------------------------------------------------
                    # Collecting daily list of 10-day stack cross-correlations
                    # --------------------------------------------------------

                    sta1_sta2_asdf_file_10_day = ASDFDataSet(glob.glob(ASDF_FILES+'CROSS_CORR_10_DAYS_STACKED_FILES/'+name_sta1[k]+'.'+name_sta2[k]+'/*')[0], mode='r')
                    stacked_10_day_data = sta1_sta2_asdf_file_10_day.auxiliary_data.CrossCorrelationStacked[name_sta2[k]][name_sta1[k]].data[::]+sta1_sta2_asdf_file_10_day.auxiliary_data.CrossCorrelationStacked[name_sta1[k]][name_sta2[k]].data[::]
                    stacked_10_day_time = sta1_sta2_asdf_file_10_day.auxiliary_data.CrossCorrelationStacked[name_sta2[k]][name_sta1[k]].parameters['crosscorr_stack_time'] +sta1_sta2_asdf_file_10_day.auxiliary_data.CrossCorrelationStacked[name_sta1[k]][name_sta2[k]].parameters['crosscorr_stack_time']
                    stacked_10_day_data_normalized = np.array(Normalize(stacked_10_day_data))

                    shift_clock_drift, coefficient_clock_drift = obscorr_window(stacked_10_day_data_normalized,stacked_10_day_time,data_normalized,time_acausal_causal,dist_pair_norm,SIGNAL_WINDOW_VMIN,SIGNAL_WINDOW_VMAX)

                    # --------------------------------------------------------------------------------------
                    date_to_plot_clock.append(crosscorr_pair_date[k])
                    data_to_plot_coefficient_clock_drift.append(coefficient_clock_drift)
                    data_to_plot_shift_clock_drift.append(shift_clock_drift)

                # --------------------------------------------------------------------------------------

                data_drift_lst.append(date_to_plot_clock)
                data_drift_lst.append(data_to_plot_coefficient_clock_drift)
                data_drift_lst.append(data_to_plot_shift_clock_drift)

            else:

                data_drift_lst.append([])
                data_drift_lst.append([])
                data_drift_lst.append([])

        # ----------------------------------------------------------------------------------------------------
        # Creating a Pandas DataFrame:
        column_info = [dist_pair[0],loc_sta1[0],loc_sta2[0]]
        for i in column_info:
            data_drift_lst.append(i)

        columns_headers1 = ['distance','loc_sta1','loc_sta2']
        for i in columns_headers1:
            columns_headers.append(i)


        clock_drift_df = pd.DataFrame(data_drift_lst, index=columns_headers).T
        # ----------------------------------------------------------------------------------------------------
        # Convert from pandas to Arrow and saving in feather formart file
        os.makedirs(FEATHER_FILES,exist_ok=True)
        file_feather_name = FEATHER_FILES+pair_sta_1+'.'+pair_sta_2+'_clock_drift_data.feather'
        feather.write_feather(clock_drift_df, file_feather_name)
        # ----------------------------------------------------------------------------------------------------

        if VERBOSE == True:

            # --------------------------------------------
            # Creating the figure and plotting Clock-drift
            # --------------------------------------------

            fig = plt.figure(figsize=(20, 15))
            fig.suptitle('Clock-drift: '+pair_sta_1+'-'+pair_sta_2+'('+str(round(dist_pair))+' km)',fontsize=20)
            fig.autofmt_xdate()
            # ----------------------------------------------------------------------------------------------------

            gs = gridspec.GridSpec(9, 2,wspace=0.5, hspace=0.8)
            map_loc = fig.add_subplot(gs[:,0],projection=ccrs.PlateCarree())

            LLCRNRLON_LARGE = -52
            URCRNRLON_LARGE = -38
            LLCRNRLAT_LARGE = -30
            URCRNRLAT_LARGE = -12

            map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
            map_loc.yaxis.set_ticks_position('both')
            map_loc.xaxis.set_ticks_position('both')

            map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+3,3), crs=ccrs.PlateCarree())
            map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+3,3), crs=ccrs.PlateCarree())
            map_loc.tick_params(labelbottom=True, labeltop=True, labelleft=True, labelright=True, labelsize=12)
            map_loc.grid(True,which='major',color='gray',linewidth=0.5,linestyle='--')

            reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
            shape_1_SHP = list(reader_1_SHP.geometries())
            plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
            map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)
            # Use the cartopy interface to create a matplotlib transform object
            # for the Geodetic coordinate system. We will use this along with
            # matplotlib's offset_copy function to define a coordinate system which
            # translates the text by 25 pixels to the left.
            geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
            text_transform = offset_copy(geodetic_transform, units='dots', y=50,x=100)

            map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k',alpha=0.5,transform=ccrs.PlateCarree())
            map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
            map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

            map_loc.text(loc_sta1[1],loc_sta1[0], pair_sta_1,fontsize=12,verticalalignment='center', horizontalalignment='right',transform=text_transform)
            map_loc.text(loc_sta2[1],loc_sta2[0], pair_sta_2,fontsize=12,verticalalignment='center', horizontalalignment='right',transform=text_transform)

            # ----------------------------------------------------------------------------------------------------
            days_major = DayLocator(interval=5)   # every 5 day
            days_minor = DayLocator(interval=1)   # every day
            months = MonthLocator(interval=3)  # every month
            yearsFmt = DateFormatter('%b-%Y')

            for z,x in enumerate(data_to_plot_coefficient_clock_drift_chan):
                if len(x) > 1:
                    # ----------------------------------------------------------------------------------------------------
                    ax0 = fig.add_subplot(gs[z,1])
                    ax0.xaxis.set_major_locator(months)
                    ax0.xaxis.set_major_formatter(yearsFmt)
                    ax0.xaxis.set_minor_locator(days_minor)
                    ax0.yaxis.set_major_locator(MultipleLocator(100))
                    ax0.yaxis.set_minor_locator(MultipleLocator(25))
                    ax0.set_ylabel('Drift ('+str(1/NEW_SAMPLING_RATE)+'s)')
                    ax0.set_title(chan_lst[z])
                    ax0.set_ylim(-200,200)

                    # -------------------------------------------------------------------------------------------------------------
                    poly_reg = PolynomialFeatures(degree=4)
                    X_poly = poly_reg.fit_transform(np.array(range(len(data_to_plot_shift_clock_drift_chan[z]))).reshape(-1, 1))
                    pol_reg = LinearRegression()
                    pol_reg.fit(X_poly, data_to_plot_shift_clock_drift_chan[z])
                    # -------------------------------------------------------------------------------------------------------------
                    for y,u in enumerate(data_to_plot_shift_clock_drift_chan[z]):
                        if data_to_plot_coefficient_clock_drift_chan[z][y] > 0.3:
                            im = ax0.scatter(date_to_plot_clock_chan[z][y],data_to_plot_shift_clock_drift_chan[z][y],c=data_to_plot_coefficient_clock_drift_chan[z][y],marker='o',edgecolors=None,cmap='magma',s=20,vmin=0,vmax=1,alpha=0.9)
                        else:
                            im = ax0.scatter(date_to_plot_clock_chan[z][y],data_to_plot_shift_clock_drift_chan[z][y],c=data_to_plot_coefficient_clock_drift_chan[z][y],marker='o',edgecolors=None,cmap='magma',s=10,vmin=0,vmax=1,alpha=0.7)

                    ax0.plot(date_to_plot_clock_chan[z], pol_reg.predict(poly_reg.fit_transform(np.array(range(len(data_to_plot_shift_clock_drift_chan[z]))).reshape(-1, 1))),'--b')

                    if z == 0:
                        axins = inset_axes(ax0,
                                               width="30%",  # width = 10% of parent_bbox width
                                               height="10%",  # height : 5%
                                               loc='upper left',
                                               bbox_to_anchor=(0.65,0.1, 1, 1),
                                               bbox_transform=ax0.transAxes,
                                               borderpad=0,
                                               )
                        plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

                else:

                    ax0 = fig.add_subplot(gs[z,1])
                    ax0.xaxis.set_major_locator(months)
                    ax0.xaxis.set_major_formatter(yearsFmt)
                    ax0.xaxis.set_minor_locator(days_minor)
                    ax0.yaxis.set_major_locator(MultipleLocator(100))
                    ax0.yaxis.set_minor_locator(MultipleLocator(25))
                    ax0.set_ylabel('Drift ('+str(1/NEW_SAMPLING_RATE)+'s)')
                    ax0.set_title(chan_lst[z])
                    ax0.set_ylim(-200,200)

            # -------------------------------------------------------------------------------------------------------------
            output_figure_CLOCK_DRIFT = CLOCK_DRIFT_OUTPUT+'CLOCK_DRIFT_FIGURES/'
            os.makedirs(output_figure_CLOCK_DRIFT,exist_ok=True)
            fig.savefig(output_figure_CLOCK_DRIFT+'CLOCK_DRIFT_BETWEEN_'+pair_sta_1+'_'+pair_sta_2+'.png',dpi=300)
            plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_stacked_cc_interstation_distance_per_obs_short(folder_name):
    '''
    Plotting Stacked Cross-correlations according to interstation distance
    @type folder_name: folder of the cross-correlations files (str)
    '''

    #Collecting daily list of cross-correlations
    crosscorr_days_lst = sorted(glob.glob(ASDF_FILES+folder_name+'/*'))

    crosscorr_pairs_lst = []
    for i,j in enumerate(crosscorr_days_lst):
    	crosscorr_file = sorted(glob.glob(j+'/*'))
    	crosscorr_pairs_lst.append(crosscorr_file)

    #Make a list of list flat
    crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

    HHE_HHE_lst = []
    HHN_HHN_lst = []
    HHZ_HHZ_lst = []

    HHE_HHN_lst = []
    HHE_HHZ_lst = []
    HHN_HHZ_lst = []

    for i in crosscorr_pairs:
        # splitting subdir/basename
        subdir, filename = os.path.split(i)
        name_pair1 = filename.split('.h5')[0].split('_')[-1]
        name_pair2 = filename.split('.h5')[0].split('_')[-2]

        if 'HHE' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HH2' in name_pair2 or 'HH2' in name_pair1 and 'HHE' in name_pair2  or 'HHE' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHE_lst.append(i)

        if 'HHN' in name_pair1 and 'HHN' in name_pair2 or 'HH1' in name_pair1 and 'HH1' in name_pair2 or 'HHN' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHN' in name_pair2:
            HHN_HHN_lst.append(i)

        if 'HHZ' in name_pair1 and 'HHZ' in name_pair2:
            HHZ_HHZ_lst.append(i)

        if 'HHE' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HH2' in name_pair2 or 'HHE' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHN_lst.append(i)

        if 'HHE' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HHE' in name_pair2 or 'HH2' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HH2' in name_pair2:
            HHE_HHZ_lst.append(i)

        if 'HHZ' in name_pair1 and 'HHN' in name_pair2 or 'HHN' in name_pair1 and 'HHZ' in name_pair2 or 'HHZ' in name_pair1 and 'HH1' in name_pair2 or 'HH1' in name_pair1 and 'HHZ' in name_pair2:
            HHN_HHZ_lst.append(i)

    CHANNEL_fig_lst = [HHE_HHE_lst,HHN_HHN_lst,HHZ_HHZ_lst,HHE_HHN_lst,HHE_HHZ_lst,HHN_HHZ_lst]
    chan_lst = ['HHE-HHE','HHN-HHN','HHZ-HHZ','HHE-HHN','HHE-HHZ','HHN-HHZ']

    for iOBS in OBS_LST:

        for ipairs,crosscorr_pairs in enumerate(CHANNEL_fig_lst):
            #Creating the figure
            fig = plt.figure(figsize=(10, 30))
            fig.suptitle(iOBS+': CrossCorrelation '+chan_lst[ipairs],fontsize=20)

            gs = gridspec.GridSpec(5, 1,wspace=0.2, hspace=0.5)

            #-------------------------------------------

            map_loc = fig.add_subplot(gs[:2],projection=ccrs.PlateCarree())

            LLCRNRLON_LARGE = -52
            URCRNRLON_LARGE = -38
            LLCRNRLAT_LARGE = -30
            URCRNRLAT_LARGE = -12

            map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
            map_loc.yaxis.set_ticks_position('both')
            map_loc.xaxis.set_ticks_position('both')

            map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE,3), crs=ccrs.PlateCarree())
            map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE,3), crs=ccrs.PlateCarree())
            map_loc.tick_params(labelbottom=True,labeltop=True,labelleft=True,labelright=True, labelsize=15)

            map_loc.grid(True,which='major',color='gray',linewidth=1,linestyle='--')

            reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
            shape_1_SHP = list(reader_1_SHP.geometries())
            plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
            map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)

            crosscorr_stack_style_lst = []
            crosscorr_stack_name_lst = []
            crosscorr_stack_data_normalized_lst = []
            crosscorr_stack_data_normalized_dist_lst = []
            crosscorr_stack_data_normalized_vmin_lst = []
            crosscorr_stack_data_normalized_vmax_lst = []
            time_to_plot = []

            for i in tqdm(crosscorr_pairs, desc='Reading '+chan_lst[ipairs]+' data'):

                #Reading data
                sta1_sta2_asdf_file = ASDFDataSet(i, mode='r')

                name_sta1 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked.list()[0]
                name_sta2 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked.list()[1]

                if iOBS in name_sta1 or iOBS in name_sta2:
                    crosscorr_stack_style_lst.append('k')
                    crosscorr_stack_name_lst.append(name_sta1+'-'+name_sta2)
                    dist_pair = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['dist']

                    loc_sta1 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['sta1_loc']
                    loc_sta2 = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['sta2_loc']

            	    #Stacked data
                    data_causal = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].data[::]
                    causal_time = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta1][name_sta2].parameters['crosscorr_stack_time']

                    data_acausal = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta2][name_sta1].data[::]
                    acausal_time = sta1_sta2_asdf_file.auxiliary_data.CrossCorrelationStacked[name_sta2][name_sta1].parameters['crosscorr_stack_time']

                    crosscorr_stack_data = data_acausal + data_causal
                    crosscorr_stack_time = acausal_time + causal_time

                    try:
                        OBS_SNR = SNR(crosscorr_stack_data,crosscorr_stack_time,dist_pair,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)
                        if OBS_SNR >= 1:

                            crosscorr_stack_data_normalized_lst.append(crosscorr_stack_data)
                            time_to_plot.append(crosscorr_stack_time)

                            crosscorr_stack_data_normalized_dist_lst.append(dist_pair)
                            crosscorr_stack_data_normalized_vmin_lst.append(dist_pair/SIGNAL_WINDOW_VMIN)
                            crosscorr_stack_data_normalized_vmax_lst.append(dist_pair/SIGNAL_WINDOW_VMAX)

                            # Use the cartopy interface to create a matplotlib transform object
                            # for the Geodetic coordinate system. We will use this along with
                            # matplotlib's offset_copy function to define a coordinate system which
                            # translates the text by 25 pixels to the left.
                            geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
                            text_transform = offset_copy(geodetic_transform, units='dots', y=0,x=80)
                            text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-15,x=15)

                            map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k',alpha=0.5, transform=ccrs.PlateCarree())
                            map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
                            map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

                    except:
                        pass


            #-------------------------------------------
            orglst = np.argsort(crosscorr_stack_data_normalized_dist_lst)
            crosscorr_stack_name_org_lst = [crosscorr_stack_name_lst[i] for i in orglst]
            crosscorr_stack_style_org_lst = [crosscorr_stack_style_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_dist_org_lst = [crosscorr_stack_data_normalized_dist_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_vmax_org_lst = [crosscorr_stack_data_normalized_vmax_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_vmin_org_lst = [crosscorr_stack_data_normalized_vmin_lst[i] for i in orglst]
            crosscorr_stack_data_normalized_org_lst = [crosscorr_stack_data_normalized_lst[i] for i in orglst]

            #--------------------------------------------------------------------------------------------------------------------
            #plot parameters
            sav_f = 11
            cmap_interpolation = 'RdBu'
            type_interpolation = 'sinc'
            #--------------------------------------------------------------------------------------------------------------------

            ax3 = fig.add_subplot(gs[2])

            vector_plot = np.array(crosscorr_stack_data_normalized_org_lst)

            extent = [-SHIFT_LEN,SHIFT_LEN,min(crosscorr_stack_data_normalized_dist_org_lst),max(crosscorr_stack_data_normalized_dist_org_lst)]
            im = ax3.imshow(vector_plot,extent=extent,origin='lower', interpolation=type_interpolation,aspect='auto',cmap=cmap_interpolation,vmax=1,vmin=-1)

            ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax3.set_xlabel('Lapse time (s)',fontsize=14)
            ax3.set_ylabel('Distance (km)',fontsize=14)
            ax3.yaxis.set_major_locator(MultipleLocator(100))
            ax3.yaxis.set_minor_locator(MultipleLocator(50))
            ax3.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax3.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            axins = inset_axes(ax3,
                               width="30%",  # width = 10% of parent_bbox width
                               height="2%",  # height : 5%
                               loc='upper left',
                               bbox_to_anchor=(0.65, 0.03, 1, 1),
                               bbox_transform=ax3.transAxes,
                               borderpad=0,
                               )
            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            #-----------------------------------------------------------------------------------------------
            crosscorr_stack_data_normalized_org_lst_20_50s = [bandpass(data_2_plot, 1.0/PERIOD_BANDS[1][1], 1.0/PERIOD_BANDS[1][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
            crosscorr_stack_data_normalized_org_lst_20_50 = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lst_20_50s]

            ax5 = fig.add_subplot(gs[3])

            vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_20_50)

            extent = [-SHIFT_LEN,SHIFT_LEN,min(crosscorr_stack_data_normalized_dist_org_lst),max(crosscorr_stack_data_normalized_dist_org_lst)]
            im = ax5.imshow(vector_plot,extent=extent,origin='lower', interpolation=type_interpolation,aspect='auto',cmap=cmap_interpolation,vmax=1,vmin=-1)

            ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax5.set_xlabel('Lapse time (s)',fontsize=14)
            ax5.set_ylabel('Distance (km)',fontsize=14)
            ax5.yaxis.set_major_locator(MultipleLocator(100))
            ax5.yaxis.set_minor_locator(MultipleLocator(50))
            ax5.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax5.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            axins = inset_axes(ax5,
                               width="30%",  # width = 10% of parent_bbox width
                               height="2%",  # height : 5%
                               loc='upper left',
                               bbox_to_anchor=(0.65, 0.03, 1, 1),
                               bbox_transform=ax5.transAxes,
                               borderpad=0,
                               )
            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            #-----------------------------------------------------------------------------------------------
            crosscorr_stack_data_normalized_org_lst_50_100s = [bandpass(data_2_plot, 1.0/PERIOD_BANDS[2][1], 1.0/PERIOD_BANDS[2][0], NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
            crosscorr_stack_data_normalized_org_lst_50_100 = [Normalize(a) for a in crosscorr_stack_data_normalized_org_lst_50_100s]

            ax7 = fig.add_subplot(gs[4])

            vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_50_100)

            extent = [-SHIFT_LEN,SHIFT_LEN,min(crosscorr_stack_data_normalized_dist_org_lst),max(crosscorr_stack_data_normalized_dist_org_lst)]
            im = ax7.imshow(vector_plot,extent=extent,origin='lower', interpolation=type_interpolation,aspect='auto',cmap=cmap_interpolation,vmax=1,vmin=-1)

            ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')
            ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],sav_f,1),savgol_filter(crosscorr_stack_data_normalized_dist_org_lst,sav_f,1),ls='--',lw=0.5,c='r')

            # adding labels
            ax7.set_xlabel('Lapse time (s)',fontsize=14)
            ax7.set_ylabel('Distance (km)',fontsize=14)
            ax7.yaxis.set_major_locator(MultipleLocator(100))
            ax7.yaxis.set_minor_locator(MultipleLocator(50))
            ax7.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
            ax7.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

            axins = inset_axes(ax7,
                               width="30%",  # width = 10% of parent_bbox width
                               height="2%",  # height : 5%
                               loc='upper left',
                               bbox_to_anchor=(0.65, 0.03, 1, 1),
                               bbox_transform=ax7.transAxes,
                               borderpad=0,
                               )
            plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            output_figure_CrossCorrelation_DAY = CLOCK_DRIFT_OUTPUT+'CROSS_CORR_STACK_INTERSTATION_DISTANCE_FIGURES/'
            os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
            fig.savefig(output_figure_CrossCorrelation_DAY+folder_name+'_INTERSTATION_DISTANCE_FIG_'+chan_lst[ipairs]+'_'+iOBS+'.png',dpi=300)
            plt.close()

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

# =======
# Classes
# =======

class CrossCorrelation:
    """
    Cross-correlation class, which contains:
    - a pair of sets of names
    - year and julian day
    - distance between stations
    - a cross-correlation data list
    """

    def __init__(self, name1, name2, lat1, lon1, lat2, lon2, pair_time_day):
        """
        @type name1: str
        @type name2: str
        @type lat1: float
        @type lon1: float
        @type lat2: float
        @type lon2: float
        @type pair_time_day: str
        """

        # names of stations
        self.name1 = name1
        self.name2 = name2

        # loc of stations
        self.lat1 = lat1
        self.lon1 = lon1
        self.lat2 = lat2
        self.lon2 = lon2

        # initializing stats
        self.croscorr_day = pair_time_day

    def __repr__(self):
        s = '<cross-correlation between stations {0}-{1}: date: {2}>'
        return s.format(self.name1, self.name2, self.croscorr_day)

    def __str__(self):
        """
        E.g., 'Cross-correlation between stations SPB - ITAB:
               365 days from 2002-01-01 to 2002-12-01'
        """
        s = ('Cross-correlation between stations '
             '{sta1}-{sta2}: '
             'Date: {crossday}')
        return s.format(sta1=self.name1,sta2=self.name2,crossday=self.croscorr_day)

    def dist(self):
        """
        Geodesic distance (in km) between stations, using the
        WGS-84 ellipsoidal model of the Earth - obspy.geodetics.base.gps2dist_azimuth
        """

        d,_,_ = gps2dist_azimuth(self.lat1, self.lon1, self.lat2, self.lon2, a=6378137.0, f=0.0033528106647474805)

        return np.array(d) / 1000.

    def add(self, tr1, tr2,sta1_hour_lst,sta2_hour_lst,xcorr=None,shift_len=SHIFT_LEN):
        """
        Stacks cross-correlation between 2 traces
        @type sta1_hour_lst: List of obspy.core.trace.Trace.hour
        @type sta2_hour_lst: List of obspy.core.trace.Trace.hour
        @type tr1: List of obspy.core.trace.Trace.data
        @type tr2: List of obspy.core.trace.Trace.data
        """

        # cross-correlation
        if xcorr is None:
            lst_day_hours = ['0:0', '0:30', '1:0', '1:30', '2:0', '2:30', '3:0', '3:30', '4:0', '4:30', '5:0', '5:30',
                             '6:0', '6:30', '7:0', '7:30', '8:0', '8:30', '9:0', '9:30', '10:0', '10:30', '11:0', '11:30',
                             '12:0', '12:30', '13:0', '13:30', '14:0', '14:30', '15:0', '15:30', '16:0', '16:30', '17:0', '17:30',
                             '18:0', '18:30', '19:0', '19:30', '20:0', '20:30', '21:0', '21:30', '22:0', '22:30', '23:0']

            # calculating cross-corr using obspy, if not already provided
            xcorr_hours = []
            for hour in lst_day_hours:
                if hour in sta1_hour_lst and hour in sta2_hour_lst:
                    xcorr_hours.append(obscorr(a=tr1[sta1_hour_lst.index(hour)], b=tr2[sta2_hour_lst.index(hour)], shift=int(round(shift_len*NEW_SAMPLING_RATE)), demean=True)[:2*shift_len*NEW_SAMPLING_RATE])
                else:
                    pass
            if len(xcorr_hours) > 0:
                xcorr = sum(xcorr_hours)/len(xcorr_hours)
                xcorr = xcorr / float((np.abs(xcorr).max()))
                xcorr_timearray = np.arange(0,shift_len,1/(2*NEW_SAMPLING_RATE))

                # normalizing cross-corr
                self.dataarray = xcorr

                # time arrya cross-corr
                self.timearray = xcorr_timearray

            else:
                self.dataarray = [0]
                self.timearray = [0]

# ============
# Main program
# ============
'''

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')
start_time = time.time()

files_STA = filelist(basedir=MSEED_DIR_STA,interval_period_date=INTERVAL_PERIOD_DATE)
files_OBS = filelist(basedir=MSEED_DIR_OBS,interval_period_date=INTERVAL_PERIOD_DATE)
files1 = files_STA+files_OBS

files_final_1 = []
for i in files1:
        files_final_1.append([i for sta in STATIONS_LST if sta in i])

files_final = [item for sublist in files_final_1 for item in sublist]

files = []
for s in files_final:
        if any(day_s in s for day_s in CHANNEL_LST):
                files.append(s)

# Total of files by station:
files_per_station = [[]]*len(STATIONS_LST)

for i,j in enumerate(STATIONS_LST):
    files_per_station[i] = len(list(filter(lambda x: j in x, files)))

print('Total of files by station:')
for i,j in enumerate(files_per_station):
    print(STATIONS_LST[i]+': '+str(j)+' files')

print('\n')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

print('==================================')
print('Opening miniseed files of each day')
print('==================================')
print('\n')

start_time = time.time()

with Pool(processes=num_processes) as p:
	max_ = len(files)
	with tqdm(total=max_) as pbar:
		for i, _ in enumerate(p.imap_unordered(get_stations_data, files)):
			pbar.update()
print('\n')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

print('====================================')
print('Calculating daily Cross-correlations:')
print('====================================')
print('\n')

days_crosscor = sorted(glob.glob(ASDF_FILES+'DATA_DAY_FILES/*'))

stationtrace_pairs_lst = []
for i,j in enumerate(days_crosscor):
	stations_file = sorted(glob.glob(j+'/*'))
	stationtrace_pairs_lst.append(list(combinations(stations_file, 2)))

stationtrace_pairs = [item for sublist in stationtrace_pairs_lst for item in sublist]

start_time = time.time()
pool = Pool(processes=num_processes)
CrossCorrelation_days_lst = []
for result in tqdm(pool.imap(func=crosscorr_func, iterable=stationtrace_pairs), total=len(stationtrace_pairs)):
	CrossCorrelation_days_lst.append(result)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('\n')
print('============================')
print('Stacking Cross-correlations:')
print('============================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(ASDF_FILES+'CROSS_CORR_DAY_FILES/*'))

crosscorr_pairs_lst = []
for i,j in enumerate(crosscorr_days_lst):
    crosscorr_file = sorted(glob.glob(j+'/*'))
    crosscorr_pairs_lst.append(crosscorr_file)

#Make a list of list flat
crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

#Separating according to pairs name
crosscorr_pairs_name_lst = []
for i in crosscorr_pairs:
    # splitting subdir/basename
    subdir, filename = os.path.split(i)
    crosscorr_pairs_name_lst.append(filename.split("_20")[0])

crosscorr_pairs_names = sorted(list(set(crosscorr_pairs_name_lst)))

crosscorr_pairs_data = [[]]*len(crosscorr_pairs_names)

for l,k in enumerate(crosscorr_pairs_names):
    crosscorr_pairs_data[l] = [j for i,j in enumerate(crosscorr_pairs) if k in j]

input_lst_crosscorr_pairs_names = []
for l,k in enumerate(crosscorr_pairs_data):
    input_lst_crosscorr_pairs_names.append([k,'CROSS_CORR_STACKED_FILES'])

#Stacking data

start_time = time.time()

pool = Pool(processes=num_processes)
CrossCorrelation_stations_lst = []
for result in tqdm(pool.imap(func=crosscorr_stack_asdf, iterable=input_lst_crosscorr_pairs_names), total=len(input_lst_crosscorr_pairs_names)):
    CrossCorrelation_stations_lst.append(result)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('\n')
print('============================================================')
print('Plotting Staked cross-correlations by interstation distance:')
print('============================================================')
print('\n')

start_time = time.time()
plot_stacked_cc_interstation_distance('CROSS_CORR_STACKED_FILES')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('\n')
print('=========================================')
print('10-day stacking daily cross-correlations:')
print('=========================================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(ASDF_FILES+'CROSS_CORR_DAY_FILES/*'))

crosscorr_pairs_lst = []
for i,j in enumerate(crosscorr_days_lst):
    crosscorr_file = sorted(glob.glob(j+'/*'))
    crosscorr_pairs_lst.append(crosscorr_file)

#Make a list of list flat
crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

#Separating according to pairs name
crosscorr_pairs_name_lst = []
for i in crosscorr_pairs:
    # splitting subdir/basename
    subdir, filename = os.path.split(i)
    crosscorr_pairs_name_lst.append(filename.split("_20")[0])

crosscorr_pairs_names = sorted(list(set(crosscorr_pairs_name_lst)))

crosscorr_pairs_data = [[]]*len(crosscorr_pairs_names)
for l,k in enumerate(crosscorr_pairs_names):
    crosscorr_pairs_data[l] = [j for i,j in enumerate(crosscorr_pairs) if k in j]

# ----------------------------------------------------------------------------------------------------------------------------------------------
for input_crosscorr_pairs_data in tqdm(crosscorr_pairs_data,desc='10 days stacking c-c'):
    crosscorr_pair_date_10day = []
    crosscorr_pairs_10day_data = []
    for crosscorr_pairs_data1 in input_crosscorr_pairs_data:
        subdir, filename = os.path.split(crosscorr_pairs_data1)
        crosscorr_pair_date = datetime.datetime.strptime(filename.split('/')[-1].split('_')[-2]+'.'+filename.split('/')[-1].split('_')[-1].split('.')[0], '%Y.%j')
        crosscorr_pair_date_10day.append(filename.split('/')[-1].split('_')[-2]+'.'+filename.split('/')[-1].split('_')[-1].split('.')[0])
        crosscorr_pairs_10day_data.append([file for file in input_crosscorr_pairs_data if datetime.datetime.strptime(file.split('/')[-1].split('_')[-2]+'.'+file.split('/')[-1].split('_')[-1].split('.')[0], '%Y.%j') >= crosscorr_pair_date-datetime.timedelta(days=5) and datetime.datetime.strptime(file.split('/')[-1].split('_')[-2]+'.'+file.split('/')[-1].split('_')[-1].split('.')[0], '%Y.%j') < crosscorr_pair_date+datetime.timedelta(days=5)])

    name_suffix_lst = ['CROSS_CORR_10_DAYS_FILES']*len(crosscorr_pair_date_10day)
    input_lst_crosscorr_pairs_names = [[crosscorr_pairs_10day_data[i],crosscorr_pair_date_10day[i],name_suffix_lst[i]] for i in range(len(crosscorr_pair_date_10day))]

    #Stacking data:
    pool = Pool(processes=num_processes)
    pool.imap(crosscorr_10_days_stack_func,input_lst_crosscorr_pairs_names)
    pool.close()
    pool.join()

print('\n')
print('============================')
print('Stacking Cross-correlations:')
print('============================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(ASDF_FILES+'CROSS_CORR_10_DAYS_FILES/*'))

crosscorr_pairs_lst = []
for i,j in enumerate(crosscorr_days_lst):
	crosscorr_file = sorted(glob.glob(j+'/*'))
	crosscorr_pairs_lst.append(crosscorr_file)

#Make a list of list flat
crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

#Separating according to pairs name
crosscorr_pairs_name_lst = []
for i in crosscorr_pairs:
	# splitting subdir/basename
	subdir, filename = os.path.split(i)
	crosscorr_pairs_name_lst.append(filename.split("_20")[0])

crosscorr_pairs_names = sorted(list(set(crosscorr_pairs_name_lst)))

crosscorr_pairs_data = [[]]*len(crosscorr_pairs_names)

for l,k in enumerate(crosscorr_pairs_names):
	crosscorr_pairs_data[l] = [j for i,j in enumerate(crosscorr_pairs) if k in j]

input_lst_crosscorr_pairs_names = []
for l,k in enumerate(crosscorr_pairs_data):
    input_lst_crosscorr_pairs_names.append([k,'CROSS_CORR_10_DAYS_STACKED_FILES'])

#Stacking data
start_time = time.time()

pool = Pool(processes=num_processes)
CrossCorrelation_stations_lst = []
for result in tqdm(pool.imap(func=crosscorr_stack_asdf, iterable=input_lst_crosscorr_pairs_names), total=len(input_lst_crosscorr_pairs_names)):
	CrossCorrelation_stations_lst.append(result)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

print('\n')
print('============================================================')
print('Plotting Staked cross-correlations by interstation distance:')
print('============================================================')
print('\n')

start_time = time.time()
plot_stacked_cc_interstation_distance('CROSS_CORR_10_DAYS_STACKED_FILES')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

#Plotting by OBS
start_time = time.time()
plot_stacked_cc_interstation_distance_per_obs_short('CROSS_CORR_10_DAYS_STACKED_FILES')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
'''
print('\n')
print('========================')
print('Clock Drift Calculating:')
print('========================')
print('\n')

start_time = time.time()
with Pool(processes=num_processes) as p:
    max_ = len(OBS_LST)
    with tqdm(total=max_,desc='Processing') as pbar:
        for i, _ in enumerate(p.imap_unordered(Calculating_clock_drift_func, OBS_LST)):
            pbar.update()
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

'''
print('\n')
print('===============================')
print('Total Clock Drift for each OBS:')
print('===============================')
print('\n')


        data_drift_OBS_dic = {
                             'name_sta1': pair_sta_1,
                             'name_sta2': pair_sta_2,
                             'dist_pair': dist_pair,
                             'loc_sta1': loc_sta1,
                             'loc_sta2': loc_sta2,
                             'date_to_plot':date_to_plot_clock_chan,
                             'coefficient_clock_drift':data_to_plot_coefficient_clock_drift_chan,
                             'shift_clock_drift':data_to_plot_shift_clock_drift_chan,
                             'chan_lst':chan_lst
                             }


clock_drift_files_lst = sorted(glob.glob(FEATHER_FILES+'/*'))

for iOBS in OBS_LST:
    clock_drift_files = [j for i,j in enumerate(clock_drift_files_lst) if iOBS in j]
    # ----------------------------------------------------------------------------------------------------

    clock_drift_df_lst = [pd.read_pickle(i) for i in clock_drift_files]
    clock_drift_df= pd.DataFrame.from_dict(clock_drift_df_lst, )


    for u,i in enumerate(clock_drift_df['coefficient_clock_drift'].values):
        print(len(i))
        for k,l in enumerate(i):
            print(clock_drift_df['chan_lst'][k][u],':',len(l))
    #print(len(clock_drift_df[['coefficient_clock_drift']].values[1]))

    print('=====')
    print('=====')
    print('=====')


    # ----------------------------------------------------------------------------------------------------
    #Creating the figure and plotting Clock-drift
    chan_lst = ['HHE-HHE','HHE-HHN','HHE-HHZ','HHN-HHN','HHN-HHE','HHN-HHZ','HHZ-HHE','HHZ-HHN','HHZ-HHZ']

    clock_drift_files_date_to_plot_ch = []
    clock_drift_files_data_to_plot_coefficient_ch = []
    clock_drift_files_data_to_plot_shift_ch = []

    for i,j in enumerate(clock_drift_files_name_sta2):
        if 'OBS' not in clock_drift_files_name_sta1[i] and clock_drift_files_name_sta2[i]:
            clock_drift_files_date_to_plot_ch.append([clock_drift_files_date_to_plot[i][k] for k,l in enumerate(chan_lst)])
            clock_drift_files_data_to_plot_coefficient_ch.append([clock_drift_files_data_to_plot_coefficient[i][k] for k,l in enumerate(chan_lst)])
            clock_drift_files_data_to_plot_shift_ch.append([clock_drift_files_data_to_plot_shift[i][k] for k,l in enumerate(chan_lst)])

    clock_drift_files_date_to_plot_ch1 = []
    clock_drift_files_data_to_plot_coefficient_ch1 = []
    clock_drift_files_data_to_plot_shift_ch1 = []
    for k,l in enumerate(clock_drift_files_date_to_plot_ch):
        clock_drift_files_date_to_plot_ch1 = []
        clock_drift_files_data_to_plot_coefficient_ch1 = []
        clock_drift_files_data_to_plot_shift_ch1 = []
        print(len(l[0]))

    print('=====')
    print('=====')
    print('=====')



    if OBS is not in clock_drift_files_name_sta1
    OBS_coefficient_clock_drift_ch = [[]]*len(chan_lst)
    OBS_shift_clock_drift_ch = [[]]*len(chan_lst)
    for k,l in enumerate(chan_lst):
        OBS_coefficient_clock_drift_ch[k].append(clock_drift_files_data_to_plot_coefficient[i][k])
        OBS_shift_clock_drift_ch[k].append(clock_drift_files_data_to_plot_shift[i][k])

    for w,e in enumerate(OBS_coefficient_clock_drift_ch):
        print(chan_lst[w],': ',len(OBS_coefficient_clock_drift_ch[w]))



    OBS_coefficient_clock_drift_ch1 = [[]]*len(chan_lst)
    OBS_shift_clock_drift_ch1 = [[]]*len(chan_lst)
    for k,l in enumerate(chan_lst):
        OBS_coefficient_clock_drift_ch1.append([item for sublist in OBS_coefficient_clock_drift_ch[k] for item in sublist])
        OBS_shift_clock_drift_ch1.append([item for sublist in OBS_shift_clock_drift_ch[k] for item in sublist])

    for w,e in enumerate(OBS_shift_clock_drift_ch1):
            print(chan_lst[w],': ',e)


        # --------------------------------------------
        # Creating the figure and plotting Clock-drift
        # --------------------------------------------

        fig = plt.figure(figsize=(20, 15))
        fig.suptitle('Clock-drift: '+pair_sta_1+'-'+pair_sta_2+'('+str(round(dist_pair))+' km)',fontsize=20)
        fig.autofmt_xdate()
        # ----------------------------------------------------------------------------------------------------

        gs = gridspec.GridSpec(9, 2,wspace=0.5, hspace=0.8)
        map_loc = fig.add_subplot(gs[:,0],projection=ccrs.PlateCarree())

        LLCRNRLON_LARGE = -52
        URCRNRLON_LARGE = -38
        LLCRNRLAT_LARGE = -30
        URCRNRLAT_LARGE = -12

        map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
        map_loc.yaxis.set_ticks_position('both')
        map_loc.xaxis.set_ticks_position('both')

        map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+3,3), crs=ccrs.PlateCarree())
        map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+3,3), crs=ccrs.PlateCarree())
        map_loc.tick_params(labelbottom=True, labeltop=True, labelleft=True, labelright=True, labelsize=12)
        map_loc.grid(True,which='major',color='gray',linewidth=0.5,linestyle='--')

        reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
        shape_1_SHP = list(reader_1_SHP.geometries())
        plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
        map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)
        # Use the cartopy interface to create a matplotlib transform object
        # for the Geodetic coordinate system. We will use this along with
        # matplotlib's offset_copy function to define a coordinate system which
        # translates the text by 25 pixels to the left.
        geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
        text_transform = offset_copy(geodetic_transform, units='dots', y=50,x=100)

        map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k',alpha=0.5,transform=ccrs.PlateCarree())
        map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
        map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

        map_loc.text(loc_sta1[1],loc_sta1[0], pair_sta_1,fontsize=12,verticalalignment='center', horizontalalignment='right',transform=text_transform)
        map_loc.text(loc_sta2[1],loc_sta2[0], pair_sta_2,fontsize=12,verticalalignment='center', horizontalalignment='right',transform=text_transform)

        # ----------------------------------------------------------------------------------------------------
        days_major = DayLocator(interval=5)   # every 5 day
        days_minor = DayLocator(interval=1)   # every day
        months = MonthLocator(interval=3)  # every month
        yearsFmt = DateFormatter('%b-%Y')

        for z,x in enumerate(data_to_plot_coefficient_clock_drift_chan):
            if len(x) > 1:
                # ----------------------------------------------------------------------------------------------------
                ax0 = fig.add_subplot(gs[z,1])
                ax0.xaxis.set_major_locator(months)
                ax0.xaxis.set_major_formatter(yearsFmt)
                ax0.xaxis.set_minor_locator(days_minor)
                ax0.yaxis.set_major_locator(MultipleLocator(100))
                ax0.yaxis.set_minor_locator(MultipleLocator(25))
                ax0.set_ylabel('Drift ('+str(1/NEW_SAMPLING_RATE)+'s)')
                ax0.set_title(chan_lst[z])
                ax0.set_ylim(-200,200)

                # -------------------------------------------------------------------------------------------------------------
                poly_reg = PolynomialFeatures(degree=4)
                X_poly = poly_reg.fit_transform(np.array(range(len(data_to_plot_shift_clock_drift_chan[z]))).reshape(-1, 1))
                pol_reg = LinearRegression()
                pol_reg.fit(X_poly, data_to_plot_shift_clock_drift_chan[z])
                # -------------------------------------------------------------------------------------------------------------
                for y,u in enumerate(data_to_plot_shift_clock_drift_chan[z]):
                    if data_to_plot_coefficient_clock_drift_chan[z][y] > 0.3:
                        im = ax0.scatter(date_to_plot_clock_chan[z][y],data_to_plot_shift_clock_drift_chan[z][y],c=data_to_plot_coefficient_clock_drift_chan[z][y],marker='o',edgecolors=None,cmap='magma',s=20,vmin=0,vmax=1,alpha=0.9)
                    else:
                        im = ax0.scatter(date_to_plot_clock_chan[z][y],data_to_plot_shift_clock_drift_chan[z][y],c=data_to_plot_coefficient_clock_drift_chan[z][y],marker='o',edgecolors=None,cmap='magma',s=10,vmin=0,vmax=1,alpha=0.7)

                ax0.plot(date_to_plot_clock_chan[z], pol_reg.predict(poly_reg.fit_transform(np.array(range(len(data_to_plot_shift_clock_drift_chan[z]))).reshape(-1, 1))),'--b')

                if z == 0:
                    axins = inset_axes(ax0,
                                           width="30%",  # width = 10% of parent_bbox width
                                           height="10%",  # height : 5%
                                           loc='upper left',
                                           bbox_to_anchor=(0.65,0.1, 1, 1),
                                           bbox_transform=ax0.transAxes,
                                           borderpad=0,
                                           )
                    plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

            else:

                ax0 = fig.add_subplot(gs[z,1])
                ax0.xaxis.set_major_locator(months)
                ax0.xaxis.set_major_formatter(yearsFmt)
                ax0.xaxis.set_minor_locator(days_minor)
                ax0.yaxis.set_major_locator(MultipleLocator(100))
                ax0.yaxis.set_minor_locator(MultipleLocator(25))
                ax0.set_ylabel('Drift ('+str(1/NEW_SAMPLING_RATE)+'s)')
                ax0.set_title(chan_lst[z])
                ax0.set_ylim(-200,200)

        # -------------------------------------------------------------------------------------------------------------

        output_figure_CLOCK_DRIFT = CLOCK_DRIFT_OUTPUT+'CLOCK_DRIFT_TOTAL_FIGURES/'
        os.makedirs(output_figure_CLOCK_DRIFT,exist_ok=True)
        fig.savefig(output_figure_CLOCK_DRIFT+'CLOCK_DRIFT_TOTAL_'+OBS_LST[i]+'.png',dpi=300)
        plt.close()
'''
