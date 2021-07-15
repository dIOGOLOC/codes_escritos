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

import json
import glob
import os
import numpy as np
from numpy.fft import rfft, irfft, fft, ifft, fftfreq
from itertools import combinations
from numpy.lib.stride_tricks import as_strided
import pandas as pd
from scipy.signal import spectrogram, detrend, resample,savgol_filter,decimate
from scipy.linalg import norm

import random
import collections
from copy import copy
import datetime
from itertools import compress

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

import pyasdf


#Configuration file

MSEED_DIR_OBS = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/obs_data_MSEED/'

MSEED_DIR_STA = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/data/'

STATIONS_LST = ['ABR01','DUB01','MAN01','OBS20','TER01','ALF01','GDU01','NAN01','TIJ01','CAJ01','GUA01',
'OBS17','PET01','TRI01','CAM01','JAC01','OBS18','RIB01','VAS01','CMC01','MAJ01','SLP01']

CHANNEL_LST = ['HHZ.D','HHN.D','HHE.D','HH1.D','HH2.D']

STATIONXML_DIR = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'
#STATIONXML_DIR = '/run/user/1000/gvfs/smb-share:server=hatabackup.local,share=dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

ORIENTATION_OUTPUT = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/ORIENTATION_OUTPUT/FIGURAS/'
#ORIENTATION_OUTPUT = '/run/user/1000/gvfs/smb-share:server=hatabackup.local,share=dados_posdoc/ON_MAR/ORIENTATION_OUTPUT/FIGURAS/'

JSON_FILES = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/ORIENTATION_OUTPUT/JSON_FILES/'
#JSON_FILES = '/run/user/1000/gvfs/smb-share:server=hatabackup.local,share=dados_posdoc/ON_MAR/ORIENTATION_OUTPUT/JSON_FILES/'

#Shapefile  boundary states
BOUNDARY_STATES_SHP = '/media/diogoloc/Backup/dados_posdoc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/UFEBRASIL.shp'
#BOUNDARY_STATES_SHP = '/run/user/1000/gvfs/smb-share:server=hatabackup.local,share=dados_posdoc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/UFEBRASIL.shp'

FIRSTDAY = '2019-08-01'
LASTDAY = '2020-06-01'

#Each hour-long seismogram is amplitude clipped at twice its standard deviation of that hour-long time window.
CLIP_FACTOR = 2

MIN_WINDOWS = 30

WINDOW_LENGTH = 3600

#max time window (s) for cross-correlation
SHIFT_LEN = 1800

PERIOD_BANDS = [[2, 5], [7, 25], [20, 50], [50, 100]]
# (these bands focus on periods ~7, 15, 25 seconds)

FREQUENCY_BANDS = [[0.5, 1], [1, 2], [2, 3], [3, 4]]

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

# ========================
# Constants and parameters
# ========================

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# ================
# MULTIPROCESSING
# ================

num_processes = 4

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

    file_lsts = sorted([i for i in files_lst if any(channel_s in i for channel_s in CHANNEL_LST)])

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
def SNR(data,time_data,dist,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE):
	"""
    Signal-to-noise ratio calculated as the peak of the absolute amplitude in the signal window divided by the standard deviation in the noise window.

    The signal window is defined by *vmin* and *vmax*:
    	dist/*vmax* < t < dist/*vmin*

    The noise window starts *signal2noise_trail* after the
    signal window and has a size of *noise_window_size*:

    	t > dist/*vmin* + *signal2noise_trail*
    	t < dist/*vmin* + *signal2noise_trail* + *noise_window_size*

    @type data: numpy array
    @type time_data: numpy array
    @type vmin: float
    @type vmax: float
    @type signal2noise_trail: float
    @type noise_window_size: float
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

	normalized_data = [2*(i-min(data)/max(data)-min(data))-1 for i in data]

	return normalized_data

#-------------------------------------------------------------------------------

def get_stations_data(f,Amp_clip=True,onebit_norm=False,white_spectral=True):
    """
    Gets stations daily data from miniseed file

    @type f: paht of the minissed file (str)
    @rtype: list of L{StationDayData}
    """

    # splitting subdir/basename
    subdir, filename = os.path.split(f)

    # network, station name and station channel in basename,
    # e.g., ON.TIJ01..HHZ.D.2020.002

    network, name = filename.split('.')[0:2]
    sta_channel_id = filename.split('.D.')[0]
    time_day = filename.split('.D.')[-1]
    year_day = time_day.split('.')[0]
    julday_day = time_day.split('.')[1]

    output_DATA_DAY = JSON_FILES+'DATA_DAY_FILES/'+year_day+'.'+julday_day+'/'

    if sta_channel_id == 'HH1':
        sta_channel_id = 'HHN'
    if sta_channel_id == 'HH2':
        sta_channel_id = 'HHE'

    if os.path.isfile(output_DATA_DAY+'DATA_DAY_'+sta_channel_id+'_'+year_day+'_'+julday_day+'.json'):
        pass

    else:

        try:
            st = read(f)

            if len(st[0].data) > 1:

        	    st_traces = [k for k in st.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH/2)]

        	    st_traces_check = []
        	    st_hours = []
        	    for k in st_traces:
        	    	if len(k[0].data) >= 3600:
        	    		st_traces_check.append(k)
        	    		st_hours.append(str(k[0].stats.starttime.hour)+':'+str(k[0].stats.starttime.minute))

        	    if len(st_hours) > MIN_WINDOWS:

        		    inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))

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

        		    traces_data_day = [ton[0].data.tolist() for ton in traces_white_spectral]


        		    lon = inv.get_coordinates(sta_channel_id)['longitude']
        		    lat = inv.get_coordinates(sta_channel_id)['latitude']

        		    # appending new channel day data
        		    station = StationDayData(name=name, network=network,fileID=sta_channel_id,lon=lon,lat=lat,time_day=time_day,hour_lst=st_hours,data_day=traces_data_day)

        		    data_dic = {
        						'name': station.name,
        						'network': station.network,
        						'fileID': station.fileID,
        						'lon': station.lon,
        						'lat': station.lat,
        						'time_day': station.time_day,
        						'hours_day': station.hour_lst,
        						'data_day': station.data_day
        			    		}

        		    os.makedirs(output_DATA_DAY,exist_ok=True)
        		    with open(output_DATA_DAY+'DATA_DAY_'+sta_channel_id+'_'+year_day+'_'+julday_day+'.json', 'w') as fp:
        		    	json.dump(data_dic, fp)
            else:
                pass
        except:
            print('Problem with the file: '+f)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Creating Dictionaries to allocate results ###
def nested_dict():
    return collections.defaultdict(nested_dict)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

def crosscorr_func(stationtrace_pairs):

    CrossCorrelation_dic = nested_dict()
    sta1 = json.load(open(stationtrace_pairs[0]))
    sta2 = json.load(open(stationtrace_pairs[1]))

    year_day = sta1['time_day'].split('.')[0]
    julday_day = sta1['time_day'].split('.')[1]

    sta1_CHANNEL = sta1['fileID'].split('..')[1]
    sta2_CHANNEL = sta2['fileID'].split('..')[1]

    day_crosscor_causal = CrossCorrelation(name1=sta1['name'],name2=sta2['name'],channel_1_name=sta1_CHANNEL,channel_2_name=sta2_CHANNEL,lat1=sta1['lat'],lon1=sta1['lon'],lat2=sta2['lat'],lon2=sta2['lon'],pair_time_day=sta1['time_day'])
    day_crosscor_acausal = CrossCorrelation(name1=sta2['name'],name2=sta1['name'],channel_1_name=sta1_CHANNEL,channel_2_name=sta2_CHANNEL,lat1=sta2['lat'],lon1=sta2['lon'],lat2=sta1['lat'],lon2=sta1['lon'],pair_time_day=sta1['time_day'])

    day_crosscor_causal.add(sta1['data_day'],sta2['data_day'],sta1['hours_day'],sta2['hours_day'])
    day_crosscor_acausal.add(sta2['data_day'],sta1['data_day'],sta2['hours_day'],sta1['hours_day'])

    if len(day_crosscor_acausal.dataarray) > 1:

        day_time_crosscor_all = day_crosscor_acausal.timearray+day_crosscor_causal.timearray
        day_data_crosscor_all = day_crosscor_acausal.dataarray+day_crosscor_causal.dataarray
        raw_SNR = SNR(day_data_crosscor_all,day_time_crosscor_all,day_crosscor_causal.dist(),vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)

        if raw_SNR > minspectSNR:

            CrossCorrelation_dic['dist'] = round(day_crosscor_causal.dist())
            CrossCorrelation_dic['date'] = sta1['time_day']
            CrossCorrelation_dic['sta1_loc'] = [sta1['lat'],sta1['lon']]
            CrossCorrelation_dic['sta1_name'] = sta1['name']
            CrossCorrelation_dic['sta2_loc'] = [sta2['lat'],sta2['lon']]
            CrossCorrelation_dic['sta2_name'] = sta2['name']
            CrossCorrelation_dic['sta1_channel'] = sta1_CHANNEL
            CrossCorrelation_dic['sta2_channel'] = sta2_CHANNEL
            CrossCorrelation_dic['crosscorr_daily_causal_time'] = day_crosscor_causal.timearray.tolist()
            CrossCorrelation_dic['crosscorr_daily_causal'] = day_crosscor_causal.dataarray.tolist()
            CrossCorrelation_dic['crosscorr_daily_acausal'] = day_crosscor_acausal.dataarray.tolist()
            CrossCorrelation_dic['crosscorr_daily_acausal_time'] = day_crosscor_acausal.timearray.tolist()

            output_CrossCorrelation_DAY = JSON_FILES+'CROSS_CORR_ORIENTATION_DAY_FILES/'+year_day+'.'+julday_day+'/'
            os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)
            with open(output_CrossCorrelation_DAY+'CROSS_CORR_ORIENTATION_DAY_FILES_'+sta1['name']+'_'+sta1_CHANNEL+'_'+sta2['name']+'_'+sta2_CHANNEL+'_'+year_day+'_'+julday_day+'.json', 'w') as fp:
                json.dump(CrossCorrelation_dic, fp)

		    # ============================
		    # Plot: map and pair crosscorr
		    # ============================

            fig = plt.figure(figsize=(15, 15))
            fig.suptitle(sta1['name']+'.'+sta1_CHANNEL+'-'+sta2['name']+'.'+sta2_CHANNEL+' - Day - '+UTCDateTime(year=int(year_day),julday=int(julday_day)).strftime('%d/%m/%Y'),fontsize=20)

            gs = gridspec.GridSpec(2, 1,wspace=0.2, hspace=0.5)

            #-------------------------------------------

            map_loc = fig.add_subplot(gs[0],projection=ccrs.PlateCarree())

            LLCRNRLON_LARGE = -52
            URCRNRLON_LARGE = -28
            LLCRNRLAT_LARGE = -30
            URCRNRLAT_LARGE = -10

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

            map_loc.scatter(sta1['lon'],sta1['lat'], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
            map_loc.scatter(sta2['lon'],sta2['lat'], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
            map_loc.plot([sta1['lon'],sta2['lon']],[sta1['lat'],sta2['lat']],c='k', transform=ccrs.PlateCarree())

            map_loc.text(sta1['lon'],sta1['lat'], sta1['name'],fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
            map_loc.text(sta2['lon'],sta2['lat'], sta2['name'],fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
			#-------------------------------------------

            ax = fig.add_subplot(gs[1])
            data_to_plot = CrossCorrelation_dic['crosscorr_daily_acausal'][::-1]+CrossCorrelation_dic['crosscorr_daily_causal']
            time_to_plot = [-1*i for i in CrossCorrelation_dic['crosscorr_daily_acausal_time'][::-1]] + CrossCorrelation_dic['crosscorr_daily_causal_time']
            ax.plot(time_to_plot,data_to_plot,color='k')
            ax.set_xlabel('time (s)',fontsize=14)
            ax.set_title('Dist = '+str(round(day_crosscor_causal.dist()))+' km',fontsize=14)

            output_figure_CrossCorrelation_DAY = ORIENTATION_OUTPUT+'CROSS_CORR_DAY_FIGURES/'+year_day+'.'+julday_day+'/'
            os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
            fig.savefig(output_figure_CrossCorrelation_DAY+'CROSS_CORR_DAY_FIG_'+sta1['name']+'_'+sta1_CHANNEL+'_'+sta2['name']+'_'+sta2_CHANNEL+'_'+year_day+'_'+julday_day+'.png')
            plt.close('all')

        return sta1['time_day']

    else:
        print("Problem: CrossCorrelation between "+sta1['name']+'.'+sta1_CHANNEL+'-'+sta2['name']+'.'+sta2_CHANNEL+" in "+sta2['time_day'])

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Calculating signal-to-noise ratio
def crosscorr_stack(crosscorr_pairs_data):
	"""
    Stacking crosscorrelation data

    @type crosscorr_pairs_data: list with JSON files
    """
    #Reading data
	name_sta1 = list(set([json.load(open(a))['sta1_name'] for a in crosscorr_pairs_data]))[0]
	name_sta2 = list(set([json.load(open(a))['sta2_name'] for a in crosscorr_pairs_data]))[0]
	channel_sta1 = list(set([json.load(open(a))['sta1_channel'] for a in crosscorr_pairs_data]))[0]
	channel_sta2 = list(set([json.load(open(a))['sta2_channel'] for a in crosscorr_pairs_data]))[0]
	dist_pair = json.load(open(crosscorr_pairs_data[0]))['dist']

	loc_sta1 = json.load(open(crosscorr_pairs_data[0]))['sta1_loc']
	loc_sta2 = json.load(open(crosscorr_pairs_data[0]))['sta2_loc']

    #Stacking data
	causal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_causal'] for a in crosscorr_pairs_data]),axis=0)
	acausal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_acausal'] for a in crosscorr_pairs_data]),axis=0)

	causal_time = np.array(json.load(open(crosscorr_pairs_data[0]))['crosscorr_daily_causal_time'])
	acausal_time = np.array(json.load(open(crosscorr_pairs_data[0]))['crosscorr_daily_acausal_time'])

	data_to_plot = acausal_lst[::-1] + causal_lst
	time_to_plot = [-1*i for i in acausal_time[::-1]] + causal_time

	#Saving CrossCorrelations stacked data:
	CrossCorrelation_stack_dic = nested_dict()

	CrossCorrelation_stack_dic['dist'] = dist_pair
	CrossCorrelation_stack_dic['sta1_name'] = name_sta1
	CrossCorrelation_stack_dic['sta2_name'] = name_sta2
	CrossCorrelation_stack_dic['sta1_loc'] = loc_sta1
	CrossCorrelation_stack_dic['sta2_loc'] = loc_sta2
	CrossCorrelation_stack_dic['sta1_channel'] = channel_sta1
	CrossCorrelation_stack_dic['sta2_channel'] = channel_sta2
	CrossCorrelation_stack_dic['crosscorr_stack'] = data_to_plot.tolist()
	CrossCorrelation_stack_dic['crosscorr_stack_time'] = time_to_plot.tolist()

	output_CrossCorrelation_DAY = JSON_FILES+'CROSS_CORR_STACKED_FILES/'+name_sta1+'_'+channel_sta1+'.'+name_sta2+'_'+channel_sta2+'/'
	os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)
	with open(output_CrossCorrelation_DAY+'CROSS_CORR_STACKED_FILES_'+name_sta1+'_'+channel_sta1+'.'+name_sta2+'_'+channel_sta2+'.json', 'w') as fp:
		json.dump(CrossCorrelation_stack_dic, fp)


	return CrossCorrelation_stack_dic['crosscorr_stack']

# =======
# Classes
# =======

class StationDayData:
    """
    Class to save station info: name, network, fileID, files directory, coordinates, datetime and data.
    """

    def __init__(self, name, network, fileID,lon=None,lat=None,time_day=None,hour_lst=None,data_day=None):
        """
        @type name: str
        @type network: str
        @type channel: str
        @type filesdir: str or unicode
        @type lon: float
        @type lat: float
        @type time_day: year + julian day
        @type time_day: list with hour of obspy.core.trace.Trace
        @type data_day: list of obspy.core.trace.Trace
        """
        self.name = name
        self.network = network
        self.fileID = fileID
        self.lon = lon
        self.lat = lat
        self.hour_lst = hour_lst
        self.time_day = time_day
        self.data_day = data_day

    def __str__(self):
        """
        @rtype: unicode
        """
        # General infos of station
        s = [u'Name    : {0}'.format(self.name),
             u'Network : {0}'.format(self.network),
             u'FilesID: {0}'.format(self.fileID),
             u'Longitude:{0}'.format(self.lon),
             u'Latitude: {0}'.format(self.lat),
             u'Time_day: {0}'.format(self.time_day),
            ]
        return u'\n'.join(s)

#-------------------------------------------------------------------------------


class CrossCorrelation:
    """
    Cross-correlation class, which contains:
    - a pair of sets of names and channels
    - year and julian day
    - distance between stations
    - a cross-correlation data list
    """

    def __init__(self, name1, name2, channel_1_name, channel_2_name, lat1, lon1, lat2, lon2, pair_time_day):
        """
        @type name1: str
        @type name2: str
        @type channel_1_name: str
        @type channel_2_name: str
        @type lat1: float
        @type lon1: float
        @type lat2: float
        @type lon2: float
        @type pair_time_day: str
        """

        # names of the stations
        self.name1 = name1
        self.name2 = name2

        # channel names of the stations
        self.channel_1_name = channel_1_name
        self.channel_2_name = channel_2_name

        # loc of stations
        self.lat1 = lat1
        self.lon1 = lon1
        self.lat2 = lat2
        self.lon2 = lon2

        # initializing stats
        self.croscorr_day = pair_time_day

    def __repr__(self):
        s = '<cross-correlation between stations {sta1}.{ch_sta1}-{sta2}.{ch_sta2}: date: {crossday}>'
        return s.format(sta1=self.name1, ch_sta1=self.channel_1_name, sta2=self.name2, ch_sta2=self.channel_2_name, crossday=self.croscorr_day)

    def __str__(self):
        """
        E.g., 'Cross-correlation between stations SPB - ITAB:
               365 days from 2002-01-01 to 2002-12-01'
        """
        s = ('Cross-correlation between stations '
             '{sta1}.{ch_sta1}-{sta2}.{ch_sta2}: '
             'Date: {crossday}')
        return s.format(sta1=self.name1, ch_sta1=self.channel_1_name, sta2=self.name2, ch_sta2=self.channel_2_name, crossday=self.croscorr_day)

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

files_final = []
for i in files1:
        files_final.append([i for sta in STATIONS_LST if sta in i])

files = [item for sublist in files_final for item in sublist]

print('Total of files = '+str(len(files)))
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

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

print('====================================')
print('Calculating daily Cross-correlations:')
print('====================================')
print('\n')

days_crosscor = sorted(glob.glob(JSON_FILES+'DATA_DAY_FILES/*'))
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
'''
print('\n')
print('====================================')
print('Stacking daily Cross-correlations:')
print('====================================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(JSON_FILES+'CROSS_CORR_ORIENTATION_DAY_FILES/*'))

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

#Stacking data
for i in crosscorr_pairs_data:
	name_sta1 = list(set([json.load(open(a))['sta1_name'] for a in i]))[0]
	channel_sta1 = list(set([json.load(open(a))['sta1_channel'] for a in i]))[0]
	name_sta2 = list(set([json.load(open(a))['sta2_name'] for a in i]))[0]
	channel_sta2 = list(set([json.load(open(a))['sta2_channel'] for a in i]))[0]
	dist_pair = json.load(open(i[0]))['dist']

	if dist_pair > 0:

		print('Pair: '+name_sta1+'.'+channel_sta1+'-'+name_sta2+'.'+channel_sta2+' - '+'days stacked: '+str(len(i)))

		causal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_causal'] for a in i]),axis=0)
		acausal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_acausal'] for a in i]),axis=0)

		causal_time = np.array(json.load(open(i[0]))['crosscorr_daily_causal_time'])
		acausal_time = np.array(json.load(open(i[0]))['crosscorr_daily_acausal_time'])

		loc_sta1 = json.load(open(i[0]))['sta1_loc']
		loc_sta2 = json.load(open(i[0]))['sta2_loc']

        # ===============================
        # Plot: map and stacked crosscorr
        # ===============================

		fig = plt.figure(figsize=(20, 7))
		fig.suptitle(name_sta1+'.'+channel_sta1+'-'+name_sta2+'.'+channel_sta2+' - Dist = '+str(round(dist_pair))+' km'+' - Days stacked: '+str(len(i)),fontsize=20)

		gs = gridspec.GridSpec(len(PERIOD_BANDS)+1, 2,wspace=0.2, hspace=0.5)

        #-------------------------------------------

		map_loc = fig.add_subplot(gs[:,0],projection=ccrs.PlateCarree())

		LLCRNRLON_LARGE = -52
		URCRNRLON_LARGE = -28
		LLCRNRLAT_LARGE = -30
		URCRNRLAT_LARGE = -10

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
		text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-15,x=15)

		map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
		map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

		map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]], transform=ccrs.PlateCarree())

		map_loc.text(loc_sta1[1],loc_sta1[0], name_sta1,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
		map_loc.text(loc_sta2[1],loc_sta2[0], name_sta2,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
    	#-------------------------------------------

		data_to_plot = acausal_lst[::-1] + causal_lst
		time_to_plot = [-1*i for i in acausal_time[::-1]] + causal_time

		ax = fig.add_subplot(gs[0,1])
		ax.plot(time_to_plot, data_to_plot,color='k')

		ylim = (data_to_plot.min(), data_to_plot.max())
		ax.set_ylim(ylim)
		ax.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		raw_SNR = SNR(data_to_plot,time_to_plot,dist_pair,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)
		ax.text(x=0.1,y=0.8,s="SNR:"+str(round(raw_SNR)),horizontalalignment='center',transform=ax.transAxes,fontsize=8,bbox={'edgecolor':'w','facecolor': 'white'})
		ax.text(x=0.9,y=0.8,s="Raw data",horizontalalignment='center',transform=ax.transAxes,fontsize=8,bbox={'edgecolor':'w','facecolor': 'white'})

    	# plotting band-filtered cross-correlation
		for i,band in enumerate(PERIOD_BANDS):
			ax = fig.add_subplot(gs[i+1,1])
			lastplot = band is PERIOD_BANDS[-1]
			firstplot = band is PERIOD_BANDS[0]

			ax.text(x=0.9,y=0.8,s="Filter:"+str(band[0])+"-"+str(band[1])+"s",transform=ax.transAxes,horizontalalignment='center',fontsize=8,bbox={'edgecolor':'w','facecolor': 'white'})

			data_to_plot_filtered = bandpass(data_to_plot, 1.0 /band[1], 1.0 / band[0], NEW_SAMPLING_RATE, corners=2, zerophase=False)
			ax.plot(time_to_plot, data_to_plot_filtered,color='k')
			ylim = (data_to_plot_filtered.min(), data_to_plot_filtered.max())
			ax.set_ylim(ylim)
			ax.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)
			ax.set_yticklabels([])

			ax.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
			ax.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
			ax.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
			ax.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)


    		# Calculating SNR
			filtered_SNR = SNR(data_to_plot_filtered,time_to_plot,dist_pair,vmin=SIGNAL_WINDOW_VMIN,vmax=SIGNAL_WINDOW_VMAX,signal2noise_trail=SIGNAL2NOISE_TRAIL,noise_window_size=NOISE_WINDOW_SIZE)
			ax.text(x=0.1,y=0.8,s="SNR:"+str(round(filtered_SNR)),horizontalalignment='center',transform=ax.transAxes,fontsize=8,bbox={'edgecolor':'w','facecolor': 'white'})

    		# formatting labels
			if not lastplot:
				ax.set_xticklabels([])
			if lastplot:
    		# adding label to time
				ax.set_xlabel('Time (s)',fontsize=14)

		output_figure_CrossCorrelation_DAY = ORIENTATION_OUTPUT+'CROSS_CORR_STACK_FIGURES/'
		os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
		fig.savefig(output_figure_CrossCorrelation_DAY+'CROSS_CORR_STACK_FIG_'+name_sta1+'_'+channel_sta1+'_'+name_sta2+'_'+channel_sta2+'.png')
		plt.close()
'''
print('\n')
print('====================================')
print('Stacking 10-day Cross-correlations:')
print('====================================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(JSON_FILES+'CROSS_CORR_ORIENTATION_DAY_FILES/*'))

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

crosscorr_pairs_data_10_day_all = []
for i in tqdm(crosscorr_pairs_data):
	crosscorr_pairs_10day_data = []
	crosscorr_pair_date_10day = []
	for k in i:
		subdir, filename = os.path.split(k)
		crosscorr_pair_date = datetime.datetime.strptime(filename.split('.')[0].split('_')[-2]+'.'+filename.split('.')[0].split('_')[-1], '%Y.%j')
		crosscorr_pair_date_10day.append(crosscorr_pair_date)
		crosscorr_pairs_10day_data.append([file for file in i if datetime.datetime.strptime(file.split('/')[-1].split('.')[0].split('_')[-2]+'.'+file.split('/')[-1].split('.')[0].split('_')[-1], '%Y.%j') >= crosscorr_pair_date and datetime.datetime.strptime(file.split('/')[-1].split('.')[0].split('_')[-2]+'.'+file.split('/')[-1].split('.')[0].split('_')[-1], '%Y.%j') < crosscorr_pair_date+datetime.timedelta(days=10)])

	year_day = filename.split('.')[0].split('_')[-2]
	julday_day = filename.split('.')[0].split('_')[-1]

	#Stacking data
	date_to_plot = []
	data_to_plot = []
	time_to_plot = []

	for ind,data10 in enumerate(crosscorr_pairs_10day_data):

		date_to_plot.append(crosscorr_pair_date_10day[ind])
		name_sta1 = json.load(open(data10[0]))['sta1_name']
		name_sta2 = json.load(open(data10[0]))['sta2_name']
		dist_pair = json.load(open(data10[0]))['dist']

		causal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_causal'] for a in data10]),axis=0)
		acausal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_acausal'] for a in data10]),axis=0)

		causal_time = np.array(json.load(open(data10[0]))['crosscorr_daily_causal_time'])
		acausal_time = np.array(json.load(open(data10[0]))['crosscorr_daily_acausal_time'])

		data_to_plot.append(acausal_lst[::-1] + causal_lst)
		time_to_plot.append([-1*i for i in acausal_time[::-1]] + causal_time)

		loc_sta1 = json.load(open(data10[0]))['sta1_loc']
		loc_sta2 = json.load(open(data10[0]))['sta2_loc']

	#crosscorr_stack_data_normalized_org_lst = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in data_to_plot]
	#stacked_data = sum(crosscorr_stack_data_normalized_org_lst)/len(crosscorr_stack_data_normalized_org_lst)
	stacked_data = sum(data_to_plot)/len(data_to_plot)

	CrossCorrelation_10_days_dic = nested_dict()
	CrossCorrelation_10_days_dic['dist'] = round(dist_pair)
	CrossCorrelation_10_days_dic['date'] = filename.split('.')[0].split('_')[-2]+'.'+filename.split('.')[0].split('_')[-1]
	CrossCorrelation_10_days_dic['sta1_loc'] = loc_sta1
	CrossCorrelation_10_days_dic['sta1_name'] = name_sta1
	CrossCorrelation_10_days_dic['sta2_loc'] = loc_sta2
	CrossCorrelation_10_days_dic['sta2_name'] = name_sta2
	CrossCorrelation_10_days_dic['crosscorr_daily_causal_time'] = time_to_plot[0].tolist()
	CrossCorrelation_10_days_dic['crosscorr_daily_10data'] = stacked_data.tolist()

	output_CrossCorrelation_DAY = JSON_FILES+'CROSS_CORR_10_DAYS_STACKED_FILES/'+name_sta1+'.'+name_sta2+'/'
	os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)
	with open(output_CrossCorrelation_DAY+'CROSS_CORR_10_DAYS_STACKED_'+name_sta1+'_'+name_sta2+'.json', 'w') as fp:
		json.dump(CrossCorrelation_10_days_dic, fp)


	#Creating the figure and plotting Stacked Cross-correlations
	fig = plt.figure(figsize=(10, 25))
	fig.suptitle('Cross-correlations: '+name_sta1+'-'+name_sta2+' - '+'days stacked: '+str(len(crosscorr_pairs_10day_data)),fontsize=20)

	gs = gridspec.GridSpec(5, 2,wspace=0.2, hspace=0.5)

	#-------------------------------------------

	map_loc = fig.add_subplot(gs[0,:],projection=ccrs.PlateCarree())

	LLCRNRLON_LARGE = -52
	URCRNRLON_LARGE = -36
	LLCRNRLAT_LARGE = -30
	URCRNRLAT_LARGE = -10

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
	text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-15,x=15)

	map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k',alpha=0.5, transform=ccrs.PlateCarree())
	map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
	map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
	map_loc.set_title('Dist = '+str(round(dist_pair))+' km',fontsize=20)

	#-------------------------------------------

	ax1 = fig.add_subplot(gs[1,:])

	ax1.plot(time_to_plot[0],stacked_data,c='k',lw=0.5)
	ax1.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)

	ax1.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax1.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax1.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax1.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)

	ax1.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

	# adding labels
	ax1.set_xlabel('Lapse time (s)',fontsize=14)
	ax1.set_title('Stacked Data')

	#--------------------------------------------------------------------------------------------------------------------

	ax2 = fig.add_subplot(gs[2,0])
	crosscorr_stack_data_normalized_org_lsts = [bandpass(data_2_plot, 1.0/25, 1.0/7, NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in data_to_plot]
	crosscorr_stack_data_normalized_org_lst = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in crosscorr_stack_data_normalized_org_lsts]

	y_factor = 1
	for i,j in enumerate(crosscorr_stack_data_normalized_org_lst):
		ax2.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst[i]],c='k',lw=0.5)

	ax2.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))])
	ax2.set_yticklabels([i.strftime("%d/%m/%y") if i==date_to_plot[0] or i==date_to_plot[-1] else None for i in date_to_plot])

	ax2.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)

	ax2.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax2.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax2.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax2.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)

	ax2.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

	# adding labels
	ax2.set_xlabel('Lapse time (s)',fontsize=14)
	ax2.set_title('Filter: 7s-25s')

	#---------------------------------------------------------

	ax3 = fig.add_subplot(gs[2,1])

	vector_plot = np.array(crosscorr_stack_data_normalized_org_lst)

	extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst)]
	im = ax3.imshow(vector_plot,extent=extent,origin='lower', interpolation='kaiser',aspect='auto',cmap='bwr')
	ax3.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-')
	ax3.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

	ax3.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax3.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax3.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax3.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax3.set_yticks([])

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

	crosscorr_stack_data_normalized_org_lst_20_50s = [bandpass(data_2_plot, 1.0/50, 1.0/20, NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in data_to_plot]
	crosscorr_stack_data_normalized_org_lst_20_50 = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in crosscorr_stack_data_normalized_org_lst_20_50s]

	ax4 = fig.add_subplot(gs[3,0])
	y_factor1 = 0.5
	for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_20_50):
		ax4.plot(time_to_plot[i],[x+i/y_factor1 for x in crosscorr_stack_data_normalized_org_lst_20_50[i]],c='k',lw=0.5)

	ax4.set_yticks([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))])
	ax4.set_yticklabels([i.strftime("%d/%m/%y") if i==date_to_plot[0] or i==date_to_plot[-1] else None for i in date_to_plot])

	ax4.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
	ax4.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)
	ax4.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax4.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax4.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax4.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	# adding labels
	ax4.set_xlabel('Lapse time (s)',fontsize=14)
	ax4.set_title('Filter: 20s-50s')

	#---------------------------------------------------------

	ax5 = fig.add_subplot(gs[3,1])

	vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_20_50)

	extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst_20_50)]
	im = ax5.imshow(vector_plot,extent=extent,origin='lower', interpolation='kaiser',aspect='auto',cmap='bwr')
	ax5.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-')
	ax5.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

	ax5.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax5.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax5.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax5.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)

	ax5.set_yticks([])

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

	crosscorr_stack_data_normalized_org_lst_50_100s = [bandpass(data_2_plot, 1.0/100, 1.0/50, NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in data_to_plot]
	crosscorr_stack_data_normalized_org_lst_50_100 = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in crosscorr_stack_data_normalized_org_lst_20_50s]

	ax6 = fig.add_subplot(gs[4,0])
	y_factor2 = 0.5
	for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_50_100):
		ax6.plot(time_to_plot[i],[x+i/y_factor2 for x in crosscorr_stack_data_normalized_org_lst_50_100[i]],c='k',lw=0.5)

	ax6.set_yticks([i/y_factor2 for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))])
	ax6.set_yticklabels([i.strftime("%d/%m/%y") if i==date_to_plot[0] or i==date_to_plot[-1] else None for i in date_to_plot])

	ax6.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
	ax6.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

	ax6.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax6.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax6.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax6.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	# adding labels
	ax6.set_xlabel('Lapse time (s)',fontsize=14)
	ax6.set_title('Filter: 50s-100s')

	#---------------------------------------------------------

	ax7 = fig.add_subplot(gs[4,1])

	vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_50_100)

	extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst_50_100)]
	im = ax7.imshow(vector_plot,extent=extent,origin='lower', interpolation='kaiser',aspect='auto',cmap='bwr')
	ax7.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-')
	ax7.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

	ax7.axvline(x=dist_pair/SIGNAL_WINDOW_VMIN, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax7.axvline(x=dist_pair/SIGNAL_WINDOW_VMAX, ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax7.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMIN), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)
	ax7.axvline(x=-(dist_pair/SIGNAL_WINDOW_VMAX), ymin=0, ymax=1,color='gray',linestyle='--',lw=1)

	ax7.set_yticks([])

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

	output_figure_CrossCorrelation_DAY = ORIENTATION_OUTPUT+'CROSS_CORR_10_STACK_FIGURES/'
	os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
	fig.savefig(output_figure_CrossCorrelation_DAY+'CROSS_CORR_10_STACK_FIG_'+name_sta1+'_'+name_sta2+'.png')
	plt.close()

print('\n')
print('===================================================================================')
print('Stacking 10-day Cross-correlations and plotting according to interstation distance:')
print('===================================================================================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(JSON_FILES+'CROSS_CORR_ORIENTATION_DAY_FILES/*'))

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

crosscorr_pairs_data_10_day_all = []
crosscorr_pairs_distance_10_day_all = []
for j,i in enumerate(crosscorr_pairs_data):
	crosscorr_pairs_10day_data = []
	crosscorr_pair_date_10day = []
	for k in i:
		subdir, filename = os.path.split(k)
		crosscorr_pair_date = datetime.datetime.strptime(filename.split('.')[0].split('_')[-2]+'.'+filename.split('.')[0].split('_')[-1], '%Y.%j')
		crosscorr_pair_date_10day.append(crosscorr_pair_date)
		crosscorr_pairs_10day_data.append([file for file in i if datetime.datetime.strptime(file.split('/')[-1].split('.')[0].split('_')[-2]+'.'+file.split('/')[-1].split('.')[0].split('_')[-1], '%Y.%j') >= crosscorr_pair_date and datetime.datetime.strptime(file.split('/')[-1].split('.')[0].split('_')[-2]+'.'+file.split('/')[-1].split('.')[0].split('_')[-1], '%Y.%j') < crosscorr_pair_date+datetime.timedelta(days=10)])

	#Stacking data
	data_to_plot = []
	for ind,data10 in enumerate(crosscorr_pairs_10day_data):
		name_sta1 = json.load(open(data10[0]))['sta1_name']
		name_sta2 = json.load(open(data10[0]))['sta2_name']
		dist_pair = json.load(open(data10[0]))['dist']

		causal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_causal'] for a in data10]),axis=0)
		acausal_lst = np.mean(np.array([json.load(open(a))['crosscorr_daily_acausal'] for a in data10]),axis=0)

		causal_time = np.array(json.load(open(data10[0]))['crosscorr_daily_causal_time'])
		acausal_time = np.array(json.load(open(data10[0]))['crosscorr_daily_acausal_time'])

		data_to_plot.append(acausal_lst[::-1] + causal_lst)
		time_to_plot = [-1*i for i in acausal_time[::-1]] + causal_time

		loc_sta1 = json.load(open(data10[0]))['sta1_loc']
		loc_sta2 = json.load(open(data10[0]))['sta2_loc']

	print('Pair '+str(j+1)+' of '+str(len(crosscorr_pairs_data))+': '+name_sta1+'-'+name_sta2+' - '+'days stacked: '+str(len(i)))

	crosscorr_pairs_data_10_day_all.append(sum(data_to_plot)/len(data_to_plot))
	crosscorr_pairs_distance_10_day_all.append(dist_pair)


print('\n')
print('============================')
print('Stacking Cross-correlations:')
print('============================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(JSON_FILES+'CROSS_CORR_ORIENTATION_DAY_FILES/*'))

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

#Stacking data

start_time = time.time()

pool = Pool(processes=num_processes)
CrossCorrelation_days_lst = []
for result in tqdm(pool.imap(func=crosscorr_stack, iterable=crosscorr_pairs_data), total=len(crosscorr_pairs_data)):
	CrossCorrelation_days_lst.append(result)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))


print('\n')
print('=======================================================================')
print('Plotting Stacked Cross-correlations according to interstation distance:')
print('=======================================================================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(JSON_FILES+'CROSS_CORR_STACKED_FILES/*'))

crosscorr_pairs_lst = []
for i,j in enumerate(crosscorr_days_lst):
	crosscorr_file = sorted(glob.glob(j+'/*'))
	crosscorr_pairs_lst.append(crosscorr_file)

#Make a list of list flat
crosscorr_pairs = [item for sublist in crosscorr_pairs_lst for item in sublist]

#Creating the figure
fig = plt.figure(figsize=(10, 20))
fig.suptitle('Cross-correlations according to interstation distance',fontsize=20)

gs = gridspec.GridSpec(4, 2,wspace=0.2, hspace=0.5)

#-------------------------------------------

map_loc = fig.add_subplot(gs[0,:],projection=ccrs.PlateCarree())

LLCRNRLON_LARGE = -52
URCRNRLON_LARGE = -36
LLCRNRLAT_LARGE = -30
URCRNRLAT_LARGE = -10

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

for i in tqdm(crosscorr_pairs):
	#Reading data
	name_sta1 = json.load(open(i))['sta1_name']
	name_sta2 = json.load(open(i))['sta2_name']

	if 'OBS' in name_sta1 or 'OBS' in name_sta2:

		crosscorr_stack_style_lst.append('k')
		crosscorr_stack_name_lst.append(name_sta1+'-'+name_sta2)
		dist_pair = json.load(open(i))['dist']
		crosscorr_stack_data_normalized_dist_lst.append(dist_pair)
		crosscorr_stack_data_normalized_vmin_lst.append(dist_pair/SIGNAL_WINDOW_VMIN)
		crosscorr_stack_data_normalized_vmax_lst.append(dist_pair/SIGNAL_WINDOW_VMAX)

		loc_sta1 = json.load(open(i))['sta1_loc']
		loc_sta2 = json.load(open(i))['sta2_loc']

	    #Stacked data
		crosscorr_stack_data = json.load(open(i))['crosscorr_stack']
		crosscorr_stack_data_normalized_lst.append(crosscorr_stack_data)
		crosscorr_stack_time = json.load(open(i))['crosscorr_stack_time']
		time_to_plot.append(crosscorr_stack_time)

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
ax2 = fig.add_subplot(gs[1,0])
crosscorr_stack_data_normalized_org_lsts = [bandpass(data_2_plot, 1.0/25, 1.0/7, NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
crosscorr_stack_data_normalized_org_lst = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in crosscorr_stack_data_normalized_org_lsts]

y_factor = 1
for i,j in enumerate(crosscorr_stack_data_normalized_org_lst):
	ax2.plot(time_to_plot[i],[x+i/y_factor for x in crosscorr_stack_data_normalized_org_lst[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

ax2.set_yticks([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))][::5])
ax2.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::5])

#ax2.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,21,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')
#ax2.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,21,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')

#ax2.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],21,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')
#ax2.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],21,1),savgol_filter([i/y_factor for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')

ax2.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
ax2.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

# adding labels
ax2.set_xlabel('Lapse time (s)',fontsize=14)
ax2.set_ylabel('Distance (km)',fontsize=14)
ax2.set_title('Filter: 7s-25s')

#---------------------------------------------------------
ax3 = fig.add_subplot(gs[1,1])

vector_plot = np.array(crosscorr_stack_data_normalized_org_lst)

extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst)]
im = ax3.imshow(vector_plot,extent=extent,origin='lower', interpolation='kaiser',aspect='auto',cmap='bwr')
ax3.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='--')
ax3.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

ax3.set_yticks([i for i in range(len(crosscorr_stack_data_normalized_org_lst))][::5])
ax3.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::5])

#ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')
#ax3.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')

#ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')
#ax3.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst))],21,1),ls='--',lw=1,c='gray')

# adding labels
ax3.set_xlabel('Lapse time (s)',fontsize=14)
ax3.set_ylabel('Distance (km)',fontsize=14)

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
crosscorr_stack_data_normalized_org_lst_20_50s = [bandpass(data_2_plot, 1.0/50, 1.0/20, NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
crosscorr_stack_data_normalized_org_lst_20_50 = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in crosscorr_stack_data_normalized_org_lst_20_50s]

ax4 = fig.add_subplot(gs[2,0])
y_factor1 = 0.5
for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_20_50):
	ax4.plot(time_to_plot[i],[x+i/y_factor1 for x in crosscorr_stack_data_normalized_org_lst_20_50[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

ax4.set_yticks([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))][::5])
ax4.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::5])

#ax4.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')
#ax4.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')

#ax4.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')
#ax4.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')

ax4.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
ax4.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

# adding labels
ax4.set_xlabel('Lapse time (s)',fontsize=14)
ax4.set_ylabel('Distance (km)',fontsize=14)
ax4.set_title('Filter: 20s-50s')

#---------------------------------------------------------
ax5 = fig.add_subplot(gs[2,1])

vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_20_50)

extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst_20_50)]
im = ax5.imshow(vector_plot,extent=extent,origin='lower', interpolation='kaiser',aspect='auto',cmap='bwr')
ax5.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='--')
ax5.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

ax5.set_yticks([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))][::5])
ax5.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::5])

#ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')
#ax5.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')

#ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')
#ax5.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_20_50))],21,1),ls='--',lw=1,c='gray')

# adding labels
ax5.set_xlabel('Lapse time (s)',fontsize=14)
ax5.set_ylabel('Distance (km)',fontsize=14)

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
crosscorr_stack_data_normalized_org_lst_50_100s = [bandpass(data_2_plot, 1.0/100, 1.0/50, NEW_SAMPLING_RATE, corners=2, zerophase=False) for data_2_plot in crosscorr_stack_data_normalized_org_lst]
crosscorr_stack_data_normalized_org_lst_50_100 = [(2*(a-a.min())/(a.max()-a.min()))-1 for a in crosscorr_stack_data_normalized_org_lst_50_100s]

ax6 = fig.add_subplot(gs[3,0])
y_factor1 = 0.5
for i,j in enumerate(crosscorr_stack_data_normalized_org_lst_50_100):
	ax6.plot(time_to_plot[i],[x+i/y_factor1 for x in crosscorr_stack_data_normalized_org_lst_50_100[i]],c=crosscorr_stack_style_org_lst[i],lw=0.5)

ax6.set_yticks([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))][::5])
ax6.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::5])

#ax6.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')
#ax6.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')

#ax6.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')
#ax6.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],21,1),savgol_filter([i/y_factor1 for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')

ax6.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='-',lw=1)
ax6.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

# adding labels
ax6.set_xlabel('Lapse time (s)',fontsize=14)
ax6.set_ylabel('Distance (km)',fontsize=14)
ax6.set_title('Filter: 50s-100s')

#---------------------------------------------------------
ax7 = fig.add_subplot(gs[3,1])

vector_plot = np.array(crosscorr_stack_data_normalized_org_lst_50_100)

extent = [-SHIFT_LEN,SHIFT_LEN,0,len(crosscorr_stack_data_normalized_org_lst_50_100)]
im = ax7.imshow(vector_plot,extent=extent,origin='lower', interpolation='kaiser',aspect='auto',cmap='bwr')
ax7.axvline(x=0, ymin=0, ymax=1,color='k',linestyle='--')
ax7.set_xlim(-SIGNAL2NOISE_TRAIL,SIGNAL2NOISE_TRAIL)

ax7.set_yticks([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))][::5])
ax7.set_yticklabels([str(int(i)) for i in crosscorr_stack_data_normalized_dist_org_lst][::5])


#ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmax_org_lst,21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')
#ax7.plot(savgol_filter(crosscorr_stack_data_normalized_vmin_org_lst,21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')

#ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmax_org_lst],21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')
#ax7.plot(savgol_filter([-i for i in  crosscorr_stack_data_normalized_vmin_org_lst],21,1),savgol_filter([i for i in range(len(crosscorr_stack_data_normalized_org_lst_50_100))],21,1),ls='--',lw=1,c='gray')


# adding labels
ax7.set_xlabel('Lapse time (s)',fontsize=14)
ax7.set_ylabel('Distance (km)',fontsize=14)

axins = inset_axes(ax7,
                   width="30%",  # width = 10% of parent_bbox width
                   height="2%",  # height : 5%
                   loc='upper left',
                   bbox_to_anchor=(0.65, 0.03, 1, 1),
                   bbox_transform=ax7.transAxes,
                   borderpad=0,
                   )
plt.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top')

output_figure_CrossCorrelation_DAY = ORIENTATION_OUTPUT+'CROSS_CORR_STACK_INTERSTATION_DISTANCE_FIGURES/'
os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)
fig.savefig(output_figure_CrossCorrelation_DAY+'CROSS_CORR_STACK_INTERSTATION_DISTANCE_FIG.png',dpi=300)
plt.close()

print('\n')
print('========================')
print('Clock Drift Calculating:')
print('========================')
print('\n')

#Collecting daily list of cross-correlations
crosscorr_days_lst = sorted(glob.glob(JSON_FILES+'CROSS_CORR_ORIENTATION_DAY_FILES/*'))

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

#Stacking data
for i in tqdm(crosscorr_pairs_data):
	crosscorr_pair_date_filename = [filename.split('/')[-1] for filename in i]
	crosscorr_pair_date = [datetime.datetime.strptime(filename.split('.')[0].split('_')[-2]+'.'+filename.split('.')[0].split('_')[-1], '%Y.%j') for filename in crosscorr_pair_date_filename]

	name_sta1 = json.load(open(i[0]))['sta1_name']
	name_sta2 = json.load(open(i[0]))['sta2_name']
	dist_pair = json.load(open(i[0]))['dist']

	if 'OBS' in name_sta1 or 'OBS' in name_sta2:

		causal_lst = np.array([json.load(open(a))['crosscorr_daily_causal'] for a in i])
		acausal_lst = np.array([json.load(open(a))['crosscorr_daily_acausal'] for a in i])

		loc_sta1 = json.load(open(i[0]))['sta1_loc']
		loc_sta2 = json.load(open(i[0]))['sta2_loc']

		#Creating the figure and plotting Clock-drift
		fig = plt.figure(figsize=(8, 15))
		fig.suptitle('Clock-drift between: '+name_sta1+'-'+name_sta2,fontsize=20)

		gs = gridspec.GridSpec(5, 1,wspace=0.2, hspace=0.5)
		map_loc = fig.add_subplot(gs[0:2],projection=ccrs.PlateCarree())

		LLCRNRLON_LARGE = -52
		URCRNRLON_LARGE = -36
		LLCRNRLAT_LARGE = -30
		URCRNRLAT_LARGE = -10

		map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
		map_loc.yaxis.set_ticks_position('both')
		map_loc.xaxis.set_ticks_position('both')

		map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+2,2), crs=ccrs.PlateCarree())
		map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+2,2), crs=ccrs.PlateCarree())
		map_loc.tick_params(labelbottom=True, labeltop=True, labelleft=True, labelright=True, labelsize=12)
		map_loc.grid(True,which='major',color='k',linewidth=1,linestyle='-')

		reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
		shape_1_SHP = list(reader_1_SHP.geometries())
		plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
		map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)
		# Use the cartopy interface to create a matplotlib transform object
		# for the Geodetic coordinate system. We will use this along with
		# matplotlib's offset_copy function to define a coordinate system which
		# translates the text by 25 pixels to the left.
		geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
		text_transform = offset_copy(geodetic_transform, units='dots', y=-5,x=80)

		map_loc.plot([loc_sta1[1],loc_sta2[1]],[loc_sta1[0],loc_sta2[0]],c='k', transform=ccrs.PlateCarree())
		map_loc.scatter(loc_sta1[1],loc_sta1[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
		map_loc.scatter(loc_sta2[1],loc_sta2[0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

		map_loc.text(loc_sta1[1],loc_sta1[0], name_sta1,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
		map_loc.text(loc_sta2[1],loc_sta2[0], name_sta2,fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)

		#-------------------------------------------

		days_major = DayLocator(interval=3)   # every day
		days_minor = DayLocator(interval=1)   # every day
		months = MonthLocator()  # every month
		yearsFmt = DateFormatter('%b-%Y')

		ax0 = fig.add_subplot(gs[2])
		ax0.xaxis.set_major_locator(months)
		ax0.xaxis.set_major_formatter(yearsFmt)
		ax0.xaxis.set_minor_locator(days_minor)
		ax0.yaxis.set_major_locator(MultipleLocator(0.1))
		ax0.yaxis.set_minor_locator(MultipleLocator(0.01))
		ax0.set_xlabel('Time (days)')
		ax0.set_ylabel('Static drift (s)')
		ax0.set_ylim(-.2,.2)

		ax1 = fig.add_subplot(gs[3])
		ax1.xaxis.set_major_locator(months)
		ax1.xaxis.set_major_formatter(yearsFmt)
		ax1.xaxis.set_minor_locator(days_minor)
		ax1.yaxis.set_major_locator(MultipleLocator(0.1))
		ax1.yaxis.set_minor_locator(MultipleLocator(0.01))
		ax1.set_xlabel('Time (days)')
		ax1.set_ylabel('Dynamic drift (s)')
		ax1.set_ylim(-.2,.2)


		ax2 = fig.add_subplot(gs[4])
		ax2.xaxis.set_major_locator(months)
		ax2.xaxis.set_major_formatter(yearsFmt)
		ax2.xaxis.set_minor_locator(days_minor)
		ax2.yaxis.set_major_locator(MultipleLocator(0.1))
		ax2.yaxis.set_minor_locator(MultipleLocator(0.01))
		ax2.set_xlabel('Time (days)')
		ax2.set_ylabel('Absolute drift (s)')
		ax2.set_ylim(-.2,.2)

		data_to_plot_clock_dynamic = []
		data_to_plot_clock_static = []
		data_to_plot_clock_absolute = []
		date_to_plot_clock = []

		sigma = 1 #70% of the data

		for k in range(len(causal_lst)):

				data_acausal_causal = np.array(acausal_lst[k][::-1] + causal_lst[k])
				data_normalized = (2*(data_acausal_causal-data_acausal_causal.min())/(data_acausal_causal.max()-data_acausal_causal.min()))-1

				#Collecting daily list of 10-day stack cross-correlations
				stacked_10_day_data = np.array(json.load(open(glob.glob(JSON_FILES+'CROSS_CORR_10_DAYS_STACKED_FILES/'+name_sta1+'.'+name_sta2+'/*')[0]))['crosscorr_daily_10data'])
				stacked_10_day_data_normalized = (2*(stacked_10_day_data-stacked_10_day_data.min())/(stacked_10_day_data.max()-stacked_10_day_data.min()))-1

				cc = obscorr(stacked_10_day_data_normalized,data_normalized,SHIFT_LEN)
				shift, max_value = xcorr_max(cc)

				cc_negative = obscorr(stacked_10_day_data_normalized[:len(stacked_10_day_data_normalized)//2], data_normalized[:len(stacked_10_day_data_normalized)//2],SHIFT_LEN)
				shift_negative, value_negative = xcorr_max(cc_negative)

				cc_positive = obscorr(stacked_10_day_data_normalized[len(stacked_10_day_data_normalized)//2:], data_normalized[len(stacked_10_day_data_normalized)//2:],SHIFT_LEN)
				shift_positive, value_positive = xcorr_max(cc_positive)

				cc_static_clock_drift = obscorr(stacked_10_day_data_normalized[:len(stacked_10_day_data_normalized)//2],stacked_10_day_data_normalized[len(stacked_10_day_data_normalized)//2:],SHIFT_LEN)
				shift_static_clock_drift, value_static_clock_drift = xcorr_max(cc_static_clock_drift)

				static_clock_drift = value_static_clock_drift

				dynamic_clock_drift = (value_positive-value_negative)/2

				absolute_clock_drift = dynamic_clock_drift + value_static_clock_drift

				#-------------------------------------------
				date_to_plot_clock.append(crosscorr_pair_date[k])
				data_to_plot_clock_static.append(static_clock_drift)
				data_to_plot_clock_dynamic.append(dynamic_clock_drift)
				data_to_plot_clock_absolute.append(absolute_clock_drift)

		date_to_plot_clock = np.array(date_to_plot_clock)
		data_to_plot_clock_static = np.array(data_to_plot_clock_static)
		data_to_plot_clock_dynamic = np.array(data_to_plot_clock_dynamic)
		data_to_plot_clock_absolute = np.array(data_to_plot_clock_absolute)

		True_False_lst_static = [True if np.mean(data_to_plot_clock_static)-sigma*np.std(data_to_plot_clock_static) <= j <= np.mean(data_to_plot_clock_static)+sigma*np.std(data_to_plot_clock_static) else False for j in data_to_plot_clock_static]
		True_False_lst_dynamic = [True if np.mean(data_to_plot_clock_dynamic)-sigma*np.std(data_to_plot_clock_dynamic) <= j <= np.mean(data_to_plot_clock_dynamic)+sigma*np.std(data_to_plot_clock_dynamic) else False for j in data_to_plot_clock_dynamic]
		True_False_lst_absolute = [True if np.mean(data_to_plot_clock_absolute)-sigma*np.std(data_to_plot_clock_absolute) <= j <= np.mean(data_to_plot_clock_absolute)+sigma*np.std(data_to_plot_clock_absolute) else False for j in data_to_plot_clock_absolute]

		date_to_plot_clock_True_static = date_to_plot_clock[True_False_lst_static]
		data_to_plot_clock_True_static = data_to_plot_clock_static[True_False_lst_static]

		date_to_plot_clock_True_dynamic = date_to_plot_clock[True_False_lst_dynamic]
		data_to_plot_clock_True_dynamic = data_to_plot_clock_dynamic[True_False_lst_dynamic]

		date_to_plot_clock_True_absolute = date_to_plot_clock[True_False_lst_absolute]
		data_to_plot_clock_True_absolute = data_to_plot_clock_absolute[True_False_lst_absolute]

		#--------------------------------------------------

		poly_reg = PolynomialFeatures(degree=4)
		X_poly = poly_reg.fit_transform(np.array(range(len(date_to_plot_clock_True_static))).reshape(-1, 1))
		pol_reg = LinearRegression()
		pol_reg.fit(X_poly, data_to_plot_clock_True_static)

		for i,j in enumerate(data_to_plot_clock_static):
			if np.mean(data_to_plot_clock_static)-sigma*np.std(data_to_plot_clock_static) <= j <= np.mean(data_to_plot_clock_static)+sigma*np.std(data_to_plot_clock_static):
				ax0.plot(date_to_plot_clock[i],data_to_plot_clock_static[i],'ok',ms=3)
			else:
				ax0.plot(date_to_plot_clock[i],data_to_plot_clock_static[i],'or',ms=3)

		ax0.plot(date_to_plot_clock_True_static, pol_reg.predict(poly_reg.fit_transform(np.array(range(len(data_to_plot_clock_True_static))).reshape(-1, 1))), color='blue')

		#--------------------------------------------------

		poly_reg = PolynomialFeatures(degree=4)
		X_poly = poly_reg.fit_transform(np.array(range(len(date_to_plot_clock_True_dynamic))).reshape(-1, 1))
		pol_reg = LinearRegression()
		pol_reg.fit(X_poly, data_to_plot_clock_True_dynamic)

		for i,j in enumerate(data_to_plot_clock_dynamic):
			if np.mean(data_to_plot_clock_dynamic)-sigma*np.std(data_to_plot_clock_dynamic) <= j <= np.mean(data_to_plot_clock_dynamic)+sigma*np.std(data_to_plot_clock_dynamic):
				l1, = ax1.plot(date_to_plot_clock[i],data_to_plot_clock_dynamic[i],'ok',ms=3)
			else:
				l2, = ax1.plot(date_to_plot_clock[i],data_to_plot_clock_dynamic[i],'or',ms=3)

		ax1.plot(date_to_plot_clock_True_dynamic, pol_reg.predict(poly_reg.fit_transform(np.array(range(len(data_to_plot_clock_True_dynamic))).reshape(-1, 1))), color='blue')
		ax1.legend((l1,l2),('%70 data','%30 data'),loc='upper right')

		#--------------------------------------------------

		poly_reg = PolynomialFeatures(degree=4)
		X_poly = poly_reg.fit_transform(np.array(range(len(date_to_plot_clock_True_absolute))).reshape(-1, 1))
		pol_reg = LinearRegression()
		pol_reg.fit(X_poly, data_to_plot_clock_True_absolute)

		for i,j in enumerate(data_to_plot_clock_absolute):
			if np.mean(data_to_plot_clock_absolute)-sigma*np.std(data_to_plot_clock_absolute) <= j <= np.mean(data_to_plot_clock_absolute)+sigma*np.std(data_to_plot_clock_absolute):
				l1, = ax2.plot(date_to_plot_clock[i],data_to_plot_clock_absolute[i],'ok',ms=3)
			else:
				l2, = ax2.plot(date_to_plot_clock[i],data_to_plot_clock_absolute[i],'or',ms=3)

		ax2.plot(date_to_plot_clock_True_absolute, pol_reg.predict(poly_reg.fit_transform(np.array(range(len(data_to_plot_clock_True_absolute))).reshape(-1, 1))), color='blue')
		ax2.legend((l1,l2),('%70 data','%30 data'),loc='upper right')

		fig.autofmt_xdate()

		output_figure_CLOCK_DRIFT = ORIENTATION_OUTPUT+'CLOCK_DRIFT_FIGURES/'
		os.makedirs(output_figure_CLOCK_DRIFT,exist_ok=True)
		fig.savefig(output_figure_CLOCK_DRIFT+'CLOCK_DRIFT_BETWEEN_'+name_sta1+'_'+name_sta2+'.png',dpi=300)
		plt.close()
'''
