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
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream, Trace
from obspy.io.xseed import Parser
from obspy.signal.filter import bandpass,lowpass
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.util import prev_pow_2
from obspy.signal.cross_correlation import correlate as obscorr

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


#Configuration file

MSEED_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/cross_cor_data/' 

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

CLOCK_DRIFT_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/FIGURAS/'    

JSON_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/JSON_FILES/'

#Shapefile  boundary states
BOUNDARY_STATES_SHP = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/brasil_estados.shp'

FIRSTDAY = '2019-08-01'
LASTDAY = '2019-08-02'

MIN_WINDOWS = 40

WINDOW_LENGTH = 3600

#freq window (Hz) to smooth ampl spectrum
WINDOW_FREQ = 0.0002 

#max time window (s) for cross-correlation
SHIFT_LEN = 3600

NEW_SAMPLING_RATE = 5

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
    files_lst = []
    for root, dirs, files in os.walk(basedir):
    	for file in files:
    		files_path  = os.path.join(root, file)
    		if any(day_s in files_path for day_s in interval_period_date):
    			files_lst.append(files_path)

    file_lsts = sorted([i for i in files_lst if 'HHZ' in i])
	    
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

def moving_avg(a, halfwindow, mask=None):
    """
    Performs a fast n-point moving average of (the last
    dimension of) array *a*, by using stride tricks to roll
    a window on *a*.
    Note that *halfwindow* gives the nb of points on each side,
    so that n = 2*halfwindow + 1.
    If *mask* is provided, values of *a* where mask = False are
    skipped.
    Returns an array of same size as *a* (which means that near
    the edges, the averaging window is actually < *npt*).
    """
    # padding array with zeros on the left and on the right:
    # e.g., if halfwindow = 2:
    # a_padded    = [0 0 a0 a1 ... aN 0 0]
    # mask_padded = [F F ?  ?      ?  F F]

    if mask is None:
        mask = np.ones_like(a, dtype='bool')

    zeros = np.zeros(a.shape[:-1] + (halfwindow,))
    falses = zeros.astype('bool')

    a_padded = np.concatenate((zeros, np.where(mask, a, 0), zeros), axis=-1)
    mask_padded = np.concatenate((falses, mask, falses), axis=-1)

    # rolling window on padded array using stride trick
    #
    # E.g., if halfwindow=2:
    # rolling_a[:, 0] = [0   0 a0 a1 ...    aN]
    # rolling_a[:, 1] = [0  a0 a1 a2 ... aN 0 ]
    # ...
    # rolling_a[:, 4] = [a2 a3 ...    aN  0  0]

    npt = 2 * halfwindow + 1  # total size of the averaging window
    rolling_a = as_strided(a_padded,
                           shape=a.shape + (npt,),
                           strides=a_padded.strides + (a.strides[-1],))
    rolling_mask = as_strided(mask_padded,
                              shape=mask.shape + (npt,),
                              strides=mask_padded.strides + (mask.strides[-1],))

    # moving average
    n = rolling_mask.sum(axis=-1)
    return np.where(n > 0, rolling_a.sum(axis=-1).astype('float') / n, np.nan)

#-------------------------------------------------------------------------------

def get_stations_data(f,onebit_norm=True,white_spectral=False):
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

    st = read(f)

    st_traces = [k for k in st.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH/2)]
    st_hours = [str(k[0].stats.starttime.hour)+':'+str(k[0].stats.starttime.minute) for k in st_traces]

    if len(st_hours) > MIN_WINDOWS: 

	    inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))

	    traces_resp = [tr.remove_response(inventory=inv,output="DISP",water_level=60) for tr in st_traces]
	    traces_demean = [tr.detrend('demean') for tr in traces_resp]
	    traces_detrend = [tr.detrend('linear') for tr in traces_demean]
	    traces_filter = [tr.filter('bandpass', freqmin=0.01,freqmax=10,corners=2, zerophase=True) for tr in traces_detrend]
	    traces_resample = [tr.resample(NEW_SAMPLING_RATE) for tr in traces_filter]

	    # ======================
	    # One-bit normalization
	    # ======================

	    if onebit_norm:
	    	traces_onebit = [np.sign(tr[0].data) for tr in traces_resample]
	    	traces_data_day = [tr.data.tolist() for tr in traces_onebit]

	    else:
	    	traces_onebit = traces_resample
	    	traces_data_day = [tr[0].data.tolist() for tr in traces_onebit]


	    # ==================
	    # Spectral whitening
	    # ==================

	    if white_spectral:
	    	window_freq = WINDOW_FREQ
	    	fft = rfft(trace.data)  # real FFT
	    	deltaf = trace.stats.sampling_rate / trace.stats.npts  # frequency step

	    	# smoothing amplitude spectrum
	    	halfwindow = int(round(window_freq / deltaf / 2.0))

	    	# normalizing spectrum and back to time domain
	    	trace.data = irfft(fft / weight, n=len(trace.data))

	    	# re bandpass to avoid low/high freq noise
	    	trace.filter(type="bandpass",freqmin=freqmin,freqmax=freqmax,corners=corners,zerophase=zerophase)

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

	    output_DATA_DAY = JSON_FILES+'DATA_DAY_FILES/'+year_day+'.'+julday_day+'/'
	    os.makedirs(output_DATA_DAY,exist_ok=True)
	    with open(output_DATA_DAY+'DATA_DAY_'+sta_channel_id+'_'+year_day+'_'+julday_day+'.json', 'w') as fp:
	    	json.dump(data_dic, fp)

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

    day_crosscor_causal = CrossCorrelation(sta1['name'],sta2['name'],sta1['lat'],sta1['lon'],sta2['lat'],sta2['lon'],sta1['time_day'])
    day_crosscor_acausal = CrossCorrelation(sta2['name'],sta1['name'],sta2['lat'],sta2['lon'],sta1['lat'],sta1['lon'],sta1['time_day'])

    day_crosscor_causal.add(sta1['data_day'],sta2['data_day'],sta1['hours_day'],sta2['hours_day'])
    day_crosscor_acausal.add(sta2['data_day'],sta1['data_day'],sta2['hours_day'],sta1['hours_day'])

    CrossCorrelation_dic['dist'] = round(day_crosscor_causal.dist())
    CrossCorrelation_dic['sta1_loc'] = [sta1['lat'],sta1['lon']]
    CrossCorrelation_dic['sta1_name'] = sta1['name']
    CrossCorrelation_dic['sta2_loc'] = [sta2['lat'],sta2['lon']]
    CrossCorrelation_dic['sta2_name'] = sta2['name']
    CrossCorrelation_dic['crosscorr_daily_causal_time'] = day_crosscor_causal.timearray.tolist()
    CrossCorrelation_dic['crosscorr_daily_causal'] = day_crosscor_causal.dataarray.tolist()
    CrossCorrelation_dic['crosscorr_daily_acausal'] = day_crosscor_acausal.dataarray.tolist()
    CrossCorrelation_dic['crosscorr_daily_acausal_time'] = day_crosscor_acausal.timearray.tolist()

    output_CrossCorrelation_DAY = JSON_FILES+'CROSS_CORR_DAY_FILES/'+year_day+'.'+julday_day+'/'
    os.makedirs(output_CrossCorrelation_DAY,exist_ok=True)
    with open(output_CrossCorrelation_DAY+'CROSS_CORR_DAY_FILES_'+sta1['name']+'_'+sta2['name']+'_'+year_day+'_'+julday_day+'.json', 'w') as fp:
        json.dump(CrossCorrelation_dic, fp)

    # ============================
    # Plot: map and pair crosscorr
    # ============================
	
    fig = plt.figure(figsize=(15, 15))
    #fig.suptitle('Dia do Evento - '+UTCDateTime(year=int(event_date[0]),julday=int(event_date[1])).strftime('%d/%m/%Y')+' - Magnitude:'+str(stZ[0].stats.sac.mag),fontsize=20)

    gs = gridspec.GridSpec(2, 1,wspace=0.2, hspace=0.5)

    #-------------------------------------------

    map_loc = fig.add_subplot(gs[0],projection=ccrs.PlateCarree())
		
    LLCRNRLON_LARGE = -52
    URCRNRLON_LARGE = -28
    LLCRNRLAT_LARGE = -28
    URCRNRLAT_LARGE = -16

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
    text_transform = offset_copy(geodetic_transform, units='dots', y=0,x=60)
    text_transform_mag = offset_copy(geodetic_transform, units='dots', y=-15,x=15)
    
    map_loc.scatter(sta1['lon'],sta1['lat'], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())    
    map_loc.scatter(sta2['lon'],sta2['lat'], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

    map_loc.plot([sta1['lon'],sta1['lat']],[sta2['lon'],sta2['lat']], transform=ccrs.PlateCarree())
    
    map_loc.text(sta1['lon'],sta1['lat'], sta1['name'],fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
    map_loc.text(sta2['lon'],sta2['lat'], sta2['name'],fontsize=15,verticalalignment='center', horizontalalignment='right',transform=text_transform)
	#-------------------------------------------
    
    ax = fig.add_subplot(gs[1])
    data_to_plot = CrossCorrelation_dic['crosscorr_daily_acausal'][::-1]+CrossCorrelation_dic['crosscorr_daily_causal']  
    time_to_plot = [-1*i for i in CrossCorrelation_dic['crosscorr_daily_acausal_time'][::-1]] + CrossCorrelation_dic['crosscorr_daily_causal_time']
    ax.plot(time_to_plot,data_to_plot,color='k')
    ax.set_xlabel('time (s)',fontsize=14)

    output_figure_CrossCorrelation_DAY = CLOCK_DRIFT_OUTPUT+'CROSS_CORR_DAY_FIGURES/'+year_day+'.'+julday_day+'/'
    os.makedirs(output_figure_CrossCorrelation_DAY,exist_ok=True)    
    fig.savefig(output_figure_CrossCorrelation_DAY+'CROSS_CORR_DAY_FIG_'+sta1['name']+'_'+sta2['name']+'_'+year_day+'_'+julday_day+'.png')
	
    return sta1['time_day']
    
# =======
# Classes
# =======

class StationDayData:
    """
    Class to save station info: name, network, channel, files directory, coordinates, datetime and data.
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
        self.lat1 = lat2
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

        	xcorr = sum(xcorr_hours)/len(xcorr_hours)
        	xcorr_timearray = np.arange(0,shift_len,1/(2*NEW_SAMPLING_RATE))
        # normalizing cross-corr
        self.dataarray = xcorr

        # time arrya cross-corr
        self.timearray = xcorr_timearray

# ============
# Main program
# ============

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

files = filelist(basedir=MSEED_DIR,interval_period_date=INTERVAL_PERIOD_DATE)

print('Total of miniseed files = '+str(len(files)))
print('\n')

print('==================================')
print('Opening miniseed files of each day')
print('==================================')
print('\n')
'''
start_time = time.time()

with Pool(processes=num_processes) as p:
	max_ = len(files)
	with tqdm(total=max_) as pbar:
		for i, _ in enumerate(p.imap_unordered(get_stations_data, files)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

'''
print('====================================')
print('Calculating daily Cross-correlation:')
print('====================================')
print('\n')

days_crosscor = sorted(glob.glob(JSON_FILES+'DATA_DAY_FILES/*'))

start_time = time.time()

for i,j in enumerate(days_crosscor):

	print('Day '+str(i+1)+' of '+str(len(days_crosscor)))
	stations_file = sorted(glob.glob(j+'/*'))

	stationtrace_pairs = list(combinations(stations_file, 2))

	pool = Pool(processes=num_processes)
	CrossCorrelation_days_lst = []
	for result in tqdm(pool.imap(func=crosscorr_func, iterable=stationtrace_pairs), total=len(stationtrace_pairs)):
	    CrossCorrelation_days_lst.append(result)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
