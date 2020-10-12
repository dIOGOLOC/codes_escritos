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


#Configuration file

MSEED_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/cross_cor_data/' 

STATIONXML_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

CLOCK_DRIFT_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/FIGURAS/'    

JSON_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/JSON_FILES/'

FIRSTDAY = '2019-08-02'
LASTDAY = '2019-08-05'

MIN_WINDOWS = 24

WINDOW_LENGTH = 3600

FILTER_DATA = [0.004,0.2]

NEW_SAMPLING_RATE = 50

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

def rotate_dir(tr1, tr2, direc):

    d = -direc*np.pi/180.+np.pi/2.
    rot_mat = np.array([[np.cos(d), -np.sin(d)],
                        [np.sin(d), np.cos(d)]])

    v12 = np.array([tr2, tr1])
    vxy = np.tensordot(rot_mat, v12, axes=1)
    tr_2 = vxy[0, :]
    tr_1 = vxy[1, :]

    return tr_1

# =======
# Classes
# =======

class StationDayData:
    """
    Class to save station info: name, network, channel, files directory, coordinates, datetime and data.
    """

    def __init__(self, name, network, fileID,lon=None,lat=None,time_day=None,data_day=None):
        """
        @type name: str
        @type network: str
        @type channel: str
        @type filesdir: str or unicode
        @type lon: float
        @type lat: float
        @type time_day: L{obspy.core.utcdatetime.UTCDateTime}
        @type data_day: L{obspy.core.trace.Trace}
        """
        self.name = name
        self.network = network
        self.fileID = fileID
        self.lon = lon
        self.lat = lat
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

def get_stations_data(f):
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
    tr = st[0]
    st.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH/2)
        
    inv = read_inventory(STATIONXML_DIR+'.'.join([network,name,'xml']))
    tr.remove_response(inventory=inv,output="DISP",water_level=60)
    
    st.detrend('demean')
    st.detrend('linear')
    st.filter('bandpass', freqmin=0.01,freqmax=10,corners=2, zerophase=True)    
    st.resample(NEW_SAMPLING_RATE)

    lon = inv.get_coordinates(sta_channel_id)['longitude']
    lat = inv.get_coordinates(sta_channel_id)['latitude']
            
    # appending new channel day data
    station = StationDayData(name=name, network=network,fileID=sta_channel_id,lon=lon,lat=lat,time_day=time_day,data_day=st[0])

    data_dic = {
    			'name': station.name,
    			'network': station.network,
    			'time_day': station.time_day,
    			'fileID': station.fileID,
    			'data_day': station.data_day.data.tolist()
    			}

    output_DATA_DAY = JSON_FILES+'DATA_DAY_FILES/'+station.network+'.'+station.name+'/'
    os.makedirs(output_DATA_DAY,exist_ok=True)
    with open(output_DATA_DAY+'DATA_DAY_'+station.network+'_'+station.name+'_'+sta_channel_id+'_'+year_day+'_'+julday_day+'.json', 'w') as fp:
    	json.dump(data_dic, fp)

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

print('==================================================================')
print('Opening miniseed files and calculating the spectrogram of each day')
print('==================================================================')
print('\n')

start_time = time.time()

with Pool(processes=num_processes) as p:
	max_ = len(files)
	with tqdm(total=max_) as pbar:
		for i, _ in enumerate(p.imap_unordered(get_stations_data, files)):
			pbar.update()

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')


#-------------------------------------------------------------------------------

print('====================')
print('Calculating the tilt')
print('====================')
print('\n')

last_daily_lst_spec = list(compress(second_daily_lst_spec, good_days.tolist()))
last_daily_lst_data = list(compress(second_daily_lst_data, good_days.tolist()))

last_good_windows_lst = list(compress(second_good_windows_lst, good_days.tolist()))

tilt_max_coh_lst = []
max_coh_lst = []
phase_max_coh_lst = []
day_time_tilt = []
for i,j in enumerate(tqdm(last_daily_lst_data)):
    
    data_HHZ = json.load(open([k for k in j if 'HHZ' in k][0]))
    data_HHE = json.load(open([k for k in j if 'HHE' in k][0]))
    data_HHN = json.load(open([k for k in j if 'HHN' in k][0]))

    year_day = data_HHZ["time_day"].split('.')[0]
    julday_day = data_HHZ["time_day"].split('.')[1]

    tr_HHZ = Trace()
    tr_HHZ.data = np.array(data_HHZ['data_day'])
    tr_HHZ.stats.starttime = UTCDateTime(year=int(year_day),julday=int(julday_day))
    tr_HHZ.stats.sampling_rate = NEW_SAMPLING_RATE
    tr_HHZ.stats.channel = 'HHZ'
    tr_HHZ.stats.station = data_HHZ['name']
    tr_HHZ.stats.network = data_HHZ['network']

    tr_HHN = Trace()
    tr_HHN.data = np.array(data_HHN['data_day'])
    tr_HHN.stats.starttime = UTCDateTime(year=int(year_day),julday=int(julday_day))
    tr_HHN.stats.sampling_rate = NEW_SAMPLING_RATE
    tr_HHN.stats.channel = 'HHN'
    tr_HHN.stats.station = data_HHZ['name']
    tr_HHN.stats.network = data_HHZ['network']
    
    tr_HHE = Trace()
    tr_HHE.data = np.array(data_HHE['data_day'])
    tr_HHE.stats.starttime = UTCDateTime(year=int(year_day),julday=int(julday_day))
    tr_HHE.stats.sampling_rate = NEW_SAMPLING_RATE
    tr_HHE.stats.channel = 'HHN'
    tr_HHE.stats.station = data_HHZ['name']
    tr_HHE.stats.network = data_HHZ['network']
    
    slide_HHZ = np.array([k.data for k in tr_HHZ.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])[last_good_windows_lst[i]]
    slide_HHE = np.array([k.data for k in tr_HHN.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])[last_good_windows_lst[i]]
    slide_HHN = np.array([k.data for k in tr_HHE.slide(window_length=WINDOW_LENGTH, step=WINDOW_LENGTH)])[last_good_windows_lst[i]]

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
    figA.savefig(daily_A_output+'Admittance.Coherence.Phase.'+data_HHZ["time_day"]+'.png', dpi=300, facecolor='w', edgecolor='w')
'''
    #-----------------------------------------------------------------------------------------------------
    
    #Calculate TILT
    
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

    print('Tilf of Maximum coherence = ', tilt)
    print('Maximum coherence = ', coh_value)
    print('Phase of Maximum coherence = ', phase_value)
    
    tilt_max_coh_lst.append(tilt)
    day_time_tilt.append([int(year_day),int(julday_day)])
    max_coh_lst.append(coh_value)
    phase_max_coh_lst.append(phase_value)

    # Now calculate spectra at tilt direction
    ftH = rotate_dir(np.array(N_lst), np.array(E_lst), tilt)

    # Get transfer functions
    cHH = np.abs(np.mean(np.conj(ftH)*ftH, axis=0))
    cHZ = np.mean(np.conj(ftH)*np.array(Z_lst), axis=0)
    
    #-----------------------------------------------------------------------------------------------------
    
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
    figureTILT.savefig(daily_A_output+'Tilt.Coherence.'+data_HHZ["time_day"]+'.png', dpi=300, facecolor='w', edgecolor='w')
   
#-------------------------------------------------------------------------------

daily_data = [mdates.num2date(UTCDateTime(year=i[0],julday=i[1]).matplotlib_date) for i in day_time_tilt]

figureTILTfinal, (ax1, ax2) = plt.subplots(2, 1,sharex=True,figsize=(10,5))
ax1.plot(daily_data, max_coh_lst, 'ok')
ax1.set_ylabel('Coherence')
ax1.set_ylim(0, 1)
ax2.set_xlim(mdates.num2date(INTERVAL_PERIOD[0].matplotlib_date),mdates.num2date(INTERVAL_PERIOD[1].matplotlib_date))


ax2.plot(daily_data, tilt_max_coh_lst, 'ok')
ax2.set_ylabel('Tilt')
ax2.set_xlabel('Day')
ax2.set_ylim(0, 360)
ax2.yaxis.set_major_locator(MultipleLocator(90))
ax2.yaxis.set_minor_locator(MultipleLocator(10))
ax2.set_xlim(mdates.num2date(INTERVAL_PERIOD[0].matplotlib_date),mdates.num2date(INTERVAL_PERIOD[-1].matplotlib_date))
plt.gcf().autofmt_xdate()


figureTILTfinal.suptitle('Station = '+STATION+' - Period = '+UTCDateTime(year=int(INTERVAL_PERIOD[0].year),julday=int(INTERVAL_PERIOD[0].julday)).strftime('%d/%m/%Y')+'--'+UTCDateTime(year=int(INTERVAL_PERIOD[-1].year),julday=int(INTERVAL_PERIOD[-1].julday)).strftime('%d/%m/%Y'),fontsize=18)
figureTILTfinal.savefig(daily_A_output+'Tilt.Coherence.total.png', dpi=300, facecolor='w', edgecolor='w')
'''