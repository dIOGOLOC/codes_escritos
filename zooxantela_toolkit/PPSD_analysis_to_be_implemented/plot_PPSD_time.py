#!/usr/bin/python -u
"""
Scritp to estimate probabilistic power spectral densities for each .SAC file
"""

import time
import numpy as np
import os
from tqdm import tqdm
from multiprocessing import Pool
from obspy import read,read_inventory, UTCDateTime, Stream
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter
import matplotlib.dates as mdates
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import glob
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx,viridis
from obspy.core.util import AttribDict

# ======
# [paths]
# ======

#STATION
STATIONS = ['OBS17','OBS18','OBS19','OBS20','OBS22']

#Directory to save PSD
OUTPUT_PSD_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_PSD_MSD/'

#Directory to save Figures
OUTPUT_FIGURE_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/Figuras/'

#Directory witg noise models
NOISE_MODEL_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/noise_models.npz'

#PPSD start date
INITIAL_DATE_1 = '2019,8,1'

#PPSD final date
FINAL_DATE_1 = '2019,8,2'

#PPSD start date
INITIAL_DATE_2 = '2019,8,1'

#PPSD final date
FINAL_DATE_2 = '2020,6,15'

# =====
# [ppsd]
# =====

#Number worker processes
NUM_PROCESS = 8

#Percentage fo the days to process and plot the PPSD?
DAY_PERCENTAGE = 2

#Restricts the data that is included in the stack by time of day and weekday.
#Monday is 1, Sunday is 7, -1 for any day of week.
#For example, using time_of_weekday=[(-1, 22, 24)]
#only individual spectra that have a starttime in between 10pm and 12am are used in the stack for all days of week
#time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
TIME_OF_WEEKDAY_DAY = -1
TIME_OF_WEEKDAY_START_HOUR = 0
TIME_OF_WEEKDAY_FINAL_HOUR = 24

#Extracting PSD values for given period in seconds.
#Selects the period bin whose center period is closest to the specified period.
PERIOD_PSD = 10

#Maximum and minimum amplitude of the PSD (sensor).
AMP_PSD_MIN = -200
AMP_PSD_MAX = -65

#Maximum and minimum amplitude of the PSD (hydrophones).
AMP_PSD_HYDROPHONE_MIN = -60
AMP_PSD_HYDROPHONE_MAX = 60

# ========================
# Constants and parameters
# ========================

DTINY = np.finfo(0.0).tiny

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)


# =========
# Functions
# =========

def filelist(basedir,interval_period_date,channel):
    """
    Returns the list of files in *basedir* whose are in the specified period
    """
    files = []
    files_list = glob.glob(basedir+'/*')

    for s in files_list:
    	if any(day_s in s for day_s in interval_period_date):
    		files.append(s)

    files = [i for i in files if channel in i]

    return sorted(files)

#------------------------------------------------------------------------------

def plot_PSD(files1,init_date1,fin_date1,files2,init_date2,fin_date2):

    ppsd1 = PPSD.load_npz(files1[0])
    [ppsd1.add_npz(i) for i in files1[1:]]
    periods1, mean_1 = ppsd1.get_mean()
    xdata1 = periods1
    ppsd1.calculate_histogram(starttime=UTCDateTime(init_date1),endtime=UTCDateTime(fin_date1),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ppsd2 = PPSD.load_npz(files2[0])
    [ppsd2.add_npz(i) for i in files2[1:]]
    periods2, mean_2 = ppsd2.get_mean()
    xdata2 = periods2
    ppsd2.calculate_histogram(starttime=UTCDateTime(init_date2),endtime=UTCDateTime(fin_date2),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    data_model_noise = np.load(NOISE_MODEL_FILE)
    periods_model_noise_ln = data_model_noise['model_periods']
    nlnm_ln = data_model_noise['low_noise']

    periods_model_noise_hn = data_model_noise['model_periods']
    nlnm_hn = data_model_noise['high_noise']
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    folder_output = OUTPUT_FIGURE_DIR+'TIME_ANALYSIS_WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd1.station+'/'
    os.makedirs(folder_output,exist_ok=True)
    filename = folder_output+ppsd1.network+'.'+ppsd1.station+'.'+ppsd1.channel+'.'+str(ppsd1.times_processed[0].year)+'.pdf'

    percentiles=[0, 25, 50, 75, 100]
    period_lim=(0.01, 179)
    cumulative_number_of_colors=20

    fig = plt.figure()
    fig.ppsd = AttribDict()
    ax = fig.add_subplot(111)

    period_lim=(0.01, 179)
    ax.plot(xdata1, mean_1, '--k', zorder=9, label='Early days')
    ax.plot(xdata2, mean_2, '-k', zorder=10, label='Late days')

    ax.plot(periods_model_noise_ln,nlnm_ln, color='gray',ls='-', linewidth=2, zorder=10, label='Low noise model')
    ax.plot(periods_model_noise_hn,nlnm_hn, color='gray',ls='-', linewidth=2, zorder=10, label='High noise model')
    ax.legend(loc='lower left')

    fig.ppsd.cmap = 'viridis'
    fig.ppsd.label = "[%]"
    fig.ppsd.max_percentage = 30
    fig.ppsd.color_limits = (0, 30)

    title = "%s   %s -- %s  (%i/%i segments)"
    title = title % (ppsd1.id,UTCDateTime(ns=ppsd1._times_processed[0]).date,
                             UTCDateTime(ns=ppsd1._times_processed[-1]).date,
                             ppsd1.current_histogram_count,len(ppsd1._times_processed))

    ax.set_title(title)
    ax.semilogx()
    ax.set_xlabel('Period [s]')
    ax.set_xlim(period_lim)
    ax.set_ylim(AMP_PSD_MIN, AMP_PSD_MAX)
    ax.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))

    ax.grid(b=True, which="major")
    ax.grid(b=True, which="minor")
    #plt.show()
    #plt.savefig(filename)

    return [xdata1, mean_1,xdata2, mean_2]

#------------------------------------------------------------------------------

def plot_PSD_PQLX(files1,init_date1,fin_date1,files2,init_date2,fin_date2):

    ppsd1 = PPSD.load_npz(files1[0])
    [ppsd1.add_npz(i) for i in files1[1:]]
    periods1, mean_1 = ppsd1.get_mean()
    xdata1 = periods1
    ppsd1.calculate_histogram(starttime=UTCDateTime(init_date1),endtime=UTCDateTime(fin_date1),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ppsd2 = PPSD.load_npz(files2[0])
    [ppsd2.add_npz(i) for i in files2[1:]]
    periods2, mean_2 = ppsd2.get_mean()
    xdata2 = periods2
    ppsd2.calculate_histogram(starttime=UTCDateTime(init_date2),endtime=UTCDateTime(fin_date2),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    data_model_noise = np.load(NOISE_MODEL_FILE)
    periods_model_noise_ln = data_model_noise['model_periods']
    nlnm_ln = data_model_noise['low_noise']

    periods_model_noise_hn = data_model_noise['model_periods']
    nlnm_hn = data_model_noise['high_noise']
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    folder_output = OUTPUT_FIGURE_DIR+'TIME_ANALYSIS_WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd1.station+'/'
    os.makedirs(folder_output,exist_ok=True)
    filename = folder_output+ppsd1.network+'.'+ppsd1.station+'.'+ppsd1.channel+'.pdf'

    percentiles=[0, 25, 50, 75, 100]
    period_lim=(0.02, 179)
    cumulative_number_of_colors=20

    fig = plt.figure()
    fig.ppsd = AttribDict()
    ax = fig.add_subplot(111)

    period_lim=(0.01, 179)
    ax.plot(xdata1, mean_1, ls='--',c='k', zorder=10, label='Early days')
    ax.plot(xdata2, mean_2, ls='-',c='k', zorder=10, label='Late days')

    ax.plot(periods_model_noise_ln,nlnm_ln, color='gray',ls='-', linewidth=2, zorder=10, label='Low noise model')
    ax.plot(periods_model_noise_hn,nlnm_hn, color='gray',ls='-', linewidth=2, zorder=10, label='High noise model')
    ax.legend(loc='lower left')

    fig.ppsd.cmap = viridis
    fig.ppsd.label = "[%]"
    fig.ppsd.max_percentage = 30
    fig.ppsd.color_limits = (0, 30)

    ax.set_title('Noise Spectra')
    ax.semilogx()
    ax.set_xlabel('Period [s]')
    ax.set_xlim(period_lim)
    ax.set_ylim(AMP_PSD_MIN, AMP_PSD_MAX)
    ax.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))

    data = (ppsd2.current_histogram*100.0/(ppsd2.current_histogram_count or 1))
    xedges = ppsd2.period_xedges

    fig.ppsd.meshgrid = np.meshgrid(xedges, ppsd2.db_bin_edges)
    ppsd = ax.pcolormesh(fig.ppsd.meshgrid[0], fig.ppsd.meshgrid[1], data.T,cmap=fig.ppsd.cmap, zorder=-1)
    fig.ppsd.quadmesh = ppsd

    cb = plt.colorbar(ppsd, ax=ax)
    #cb.set_clim(*fig.ppsd.color_limits)
    cb.set_label(fig.ppsd.label)
    fig.ppsd.colorbar = cb

    ppsd.set_clim(*fig.ppsd.color_limits)

    ax.grid(b=True, which="major")
    ax.grid(b=True, which="minor")
    plt.show()
    #plt.savefig(filename)

#------------------------------------------------------------------------------

def plot_PSD_hydrophone(files1,init_date1,fin_date1,files2,init_date2,fin_date2):

    ppsd1 = PPSD.load_npz(files1[0])
    [ppsd1.add_npz(i) for i in files1[1:]]
    periods1, mean_1 = ppsd1.get_mean()
    xdata1 = periods1
    ppsd1.calculate_histogram(starttime=UTCDateTime(init_date1),endtime=UTCDateTime(fin_date1),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ppsd2 = PPSD.load_npz(files2[0])
    [ppsd2.add_npz(i) for i in files2[1:]]
    periods2, mean_2 = ppsd2.get_mean()
    xdata2 = periods2
    ppsd2.calculate_histogram(starttime=UTCDateTime(init_date2),endtime=UTCDateTime(fin_date2),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    data_model_noise = np.load(NOISE_MODEL_FILE)
    periods_model_noise_ln = data_model_noise['model_periods']
    nlnm_ln = data_model_noise['low_noise']

    periods_model_noise_hn = data_model_noise['model_periods']
    nlnm_hn = data_model_noise['high_noise']
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    folder_output = OUTPUT_FIGURE_DIR+'TIME_ANALYSIS_WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd1.station+'/'
    os.makedirs(folder_output,exist_ok=True)
    filename = folder_output+ppsd1.network+'.'+ppsd1.station+'.'+ppsd1.channel+'.'+str(ppsd1.times_processed[0].year)+'.pdf'

    percentiles=[0, 25, 50, 75, 100]
    period_lim=(0.01, 10)
    cumulative_number_of_colors=20

    fig = plt.figure()
    fig.ppsd = AttribDict()

    ax = fig.add_subplot(111)

    period_lim=(0.01, 179)
    ax.plot(xdata1, mean_1, '-k', zorder=9)
    ax.plot(xdata2, mean_2, '--g', zorder=10)

    fig.ppsd.cmap = viridis
    fig.ppsd.label = "[%]"
    fig.ppsd.max_percentage = 30
    fig.ppsd.color_limits = (0, 30)

    title = "%s   %s -- %s  (%i/%i segments)"
    title = title % (ppsd1.id,UTCDateTime(ns=ppsd1._times_processed[0]).date,
                             UTCDateTime(ns=ppsd1._times_processed[-1]).date,
                             ppsd1.current_histogram_count,len(ppsd1._times_processed))

    ax.set_title(title)
    ax.semilogx()
    ax.set_xlabel('Period [s]')
    ax.set_xlim(period_lim)
    ax.set_ylim(AMP_PSD_HYDROPHONE_MIN, AMP_PSD_HYDROPHONE_MAX)
    ax.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))

    ax.grid(b=True, which="major")
    ax.grid(b=True, which="minor")
    #plt.show()
    #plt.savefig(filename)

    return [xdata1, mean_1,xdata2, mean_2]

# =====================
# Retrieving .NPZ files
# =====================

print('\n')
print('Retrieving NPZ files')
print('\n')

HHZ_dir = []
HHE_dir = []
HHN_dir = []
HHX_dir = []

for station in STATIONS:
	for root, dirs, files in os.walk(OUTPUT_PSD_DIR+station+'/'):
	    for directory in dirs:
	        if "HHZ.PPSD" in directory:
	            HHZ_dir.append(os.path.join(root, directory))
	        if "HHE.PPSD" in directory:
	        	HHE_dir.append(os.path.join(root, directory))
	        if "HHN.PPSD" in directory:
	          	HHN_dir.append(os.path.join(root, directory))
	        if "HHX.PPSD" in directory:
	            HHX_dir.append(os.path.join(root, directory))

# =================
# Filtering by date
# =================

fday1 = UTCDateTime(INITIAL_DATE_1)
lday1 = UTCDateTime(FINAL_DATE_1)
INTERVAL_PERIOD1 = [UTCDateTime(x.astype(str)) for x in np.arange(fday1.datetime,lday1.datetime+ONEDAY,ONEDAY)]
INTERVAL_PERIOD_DATE1 = [str(x.year)+'.'+"%03d" % x.julday for x in INTERVAL_PERIOD1]

fday2 = UTCDateTime(INITIAL_DATE_2)
lday2 = UTCDateTime(FINAL_DATE_2)
INTERVAL_PERIOD2 = [UTCDateTime(x.astype(str)) for x in np.arange(fday2.datetime,lday2.datetime+ONEDAY,ONEDAY)]
INTERVAL_PERIOD_DATE2 = [str(x.year)+'.'+"%03d" % x.julday for x in INTERVAL_PERIOD2]

print('===============================')
print('Scanning name of miniseed files')
print('===============================')
print('\n')

files_HHZ1 = []
files_HHE1 = []
files_HHN1 = []
files_HHX1 = []

files_HHZ2 = []
files_HHE2 = []
files_HHN2 = []
files_HHX2 = []

for i,j in enumerate(HHX_dir):
	# initializing list of stations by scanning name of miniseed files
	files_HHZ1.append(filelist(basedir=HHZ_dir[i],interval_period_date=INTERVAL_PERIOD_DATE1,channel='HHZ'))
	files_HHE1.append(filelist(basedir=HHE_dir[i],interval_period_date=INTERVAL_PERIOD_DATE1,channel='HHE'))
	files_HHN1.append(filelist(basedir=HHN_dir[i],interval_period_date=INTERVAL_PERIOD_DATE1,channel='HHN'))
	files_HHX1.append(filelist(basedir=HHX_dir[i],interval_period_date=INTERVAL_PERIOD_DATE1,channel='HHX'))

	files_HHZ2.append(filelist(basedir=HHZ_dir[i],interval_period_date=INTERVAL_PERIOD_DATE2,channel='HHZ'))
	files_HHE2.append(filelist(basedir=HHE_dir[i],interval_period_date=INTERVAL_PERIOD_DATE2,channel='HHE'))
	files_HHN2.append(filelist(basedir=HHN_dir[i],interval_period_date=INTERVAL_PERIOD_DATE2,channel='HHN'))
	files_HHX2.append(filelist(basedir=HHX_dir[i],interval_period_date=INTERVAL_PERIOD_DATE2,channel='HHX'))

# ===================
# Ploting .NPZ files
# ===================

print('Ploting .NPZ files')
start_time = time.time()
print('\n')

data_HHZ = []
data_HHE = []
data_HHN = []
data_HHX = []

print('Channel: HHZ')
for i,j in enumerate(files_HHZ1):
	plot_PSD_PQLX(files_HHZ1[i],INITIAL_DATE_1,FINAL_DATE_1,files_HHZ2[i],INITIAL_DATE_2,FINAL_DATE_2)
print('\n')


#--------------------------------------------------------------------------------------------------------------------
print('Channel: HHZ')
for i,j in enumerate(files_HHZ1):
	data_HHZ.append(plot_PSD(files_HHZ1[i],INITIAL_DATE_1,FINAL_DATE_1,files_HHZ2[i],INITIAL_DATE_2,FINAL_DATE_2))
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHE')
for i,j in enumerate(files_HHE1):
	data_HHE.append(plot_PSD(files_HHE1[i],INITIAL_DATE_1,FINAL_DATE_1,files_HHE2[i],INITIAL_DATE_2,FINAL_DATE_2))
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHN')
for i,j in enumerate(files_HHN1):
	data_HHN.append(plot_PSD(files_HHN1[i],INITIAL_DATE_1,FINAL_DATE_1,files_HHN2[i],INITIAL_DATE_2,FINAL_DATE_2))
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHX')
for i,j in enumerate(files_HHX1):
	data_HHX.append(plot_PSD_hydrophone(files_HHX1[i],INITIAL_DATE_1,FINAL_DATE_1,files_HHX2[i],INITIAL_DATE_2,FINAL_DATE_2))


# =================
# Ploting all files
# =================

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_model_noise = np.load(NOISE_MODEL_FILE)
periods_model_noise_ln = data_model_noise['model_periods']
nlnm_ln = data_model_noise['low_noise']

periods_model_noise_hn = data_model_noise['model_periods']
nlnm_hn = data_model_noise['high_noise']
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#folder_output = OUTPUT_FIGURE_DIR+'TIME_ANALYSIS_WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd1.station+'/'
#os.makedirs(folder_output,exist_ok=True)
#filename = folder_output+'.network+'.'+ppsd1.station+'.'+ppsd1.channel+'.pdf'

percentiles=[0, 25, 50, 75, 100]
period_lim=(0.01, 179)
cumulative_number_of_colors=20

fig = plt.figure()
fig.ppsd = AttribDict()
ax = fig.add_subplot(111)

period_lim=(0.01, 179)

STATIONS_colors = ['k','r','g','b','y']

for i,j in enumerate(data_HHZ):
    ax.plot(j[2], j[3],c=STATIONS_colors[i],ls='--', zorder=10, label=STATIONS[i])

ax.plot(periods_model_noise_ln,nlnm_ln, color='gray',ls='-', linewidth=2, zorder=10, label='NLNM')
ax.plot(periods_model_noise_hn,nlnm_hn, color='gray',ls='-', linewidth=2, zorder=10, label='HLNM')
ax.legend(loc='lower left')

fig.ppsd.cmap = viridis
fig.ppsd.label = "[%]"
fig.ppsd.max_percentage = 30
fig.ppsd.color_limits = (0, 30)

ax.set_title('Probabilistic Power Spectral Densities Mean')
ax.semilogx()
ax.set_xlabel('Period [s]')
ax.set_xlim(period_lim)
ax.set_ylim(AMP_PSD_MIN, AMP_PSD_MAX)
ax.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))

ax.grid(b=True, which="major")
ax.grid(b=True, which="minor")
plt.show()
    #plt.savefig(filename)

print('\n')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
