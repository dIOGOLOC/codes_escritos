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
import matplotlib.gridspec as gridspec

# ======
# [paths]
# ======

#STATION
STATIONS = ['OBS17','OBS18','OBS19','OBS20','OBS22']

#Directory to save PSD
OUTPUT_PSD_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_PSD_MSD/'

#Directory to save Figures
OUTPUT_FIGURE_DIR = '/home/diogoloc/dados_posdoc/ON_MAR/Figuras/PQLX/'

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

def plot_PSD_PQLX(files2_lst,init_date2,fin_date2,channel):

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	data_model_noise = np.load(NOISE_MODEL_FILE)
	periods_model_noise_ln = data_model_noise['model_periods']
	nlnm_ln = data_model_noise['low_noise']

	periods_model_noise_hn = data_model_noise['model_periods']
	nlnm_hn = data_model_noise['high_noise']
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	folder_output = OUTPUT_FIGURE_DIR+'NOISE_ANALYSIS/'
	os.makedirs(folder_output,exist_ok=True)
	filename = folder_output+channel+'_all_stations.png'

	percentiles=[0, 25, 50, 75, 100]
	period_lim=(0.02, 200)
	cumulative_number_of_colors=20

	fig = plt.figure(figsize=(10,5))
	fig.ppsd = AttribDict()
	ax = fig.add_subplot(111)

	stations_line_style = [(0, (5, 1)),(0, (3, 1, 1, 1, 1, 1)),'dotted', 'solid', (0, (3, 1, 1, 1))]

	for i, files2 in enumerate(files2_lst):
		ppsd2 = PPSD.load_npz(files2[0])
		[ppsd2.add_npz(i) for i in files2[1:]]
		periods2, mean_2 = ppsd2.get_mean()
		xdata2 = periods2
		ppsd2.calculate_histogram(starttime=UTCDateTime(init_date2),endtime=UTCDateTime(fin_date2),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])

		ax.plot(xdata2, mean_2, ls=stations_line_style[i],c='k', zorder=10, label=STATIONS[i])
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	ax.plot(periods_model_noise_ln,nlnm_ln, color='gray',ls='-', linewidth=2, zorder=10, label='Noise model limits')
	ax.plot(periods_model_noise_hn,nlnm_hn, color='gray',ls='-', linewidth=2, zorder=10)
	ax.legend(loc='lower left')

	fig.ppsd.cmap = 'Blues_r'
	fig.ppsd.label = "[%]"
	fig.ppsd.max_percentage = 30
	fig.ppsd.color_limits = (0, 30)

	ax.set_title('Noise Spectra ('+channel+')')
	ax.semilogx()
	ax.set_xlabel('Period [s]')
	ax.set_xlim(period_lim)
	ax.set_ylim(AMP_PSD_MIN, AMP_PSD_MAX)
	ax.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
	ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))

	files2_obs17 = files2_lst[0]
	ppsd_obs17 = PPSD.load_npz(files2_obs17[0])
	[ppsd_obs17.add_npz(i) for i in files2_obs17[1:]]
	data = (ppsd_obs17.current_histogram*100.0/(ppsd_obs17.current_histogram_count or 1))
	xedges = ppsd_obs17.period_xedges

	fig.ppsd.meshgrid = np.meshgrid(xedges, ppsd_obs17.db_bin_edges)
	ppsd = ax.pcolormesh(fig.ppsd.meshgrid[0], fig.ppsd.meshgrid[1], data.T,cmap=fig.ppsd.cmap, zorder=-1)
	fig.ppsd.quadmesh = ppsd

	cb = plt.colorbar(ppsd, ax=ax)
	#cb.set_clim(*fig.ppsd.color_limits)
	cb.set_label(fig.ppsd.label)
	cb.ax.tick_params(labelcolor='white')
	cb.ax.yaxis.label.set_color('white')

	fig.ppsd.colorbar = cb

	ppsd.set_clim(*fig.ppsd.color_limits)

	ax.tick_params(labelcolor='white')
	ax.xaxis.label.set_color('white')
	ax.yaxis.label.set_color('white')
	ax.title.set_color('white')

	ax.grid(b=True, which="major")
	#ax.grid(b=True, which="minor")
	plt.savefig(filename,facecolor='None',dpi=300,bbox_inches='tight',pad_inches=0.3)

#------------------------------------------------------------------------------

def plot_PSD_PQLX_all(files2_lst,init_date2,fin_date2,channel):

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    data_model_noise = np.load(NOISE_MODEL_FILE)
    periods_model_noise_ln = data_model_noise['model_periods']
    nlnm_ln = data_model_noise['low_noise']

    periods_model_noise_hn = data_model_noise['model_periods']
    nlnm_hn = data_model_noise['high_noise']
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    folder_output = OUTPUT_FIGURE_DIR+'NOISE_ANALYSIS/'
    os.makedirs(folder_output,exist_ok=True)
    filename = folder_output+channel+'_PQLX_all_stations.png'

    percentiles=[0, 25, 50, 75, 100]
    period_lim=(0.02, 200)
    cumulative_number_of_colors=20

    fig = plt.figure(figsize=(15,10))
    fig.ppsd = AttribDict()

    gs = gridspec.GridSpec(2, 6)
    gs.update(wspace = 0.2, hspace = 0.2)

    ax1 = plt.subplot(gs[0, :2])
    ax2 = plt.subplot(gs[0, 2:4],sharey=ax1)
    ax3 = plt.subplot(gs[0, 4:6],sharey=ax1)
    ax4 = plt.subplot(gs[1, 1:3],sharey=ax1)
    ax5 = plt.subplot(gs[1, 3:5],sharey=ax1)

    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)

    ax = [ax1,ax2,ax3,ax4,ax5]

    for i, files2 in enumerate(files2_lst):

        ppsd2 = PPSD.load_npz(files2[0])
        [ppsd2.add_npz(i) for i in files2[1:]]
        periods2, mean_2 = ppsd2.get_mean()
        xdata2 = periods2
        ppsd2.calculate_histogram(starttime=UTCDateTime(init_date2),endtime=UTCDateTime(fin_date2),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])

        ax[i].plot(xdata2, mean_2, ls='--',c='k', zorder=10, label='mean')
        #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        ax[i].plot(periods_model_noise_ln,nlnm_ln, color='gray',ls='-', linewidth=2, zorder=10, label='Noise model limits')
        ax[i].plot(periods_model_noise_hn,nlnm_hn, color='gray',ls='-', linewidth=2, zorder=10)
        ax[i].legend(loc='lower left')

        fig.ppsd.cmap = 'Blues_r'
        fig.ppsd.label = "[%]"
        fig.ppsd.max_percentage = 30
        fig.ppsd.color_limits = (0, 30)

        title = "%s   %s--%s"
        title = title % (ppsd2.id,UTCDateTime(ns=ppsd2._times_processed[0]).date,
                             UTCDateTime(ns=ppsd2._times_processed[-1]).date)

        ax[i].set_title(title)
        ax[i].semilogx()
        ax[i].set_xlim(period_lim)
        ax[i].set_ylim(AMP_PSD_MIN, AMP_PSD_MAX)
        ax[i].xaxis.set_major_formatter(FormatStrFormatter("%g"))

        if ax[i] not in [ax2,ax3,ax5]:
            ax[i].set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
        else:
            pass

        if ax[i] not in [ax1,ax2,ax3]:
            ax[i].set_xlabel('Period [s]')
        else:
            pass

        if ax[i] not in [ax1,ax2,ax4,ax5]:
            ax[i].tick_params(labelbottom=True, labeltop=False, labelright=True)
            ax[i].set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
            ax[i].yaxis.set_label_position("right")
        else:
            pass

        files2_obs17 = files2_lst[i]
        ppsd_obs17 = PPSD.load_npz(files2_obs17[0])
        [ppsd_obs17.add_npz(i) for i in files2_obs17[1:]]
        data = (ppsd_obs17.current_histogram*100.0/(ppsd_obs17.current_histogram_count or 1))
        xedges = ppsd_obs17.period_xedges

        fig.ppsd.meshgrid = np.meshgrid(xedges, ppsd_obs17.db_bin_edges)
        ppsd = ax[i].pcolormesh(fig.ppsd.meshgrid[0], fig.ppsd.meshgrid[1], data.T,cmap=fig.ppsd.cmap, zorder=-1)
        fig.ppsd.quadmesh = ppsd

        ax[i].tick_params(labelcolor='white')
        ax[i].xaxis.label.set_color('white')
        ax[i].yaxis.label.set_color('white')
        ax[i].title.set_color('white')
        ax[i].tick_params(bottom=True, top=True, left=True, right=True)
        ax[i].grid(b=True, which="major")
        #ax[i].grid(b=True, which="minor")

    cbar_ax = fig.add_axes([0.8, 0.15, 0.025, 0.3])
    cb = plt.colorbar(ppsd, cax=cbar_ax)
    cb.set_label(fig.ppsd.label)
    cb.ax.tick_params(labelcolor='white')
    cb.ax.yaxis.label.set_color('white')
    fig.ppsd.colorbar = cb

    ppsd.set_clim(*fig.ppsd.color_limits)
    plt.savefig(filename,facecolor='None',dpi=300,bbox_inches='tight',pad_inches=0.3)

#------------------------------------------------------------------------------

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

print('Channel: HHZ')
plot_PSD_PQLX_all(files_HHZ2,INITIAL_DATE_2,FINAL_DATE_2,'HHZ')
#--------------------------------------------------------------------------------------------------------------------
print('Channel: HHE')
plot_PSD_PQLX_all(files_HHE2,INITIAL_DATE_2,FINAL_DATE_2,'HHE')
#--------------------------------------------------------------------------------------------------------------------
print('Channel: HHN')
plot_PSD_PQLX_all(files_HHN2,INITIAL_DATE_2,FINAL_DATE_2,'HHN')

#--------------------------------------------------------------------------------------------------------------------
print('Channel: HHZ')
plot_PSD_PQLX(files_HHZ2,INITIAL_DATE_2,FINAL_DATE_2,'HHZ')

#--------------------------------------------------------------------------------------------------------------------
print('Channel: HHE')
plot_PSD_PQLX(files_HHE2,INITIAL_DATE_2,FINAL_DATE_2,'HHE')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHN')
plot_PSD_PQLX(files_HHN2,INITIAL_DATE_2,FINAL_DATE_2,'HHN')
#--------------------------------------------------------------------------------------------------------------------

print('\n')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
