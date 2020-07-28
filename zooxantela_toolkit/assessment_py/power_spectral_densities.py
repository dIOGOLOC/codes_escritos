'''
Script to estimate probabilistic power spectral densities for 
one combination of network/station/location/channel/sampling_rate.
(https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html)

Calculations are based on the routine used by [McNamara2004]:
McNamara, D. E. and Buland, R. P. (2004),
Ambient Noise Levels in the Continental United States,
Bulletin of the Seismological Society of America, 94 (4), 1517-1527.
http://www.bssaonline.org/content/94/4/1517.abstract. 


For information on New High/Low Noise Model see [Peterson1993]:
Peterson, J. (1993),
Observations and Modeling of Seismic Background Noise,
U.S. Geological Survey open-file report 93-322, Albuquerque, N.M.
http://ehp3-earthquake.wr.usgs.gov/regional/asl/pubs/files/ofr93-322.pdf
'''


import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from obspy import UTCDateTime
import obspy
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,XML_FILE,OUTPUT_PSD_DIR,INITIAL_DATE,FINAL_DATE,
					TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR
				   )

# ====================================
# Function to calculate PSD from file
# ====================================

def calc_PSD(file):
    st = obspy.read(file)
    l = st[0]

    sta_name = l.stats.station
    NETWORK_CODE = l.stats.network
    sta_channel = l.stats.channel
        
    time_data = l.stats.starttime
    time_data_year = '{:04}'.format(time_data.year)
    time_data_julday = '{:03}'.format(time_data.julday)
    time_data_hour = '{:02}'.format(time_data.hour)
    time_data_minute = '{:02}'.format(time_data.minute)

    inv = obspy.read_inventory(XML_FILE)

    ppsd = PPSD(l.stats, inv)
    ppsd.add(st) 
    os.makedirs(OUTPUT_PSD_DIR+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/',exist_ok=True)
    ppsd.save_npz(OUTPUT_PSD_DIR+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+sta_channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
    
    return 0


# =======================
# Ploting TOTAL PPSD DATA
# =======================

def plot_PSD(directory):
	files = sorted(glob.glob(directory+'/*'))

	ppsd = PPSD.load_npz(files[0])

	[ppsd.add_npz(i) for i in files[1:]]

	ppsd.calculate_histogram(starttime=UTCDateTime(INITIAL_DATE),endtime=UTCDateTime(FINAL_DATE),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])    
	folder_output = OUTPUT_FIGURE_DIR+'WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd.station+'/'
	os.makedirs(folder_output,exist_ok=True)
	ppsd.plot(cmap=pqlx,filename=folder_output+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')