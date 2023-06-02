'''
Script to get information about the header of the raw data
(https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#obspy.core.stream.read)
and plot a mosaic the the Data availability.
'''

import matplotlib.pyplot as plt
import obspy
import os
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


from parameters_py.config import (
					OUTPUT_FIGURE_DIR
				   )

# ==============================
# Function to estimate the date
# ==============================

def get_date_file(input_list):
    #print('Processing day: '+input_list)
    channel_lst = obspy.read(input_list,headonly=True)
    endtime = channel_lst[0].stats.endtime
        
    return  {'input_list':input_list,'endtime':str(endtime)}
    
# ====================================
# Function to plot DATA availability
# ====================================

def plot_data_mosaic(date_lst,kstnm):
    #days = DayLocator()   # every day
    days = MonthLocator(1)   # every day
    months = MonthLocator(6)  # every month
    yearsFmt = DateFormatter('%Y-%b')


    fig, ax = plt.subplots(nrows=len(date_lst), ncols=1,figsize=(20,10),sharex=True)
    plt.subplots_adjust(hspace=.1)

    for i,j in enumerate(date_lst):
        data_y = np.ones_like(j).tolist()
        #ax[i].plot(j,data_y,"|",color='k',markersize=50)
        ax[i].plot(j,data_y,"s",color='k',markersize=50)
        ax[i].xaxis.set_major_locator(months)
        ax[i].xaxis.set_major_formatter(yearsFmt)
        ax[i].xaxis.set_minor_locator(days)
        ax[i].set_ylim(0.99,1.01)
        ax[i].set_yticks([])
        ax[i].grid('on')
        ax[i].set_ylabel(kstnm[i],rotation=0, fontsize=20, labelpad=50)
        plt.setp(ax[i].xaxis.get_majorticklabels(), rotation=30 )
    fig.suptitle('Data Availability')
    os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
    fig.savefig(OUTPUT_FIGURE_DIR+'Network_Data_Availability.pdf',dpi=300)

    plt.show()