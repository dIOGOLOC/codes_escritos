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
					OUTPUT_FIGURE_DIR,INITIAL_DATE,FINAL_DATE
				   )

# ==============================
# Function to estimate the date
# ==============================

def get_date_file(directory_data):
    print('Processing folder: '+directory_data)
    os.chdir(directory_data)
    files = glob.glob('*')
    files = sorted(files)
    input_list = []
    starttime = []
    for i, j in enumerate(files):
        print("Calculating "+str(i+1)+" of "+str(len(files)))
        print("File = "+directory_data+j)
        channel_lst = obspy.read(directory_data+j,headonly=True)
        starttime.append(str(channel_lst[0].stats.starttime))
        input_list.append(directory_data+j)
    return  {'input_list':input_list,'starttime':starttime}
    
# ====================================
# Function to plot DATA availability
# ====================================

def plot_data_mosaic(date_lst,kstnm):
    days1 = DayLocator(interval=1)   # every day
    days5 = DayLocator(interval=5)   # every day
    months = MonthLocator()  # every month
    yearsFmt = DateFormatter('%Y-%m-%d')

    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,5))
    for i,j in enumerate(date_lst):
        data_y = np.ones_like(j).tolist()
        ax.bar(j,data_y,width=1,color='k')
        ax.xaxis.set_major_locator(days5)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(days1)
        ax.set_ylim(0,1.)
        ax.set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
        ax.set_yticks([])
        ax.grid(True)
        ax.set_ylabel(kstnm,rotation=0, fontsize=20, labelpad=50)
        plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)
        ax.tick_params(which='minor', length=4)
        ax.tick_params(which='major', length=10)

    fig.suptitle('Data Availability')
    os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
    fig.savefig(OUTPUT_FIGURE_DIR+'Network_Data_Availability.pdf',dpi=300)

    plt.show()