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
    print('Processing file: '+input_list)
    try:
        channel_lst = obspy.read(input_list,headonly=True)
        for i,j in enumerate(channel_lst):
            endtime = j.stats.endtime
            date_time = datetime.date(endtime.year,endtime.month,endtime.day)
        
        return  {'input_list':input_list,'endtime':str(endtime),'date_time':str(date_time)}
    except:
        return 1

# ====================================
# Function to plot DATA availability
# ====================================

def plot_data_mosaic(date_lst,kstnm):
    days = DayLocator()   # every day
    months = MonthLocator()  # every month
    yearsFmt = DateFormatter('%Y-%m-%d')


    fig, ax = plt.subplots(nrows=len(date_lst), ncols=1,figsize=(20,10),sharex=True)
    plt.subplots_adjust(hspace=.1)

    

    for i,j in enumerate(date_lst):
        data_y = np.ones_like(j).tolist()
        ax[i].plot(j,data_y,"s",color='k',markersize=30)
        ax[i].xaxis.set_major_locator(days)
        ax[i].xaxis.set_major_formatter(yearsFmt)
        ax[i].xaxis.set_minor_locator(days)
        ax[i].set_ylim(0.99,1.01)
        ax[i].set_yticks([])
        ax[i].grid('on')
        ax[i].set_ylabel(kstnm[i],rotation=0, fontsize=20, labelpad=50)
        plt.setp(ax[i].xaxis.get_majorticklabels(), rotation=30 )
    fig.suptitle('Network Data Availability')
    os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
    fig.savefig(OUTPUT_FIGURE_DIR+'Network_Data_Availability.pdf',dpi=300)

    plt.show()