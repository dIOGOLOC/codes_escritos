'''
Script to get information about the status header of the raw data
(https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#obspy.core.stream.read)
'''

import matplotlib.pyplot as plt
import pandas as pd
import obspy
from obspy import read,read_inventory, UTCDateTime, Stream
import os
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.signal import PPSD
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# ============================
# Function to find status data
# ============================

def filelist(basedir):
    """
    Returns the list of status files in *basedir*
    """
    files = []
    for folder in os.listdir(basedir):
        if folder.endswith("m8")  or folder.endswith("m9") or folder.endswith("ma") or folder.endswith("me"):
            folder_path = os.path.join(basedir, folder)
            s = glob.glob(folder_path+'/*')
            files.append(s)
    files = sorted(files)

    flat_list_files = [item for sublist in files for item in sublist]

    return flat_list_files

# =======================================
# Function to separate status data by day
# =======================================


def list_split_day(lst_status_files):
	"""
	Returns the list of status files filtered by day
	"""
	day_lst = []
	for i in lst_status_files:
		day_lst.append(i.split('/')[-1].split('_')[0])
	day_lst = sorted(list(set(day_lst)))

	daily_lst = [[]]*len(day_lst)
	for i,j in enumerate(day_lst):
		daily_lst[i] = [k for k in lst_status_files if j == k.split('/')[-1].split('_')[0]]

	return daily_lst

# ===============================
# Function to merge and plot data
# ===============================


def get_status_file_GURALP(lst):
	"""
	Extract data from GCF files
	"""

	try: 
		if len(lst) > 1: 
			st = Stream()
			for i in lst:
				st += read(i)
			data_lst = [k.data for k in st]
			
			flat_data_lst = [item for sublist in data_lst for item in sublist]

			s = str(st[0].stats.starttime.year)+'-'+"%03d" %st[0].stats.starttime.julday
			dataframe_lst = [st[0].stats.station,st[0].stats.channel,s,np.mean(flat_data_lst)]

		else:
			st = read(lst[0])
			s = str(st[0].stats.starttime.year)+'-'+"%03d" %st[0].stats.starttime.julday
			dataframe_lst = [st[0].stats.station,st[0].stats.channel,s,np.mean(st[0].data)]

		return dataframe_lst

	except:
		pass
	
