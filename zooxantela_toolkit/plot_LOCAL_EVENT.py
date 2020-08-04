#!/usr/bin/python -u
"""
Scritp to estimate probabilistic power spectral densities for each .SAC file
"""

import time
import os
from tqdm import tqdm
from multiprocessing import Pool
from obspy import read,read_inventory, UTCDateTime, Stream
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter
import matplotlib.dates as mdates
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt


from visual_py.event_plot import plot_event_data

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,OUTPUT_EV_DIR,DIR_STATUS,NUM_PROCESS,OUTPUT_PSD_DIR,XML_FILE
					)

# ===================
# Retrieving XML file
# ===================

inv = read_inventory(XML_FILE)

# =====================
# Retrieving .SAC files
# =====================

print('\n')
print('Retrieving Events files')
print('\n')

EVENT_dir = []

for root, dirs, files in os.walk(OUTPUT_EV_DIR+'Local/'):
    for directory in dirs:
        if len(directory.split('.')) > 5:
            EVENT_dir.append(os.path.join(directory))

EVENT_dir = sorted(list(set(EVENT_dir)))

#--------------------------------------------------------------------------------------------------------------------

HHZ_files = []

for root, dirs, files in os.walk(OUTPUT_EV_DIR+'Local/'):
    for file in files:
        if '.Z' in file:
            HHZ_files.append(os.path.join(root,file))

#--------------------------------------------------------------------------------------------------------------------


lst_events = []

for i in EVENT_dir:
  lst_events.append([k for k in HHZ_files if i in k])

# ===================
# Ploting EVENT files
# ===================
print('Ploting EVENT files')
start_time = time.time()

for i,j in enumerate(lst_events):
  plot_event_data(j,inv,EVENT_dir[i],'Local_HHZ/',2.0,10.0)

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')