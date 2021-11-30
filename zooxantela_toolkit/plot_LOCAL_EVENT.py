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


from visual_py.event_plot import plot_event_data,plot_map_event_data,plot_map_event_data_hydrophone

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,OUTPUT_EV_DIR,DIR_STATUS,NUM_PROCESS,OUTPUT_PSD_DIR,XML_FILE
					)

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

lst_eventsZ = []

for i in EVENT_dir:
  lst_eventsZ.append([k for k in HHZ_files if i in k ])

#--------------------------------------------------------------------------------------------------------------------

HHN_files = []

for root, dirs, files in os.walk(OUTPUT_EV_DIR+'Local/'):
    for file in files:
        if '.N' in file:
            HHN_files.append(os.path.join(root,file))

lst_eventsN = []

for i in EVENT_dir:
  lst_eventsN.append([k for k in HHN_files if i in k ])

#--------------------------------------------------------------------------------------------------------------------

HHE_files = []

for root, dirs, files in os.walk(OUTPUT_EV_DIR+'Local/'):
    for file in files:
        if '.E' in file:
            HHE_files.append(os.path.join(root,file))

lst_eventsE = []

for i in EVENT_dir:
  lst_eventsE.append([k for k in HHE_files if i in k ])

#-------------------------------------------------------------------------------------------------------------------

# ===================
# Ploting EVENT files
# ===================

print('Ploting Station x Event')
start_time = time.time()

for i,j in enumerate(lst_eventsZ):
  plot_map_event_data(lst_eventsZ[i],lst_eventsN[i],lst_eventsE[i],EVENT_dir[i],'Local/',10,20)

#-------------------------------------------------------------------------------------------------------------------
'''
HHX_files = []

for root, dirs, files in os.walk(OUTPUT_EV_DIR+'Local/'):
    for file in files:
        if '.X' in file:
            HHX_files.append(os.path.join(root,file))

HHX_files_RSBR = []

lst_eventsX = []

for i in EVENT_dir:
  lst_eventsX.append([k for k in HHX_files if i in k ])

#--------------------------------------------------------------------------------------------------------------------
for i,j in enumerate(lst_eventsX):
    if len(j) != 0:
        plot_map_event_data_hydrophone(lst_eventsX[i],EVENT_dir[i],'Local/',2,20)

#--------------------------------------------------------------------------------------------------------------------
'''
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
