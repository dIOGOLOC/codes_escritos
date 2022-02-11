#!/usr/bin/python -u

'''
--------------------------------------------------------------------------------
       Function to plot local the dataset according to local events time
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 02/2022


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will plot the local dataset according to a given list of an events.


Inputs:
Event traces (format: SAC)


Outputs:
Figures (PDF)


Examples of Usage (in command line):
   >> python plot_LOCAL_EVENT.py

'''

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


if LABEL_LANG == 'br':
    print('\n')
    print('Procurando os arquivos dos eventos locais.')
    print('\n')
else:
    print('\n')
    print('Collecting local events files.')
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

#--------------------------------------------------------------------------------------------------------------------

HHX_files = []

for root, dirs, files in os.walk(OUTPUT_EV_DIR+'Local/'):
    for file in files:
        if '.X' in file:
            HHX_files.append(os.path.join(root,file))

HHX_files_RSBR = []

lst_eventsX = []

for i in EVENT_dir:
  lst_eventsX.append([k for k in HHX_files if i in k ])


# -------------------------------------
# Saving all files into a MSEED stream:
# -------------------------------------

for i,j in enumerate(lst_eventsZ):
	event_name_info = j[0].split('/')[-2]
	st_all = Stream()

	z = [read(k)[0] for k in lst_eventsZ[i]]
	n = [read(k)[0] for k in lst_eventsN[i]]
	e = [read(k)[0] for k in lst_eventsE[i]]

	[st_all.append(z1) for z1 in z]
	[st_all.append(n1) for n1 in n]
	[st_all.append(e1) for e1 in e]

	output = OUTPUT_EV_DIR+'Local/NETWORK_MSEED_FILES/'+event_name_info+'/'
	os.makedirs(output,exist_ok=True)
	st_all.write(output+'event_'+event_name_info+'.mseed')

# ===================
# Ploting EVENT files
# ===================

if LABEL_LANG == 'br':
    print('Plotando Estação x Evento')

else:
    print('Ploting Station x Event')

start_time = time.time()

for i,j in enumerate(lst_eventsZ):
  plot_map_event_data(lst_eventsZ[i],lst_eventsN[i],lst_eventsE[i],EVENT_dir[i],'Local/',10,20)


#-------------------------------------------------------------------------------------------------------------------


if LABEL_LANG == 'br':
    print('Plotando Hidrofone x Evento')

else:
    print('Ploting Hydrophone x Event')

#--------------------------------------------------------------------------------------------------------------------
for i,j in enumerate(lst_eventsX):
    if len(j) != 0:
        plot_map_event_data_hydrophone(lst_eventsX[i],EVENT_dir[i],'Local/',2,20)

#--------------------------------------------------------------------------------------------------------------------

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
