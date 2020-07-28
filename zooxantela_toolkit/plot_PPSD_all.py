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


from assessment_py.power_spectral_densities import plot_PSD

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,DIR_STATUS,NUM_PROCESS,OUTPUT_PSD_DIR
					)


# =====================
# Retrieving .NPZ files
# =====================

print('\n')
print('Retrieving SAC files')
print('\n')

HHZ_dir = []
HHE_dir = []
HHN_dir = []

for root, dirs, files in os.walk(OUTPUT_PSD_DIR):
    for directory in dirs:
        if "HHZ.PPSD" in directory:
            HHZ_dir.append(os.path.join(root, directory))
        if "HHE.PPSD" in directory:
        	HHE_dir.append(os.path.join(root, directory))
        if "HHN.PPSD" in directory:
       		HHN_dir.append(os.path.join(root, directory))


# ===================
# Ploting .NPZ files
# ===================

print('Ploting .NPZ files')
start_time = time.time()
print('\n')
#--------------------------------------------------------------------------------------------------------------------
print('Channel: HHZ')
for i,j in enumerate(HHZ_dir):
  plot_PSD(j)
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHE')
for i,j in enumerate(HHE_dir):
  plot_PSD(j)
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHN')
for i,j in enumerate(HHN_dir):
  plot_PSD(j)
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')