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


from assessment_py.power_spectral_densities import calc_PSD

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,DIR_STATUS,NUM_PROCESS,OUTPUT_PSD_DIR
					)


# =====================
# Retrieving SAC files
# =====================

print('\n')
print('Retrieving SAC files')
print('\n')

HHZ_files = []
HHE_files = []
HHN_files = []
HHX_files = []

for root, dirs, files in os.walk(DIR_DATA):
    for file in files:
        if "HHZ" in file:
            HHZ_files.append(os.path.join(root, file))
        elif "HHE" in file:
        	HHE_files.append(os.path.join(root, file))
       	elif "HHN" in file:
       		HHN_files.append(os.path.join(root, file))
        elif "HHX" in file:
          HHX_files.append(os.path.join(root, file))

HHZ_files = sorted(HHZ_files)
HHE_files = sorted(HHE_files)
HHN_files = sorted(HHN_files)
HHX_files = sorted(HHX_files)

# =========================
# Multiprocessing SAC files
# =========================

print('Multiprocessing PPSD')
start_time = time.time()
print('\n')
#--------------------------------------------------------------------------------------------------------------------
'''
print('Channel: HHZ')
result_HHZ = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=calc_PSD, iterable=HHZ_files), total=len(HHZ_files)):
	result_HHZ.append(result)
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHE')
result_HHE = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=calc_PSD, iterable=HHE_files), total=len(HHE_files)):
	result_HHE.append(result)
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print('Channel: HHN')
result_HHN = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=calc_PSD, iterable=HHN_files), total=len(HHN_files)):
	result_HHN.append(result)
print('\n')
#--------------------------------------------------------------------------------------------------------------------
'''
print('Channel: HHX')
result_HHX = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=calc_PSD, iterable=HHX_files), total=len(HHX_files)):
  result_HHX.append(result)
print('\n')
#--------------------------------------------------------------------------------------------------------------------

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
