#!/usr/bin/python -u
"""
Scritp to estimate probabilistic power spectral densities for each .SAC file
"""

import time
import os
import glob
from tqdm import tqdm
from multiprocessing import Pool
from obspy import read,read_inventory, UTCDateTime, Stream
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter
import matplotlib.dates as mdates
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt


from assessment_py.power_spectral_densities import plot_PPSD_by_period_sensor,plot_PPSD_by_period_hydrophone

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,DIR_STATUS,NUM_PROCESS,OUTPUT_PSD_DIR
					)


# =====================
# Retrieving .NPZ files
# =====================

print('\n')
print('Retrieving .NPZ files')
print('\n')

dirs_PPSD = glob.glob(OUTPUT_PSD_DIR+'*')

NPZ_dirs = sorted([i+'/' for i in dirs_PPSD])

# ===================
# Ploting .NPZ files
# ===================

start_time = time.time()
print('\n')
#--------------------------------------------------------------------------------------------------------------------
for i,j in enumerate(NPZ_dirs):
  plot_PPSD_by_period_sensor(j)
print('\n')
#--------------------------------------------------------------------------------------------------------------------
for i,j in enumerate(NPZ_dirs):
  plot_PPSD_by_period_hydrophone(j)
print('\n')
#--------------------------------------------------------------------------------------------------------------------


print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')