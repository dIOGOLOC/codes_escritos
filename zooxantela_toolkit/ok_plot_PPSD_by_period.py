#!/usr/bin/python -u

'''
--------------------------------------------------------------------------------
   Plotting probabilistic power spectral densities for each daily .SAC file
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will plot probabilistic power spectral densities amplitude filtered for each
channel according to a specific period in seconds.


More information in:
https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html

Calculations are based on the routine used by [McNamara2004]:
McNamara, D. E. and Buland, R. P. (2004), Ambient Noise Levels in the Continental
United States, Bulletin of the Seismological Society of America, 94 (4), 1517-1527.
http://www.bssaonline.org/content/94/4/1517.abstract.

For information on New High/Low Noise Model see [Peterson1993]:
Peterson, J. (1993), Observations and Modeling of Seismic Background Noise,
U.S. Geological Survey open-file report 93-322, Albuquerque, N.M.
http://ehp3-earthquake.wr.usgs.gov/regional/asl/pubs/files/ofr93-322.pdf


Inputs:
Daily PPSD file

Outputs:
Images of the PSD amplitude according to the period (format: PDF)

Examples of Usage (in command line):
   >> python plot_PPSD_by_period.py

'''

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
