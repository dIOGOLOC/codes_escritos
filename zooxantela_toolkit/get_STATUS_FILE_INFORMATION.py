#!/usr/bin/python -u

'''
--------------------------------------------------------------------------------
      Getting and Plotting information from status files (Guralp format)
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will retrieve and plot information from status files recorded by a
Güralp digitiser.


More information in:
https://www.guralp.com/documents/DM24.pdf


Inputs:
Status digitiser files (m8,m9,m0 e me folders).
Description:
• 00 is the digitiser status stream (notice no sample rate),
• M8, M9, MA are sensor mass positions for Z, N, E channels
• MB, ME, MF are three of the optional eight 16bit channels available

Outputs:
Images of the status file according to the date (format: PDF)

Examples of Usage (in command line):
   >> python get_STATUS_FILE_INFORMATION.py

'''

import time
from tqdm import tqdm
from multiprocessing import Pool
from obspy import read,read_inventory, UTCDateTime, Stream
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter
import matplotlib.dates as mdates
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt


from assessment_py.status_assessment import filelist,get_status_file_GURALP,list_split_day


from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,DIR_STATUS,NUM_PROCESS,LABEL_LANG
					)


# ============================
# Retrieving GCF status files
# ============================
if LABEL_LANG == 'br':
    print('Obtendo arquivos de funcionamento GCF.')
    print('\n')

else:
    print('\n')
    print('Getting GCF status files.')
    print('\n')


files = filelist(DIR_STATUS)

m8_folder = [i for i in files if "m8" in i]
m9_folder = [i for i in files if "m9" in i]
ma_folder = [i for i in files if "ma" in i]
me_folder = [i for i in files if "me" in i]


# ==================================
# Separating GCF status files by day
# ==================================
if LABEL_LANG == 'br':
    print('Separando os arquivos por dia.')
    print('\n')

else:
    print('Separating GCF status files by day.')
    print('\n')


daily_lst_m8 = list_split_day(m8_folder)
daily_lst_m9 = list_split_day(m9_folder)
daily_lst_ma = list_split_day(ma_folder)
daily_lst_me = list_split_day(me_folder)


# ================================
# Multiprocessing GCF status files
# ================================
if LABEL_LANG == 'br':
    print('Multiprocessamento dos arquivos de funcionamento GCF.')
    print('\n')

else:
    print('Multiprocessing GCF status files.')
    print('\n')

start_time = time.time()
#--------------------------------------------------------------------------------------------------------------------
if LABEL_LANG == 'br':
    print('Canal: m8')

else:
    print('Channel: m8')
    print('\n')

result_m8 = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_m8), total=len(daily_lst_m8)):
    result_m8.append(result)
#--------------------------------------------------------------------------------------------------------------------

if LABEL_LANG == 'br':
    print('Canal: m9')

else:
    print('Channel: m9')
    print('\n')

result_m9 = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_m9), total=len(daily_lst_m9)):
    result_m9.append(result)
#--------------------------------------------------------------------------------------------------------------------

if LABEL_LANG == 'br':
    print('Canal: ma')

else:
    print('Channel: ma')
    print('\n')

result_ma = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_ma), total=len(daily_lst_ma)):
    result_ma.append(result)
#--------------------------------------------------------------------------------------------------------------------
if LABEL_LANG == 'br':
    print('Canal: me')

else:
    print('Channel: me')
    print('\n')

result_me = []
pool = Pool(processes=NUM_PROCESS)
for result in tqdm(pool.imap_unordered(func=get_status_file_GURALP, iterable=daily_lst_me), total=len(daily_lst_me)):
    result_me.append(result)
#--------------------------------------------------------------------------------------------------------------------

print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')

# ====================
# Plotting status data
# ====================
print('Plotting status data')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_m8:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('Z mass position (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_Z_mass_position.png', dpi=300, facecolor='w')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_m9:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('N mass position (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_N_mass_position.png', dpi=300, facecolor='w')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_ma:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('E mass position (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_E_mass_position.png', dpi=300, facecolor='w')

#--------------------------------------------------------------------------------------------------------------------
days = DayLocator()   # every year
months = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y-%m-%d')

fig, ax = plt.subplots(figsize=(15,7))
for i in result_me:
	try:
		ax.plot(UTCDateTime(i[2]).matplotlib_date, i[3],'ok')
	except:
		pass
ax.set_ylabel('Temperature (microvolts)')
y_lim = ax.get_ylim()
ax.set_ylim(y_lim[0],y_lim[1])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(days)

ax.grid(True)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()
fig.savefig(OUTPUT_FIGURE_DIR+result_me[0][0]+'_temperature.png', dpi=300, facecolor='w')
