'''
Script to plot PPSD data based in 
https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html
'''

import os
import glob
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from obspy import UTCDateTime

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,INITIAL_DATE,FINAL_DATE,TIME_OF_WEEKDAY_DAY,TIME_OF_WEEKDAY_START_HOUR,TIME_OF_WEEKDAY_FINAL_HOUR
				   )

# ==================================
# Function to plot TOTAL PPSD DATA
# ==================================

def plot_PPSD_TOTAL_data(date_lst):
    os.chdir(date_lst)
    files = sorted(glob.glob('*.npz'))
    ppsd = PPSD.load_npz(files[0])
    [ppsd.add_npz(i) for i in files[1:]]
    os.makedirs(OUTPUT_FIGURE_DIR+'TOTAL/'+ppsd.station+'/',exist_ok=True)
    ppsd.plot(cmap=pqlx,filename=OUTPUT_FIGURE_DIR+'TOTAL/'+ppsd.station+'/'+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')

def plot_PPSD_WINDOWED_data(date_lst):
    os.chdir(date_lst)
    files = sorted(glob.glob('*.npz'))
    print(files)
    ppsd = PPSD.load_npz(files[0])
    print(ppsd) 

    [ppsd.add_npz(i) for i in files[1:]]
    ppsd.calculate_histogram(starttime=UTCDateTime(INITIAL_DATE),endtime=UTCDateTime(FINAL_DATE),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])    
    folder_output = OUTPUT_FIGURE_DIR+'WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd.station+'/'
    os.makedirs(folder_output,exist_ok=True)
    ppsd.plot(cmap=pqlx,filename=folder_output+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')
