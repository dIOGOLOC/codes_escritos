'''
--------------------------------------------------------------------------------
 Function to trim/plot local the dataset according to local events time
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 08/2023


Description:
This code will trim and plot the local datase according to a given an event time
and a list of stations.

More information in:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.trim.html


Inputs:
JSON file with event description:
    ev_timeUTC: event time in UTC (str)
    ev_year: year of the event
    ev_month: month of the event
    ev_day: day of the event
    ev_julday: julian day of the event
    ev_hour: hour of the event
    ev_minute: minute of the event
    ev_second: second of the event
    ev_microsecond: microsecond of the event
    evla: latitude of the event
    evlo: longitude of the event
    evdp: depth of the event
    mag: magnitude of the event

'''


import os
import glob
import obspy as op
import json
from tqdm import tqdm
from multiprocessing import Pool

from obspy import read,read_inventory, UTCDateTime, Stream
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace
from obspy.signal.trigger import classic_sta_lta, trigger_onset, coincidence_trigger,recursive_sta_lta,plot_trigger

import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter,LatitudeLocator,LongitudeLocator


# =====
# paths
# =====

# directory of raw files
DIR_DATA = '/home/sysop/dados_posdoc/MTZ_2024/DATA_MTZ/DATA_2024_NO_PP_FILTER/'

#Directory to save EVENT data
OUTPUT_EV_DIR = '//home/sysop/dados_posdoc/MTZ_2024/DATA_MTZ/'

# =====
# event
# =====

#Taup_time model to calculate travel times
TAUPY_MODEL = 'iasp91'

#Minimum event distance
EV_GCARC_MIN = 30

#Maximum event distance
EV_GCARC_MAX = 90   

#Minimum event magnitude
EV_MAGNITUDE_MB = 5.5

#Time to trim the seismogram before P-wave arrival
CUT_BEFORE_P = 10

#Time to trim the seismogram after P-wave arrival
CUT_AFTER_P = 260

# ===============================
# Function to cut and plot event:
# ===============================

def cut_data_by_event(input):
	
    ########################################################################################################################################################
    #STREAM
    #-----------------------------------
    #Component E

    ev_files = sorted(glob.glob(input+'/*'))

    if len(ev_files) == 4:

        stE = read(ev_files[0])
        
        #Calculating ray parameter
        model = TauPyModel(model=TAUPY_MODEL)
        arrivals = model.get_travel_times(source_depth_in_km=stE[0].stats.sac.evdp, distance_in_degree=stE[0].stats.sac.gcarc, phase_list=["P"])
        arr = arrivals[0]

        arrivalsPP =  model.get_travel_times(source_depth_in_km=stE[0].stats.sac.evdp, distance_in_degree=stE[0].stats.sac.gcarc, phase_list=["PP"])
        arrPP = arrivalsPP[0]	

        time_PP = arrPP.time - arr.time

        #Creating Event Directory
        event_directory = OUTPUT_EV_DIR+'DATA_2024_YES_PP_FILTER/'+'/'.join(ev_files[0].split('/DATA_2024_NO_PP_FILTER/')[-1].split('/')[:-1])+'/'
        os.makedirs(event_directory, exist_ok=True)
        
        data_stE = []
        for k,l in enumerate(stE[0].times()):
            if  l < time_PP - 10:
                data_stE.append(stE[0].data[k])
            elif time_PP - 10.0 <= l <= time_PP - 5.0:
                data_stE.append(stE[0].data[k] * np.exp(-0.2*(l - time_PP + 10.0)))
            elif time_PP - 5.0 <= l <= time_PP + 5.0: 
                data_stE.append(stE[0].data[k] * 0.0)
            elif time_PP + 5.0 <= l <= time_PP + 10.0: 
                data_stE.append(stE[0].data[k] * np.exp(-0.2*(-l + time_PP + 10.0)))
            elif l > time_PP + 10.0: 
                data_stE.append(stE[0].data[k])
                        
        headerHHE = {
                    'kstnm': stE[0].stats.sac.kstnm, 'kcmpnm': 'HHE','knetwk':'BP',
                    'stla': stE[0].stats.sac.stla, 'stlo': stE[0].stats.sac.stlo,
                    'evdp': float(stE[0].stats.sac.evdp), 'evla': float(stE[0].stats.sac.evla), 'evlo': float(stE[0].stats.sac.evlo), 'mag': float(stE[0].stats.sac.mag),
                    'nzhour': int(stE[0].stats.sac.nzhour), 'nzjday': int(stE[0].stats.sac.nzjday), 'nzmin': int(stE[0].stats.sac.nzmin),'nzmsec': int(stE[0].stats.sac.nzmsec),
                    'nzsec': int(stE[0].stats.sac.nzsec), 'nzyear': int(stE[0].stats.sac.nzyear),'cmpaz': 90.0, 'cmpinc': 90.0,
                    'dist': float(stE[0].stats.sac.dist), 'gcarc': float(stE[0].stats.sac.gcarc),'az': float(stE[0].stats.sac.az), 'baz': float(stE[0].stats.sac.baz),
                    'b':float(-CUT_BEFORE_P),'a':0,'user0':float(arr.ray_param/6371),'delta':stE[0].stats.delta
                    }
        
        sacHHE = SACTrace(data=np.array(data_stE), **headerHHE)
        sacHHE.write(event_directory+ev_files[0].split('/')[-1])

        #-----------------------------------
        #Component N

        stN = read(ev_files[1])

        data_stN = []
        for k,l in enumerate(stN[0].times()):
            if  l < time_PP - 10:
                data_stN.append(stN[0].data[k])
            elif time_PP - 10.0 <= l <= time_PP - 5.0:
                data_stN.append(stN[0].data[k] * np.exp(-0.2*(l - time_PP + 10.0)))
            elif time_PP - 5.0 <= l <= time_PP + 5.0: 
                data_stN.append(stN[0].data[k] * 0.0)
            elif time_PP + 5.0 <= l <= time_PP + 10.0: 
                data_stN.append(stN[0].data[k] * np.exp(-0.2*(-l + time_PP + 10.0)))
            elif l > time_PP + 10.0: 
                data_stN.append(stN[0].data[k])

        headerHHN = {
                    'kstnm': stN[0].stats.sac.kstnm, 'kcmpnm': 'HHN','knetwk':'BP',
                    'stla': stN[0].stats.sac.stla, 'stlo': stN[0].stats.sac.stlo,
                    'evdp': float(stN[0].stats.sac.evdp), 'evla': float(stN[0].stats.sac.evla), 'evlo': float(stN[0].stats.sac.evlo), 'mag': float(stN[0].stats.sac.mag),
                    'nzhour': int(stN[0].stats.sac.nzhour), 'nzjday': int(stN[0].stats.sac.nzjday), 'nzmin': int(stN[0].stats.sac.nzmin),'nzmsec': int(stN[0].stats.sac.nzmsec),
                    'nzsec': int(stN[0].stats.sac.nzsec), 'nzyear': int(stN[0].stats.sac.nzyear),'cmpaz': 0.0, 'cmpinc': 90.0,
                    'dist': float(stN[0].stats.sac.dist), 'gcarc': float(stN[0].stats.sac.gcarc),'az': float(stN[0].stats.sac.az), 'baz': float(stN[0].stats.sac.baz),
                    'b':float(-CUT_BEFORE_P),'a':0,'user0':float(arr.ray_param/6371),'delta':stN[0].stats.delta
                    }
       
        sacHHN = SACTrace(data=np.array(data_stN), **headerHHN)
        sacHHN.write(event_directory+ev_files[1].split('/')[-1])
        
        #-----------------------------------
        #Component Z

        stZ = read(ev_files[2])

        data_stZ = []
        for k,l in enumerate(stZ[0].times()):
            if  l < time_PP - 10:
                data_stZ.append(stZ[0].data[k])
            elif time_PP - 10.0 <= l <= time_PP - 5.0:
                data_stZ.append(stZ[0].data[k] * np.exp(-0.2*(l - time_PP + 10.0)))
            elif time_PP - 5.0 <= l <= time_PP + 5.0: 
                data_stZ.append(stZ[0].data[k] * 0.0)
            elif time_PP + 5.0 <= l <= time_PP + 10.0: 
                data_stZ.append(stZ[0].data[k] * np.exp(-0.2*(-l + time_PP + 10.0)))
            elif l > time_PP + 10.0: 
                data_stZ.append(stZ[0].data[k])
        
        headerHHZ = {
                    'kstnm': stZ[0].stats.sac.kstnm, 'kcmpnm': 'HHZ','knetwk':'BP',
                    'stla': stZ[0].stats.sac.stla, 'stlo': stZ[0].stats.sac.stlo,
                    'evdp': float(stZ[0].stats.sac.evdp), 'evla': float(stZ[0].stats.sac.evla), 'evlo': float(stZ[0].stats.sac.evlo), 'mag': float(stZ[0].stats.sac.mag),
                    'nzhour': int(stZ[0].stats.sac.nzhour), 'nzjday': int(stZ[0].stats.sac.nzjday), 'nzmin': int(stZ[0].stats.sac.nzmin),'nzmsec': int(stZ[0].stats.sac.nzmsec),
                    'nzsec': int(stZ[0].stats.sac.nzsec), 'nzyear': int(stZ[0].stats.sac.nzyear),'cmpaz': 0.0, 'cmpinc': 0.0,
                    'dist': float(stZ[0].stats.sac.dist), 'gcarc': float(stZ[0].stats.sac.gcarc),'az': float(stZ[0].stats.sac.az), 'baz': float(stZ[0].stats.sac.baz),
                    'b':float(-CUT_BEFORE_P),'a':0,'user0':float(arr.ray_param/6371),'delta':stZ[0].stats.delta
                    }
       
        sacHHZ = SACTrace(data=np.array(data_stZ), **headerHHZ)
        sacHHZ.write(event_directory+ev_files[2].split('/')[-1])

    #----------------------------------
    ########################################################################################################################################################

# ==================================================
#  Importing Local Event dictionary from JSON file
# ==================================================

print('\n')
print('Looking for Events data')
print('\n')

ev_folders = glob.glob(DIR_DATA+'*/*/*/*')

# ================
# trim event data
# ================

print('PP phase filter for each station')
print('\n')

with Pool(processes=20) as p:
	max_ = len(ev_folders)
	with tqdm(total=max_,desc='Filtering') as pbar:
		for i, _ in enumerate(p.imap_unordered(cut_data_by_event,ev_folders)):
			pbar.update()

print('Finished!')