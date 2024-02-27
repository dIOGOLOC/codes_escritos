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
DIR_DATA = '/medata01/SEISCOMP_DATA/'

# Station location file path
STA_LOC_PATH = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/MTZ_2024/STA_COORD/cood_rede_RSBR_MTZ_new.txt'

# Event list file path
EV_LOC_PATH = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/MTZ_2024/EVENTS/neic_lst_2016_2024.csv'

#Directory to save EVENT data
OUTPUT_EV_DIR = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/MTZ_2024/DATA_MTZ/REDE_RSBR/event_data_raw/'

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
	
    knetwk = input[0]
    kstnm = input[1]
    stla = input[2]
    stlo = input[3]
    ev_timeUTC = input[4]
    ev_lat = input[5]
    ev_long = input[6]
    ev_depth = input[7]
    ev_mag = input[8]

    #Calculating distance, azimuth and backazimuth
    dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
    gcarc = op.geodetics.kilometer2degrees(dist/1000)
    if EV_GCARC_MIN <= gcarc <= EV_GCARC_MAX:

        #Calculating ray parameter
        model = TauPyModel(model=TAUPY_MODEL)
        arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
        arr = arrivals[0]

        #Reference time
        starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P
        endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P

        ########################################################################################################################################################
        #STREAM
        try:
            #-----------------------------------
            #Component E

            if os.path.isfile(DIR_DATA+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/HHE.D'+'/'+knetwk+'.'+kstnm+'..HHE.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:
                #Creating Event Directory
                event_directory = OUTPUT_EV_DIR+'NO_FILTER_PP/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
                os.makedirs(event_directory, exist_ok=True)

                stE = read(DIR_DATA+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/HHE.D'+'/'+knetwk+'.'+kstnm+'..HHE.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
                stE.trim(starttime,endtime)

                headerHHE = {
                            'kstnm': kstnm, 'kcmpnm': 'HHE','knetwk':knetwk,
                            'stla': float(stla), 'stlo': float(stlo),
                            'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
                            'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 'nzyear': int(starttime.year),
                            'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc),'az': float(az), 'baz': float(baz),
                            'b':float(-CUT_BEFORE_P),'a':0,'user0':float(arr.ray_param/6371),
                            'delta':stE[0].stats.delta
                            }

                sacHHE = SACTrace(data=stE[0].data, **headerHHE)
                sacHHE.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.E')

            #-----------------------------------
            #Component N
            if os.path.isfile(DIR_DATA+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/HHN.D'+'/'+knetwk+'.'+kstnm+'..HHN.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:

                #Creating Event Directory
                event_directory = OUTPUT_EV_DIR+'NO_FILTER_PP/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
                os.makedirs(event_directory, exist_ok=True)

                stN = read(DIR_DATA+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/HHN.D'+'/'+knetwk+'.'+kstnm+'..HHN.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
                stN.trim(starttime,endtime)

                headerHHY = {
                            'kstnm': kstnm, 'kcmpnm': 'HHN','knetwk':knetwk,
                            'stla': float(stla), 'stlo': float(stlo),
                            'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
                            'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
                            'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
                            'b':float(-CUT_BEFORE_P),'a':0,'user0':float(arr.ray_param/6371),

                            'delta':stN[0].stats.delta
                            }

                sacHHY = SACTrace(data=stN[0].data, **headerHHY)
                sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.N')

            #-----------------------------------
            #Component Z

            if os.path.isfile(DIR_DATA+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/HHZ.D'+'/'+knetwk+'.'+kstnm+'..HHZ.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)) == True:

                #Creating Event Directory
                event_directory = OUTPUT_EV_DIR+'NO_FILTER_PP/'+knetwk+'/'+kstnm+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]
                os.makedirs(event_directory, exist_ok=True)

                stZ = read(DIR_DATA+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/HHZ.D'+'/'+knetwk+'.'+kstnm+'..HHZ.D.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))
                stZ.trim(starttime,endtime)

                headerHHZ = {
                            'kstnm': kstnm, 'kcmpnm': 'HHZ','knetwk':knetwk,
                            'stla': float(stla), 'stlo': float(stlo),
                            'evdp': float(ev_depth), 'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag),
                            'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
                            'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
                            'b':float(-CUT_BEFORE_P),'a':0,'user0':float(arr.ray_param/6371),
                            'delta':stZ[0].stats.delta
                            }

                sacHHZ = SACTrace(data=stZ[0].data, **headerHHZ)
                sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.Z')
        except:
            pass
        
        #----------------------------------

# ============================================
#  Importing station dictionary from JSON file
# ============================================

print('\n')
print('STATIONS data')
print('\n')

sta_loc_info = np.genfromtxt(STA_LOC_PATH,delimiter=';',skip_header=1,dtype='str')

knetwk = sta_loc_info[:,0]
kstnm = sta_loc_info[:,1]
stla = [float(i) for i in sta_loc_info[:,2]]
stlo = [float(i) for i in sta_loc_info[:,3]]
stel = [float(i) for i in sta_loc_info[:,4]]

for i in kstnm:
	print('Station = ',i)
print('\n')

# ==================================================
#  Importing Local Event dictionary from JSON file
# ==================================================

print('\n')
print('Looking for Events data')
print('\n')

filename_STA = EV_LOC_PATH

event_dic = np.genfromtxt(EV_LOC_PATH,delimiter=',',skip_header=1,usecols=(0,1,2,3,4),dtype='str')

ev_timeUTC = event_dic[:,0]
evla = [float(i) for i in event_dic[:,1]]
evlo = [float(i) for i in event_dic[:,2]]
evdp = [float(i) for i in event_dic[:,3]]
mag = [float(i) for i in event_dic[:,4]]

print('Number of events = '+str(len(mag)))
print('\n')

# ================
# trim event data
# ================

print('Trim local event data for each station')
print('\n')

for i,j in enumerate(kstnm):

	with Pool(processes=8) as p:
		max_ = len(mag)
		with tqdm(total=max_,desc='station: '+j) as pbar:
			for i, _ in enumerate(p.imap_unordered(cut_data_by_event, [[knetwk[i],kstnm[i],stla[i],stlo[i],ev_timeUTC[k],evla[k],evlo[k],evdp[k],mag[k]] for k,l in tqdm(enumerate(mag),total=len(mag))])):
				pbar.update()

print('Finished!')