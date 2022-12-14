'''
--------------------------------------------------------------------------------
  Function trim local event data 
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 11/2022

'''

import os
import glob
from tarfile import StreamError
import obspy as op
from obspy import Stream, Inventory
import pandas as pd
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.io.sac.sactrace import SACTrace

DIR_DATA = '/home/diogoloc/dados_posdoc/ON02_analysis/'

CUT_BEFORE_P_LOCAL = 300

CUT_AFTER_P_LOCAL = 600

XML_station_dir = '/home/diogoloc/dados_posdoc/ON02_analysis/XML/'

events_files_dir = '/home/diogoloc/dados_posdoc/ON02_analysis/Events_info/'

# --------- #
# Functions #
# --------- #

def get_event_data(csv_file):
    event_data = pd.read_csv(csv_file, delimiter=';', usecols=[1,2,3,4,5])
    
    return event_data

def get_station_data(xml_file):
    station_data = op.read_inventory(xml_file)
    
    return station_data


def cut_data_by_local_event(knetwk,kstnm,stla,stlo,evla,evlo,evdp,evmag,ev_timeUTC,inventory):
    data_sta = DIR_DATA+'Data/'+str(op.UTCDateTime(ev_timeUTC).year)+'/'+knetwk+'/'+kstnm+'/'

    if os.path.isdir(data_sta) == True:

        #Calculating distance, azimuth and backazimuth
        dist,az,baz = op.geodetics.gps2dist_azimuth(evla,evlo,stla,stlo)
        gcarc = op.geodetics.kilometer2degrees(dist/1000)

        #Calculating ray parameter
        model = TauPyModel(model='iasp91')
        arrivals = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=gcarc)
        arr = arrivals[0]

        #Reference time
        starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P_LOCAL
        endtime = (op.UTCDateTime(ev_timeUTC)+arr.time)+CUT_AFTER_P_LOCAL

        ########################################################################################################################################################
        #STREAM

        for i in glob.glob(data_sta+'*'):

            #Creating Event Directory
            OUTPUT_EV_DIR = DIR_DATA+'Events_data/'
            event_directory = OUTPUT_EV_DIR+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'/'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'/'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'/'
            os.makedirs(event_directory, exist_ok=True)
                  
            #-----------------------------------
            #Component E
            if 'HE.D' in i:
                channel = i.split('/')[-1]
                HE_FILE = glob.glob(data_sta+channel+'/*'+channel+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))[0]
                    
                stE = op.read(HE_FILE)
                stE.trim(starttime,endtime)
                inv = op.read_inventory(inventory)
                stE.remove_response(inventory=inv,output='DISP')

                headerHHE = {
                                'kstnm': kstnm, 'kcmpnm': channel.split('.')[0],'knetwk':knetwk,
                                'stla': float(stla), 'stlo': float(stlo),
                                'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag),
                                'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 'nzyear': int(starttime.year),
                                'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc),'az': float(az), 'baz': float(baz),
                                'o':float(CUT_BEFORE_P_LOCAL),
                                'delta':stE[0].stats.delta
                            }

                sacHHE = SACTrace(data=stE[0].data, **headerHHE)
                sacHHE.write(event_directory+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.E')

            #-----------------------------------
            #Component N
            if 'HN.D' in i:
                channel = i.split('/')[-1]
                HN_FILE = glob.glob(data_sta+channel+'/*'+channel+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))[0]

                stN = op.read(HN_FILE)
                stN.trim(starttime,endtime)
                inv = op.read_inventory(inventory)
                stN.remove_response(inventory=inv,output='DISP')  

                headerHHN = {
                                'kstnm': kstnm, 'kcmpnm': channel.split('.')[0],'knetwk':knetwk,
                                'stla': float(stla), 'stlo': float(stlo),
                                'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag),
                                'nzhour': int(starttime.hour), 'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
                                'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
                                'o':float(CUT_BEFORE_P_LOCAL),
                                'delta':stN[0].stats.delta
                            }

                sacHHN = SACTrace(data=stN[0].data, **headerHHN)
                sacHHN.write(event_directory+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.N')

            #-----------------------------------
            #Component Z
            if 'HZ.D' in i:
                channel = i.split('/')[-1]
                HZ_FILE = glob.glob(data_sta+channel+'/*'+channel+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday))[0]

                stZ = op.read(HZ_FILE)
                stZ.trim(starttime,endtime)
                inv = op.read_inventory(inventory)
                stZ.remove_response(inventory=inv,output='DISP')

                headerHHZ = {
                            'kstnm': kstnm, 'kcmpnm': channel.split('.')[0],'knetwk':knetwk,
                            'stla': float(stla), 'stlo': float(stlo),
                            'evdp': float(evdp), 'evla': float(evla), 'evlo': float(evlo), 'mag': float(evmag),
                            'nzhour': int(starttime.hour),'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),'nzyear': int(starttime.year),
                            'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 'az': float(az), 'baz': float(baz),
                            'o':float(CUT_BEFORE_P_LOCAL),
                            'delta':stZ[0].stats.delta
                            }

                sacHHZ = SACTrace(data=stZ[0].data, **headerHHZ)
                sacHHZ.write(event_directory+knetwk+'.'+kstnm+'.'+'{:04}'.format(op.UTCDateTime(ev_timeUTC).year)+'.'+'{:03}'.format(op.UTCDateTime(ev_timeUTC).julday)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).hour)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).minute)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).second)+'.'+'{:02}'.format(op.UTCDateTime(ev_timeUTC).microsecond)[:3]+'.Z')


def create_snuffler_event_file(directory_events):
    #-----------------------------------
    # Reading events files
    for i in glob.glob(directory_events+'*/*/*'):
        stream_all = Stream()
        for st_file in glob.glob(i+'/*'):
            st = op.read(st_file)
            stream_all+=st 

        dir_name = i.split('/')[-1]
        # Creating Event Snuffler Directory
        OUTPUT_EV_DIR_SNU = DIR_DATA+'Events_Snuffler/'
        event_directory_snu = OUTPUT_EV_DIR_SNU+dir_name
        os.makedirs(event_directory_snu, exist_ok=True)

        stream_all.write(event_directory_snu+'/'+dir_name+'.mseed',format="MSEED")

        evla = st[0].stats.sac.evla
        evlo = st[0].stats.sac.evlo
        evmag = st[0].stats.sac.mag
        ev_timeUTC = st[0].stats.starttime+st[0].stats.sac.o

        event_dic = {'name':dir_name,
                    'time':op.UTCDateTime(ev_timeUTC).strftime('%Y-%m-%d %H:%M:%S.ffffff'),
                    'latitude':str(evla),
                    'longitude':str(evlo),
                    'magnitude':str(evmag)}

        fout = event_directory_snu+'/'+dir_name+'.txt'
        fo = open(fout, "w")

        for k, v in event_dic.items():
            fo.write(str(k) + ' = '+ str(v) + '\n')

        fo.close()


# ------------ # 
# Main Program #
# ------------ # 

# Slicing data according to events time:

# importing stations table ...

XML_files = glob.glob(XML_station_dir+'*.xml')

for XML_file in XML_files:
    inventory = get_station_data(XML_file)
    knetwk = inventory.get_contents()['networks'][0]
    kstnm = inventory.get_contents()['channels'][0].split('.')[1]
    stlo = inventory.get_coordinates(inventory.get_contents()['channels'][0])['longitude']
    stla = inventory.get_coordinates(inventory.get_contents()['channels'][0])['latitude']
 
    # importing events table ...
    events_files_lst = glob.glob(events_files_dir+'*.csv')
    for event_file in events_files_lst:
        events_table = get_event_data(event_file)

        for i in events_table.iterrows():

            ev_timeUTC = i[1]['origin']
            evlo = i[1]['longitude']
            evla = i[1]['latitude']
            evdp = i[1]['depth']
            evmag = i[1]['magnitude']

            cut_data_by_local_event(knetwk,kstnm,stla,stlo,evla,evlo,evdp,evmag,ev_timeUTC,XML_file)

# Creating MSEED to use in SNUFFLER:

create_snuffler_event_file(DIR_DATA+'Events_data/')