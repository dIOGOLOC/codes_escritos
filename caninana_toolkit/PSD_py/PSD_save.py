'''
Script to estimate probabilistic power spectral densities for 
one combination of network/station/location/channel/sampling_rate.
(https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html)

Calculations are based on the routine used by [McNamara2004]:
McNamara, D. E. and Buland, R. P. (2004),
Ambient Noise Levels in the Continental United States,
Bulletin of the Seismological Society of America, 94 (4), 1517-1527.
http://www.bssaonline.org/content/94/4/1517.abstract. 


For information on New High/Low Noise Model see [Peterson1993]:
Peterson, J. (1993),
Observations and Modeling of Seismic Background Noise,
U.S. Geological Survey open-file report 93-322, Albuquerque, N.M.
http://ehp3-earthquake.wr.usgs.gov/regional/asl/pubs/files/ofr93-322.pdf
'''


import obspy 
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, Stream,read_inventory
import os
import glob
from obspy.io.xseed import Parser
from obspy import read_inventory
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx

from parameters_py.config import (
					DIR_DATA,SOURCE,NETWORK_CODE,NETWORK_DESCRIPTION,START_DATE,SAMPLING_RATE,LOCATION,
                    OUTPUT_XML_FILE_DIR,OUTPUT_JSON_FILE_DIR,OUTPUT_PSD_DIR,
                    TIME_OF_WEEKDAY_DAY,TIME_OF_WEEKDAY_START_HOUR,TIME_OF_WEEKDAY_FINAL_HOUR
                    
				   )

print('Importing XML file')

inv = read_inventory(OUTPUT_XML_FILE_DIR+NETWORK_CODE+".xml")
print(inv)

def calc_PSD(data,sta_name):
    st = read(data)
    st.merge()
    for k,l in enumerate(st):
        l.stats.station = sta_name
        l.stats.network = NETWORK_CODE
        l.stats.location = LOCATION
        
        time_data = l.stats.starttime
        time_data_year = '{:04}'.format(time_data.year)
        time_data_julday = '{:03}'.format(time_data.julday)
        time_data_hour = '{:02}'.format(time_data.hour)
        time_data_minute = '{:02}'.format(time_data.minute)

        sta_channel = l.stats.channel
        if sta_channel == 'HH1':
            l.stats.channel = 'HHN'

            if os.path.isfile(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz'):
                pass
            else:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
                #print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

        if sta_channel == 'HHY':
            l.stats.channel = 'HHN'

            if os.path.isfile(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz'):
                pass
            else:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
                #print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

        elif sta_channel == 'HH2':
            l.stats.channel = 'HHE'

            if os.path.isfile(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz'):
                pass
            else:            
                try:
                    ppsd = PPSD(l.stats, metadata=inv)
                    ppsd.add(st) 
                    os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
                    #print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                    ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                except:
                    pass

        elif sta_channel == 'HHX':
            l.stats.channel = 'HHE'

            if os.path.isfile(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz'):
                pass
            else:            
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
                #print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

        else:
            if os.path.isfile(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz'):
                pass
            else:            
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
                #print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')