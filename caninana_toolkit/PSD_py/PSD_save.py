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
                    OUTPUT_XML_FILE_DIR,OUTPUT_JSON_FILE_DIR,OUTPUT_PSD_DIR
                    
				   )

print('Importing XML file')

inv = read_inventory(OUTPUT_XML_FILE_DIR+NETWORK_CODE+".xml")
print(inv)

def calc_PSD(data,sta_name):
    os.chdir(data)
    st = read('*')
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
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)

            ppsd = PPSD(l.stats, metadata=inv)
            ppsd.add(st) 
            os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
            print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

        if sta_channel == 'HHY':
            l.stats.channel = 'HHN'
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)

            ppsd = PPSD(l.stats, metadata=inv)
            ppsd.add(st) 
            os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
            print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

        elif sta_channel == 'HH1j':
            l.stats.channel = 'HHN'
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)
        
            ppsd = PPSD(l.stats, metadata=inv)
            ppsd.add(st) 
            os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
            print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')


        elif sta_channel == 'HH2':
            l.stats.channel = 'HHE'
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)
        
            try:
                ppsd = PPSD(l.stats, metadata=inv)
                ppsd.add(st) 
                os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
                print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
                ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            except:
                pass

        elif sta_channel == 'HH2j':
            l.stats.channel = 'HHE'
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)
         
            ppsd = PPSD(l.stats, metadata=inv)
            ppsd.add(st) 
            os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
            print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')


        elif sta_channel == 'HHX':
            l.stats.channel = 'HHE'
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)
        
            ppsd = PPSD(l.stats, metadata=inv)
            ppsd.add(st) 
            os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
            print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')


        else:
            print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel)
        
            ppsd = PPSD(l.stats, metadata=inv)
            ppsd.add(st) 
            os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/',exist_ok=True)
            print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
            ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+l.stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+l.stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')