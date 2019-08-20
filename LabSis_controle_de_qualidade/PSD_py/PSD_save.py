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
from obspy.clients.arclink.client import Client


from parameters_py.config import (
					DIR_DATA,OUTPUT_JSON_FILE_DIR,OUTPUT_PSD_DIR,INITIAL_DATE,FINAL_DATE,XML_FILE                    
				   )

def calc_PSD(data):
	print('Importing XML file')
	inv = read_inventory(XML_FILE)
	print(inv)

	for i,j in enumerate(data):
	    st = read(j)
	    st.merge()

	    time_data = st[0].stats.starttime
	    time_data_year = '{:04}'.format(time_data.year)
	    time_data_julday = '{:03}'.format(time_data.julday)
	    time_data_hour = '{:02}'.format(time_data.hour)
	    time_data_minute = '{:02}'.format(time_data.minute)

	    sta_name = st[0].stats.station
	    sta_channel = st[0].stats.channel
	    print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel+' / '+str(i+1)+' of '+str(len(data)))

	    ppsd = PPSD(st[0].stats, metadata=inv)
	    ppsd.add(st)

	    os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/',exist_ok=True)
	    print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/'+sta_name+'..'+sta_channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
	    ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/'+sta_name+'..'+sta_channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')


def calc_PSD_client(sta_name):
    
    client = Client(user=USER,host=HOST, port=PORT, institution=INSTITUTION)

    data_lista = np.arange(obspy.UTCDateTime(INITIAL_DATE),obspy.UTCDateTime(FINAL_DATE),86400) 

    for i,j in enumerate(data_lista):
	    st = client.get_waveforms(NETWORK_CODE,sta_name,LOCATION,CHANNEL_CODE,j,j+86400)
	    st.merge()

	    st[0].stats.station = sta_name
	    st[0].stats.network = NETWORK_CODE
	    st[0].stats.location = LOCATION
		
	    time_data = st[0].stats.starttime
	    time_data_year = '{:04}'.format(time_data.year)
	    time_data_julday = '{:03}'.format(time_data.julday)
	    time_data_hour = '{:02}'.format(time_data.hour)
	    time_data_minute = '{:02}'.format(time_data.minute)

	    sta_channel = st[0].stats.channel
	    print('Calculating PPSD: station: '+sta_name+' / channel: '+sta_channel+' / '+str(i+1)+' of '+str(len(data_lista)))
		
	    ppsd = PPSD(st[0].stats, metadata=inv)
	    ppsd.add(st) 
	    os.makedirs(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+st[0].stats.channel+'.PPSD'+'/',exist_ok=True)
	    print(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+st[0].stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+st[0].stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')
	    ppsd.save_npz(OUTPUT_PSD_DIR+time_data_year+'/'+sta_name+'/'+st[0].stats.channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+st[0].stats.channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')