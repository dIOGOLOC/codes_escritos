'''
Script to create a custom StationXML file with ObsPy based in
https://docs.obspy.org/master/tutorial/code_snippets/stationxml_file_from_scratch.html.

Instrument Response can be looked up and attached to the channels from the 
IRIS DMC Library of Nominal Responses for Seismic Instruments (http://ds.iris.edu/NRL/)
using ObsPy’s NRL client
(https://docs.obspy.org/master/packages/obspy.clients.nrl.html#module-obspy.clients.nrl).
'''

import obspy 
from obspy import read, Stream,read_inventory
import os
import json
from obspy import read_inventory
from obspy.clients.nrl import NRL
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

from parameters_py.config import (
					DIR_DATA,SOURCE,NETWORK_CODE,NETWORK_DESCRIPTION,START_DATE,SAMPLING_RATE,LOCATION,
                    OUTPUT_XML_FILE_DIR,OUTPUT_JSON_FILE_DIR
                    
				   )

# Finding IRIS DMC Library of Nominal Responses to instrumental response

# By default this accesses the NRL online. Offline copies of the NRL can
# also be used instead
nrl = NRL()
# The contents of the NRL can be explored interactively in a Python prompt,
# see API documentation of NRL submodule:
# http://docs.obspy.org/packages/obspy.clients.nrl.html
# Here we assume that the end point of data logger and sensor are already
# known:


#Creating inventory:

inv = Inventory(
    # We'll add networks later.
    networks=[],
    # The source should be the id whoever create the file.
    source=SOURCE)


# Creating a list with the network:

net = Network(
    # This is the network code according to the SEED standard.
    code=NETWORK_CODE,
    # A list of stations. We'll add one later.
    stations=[],
    description=NETWORK_DESCRIPTION,
    # Start-and end dates are optional.
    start_date=obspy.UTCDateTime(START_DATE)
    )


# Importing station list
print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['KSTNM']
stla = sta_dic['STLA']
stlo = sta_dic['STLO']
stel = sta_dic['STEL']
sensor_keys = sta_dic['SENSOR_KEYS']
datalogger_keys = sta_dic['DATALOGGER_KEYS']

for i,j in enumerate(kstnm):
    print('Importing '+j+' station parameters')
    print('\n')

    sta = Station(
        # This is the station code according to the SEED standard.
        code=j,
        latitude=float(stla[i]),
        longitude=float(stla[i]),
        elevation=float(stel[i]),
        creation_date=obspy.UTCDateTime(START_DATE),
        site=Site(name=NETWORK_DESCRIPTION))
    
    cha_HHZ = Channel(
        # This is the channel code according to the SEED standard.
        code="HHZ",
        # This is the location code according to the SEED standard.
        location_code=LOCATION,
        # Note that these coordinates can differ from the station coordinates.
        latitude=float(stla[i]),
        longitude=float(stla[i]),
        elevation=float(stel[i]),
        depth=0.0,
        azimuth=0.0,
        dip=-90.0,
        sample_rate=SAMPLING_RATE)

    cha_HHE = Channel(
        # This is the channel code according to the SEED standard.
        code="HHE",
        # This is the location code according to the SEED standard.
        location_code=LOCATION,
        # Note that these coordinates can differ from the station coordinates.
        latitude=float(stla[i]),
        longitude=float(stla[i]),
        elevation=float(stel[i]),
        depth=0.0,
        azimuth=90.0,
        dip=0.0,
        sample_rate=SAMPLING_RATE)

    cha_HHN = Channel(
        # This is the channel code according to the SEED standard.
        code="HHN",
        # This is the location code according to the SEED standard.
        location_code=LOCATION,
        # Note that these coordinates can differ from the station coordinates.
        latitude=float(stla[i]),
        longitude=float(stla[i]),
        elevation=float(stel[i]),
        depth=0.0,
        azimuth=0.0,
        dip=0.0,
        sample_rate=SAMPLING_RATE)
    
    # Now tie it all together.

    response = nrl.get_response(sensor_keys = sensor_keys[i].split(','),datalogger_keys = datalogger_keys[i].split(','))

    cha_HHZ.response = response
    cha_HHN.response = response
    cha_HHE.response = response
    channel_sta = [cha_HHZ,cha_HHN,cha_HHE]
    for k in channel_sta:
        sta.channels.append(k)
    net.stations.append(sta)


inv.networks.append(net)  
os.makedirs(OUTPUT_XML_FILE_DIR,exist_ok=True)
inv.write(OUTPUT_XML_FILE_DIR+NETWORK_CODE+".xml", format="stationxml", validate=True)
print(inv)
print('\n')

print('XML file created')