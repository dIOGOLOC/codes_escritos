#!/usr/bin/env python
# coding: utf-8

import time
from tqdm import tqdm
from multiprocessing import Pool

import glob
import os
import numpy as np

import sys
from obspy import read, Trace, Stream, UTCDateTime,read_inventory
from obspy.core import AttribDict
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader, SEGYFile, SEGYTrace
from obspy.io.segy.core import _read_segy

# ==================
# Configuration file
# ==================

# StationXML folder input
XML_FILES_INPUT = '/home/diogoloc/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

# Folders input
MSEED_FILES_INPUT = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_MSEED/'

# Folders output
SEGY_FILES_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/obs_data_SEGY/'

#Number of threads
num_processes = 8

# ========
# Function
# ========

#Convert MSEED data to SEGY data:
def convert_fast(input):
    input_file = input[0]
    output_file = input[1]+input[2]
    XML_FILE_NAME = input[2].split('..')[0]+'.xml'
    stationxml_file = read_inventory(XML_FILES_INPUT+XML_FILE_NAME)

    if os.path.isfile(output_file):
        pass

    else:

        os.makedirs(input[1],exist_ok=True)
        mseed_stream = read(input_file)
        stream_starttime = mseed_stream[0].stats.starttime
        mseed_stream[0].stats.starttime = UTCDateTime(stream_starttime.year,stream_starttime.month,stream_starttime.day,stream_starttime.hour,stream_starttime.minute,stream_starttime.second)
        mseed_stream[0].data = np.float32(mseed_stream[0].data)

        mseed_stream[0].stats.network = stationxml_file[0].code
        mseed_stream[0].stats.station = stationxml_file[0][0].code

        mseed_stream[0].stats.coordinates = AttribDict()
        mseed_stream[0].stats.coordinates.elevation = stationxml_file.get_coordinates(input[2].split('..')[0]+'..HHZ')['elevation']
        mseed_stream[0].stats.coordinates.latitude = stationxml_file.get_coordinates(input[2].split('..')[0]+'..HHZ')['latitude']
        mseed_stream[0].stats.coordinates.longitude = stationxml_file.get_coordinates(input[2].split('..')[0]+'..HHZ')['longitude']

        mseed_streams = Stream()
        for k in mseed_stream.slide(window_length=300, step=300,offset=0, include_partial_windows=True, nearest_sample=True):
            mseed_streams += k
        mseed_streams.write(output_file, format='SEGY')

# =======
# Program
# =======

MSEED_FILES = sorted(glob.glob(MSEED_FILES_INPUT+'**/**/**/**/*'))

input_lst = []

for in_file in MSEED_FILES:
    # splitting subdir/basename
    name_folder1, name_file = os.path.split(in_file)

    name_folder = name_folder1.split('/')[-3:]

    input_lst.append([in_file,SEGY_FILES_OUTPUT+'/'.join(name_folder)+'/',name_file+'.segy'])

#MULTIPROCESSING

start_time = time.time()

with Pool(processes=num_processes) as p:
	max_ = len(input_lst)
	with tqdm(total=max_) as pbar:
		for i, _ in enumerate(p.imap_unordered(convert_fast, input_lst)):
			pbar.update()
print('\n')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
