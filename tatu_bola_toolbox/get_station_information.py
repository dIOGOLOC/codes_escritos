#!/usr/bin/python -u
"""
Scritp to get information from Stations

Exemple of .TXT file to retrieve station 
informations about file type, name, longitude and latitude.

TYPE		NAME	LAT			LON
GCF			CRAT	-5.2712		-40.4705
GCF		 	BDCO	-5.4517		-45.0203
REFTEK130	TUTU	-5.3542		-44.7615
REFTEK130	SAAL	-4.9732		-42.0217
MSEED		BPPF	-6.2271		-47.2518	

To see the examples of file types check the following link:
(https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#supported-formats)

"""
from parameters_py.config import (
					DIR_SAC,STA_CSV_FILE
				   )
import os


# ====================================================
# Check some informations in the configuration file  
# ====================================================

if os.path.isfile(STA_CSV_FILE) == True:
	pass
else:
	print(STA_CSV_FILE+' is not a file! Please, check your configuration file')

print('\n')
print('++++ Input file and directories are OK! ++++')
print('\n')

from pre_processing_py import get_station_information
