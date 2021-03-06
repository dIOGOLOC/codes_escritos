#!/usr/bin/python -u
"""
Scritp to get information from Stations and Events
"""
from parameters_py.config import (
					DIR_SAC,STA_CSV_FILE
				   )
import os


# ====================================================
# Check some informations in the configuration file  
# ====================================================

if os.path.isdir(DIR_SAC) == True:
	pass
else:
	print(DIR_SAC+' is not a directory! Please, check your configuration file')

if os.path.isfile(STA_CSV_FILE) == True:
	pass
else:
	print(STA_CSV_FILE+' is not a file! Please, check your configuration file')

print('\n')
print('++++ Input file and directories are OK! ++++')
print('\n')

from get_information_py import get_station_information
