#!/usr/bin/python -u
"""
Scritp to get information from Stations and Events
"""
from parameters_py.config import (
					DIR_SAC,NEIC_CSV_FILE
				   )
import os


# ====================================================
# Check some informations in the configuration file  
# ====================================================

if os.path.isdir(DIR_SAC) == True:
	pass
else:
	print(DIR_SAC+' is not a directory! Please, check your configuration file')

if os.path.isfile(NEIC_CSV_FILE) == True:
	pass
else:
	print(NEIC_CSV_FILE+' is not a file! Please, check your configuration file')

print('\n')
print('++++ Input file and directories are OK! ++++')
print('\n')

from pre_processing_py import get_events_information
