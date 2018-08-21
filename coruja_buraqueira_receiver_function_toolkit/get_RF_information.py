#!/usr/bin/python -u
"""
Scritp to get information from Stations and Events
"""
from parameters_py.config import (
					DIR_SAC
				   )
import os


# ====================================================
# Check some informations in the configuration file  
# ====================================================

if os.path.isdir(DIR_SAC) == True:
	pass
else:
	print(DIR_SAC+' is not a directory! Please, check your configuration file')

print('\n')
print('++++ Input file and directories are OK! ++++')
print('\n')

from visual_py import get_RF_information
