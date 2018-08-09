'''
Script to copy raw data before merge
'''

#Useful modules

import os
from obspy import read
import shutil

from parameters_py.config import (
					INTERPOLATE,SAMPLING_RATE,INTERPOLATE_TYPE
				   )

def copy_convert_raw_data(FOLDER_OUTPUT,DATA_RAW,FILE_NAME):
	os.makedirs(FOLDER_OUTPUT, exist_ok=True)
	print('========================= Copying Data ========================= ')
	tr = read(DATA_RAW)
	if INTERPOLATE == True:
		tr.interpolate(sampling_rate=SAMPLING_RATE, method=INTERPOLATE_TYPE)
	tr.write(FOLDER_OUTPUT+FILE_NAME)
	return 	'File = '+FOLDER_OUTPUT+FILE_NAME
