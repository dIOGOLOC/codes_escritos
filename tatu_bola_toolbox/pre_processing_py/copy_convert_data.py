'''
Script to copy raw data before merge
'''

#Useful modules

import os
from obspy import read
import shutil

from parameters_py.config import ( 
					knetwk,KCMPNM_N,KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,FILE_TYPE,
					NAME_SUFFIX_Z,INTERPOLATE,SAMPLING_RATE,FILTER,LOWPASS_FREQ,LOWPASS_CORNER,FILE_SUFFIX,
					LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,HIGHPASS_ZEROPHASE,FILTERS_TAPER,
					TAPER_TYPE,TAPER_MAX_PERCENTAGE,FILTERS_DETREND,DETREND_TYPE,FILTERS_RMEAN,RMEAN_TYPE
)

def copy_convert_raw_data(FOLDER_OUTPUT,DATA_RAW,STA_NAME):
	print('========================= Copying Data ========================= ')
	try:
		tr = read(DATA_RAW)
		for i,j in enumerate(tr):
				
			FILE_NETWORK = knetwk
			FILE_CHANNEL = j.stats.channel
			FILE_YEAR =  '{:04}'.format(j.stats.starttime.year)
			FILE_JULDAY =  '{:03}'.format(j.stats.starttime.julday)
			FILE_HOUR =  '{:02}'.format(j.stats.starttime.hour)
			FILE_MINUTE =  '{:02}'.format(j.stats.starttime.minute)
			FILE_SECOND =  '{:02}'.format(j.stats.starttime.second)
			FILE_MSECOND =  '{:03}'.format(j.stats.starttime.microsecond)
				
			FILE_NAME = FILE_NETWORK+'.'+STA_NAME+'.'+FILE_CHANNEL+'.'+FILE_YEAR+'.'+FILE_JULDAY+'.'+FILE_HOUR+'.'+FILE_MINUTE+'.'+FILE_SECOND+'.'+FILE_MSECOND+'.'+FILE_SUFFIX
			OUTPUT_FOLDER_FILE = FOLDER_OUTPUT+'/'+STA_NAME+'/'+FILE_YEAR+'/'+FILE_JULDAY+'/'
			os.makedirs(OUTPUT_FOLDER_FILE, exist_ok=True)
			if FILTERS_RMEAN == True:
				j.detrend(type=RMEAN_TYPE)  
			if FILTERS_DETREND == True:
				j.detrend(type=DETREND_TYPE) 
			if FILTERS_TAPER == True:
				j.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
			if FILTER == True:
				j.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
				j.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
			if INTERPOLATE == True:
				j.interpolate(sampling_rate=SAMPLING_RATE)
				j.write(OUTPUT_FOLDER_FILE+FILE_NAME,format=FILE_TYPE)
		return 	'File ok = '+DATA_RAW

	except:
		return 	'File error '+DATA_RAW
		pass