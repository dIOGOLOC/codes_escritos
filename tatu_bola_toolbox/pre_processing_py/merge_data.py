'''
Function to merge raw sac files
'''

import obspy as op
import os
import shutil

def merge_data_ZNE(folder_name,knetwk,kstnm,year,month,day,NAME_SUFFIX_E,NAME_SUFFIX_N,NAME_SUFFIX_Z,FILTERS,RMEAN_TYPE,DETREND_TYPE,
			TAPER_TYPE,TAPER_MAX_PERCENTAGE,LOWPASS_FREQ,LOWPASS_CORNER,LOWPASS_ZEROPHASE,HIGHPASS_FREQ,HIGHPASS_CORNER,
			HIGHPASS_ZEROPHASE,SAMPLING_RATE):
	try:
		print('========================= Merging .SAC file ========================= ')
		file_nameX = knetwk+'.'+kstnm+'.'+year+'.'+month+'.'+day+'.'+NAME_SUFFIX_E
		os.chdir(folder_name)
		HHX = op.read('*HHE*')
		HHX.merge(method=1, fill_value='interpolate')
		if FILTERS == True:
			HHX.detrend(type=RMEAN_TYPE)  
			HHX.detrend(type=DETREND_TYPE) 
			HHX.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
			HHX.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
			HHX.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
			HHX.interpolate(sampling_rate=SAMPLING_RATE)
			HHX.write(folder_name+file_nameX,'SAC')
		else:
			HHX.write(folder_name+file_nameX,'SAC')
		os.system('rm *HHE*.sac')
		
		file_nameY = knetwk+'.'+kstnm+'.'+year+'.'+month+'.'+day+'.'+NAME_SUFFIX_N
		HHY = op.read('*HHN*')
		HHY.merge(method=1, fill_value='interpolate')
		if FILTERS == True:
			HHY.detrend(type=RMEAN_TYPE)  
			HHY.detrend(type=DETREND_TYPE) 
			HHY.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
			HHY.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
			HHY.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
			HHY.interpolate(sampling_rate=SAMPLING_RATE)
			HHY.write(folder_name+file_nameY,'SAC')
		else:
			HHY.write(folder_name+file_nameY,'SAC')
		os.system('rm *HHN*.sac')

		file_nameZ = knetwk+'.'+kstnm+'.'+year+'.'+month+'.'+day+'.'+NAME_SUFFIX_Z
		HHZ = op.read('*HHZ*')
		HHZ.merge(method=1, fill_value='interpolate')
		if FILTERS == True:
			HHZ.detrend(type=RMEAN_TYPE)  
			HHZ.detrend(type=DETREND_TYPE) 
			HHZ.taper(type=TAPER_TYPE,max_percentage=TAPER_MAX_PERCENTAGE) 
			HHZ.filter('lowpass',freq=LOWPASS_FREQ,corners=LOWPASS_CORNER,zerophase=LOWPASS_ZEROPHASE) 
			HHZ.filter('highpass',freq=HIGHPASS_FREQ,corners=HIGHPASS_CORNER,zerophase=HIGHPASS_ZEROPHASE)
			HHZ.interpolate(sampling_rate=SAMPLING_RATE)
			HHZ.write(folder_name+file_nameZ,'SAC')
		else:
			HHZ.write(folder_name+file_nameZ,'SAC')
		os.system('rm *HHZ*.sac')
		
		return 'Data OK - Folder = '+folder_name

	except: 
		return 'Data error - Folder = '+folder_name
