'''
Function to trim raw sac files
'''


import os
import obspy as op
from obspy.taup import TauPyModel
import numpy as np

from parameters_py.config import (
					DIR_SAC,DIR_EVENT,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES
				   )

#Function to cut each file in each folder from every event without PP wave:


def rotate_data_by_event(kstnm,event_dir):
		try:
			print('Current Directory = '+event_dir)		
			os.chdir(DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday)
			
			#########################################################################################################################################################
			#Component E
			
			HH = op.read('*')

			headerHHX = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_E, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(starttime.hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute),'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]), 'nzsec': int(starttime.second), 
						'nzyear': int(starttime.year),'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P),'delta':HHX[0].stats.delta
						}

				sacHHX = op.io.sac.sactrace.SACTrace(data=np.array(data_HHX), **headerHHX)
				sacHHX.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_E)
				
				#########################################################################################################################################################
				#Component N
				HHY = op.read('*.'+NAME_SUFFIX_N)
				HHY.trim(starttime,endtime)

				data_HHY = []
		
				for k,l in enumerate(HHY[0].times()):
					if time_PP - 30.0 <= l <= time_PP - 5.0:
						data_HHY.append(HHY[0].data[k] * np.exp(-0.2*(l - time_PP + 30.0)))
					elif time_PP - 5.0 <= l <= time_PP + 5.0: 
						data_HHY.append(HHY[0].data[k] * 0.0)
					elif time_PP + 5.0 <= l <= time_PP + 30.0: 
						data_HHY.append(HHY[0].data[k] * np.exp(-0.2*(-l + time_PP + 30.0)))
					else:
						data_HHY.append(HHY[0].data[k])

				headerHHY = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_N, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': int(starttime.hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),
						'nzyear': int(starttime.year),	'cmpaz': 0.0, 'cmpinc': 90.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P),'delta':HHY[0].stats.delta
						}

				sacHHY = op.io.sac.sactrace.SACTrace(data=np.array(data_HHY), **headerHHY)
				sacHHY.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:03}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_N)


				#########################################################################################################################################################
				#Component Z

				HHZ = op.read('*.'+NAME_SUFFIX_Z)
				HHZ.trim(starttime,endtime)

				data_HHZ = []
					
				for k,l in enumerate(HHZ[0].times()):
					if time_PP - 30.0 <= l <= time_PP - 5.0:
						data_HHZ.append(HHZ[0].data[k] * np.exp(-0.2*(l - time_PP + 30.0)))
					elif time_PP - 5.0 <= l <= time_PP + 5.0: 
						data_HHZ.append(HHZ[0].data[k] * 0.0)
					elif time_PP + 5.0 <= l <= time_PP + 30.0: 
						data_HHZ.append(HHZ[0].data[k] * np.exp(-0.2*(-l + time_PP + 30.0)))
					else:
						data_HHZ.append(HHZ[0].data[k])

				headerHHZ = {
						'kstnm': kstnm, 'kcmpnm': KCMPNM_Z, 'stla': float(stla), 'stlo': float(stlo), 'evdp': float(ev_depth),  
						'evla': float(ev_lat), 'evlo': float(ev_long), 'mag': float(ev_mag), 'nzhour': float(ev_hour),
						'nzjday': int(starttime.julday), 'nzmin': int(starttime.minute), 'nzmsec': int('{:03}'.format(starttime.microsecond)[:3]),'nzsec': int(starttime.second),
						'nzyear': int(starttime.year),	'cmpaz': 0.0, 'cmpinc': 0.0, 'dist': float(dist/1000), 'gcarc': float(gcarc), 
						'az': float(az), 'baz': float(baz), 'user8': float(arr.ray_param/6371),'o':float(CUT_BEFORE_P),'delta':HHZ[0].stats.delta
						}

				sacHHZ = op.io.sac.sactrace.SACTrace(data=np.array(data_HHZ), **headerHHZ)
				print(sacHHZ)
				sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_Z)
				
			return 'Data OK with this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday
		except: 
			return 'Data error with this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday
		
	return 'No data for this event = '+DIR_SAC+kstnm+'/'+ev_year+'/'+ev_julday
