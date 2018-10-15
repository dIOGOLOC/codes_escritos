'''
Function to trim raw sac files
'''


import os
import obspy as op
from obspy.taup import TauPyModel
import numpy as np

from parameters_py.config import (
					DIR_SAC,DIR_EVENT,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,DIR_EVENT_NO_PP,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES
				   )
#Function to cut each file in each folder from every event without PP wave:


def cut_PP_data_by_event(kstnm,stla,stlo,ev_timeUTC,ev_julday,ev_year,ev_month,ev_day,ev_hour,ev_minute,ev_second,ev_microsecond,ev_lat,ev_long,ev_depth,ev_mag):
	#Calculating distance, azimuth and backazimuth
	dist,az,baz = op.geodetics.gps2dist_azimuth(ev_lat,ev_long,stla,stlo)
	gcarc = op.geodetics.kilometer2degrees(dist/1000)
	
	if EV_GCARC_MIN <= gcarc <= EV_GCARC_MAX and EV_MAGNITUDE_MB <= ev_mag:
		try:
			#Calculating ray parameter
			model = TauPyModel(model=TAUPY_MODEL)
			arrivals = model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["P"])
			arr = arrivals[0]
			
			#Reference time
			starttime = (op.UTCDateTime(ev_timeUTC)+arr.time)-CUT_BEFORE_P
			
			event_directory_YES_PP = DIR_EVENT+kstnm+'/'+'{:04}'.format(starttime.year)+'/'+'{:03}'.format(starttime.julday)+'/'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]
			
			if os.path.isdir(event_directory_YES_PP) == True:
				#Change to Event Directory 
				print(event_directory_YES_PP)
				os.chdir(event_directory_YES_PP)	

				#Creating Event Directory NP PP
				event_directory = DIR_EVENT_NO_PP+kstnm+'/'+'{:04}'.format(starttime.year)+'/'+'{:03}'.format(starttime.julday)+'/'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]
				os.makedirs(event_directory, exist_ok=True)


				



				
				#########################################################################################################################################################
				#Component E
				HHX = op.read('*.'+NAME_SUFFIX_E)

				arrivalsPP =  model.get_travel_times(source_depth_in_km=ev_depth, distance_in_degree=gcarc, phase_list=["PP"])
				arrPP = arrivalsPP[0]	

				time_PP = arrPP.time - arr.time + CUT_BEFORE_P
				data_HHX = []
					
				for k,l in enumerate(HHX[0].times()):
					if time_PP - 30.0 <= l <= time_PP - 5.0:
						data_HHX.append(HHX[0].data[k] * np.exp(-0.2*(l - time_PP + 30.0)))
					elif time_PP - 5.0 <= l <= time_PP + 5.0: 
						data_HHX.append(HHX[0].data[k] * 0.0)
					elif time_PP + 5.0 <= l <= time_PP + 30.0: 
						data_HHX.append(HHX[0].data[k] * np.exp(-0.2*(-l + time_PP + 30.0)))
					else:
						data_HHX.append(HHX[0].data[k])
					
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
				HHY = op.read('*.'+NAME_SUFFIX_Z)
				
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
				sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_Z)
							
				return 'Event = '+event_directory
		except:
				return 'Problem - '+event_directory_YES_PP