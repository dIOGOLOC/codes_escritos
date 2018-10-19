'''
Function to trim raw sac files
'''


import os
import obspy as op
from obspy.taup import TauPyModel
import numpy as np

from parameters_py.config import (
					DIR_EVENT,STA_CSV_FILE,OUTPUT_JSON_FILE_DIR,TAUPY_MODEL,DIR_EVENT_NO_PP,
					EV_GCARC_MIN,EV_GCARC_MAX,EV_MAGNITUDE_MB,CUT_BEFORE_P,CUT_AFTER_P,KCMPNM_N,
					KCMPNM_E,KCMPNM_Z,knetwk,NAME_SUFFIX_N,NAME_SUFFIX_E,NAME_SUFFIX_Z,MP_PROCESSES
				   )
#Function to cut each file in each folder from every event without PP wave:


def cut_PP_data_by_event(data_folder,kstnm):
	event_directory_YES_PP = data_folder
	if os.path.isdir(event_directory_YES_PP) == True:
		#Change to Event Directory 
		print(event_directory_YES_PP)
		os.chdir(event_directory_YES_PP)	

		#########################################################################################################################################################
		#Component E
		HHX = op.read('*.'+NAME_SUFFIX_E)

		#Calculating distance, azimuth and backazimuth
		dist,az,baz = op.geodetics.gps2dist_azimuth(HHX[0].stats.sac.evla,HHX[0].stats.sac.evlo,HHX[0].stats.sac.stla,HHX[0].stats.sac.stlo)
		gcarc = op.geodetics.kilometer2degrees(dist/1000)

		#Reference time
		starttime = HHX[0].stats.starttime
		
		#Creating Event Directory NP PP
		event_directory = DIR_EVENT_NO_PP+kstnm+'/'+'{:04}'.format(starttime.year)+'/'+'{:03}'.format(starttime.julday)+'/'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]
		os.makedirs(event_directory, exist_ok=True)
		
		model = TauPyModel(model=TAUPY_MODEL)

		arrivalsP =  model.get_travel_times(source_depth_in_km=HHX[0].stats.sac.evdp, distance_in_degree=HHX[0].stats.sac.gcarc, phase_list=["P"])
		arrP = arrivalsP[0]	

		arrivalsPP =  model.get_travel_times(source_depth_in_km=HHX[0].stats.sac.evdp, distance_in_degree=HHX[0].stats.sac.gcarc, phase_list=["PP"])
		arrPP = arrivalsPP[0]	

		time_PP = arrPP.time - arrP.time + CUT_BEFORE_P

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
					'kstnm': kstnm, 'kcmpnm': KCMPNM_E, 'stla': float(HHX[0].stats.sac.stla), 'stlo': float(HHX[0].stats.sac.stlo), 'evdp': float(HHX[0].stats.sac.evdp),  
					'evla': float(HHX[0].stats.sac.evla), 'evlo': float(HHX[0].stats.sac.evlo), 'mag': float(HHX[0].stats.sac.mag), 'nzhour': int(HHX[0].stats.starttime.hour),
					'nzjday': int(HHX[0].stats.starttime.julday), 'nzmin': int(HHX[0].stats.starttime.minute),'nzmsec': int('{:03}'.format(HHX[0].stats.starttime.microsecond)[:3]), 'nzsec': int(HHX[0].stats.starttime.second), 
					'nzyear': int(HHX[0].stats.starttime.year),'cmpaz': 90.0, 'cmpinc': 90.0,'dist': float(HHX[0].stats.sac.dist), 'gcarc': float(HHX[0].stats.sac.gcarc), 
					'az': float(HHX[0].stats.sac.az), 'baz': float(HHX[0].stats.sac.baz), 'user8': float(HHX[0].stats.sac.user8),'o':float(CUT_BEFORE_P),'delta':HHX[0].stats.sac.delta
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
					'kstnm': kstnm, 'kcmpnm': KCMPNM_N, 'stla': float(HHY[0].stats.sac.stla), 'stlo': float(HHY[0].stats.sac.stlo), 'evdp': float(HHY[0].stats.sac.evdp),  
					'evla': float(HHY[0].stats.sac.evla), 'evlo': float(HHY[0].stats.sac.evlo), 'mag': float(HHY[0].stats.sac.mag), 'nzhour': int(HHY[0].stats.starttime.hour),
					'nzjday': int(HHY[0].stats.starttime.julday), 'nzmin': int(HHY[0].stats.starttime.minute),'nzmsec': int('{:03}'.format(HHY[0].stats.starttime.microsecond)[:3]), 'nzsec': int(HHY[0].stats.starttime.second), 
					'nzyear': int(HHY[0].stats.starttime.year),'cmpaz': 0.0, 'cmpinc': 90.0,'dist': float(HHY[0].stats.sac.dist), 'gcarc': float(HHY[0].stats.sac.gcarc), 
					'az': float(HHY[0].stats.sac.az), 'baz': float(HHY[0].stats.sac.baz), 'user8': float(HHY[0].stats.sac.user8),'o':float(CUT_BEFORE_P),'delta':HHY[0].stats.sac.delta
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
					'kstnm': kstnm, 'kcmpnm': KCMPNM_Z, 'stla': float(HHY[0].stats.sac.stla), 'stlo': float(HHY[0].stats.sac.stlo), 'evdp': float(HHY[0].stats.sac.evdp),  
					'evla': float(HHY[0].stats.sac.evla), 'evlo': float(HHY[0].stats.sac.evlo), 'mag': float(HHY[0].stats.sac.mag), 'nzhour': int(HHY[0].stats.starttime.hour),
					'nzjday': int(HHY[0].stats.starttime.julday), 'nzmin': int(HHY[0].stats.starttime.minute),'nzmsec': int('{:03}'.format(HHY[0].stats.starttime.microsecond)[:3]), 'nzsec': int(HHY[0].stats.starttime.second), 
					'nzyear': int(HHY[0].stats.starttime.year),'cmpaz': 0.0, 'cmpinc': 0.0,'dist': float(HHY[0].stats.sac.dist), 'gcarc': float(HHY[0].stats.sac.gcarc), 
					'az': float(HHY[0].stats.sac.az), 'baz': float(HHY[0].stats.sac.baz), 'user8': float(HHY[0].stats.sac.user8),'o':float(CUT_BEFORE_P),'delta':HHY[0].stats.sac.delta
					}


		sacHHZ = op.io.sac.sactrace.SACTrace(data=np.array(data_HHZ), **headerHHZ)
		sacHHZ.write(event_directory+'/'+knetwk+'.'+kstnm+'.'+'{:04}'.format(starttime.year)+'.'+'{:03}'.format(starttime.julday)+'.'+'{:02}'.format(starttime.hour)+'.'+'{:02}'.format(starttime.minute)+'.'+'{:02}'.format(starttime.second)+'.'+'{:03}'.format(starttime.microsecond)[:3]+'.'+NAME_SUFFIX_Z)			