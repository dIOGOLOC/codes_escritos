'''
Function to merge raw sac files
'''

import obspy as op
import os
import shutil

def merge_data_ZNE(folder_name,knetwk,kstnm,NAME_SUFFIX_E,NAME_SUFFIX_N,NAME_SUFFIX_Z,FILE_FORMAT):
	try:
		print('========================= Merging .SAC file ========================= ')
		os.chdir(folder_name)
		HH = op.read('*')
		HH.merge(method=1, fill_value='interpolate')
		
		for i,j in enumerate(HH):
			year =  '{:04}'.format(j.stats.starttime.year)
			julday =  '{:03}'.format(j.stats.starttime.julday)
			if j.stats.channel in ['HHZ']:
				file_nameZ = knetwk+'.'+kstnm+'.'+year+'.'+julday+'.'+NAME_SUFFIX_Z
				j.write(folder_name+'/'+file_nameZ,FILE_FORMAT)
			if j.stats.channel in ['HHE','HH2','HHX']:
				file_nameX = knetwk+'.'+kstnm+'.'+year+'.'+julday+'.'+NAME_SUFFIX_E
				j.write(folder_name+'/'+file_nameX,FILE_FORMAT)
			if j.stats.channel in ['HHN','HH1','HHY']:
				file_nameY = knetwk+'.'+kstnm+'.'+year+'.'+julday+'.'+NAME_SUFFIX_N
				j.write(folder_name+'/'+file_nameY,FILE_FORMAT)
		
		os.system('rm *.sac')
		
		return 'Data OK - Folder = '+folder_name

	except: 
		return 'Data error - Folder = '+folder_name
