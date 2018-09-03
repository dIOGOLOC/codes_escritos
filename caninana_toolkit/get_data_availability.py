#!/usr/bin/python -u
"""
Scritp to generate DATA availability plots
"""
import os
import json
import obspy
import datetime
from multiprocessing import Pool




# ==============================
# Generating DATA availability
# ==============================

from visual_py.data_availability import get_date_file


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,OUTPUT_FIGURE_DIR,EXAMPLE_OF_FILE
				   )

# ===================================
#  Function to call get data function
# ===================================

def parallel_get_data(lst):
	result = get_date_file(input_list=lst)
	
	if isinstance(result,dict):
		return result

# Importing stations list
print('\n')
print('Looking for STATIONS data in JSON file in '+OUTPUT_JSON_FILE_DIR)
print('\n')

filename_STA = OUTPUT_JSON_FILE_DIR+'STA_dic.json'

sta_dic = json.load(open(filename_STA))

kstnm = sta_dic['KSTNM']
stla = sta_dic['STLA']
stlo = sta_dic['STLO']
stel = sta_dic['STEL']
sensor_keys = sta_dic['SENSOR_KEYS']
datalogger_keys = sta_dic['DATALOGGER_KEYS']

print('Creating input list with endtime of each raw file')
print('\n')

input_list = [[]]*len(kstnm)

for i,j in enumerate(datalogger_keys):

	if  j.split(',')[0] == 'Nanometrics':
		datafile_lst = [] 
		for root, dirs, files in os.walk(DIR_DATA):
			for directories in dirs:
				datafile_name = os.path.join(root, directories)
				
				if '/'+kstnm[i]+'/' in datafile_name and len(datafile_name.split('/')) == len(EXAMPLE_OF_FILE.split('/')):
						datafile_lst.append(datafile_name)
		datafile_lstS = sorted(datafile_lst)
		
		print(' Number of days = '+ str(len(datafile_lstS)))		
		
		# ==============================
		#  Creating stations Input lists
		# ==============================

		print('Creating stations input lists')
		print('\n')

		input_list[i] = [
					[l] for l in datafile_lstS
					]
		print('==============================')


	if  j.split(',')[0] == 'REF TEK':
		datafile_lst = [] 
		for root, dirs, files in os.walk(DIR_DATA):
			for directories in dirs:
				datafile_name = os.path.join(root, directories)
				if '/'+kstnm[i]+'/' in datafile_name and len(datafile_name.split('/')) == len(EXAMPLE_OF_FILE.split('/')):
					datafile_lst.append(datafile_name)
		datafile_lstS = sorted(datafile_lst)

		print(' Number of days = '+ str(len(datafile_lstS)))

		# ==============================
		#  Creating stations Input lists
		# ==============================

		print('Creating stations input lists')
		print('\n')

		input_list[i] = [
					[l] for l in datafile_lstS
					]
		print('==============================')
		

# ==============
#  Get time Data 
# ==============



print('Get time Data  for each station')
print('\n')

os.makedirs(OUTPUT_JSON_FILE_DIR,exist_ok=True)
pool = Pool(MP_PROCESSES)
result_dic = {'kstnm':[],'data':[]}
for x,y in enumerate(input_list):
	result_dic['data'].append(pool.starmap(parallel_get_data, y))
	result_dic['kstnm'].append(kstnm[x])

with open(OUTPUT_JSON_FILE_DIR+'TIME_dic.json', 'w') as fp:
		json.dump(result_dic, fp)