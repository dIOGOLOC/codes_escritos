#!/usr/bin/python -u
"""
Scritp to plot PPSD windowed DATA
"""

import os

# ==============================
# Generating DATA availability
# ==============================

from visual_py.plot_PSD_DATA import plot_PPSD_WINDOWED_data


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,MP_PROCESSES,NETWORK_CODE,OUTPUT_PSD_DIR,DAY_PERCENTAGE
				   )

# Importing stations list

# ===========================
# Finding stations PPSD data
# ===========================

print('\n')
print('Looking for PPSD STATIONS data in'+OUTPUT_PSD_DIR)
print('\n')

datafile_lst = [] 
for root, dirs, files in os.walk(OUTPUT_PSD_DIR):
	for directories in dirs:
		datafile_name = os.path.join(root, directories)
		if '.PPSD' in datafile_name:
			datafile_lst.append(datafile_name)
datafile_lstS = sorted(datafile_lst)

# ================
#  plot PPSD Data 
# ================

for i,j in enumerate(datafile_lstS):
	plot_PPSD_WINDOWED_data(j)