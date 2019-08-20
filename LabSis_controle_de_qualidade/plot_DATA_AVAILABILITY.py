#!/usr/bin/python -u
"""
Script to get raw DATA availability
"""


import os
import json
import obspy
import datetime
from multiprocessing import Pool



# ==============================
# Generating DATA availability
# ==============================

from visual_py.data_availability import plot_date_file


from parameters_py.config import (
					XML_FILE,DIR_DATA,OUTPUT_FIGURE_DIR
				   )
				   

# ======================
#  Plot time Data 
# ======================

print('\n')
print('Plotting Data Availability')
print('\n')

plot_date_file(FIG_FOLDER_OUTPUT=OUTPUT_FIGURE_DIR,directory_data=DIR_DATA,XML_FILE=XML_FILE)
