#!/usr/bin/python -u
'''

--------------------------------------------------------------------------------
          Getting information about the header of the raw data
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 12/2021


Project: Monitoramento Sismo-Oceanográfico
P. Number: 2015/00515-6


Description:
This code will retrieve information about the header of the raw data for each
daily file and plot a mosaic the the Data availability.


More information in:
(https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#obspy.core.stream.read)

Inputs:
Daily .SAC file

Outputs:
Images of the quantity of files per hour (format: PDF)

Examples of Usage (in command line):
   >> python plot_DATA_AVAILABILITY.py

'''


import os
import json
import obspy
import datetime
from multiprocessing import Pool



# ==============================
# Generating DATA availability
# ==============================

from assessment_py.data_availability import plot_date_file


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
