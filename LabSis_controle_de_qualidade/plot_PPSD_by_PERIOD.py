#!/usr/bin/python -u
"""
Scritp to plot PPSD windowed DATA
"""

import os

# ==============================
# Generating DATA availability
# ==============================

from visual_py.plot_PSD_DATA import plot_PPSD_by_period


from parameters_py.config import (
					OUTPUT_JSON_FILE_DIR,DIR_DATA,OUTPUT_PSD_DIR
				   )

# Importing stations list

# ===========================
# Finding stations PPSD data
# ===========================

print('\n')
print('Looking for PPSD STATIONS data in'+OUTPUT_PSD_DIR)
print('\n')

# ================
#  plot PPSD Data 
# ================

plot_PPSD_by_period(OUTPUT_PSD_DIR)