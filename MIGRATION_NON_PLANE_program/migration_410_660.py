#!/usr/bin/python -u
"""
This script reads receiver function data from a set of stations, and
calculates the ray path between all pairs station-event. The data MUST be 
located in some root folder (RF_DIR). The station-event informations 
(coordinates) MUST be in the header of the files. In the current version 
of the program, files should be organized inside their directory as you prefer, 
but they need to finish with some extension name (RF_EXT), e.g.:
*.mseed, *.sac, *.itr, *.eqr
The implemented algorithm follows the lines of Gao and Liu (2014),
"Imaging mantle discontinuities using multiply-reflected P-to-S
conversions", Earth and Planetary Science Letters 402 (2014) 99–106.
Here we utilize both the P-to-S converted phase (Pds) and the multiply reflected 
and converted phase (Ppds) at the discontinuities to simultaneously determine the 
depth of mantle discontinuities and velocity anomalies in the overlying layer. 
Note that all the parameters are defined in the configuration file.
"""

from time_scripts import time_P_Pds
from variable_scripts import mgconfig
import os
import sys
import warnings

# ====================================================
# parsing configuration file to import some parameters
# ====================================================

from variable_scripts.mgconfig import (RF_DIR,PROG_MIGRATION_DIR,MODEL_FILE_NPZ,DIST_T_DIR)

print("- dir of receiver function data: " + RF_DIR)
print("- dir of dist time: " + DIST_T_DIR)
print("- velocity model file: " + MODEL_FILE_NPZ)
print("- dir of the program: " + PROG_MIGRATION_DIR)
