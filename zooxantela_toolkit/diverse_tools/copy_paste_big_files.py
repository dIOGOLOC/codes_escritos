#!/usr/bin/env python
# coding: utf-8

import time
from tqdm import tqdm
from multiprocessing import Pool

import glob
import os
import numpy as np


# ==================
# Configuration file
# ==================

# Folders input
ASDF_FILES_INPUT = '/run/user/1000/gvfs/smb-share:server=hatabackup.local,share=on_mar/CLOCK_DRIFT_OUTPUT/ASDF_FILES/CROSS_CORR_DAY_FILES/'

# Folders output
ASDF_FILES_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/ASDF_FILES/'

#Number of threads
num_processes = 12

# ========
# Function
# ========

#Copy and paste files:
def copy_fast(input):
    input_folder = input[0]
    output_folder = input[1]

    os.system('rsync -avPr '+input_folder+' '+output_folder)


# =======
# Program
# =======

days_crosscor_folders = sorted(glob.glob(ASDF_FILES_INPUT+'/*'))

input_lst = []

for in_file in days_crosscor_folders:
    date_folder = in_file.split('/')[-1]

    input_lst.append([in_file+'/',ASDF_FILES_OUTPUT+'CROSS_CORR_DAY_FILES/'+date_folder+'/'])

#MULTIPROCESSING

start_time = time.time()

with Pool(processes=num_processes) as p:
	max_ = len(input_lst)
	with tqdm(total=max_) as pbar:
		for i, _ in enumerate(p.imap_unordered(copy_fast, input_lst)):
			pbar.update()
print('\n')
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))
print('\n')
