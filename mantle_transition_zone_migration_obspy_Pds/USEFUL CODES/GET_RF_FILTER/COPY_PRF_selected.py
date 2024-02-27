# Importando as bibliotecas em python

import numpy as np
import matplotlib.pyplot as plt
import obspy as op
import glob
import os
import shutil
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator,Locator
from matplotlib import colormaps
import pandas as pd
from tqdm import tqdm

# # Pastas com as entradas e saídas

#INPUT FOLDER (Pasta com as RFs)
INPUT_DIR = '/home/sysop/dados_posdoc/MTZ_2024/PRF_SEISPY_DATA_YES_PP_FILTER/'

#OUTPUT FOLDER (Pasta para salvar a figura)
OUTPUT_DIR = '/home/sysop/dados_posdoc/MTZ_2024/OUTPUT/'


# # Lendo as Funções do Receptor

RF = glob.glob(INPUT_DIR+'*/*P_R*')

# # Lendo a lista de Funções do Receptor selecionadas pelo SEISPY

RF_select_list = glob.glob(INPUT_DIR+'*/*finallist.*')
for k in RF_select_list:
    table_RF = np.genfromtxt(k,delimiter=' P ',dtype='str')

    name_RF_sel = []
    for i,j in enumerate(table_RF):
        name_RF_sel.append(j[0])

    final_RF_sel_lst = [ ]
    for i in name_RF_sel:
        try:
            final_RF_sel_lst.append([k for k in RF if i in k][0]) 
        except:
            pass
    final_RF_sel_lst = sorted(final_RF_sel_lst)

    # --------------------
    # saving PRFs selected

    df_lst = []
    for i in tqdm(sorted(final_RF_sel_lst)):
        dir_name = os.path.dirname(i)
        RF_sel_directory = OUTPUT_DIR+'PRF_selected_YES_PP_FILTER/'+i.split('/')[-2]+'/'
        os.makedirs(RF_sel_directory,exist_ok=True)
            
        shutil.copy2(dir_name+'/'+i.split('/')[-1].split('_')[0]+'_P_R.sac',RF_sel_directory+i.split('/')[-1].split('_')[0]+'_P_R.sac')
        shutil.copy2(dir_name+'/'+i.split('/')[-1].split('_')[0]+'_P_T.sac',RF_sel_directory+i.split('/')[-1].split('_')[0]+'_P_T.sac')
