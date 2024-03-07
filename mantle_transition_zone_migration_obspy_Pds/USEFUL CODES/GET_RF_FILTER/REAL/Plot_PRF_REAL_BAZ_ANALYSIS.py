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
INPUT_DIR = '/home/sysop/dados_posdoc/MTZ_2024/PRF_selected_YES_PP_FILTER_POST/'

#OUTPUT FOLDER (Pasta para salvar a figura)
OUTPUT_DIR = '/home/sysop/dados_posdoc/MTZ_2024/OUTPUT/'


# # Lendo as Funções do Receptor

RF_folders = glob.glob(INPUT_DIR+'*/')

for folder in tqdm(RF_folders,total=len(RF_folders),desc='Plotting'):
    final_RF_sel_lst = glob.glob(folder+'*_P_R.sac')

    # # Empilhando as Funções do Receptor

    RF_stack_lst = np.array([op.read(i)[0].data for i in final_RF_sel_lst])
    RF_stack = np.sum(RF_stack_lst, axis=0)/len(RF_stack_lst)

    # # Organizando as Funções do Receptor em função do backazimuth (°)

    gcarc_lst = ([op.read(i,headonly=True)[0].stats.sac.baz for i in final_RF_sel_lst])
    orglisl = np.argsort(gcarc_lst)

    # # Criando a figura através do  gridspec:
    # ## Cada nova função do receptor graficada é adicionada um fator (0.05) para plotar todas no mesmo eixo

    fig = plt.figure(figsize=(15,25),facecolor='white')
    # set up subplot grid
    gs = gridspec.GridSpec(2,1,wspace=1,hspace=0,height_ratios=[1,20])

    factor = 0

    ax1 = fig.add_subplot(gs[1])
    for j,i in enumerate(orglisl):
        RF_data = op.read(final_RF_sel_lst[i])[0]
        factor += 0.05

        vetor_normalizado = RF_data.data

        ax1.plot(RF_data.times()-10,factor+vetor_normalizado,'k',linewidth=1)
        ax1.fill_between(RF_data.times()-10,factor+vetor_normalizado,factor,where=(factor+vetor_normalizado)>=factor, facecolor='midnightblue',alpha=0.5, interpolate=True)
        ax1.fill_between(RF_data.times()-10,factor+vetor_normalizado,factor,where=(factor+vetor_normalizado)<=factor, facecolor='lightskyblue',alpha=0.5, interpolate=True)

        var_y = factor+vetor_normalizado[0]

        if j % 10 == 0:
            plt.text(160.1,factor+vetor_normalizado[0],"{0:.2f}".format(RF_data.stats.sac.baz)+'°', fontsize='xx-large')

    ax1.set_yticks([])
    ax1.set_ylim(0,factor+0.5)
    ax1.set_ylabel('Backazimuth (°)',labelpad=2,fontsize='xx-large')
    ax1.set_xlim(-5,160)
    ax1.text(0.93,0.98, 'n:'+str(len(final_RF_sel_lst)),transform=ax1.transAxes, fontsize='xx-large')

    # -------

    ax2 = fig.add_subplot(gs[0])

    ax2.plot(RF_data.times()-10,RF_stack,'k',linewidth=1)
    ax2.fill_between(RF_data.times()-10,RF_stack,0,where=(RF_stack)>=0, facecolor='midnightblue', interpolate=True)
    ax2.fill_between(RF_data.times()-10,RF_stack,0,where=(RF_stack)<=0, facecolor='lightskyblue', interpolate=True)

    ax2.set_yticks([])
    ax2.set_xticklabels([])
    ax2.set_ylabel('Stack',labelpad=2,fontsize='xx-large')
    ax2.set_xlim(-5,160)

    ax2.set_title(RF_data.stats.network+'.'+RF_data.stats.station, va='center', fontsize='xx-large')
    os.makedirs(OUTPUT_DIR+'FIGURES_YES_PP_FILTER_REAL',exist_ok=True)
    plt.savefig(OUTPUT_DIR+'FIGURES_YES_PP_FILTER_REAL/RF_'+RF_data.stats.network+'_'+RF_data.stats.station+'_ALL_BAZ.png',dpi=300)
    plt.close()