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


# # Pastas com as entradas e saídas

#INPUT FOLDER (Pasta com as RFs)
INPUT_DIR = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/MTZ_2024/PRF_SEISPY_DATA_YES_PP_FILTER/BP.TUTU/'

#OUTPUT FOLDER (Pasta para salvar a figura)
OUTPUT_DIR = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/MTZ_2024/OUTPUT/'


# # Lendo as Funções do Receptor

RF = glob.glob(INPUT_DIR+'*P_R*')


# # Lendo a lista de Funções do Receptor selecionadas pelo SEISPY

RF_select_list = glob.glob(INPUT_DIR+'*finallist.*')[0]


table_RF = np.genfromtxt(RF_select_list,delimiter=' P ',dtype='str')


name_RF_sel = []
for i,j in enumerate(table_RF):
    name_RF_sel.append(j[0])

final_RF_sel_lst = [ ]
for i in name_RF_sel:
    final_RF_sel_lst.append([k for k in RF if i in k][0]) 
final_RF_sel_lst = sorted(final_RF_sel_lst)

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

plt.savefig(OUTPUT_DIR+'FIGURES_YES_PP_FILTER/RF_'+RF_data.stats.network+'_'+RF_data.stats.station+'_ALL_BAZ.png',dpi=300)
