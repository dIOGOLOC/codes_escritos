'''
Script to get information about the header of the raw data
(https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#obspy.core.stream.read)
and plot a mosaic the the Data availability.
'''

import matplotlib.pyplot as plt
import pandas as pd
import obspy
import os
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.signal import PPSD
from obspy.clients.arclink.client import Client
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,XML_FILE,INITIAL_DATE,FINAL_DATE,USER,HOST,PORT,INSTITUTION
				   )

# =================================
# Function to plot all channel data
# =================================

def plot_date_file(FIG_FOLDER_OUTPUT,directory_data,XML_FILE):

    data_lista = []

    print('Looking for data in the directory = '+directory_data)

    for root, dirs, files in os.walk(directory_data):
        for name in files:
            data_lista.append(os.path.join(root, name))

    data_lista = sorted(data_lista)
    
    dataframe_lista = []
    #create a empty dataframe with pandas
    for i,j in enumerate(data_lista):
        print("Extracting data from header "+str(i+1)+" of "+str(len(data_lista)))
                
        #Reading header from data
        st = obspy.read(j)
                
        #----------------------------
        #Dataframe starting

        network = st[0].stats.network
        station = st[0].stats.station
        channel = st[0].stats.channel

        time_lst = []

        for t,trace in enumerate(st):
            starttime = trace.stats.starttime
            endtime = trace.stats.endtime
            time_lst.append(np.arange(starttime,endtime,60))

        flat_time_lst = [item for sublist in time_lst for item in sublist]

        DATETIME = str(st[0].stats.starttime.year)+','+str(st[0].stats.starttime.month)+','+str(st[0].stats.starttime.day)

        #Contador da lista de horas
        time_flat_time_lst = [[]]*24
        for g,h in enumerate(np.arange(24)):
            lst_time = []
            for x,c in enumerate(flat_time_lst):
                if c.hour == h:   
                    lst_time.append(c.hour)
            time_flat_time_lst[g] = lst_time

        NUMBER_HOUR = [[]]*24
        for q,w in enumerate(time_flat_time_lst):
            NUMBER_HOUR[q] = len(w)
          

        dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[NUMBER_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','NUMBER_HOUR']).T)
        print(pd.DataFrame([[network],[station],[channel],[DATETIME],[NUMBER_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','NUMBER_HOUR']).T)
        print()
        print()

        #Dataframe ending
        #----------------------------

    df = pd.concat(dataframe_lista, ignore_index=True)

    #Sorting according to station

    station_lista = list(set(df['STATION']))

    for i,j in enumerate(station_lista):
        df_sta = df[df['STATION'] == j]

        channel_lista = list(set(df_sta['CHANNEL']))
        channel_lista = sorted(channel_lista)

        # ==========================================================
        # Calculating datetime between INITIAL_DATE and  FINAL_DATE
        # ==========================================================

        datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)

        datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)

        datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)


        xlim_initial = mdates.date2num(datatime_initial)
        xlim_final = mdates.date2num(datatime_final)
            
        #----------------------------
        #Function to check if the dates in data set are inside the period chosen (INITIAL_DATE to FINAL_DATE)

        def check_datetime_in_period(datetime_lst,df_DATETIME,df_NUMBER_HOUR):

            array_to_plot_by_xlim = []
            for x,c in enumerate(datetime_lst):
                lista_temp = []
                for t,y in enumerate(df_DATETIME):
                        if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c:
                                lista_temp.append(df_NUMBER_HOUR[df_DATETIME[df_DATETIME == y].index[0]])
                array_to_plot_by_xlim.append(lista_temp)

            data_x_axis = []
            for x,c in enumerate(array_to_plot_by_xlim):
                if c != []:
                    data_x_axis.append(c[0][::-1])
                else:
                    data_x_axis.append(np.zeros_like(np.arange(24)))

            data_x_axis = np.array(data_x_axis).T

            return data_x_axis

        # ====================================
        # Function to plot DATA availability
        # ====================================

        #x axis parameters

        days1 = DayLocator(interval=1)   # every day
        days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every day
        months = MonthLocator()  # every month
        yearsFmt = DateFormatter('%Y-%m-%d')
        
        days1.MAXTICKS = 10000


        #Matplotlib parameters
        fig, ax = plt.subplots(nrows=len(channel_lista), ncols=1,sharex=True,sharey=True,figsize=(40,15))
        fig.suptitle(j,fontsize=25,y=0.9)

        for k,l in enumerate(channel_lista):

            df_ch = df_sta[df_sta['CHANNEL'] == l]

            data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['NUMBER_HOUR'])

            im = ax[k].imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=60)
            ax[k].set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
            ax[k].yaxis.set_major_locator(MultipleLocator(4))
            ax[k].yaxis.set_minor_locator(MultipleLocator(1))
            ax[k].xaxis.set_major_locator(days5)
            ax[k].xaxis.set_major_formatter(yearsFmt)
            ax[k].xaxis.set_minor_locator(days1)
            ax[k].tick_params(which='minor', length=4)
            ax[k].tick_params(which='major', length=10)
            ax[k].set_ylim(0,24)
            ax[k].set_ylabel(l,fontsize=15)
            ax[k].grid(b=True, which='major', color='k', linestyle='-')
            ax[k].grid(b=True, which='minor', color='k', linestyle='-')

        
        plt.setp(ax[k].xaxis.get_majorticklabels(), fontsize=10, rotation=30)
        ax[-1].set_xlabel('Time', fontsize=20)
        
        #criando a localização da barra de cores:
        axins = inset_axes(ax[0],
                           width="10%",  # width = 10% of parent_bbox width
                           height="5%",  # height : 50%
                           loc='upper left',
                           bbox_to_anchor=(0.85, 0.1, 1, 1),
                           bbox_transform=ax[0].transAxes,
                           borderpad=0,
                           )
        cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,30,60],label='Files per hour')
        cbar.ax.set_xticklabels(['0%','50%','100%'])
        
        os.makedirs(FIG_FOLDER_OUTPUT,exist_ok=True)
        fig.savefig(FIG_FOLDER_OUTPUT+j+'_'+'COMPLETENESS_'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.pdf',dpi=500)
        #plt.show()

# ==================================
# Function to plot each channel data
# ==================================

def plot_date_file_channel(FIG_FOLDER_OUTPUT,directory_data,XML_FILE):

    data_lista = []

    print('Looking for data in the directory = '+directory_data)

    for root, dirs, files in os.walk(directory_data):
        for name in files:
            data_lista.append(os.path.join(root, name))

    data_lista = sorted(data_lista)
    
    dataframe_lista = []
    #create a empty dataframe with pandas
    for i,j in enumerate(data_lista):
        print("Extracting data from header "+str(i+1)+" of "+str(len(data_lista)))
                
        #Reading header from data
        st = obspy.read(j)
                
        #----------------------------
        #Dataframe starting

        network = st[0].stats.network
        station = st[0].stats.station
        channel = st[0].stats.channel

        time_lst = []

        for t,trace in enumerate(st):
            starttime = trace.stats.starttime
            endtime = trace.stats.endtime
            time_lst.append(np.arange(starttime,endtime,60))

        flat_time_lst = [item for sublist in time_lst for item in sublist]

        DATETIME = str(st[0].stats.starttime.year)+','+str(st[0].stats.starttime.month)+','+str(st[0].stats.starttime.day)

        #Contador da lista de horas
        time_flat_time_lst = [[]]*24
        for g,h in enumerate(np.arange(24)):
            lst_time = []
            for x,c in enumerate(flat_time_lst):
                if c.hour == h:   
                    lst_time.append(c.hour)
            time_flat_time_lst[g] = lst_time

        NUMBER_HOUR = [[]]*24
        for q,w in enumerate(time_flat_time_lst):
            NUMBER_HOUR[q] = len(w)
          

        dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[NUMBER_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','NUMBER_HOUR']).T)
        print(pd.DataFrame([[network],[station],[channel],[DATETIME],[NUMBER_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','NUMBER_HOUR']).T)
        print()
        print()

        #Dataframe ending
        #----------------------------

    df = pd.concat(dataframe_lista, ignore_index=True)

    #Sorting according to station

    station_lista = list(set(df['STATION']))

    for i,j in enumerate(station_lista):
        df_sta = df[df['STATION'] == j]

        channel_lista = list(set(df_sta['CHANNEL']))
        channel_lista = sorted(channel_lista)

        # ==========================================================
        # Calculating datetime between INITIAL_DATE and  FINAL_DATE
        # ==========================================================

        datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)

        datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)

        datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)


        xlim_initial = mdates.date2num(datatime_initial)
        xlim_final = mdates.date2num(datatime_final)
            
        #----------------------------
        #Function to check if the dates in data set are inside the period chosen (INITIAL_DATE to FINAL_DATE)

        def check_datetime_in_period(datetime_lst,df_DATETIME,df_NUMBER_HOUR):

            array_to_plot_by_xlim = []
            for x,c in enumerate(datetime_lst):
                lista_temp = []
                for t,y in enumerate(df_DATETIME):
                        if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c:
                                lista_temp.append(df_NUMBER_HOUR[df_DATETIME[df_DATETIME == y].index[0]])
                array_to_plot_by_xlim.append(lista_temp)

            data_x_axis = []
            for x,c in enumerate(array_to_plot_by_xlim):
                if c != []:
                    data_x_axis.append(c[0][::-1])
                else:
                    data_x_axis.append(np.zeros_like(np.arange(24)))

            data_x_axis = np.array(data_x_axis).T

            return data_x_axis

        # ====================================
        # Function to plot DATA availability
        # ====================================

        #x axis parameters

        days1 = DayLocator(interval=1)   # every day
        days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every day
        months = MonthLocator()  # every month
        yearsFmt = DateFormatter('%Y-%m-%d')
        
        days1.MAXTICKS = 10000

        for k,l in enumerate(channel_lista):
            #Matplotlib parameters
            fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(15,5))

            df_ch = df_sta[df_sta['CHANNEL'] == l]

            data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['NUMBER_HOUR'])

            #Percentage of data
            total_length = sum(len(row) for row in data_x_axis)
            total_no_data  = np.count_nonzero(data_x_axis)

            data_percentage = int(round((total_no_data*100)/total_length))

            fig.suptitle(j+'.'+l+' - Completude: '+str(data_percentage)+'%',fontsize=25,y=0.9)

            im = ax.imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=60)
            ax.set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
            ax.yaxis.set_major_locator(MultipleLocator(4))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            ax.xaxis.set_major_locator(days5)
            ax.xaxis.set_major_formatter(yearsFmt)
            ax.xaxis.set_minor_locator(days1)
            ax.tick_params(which='minor', length=4)
            ax.tick_params(which='major', length=10)
            ax.set_ylim(0,24)
            ax.set_ylabel(l,fontsize=15)
            ax.grid(b=True, which='major', color='k', linestyle='-')
            ax.grid(b=True, which='minor', color='k', linestyle='-')
        
            plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)
            
            #criando a localização da barra de cores:
            axins = inset_axes(ax,
	                           width="10%",  # width = 10% of parent_bbox width
	                           height="5%",  # height : 50%
	                           loc='upper left',
	                           bbox_to_anchor=(0.85, 0.1, 1, 1),
	                           bbox_transform=ax.transAxes,
	                           borderpad=0,
	                           )
            cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,30,60],label='Arquivos por hora')
            cbar.ax.set_xticklabels(['0%','50%','100%'])
	        
            os.makedirs(FIG_FOLDER_OUTPUT,exist_ok=True)
            fig.savefig(FIG_FOLDER_OUTPUT+j+'_'+'COMPLETENESS_'+l+'_'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.pdf',dpi=500)

def get_date_file_via_client(FIG_FOLDER_OUTPUT,STATION_NAME):
    
    #Conectando ao arclink
    client = Client(user=USER,host=HOST, port=PORT, institution=INSTITUTION)

    #Crio a lista com dias de acordo tempo inicial e final
    data_lista = np.arange(obspy.UTCDateTime(INITIAL_DATE),obspy.UTCDateTime(FINAL_DATE),86400) 

    #Criando um dicionário em branco
    time_dic = nested_dict()

    #Varredura na lista com dias e lendo o cabeçalho dos arquivos
    for i,j in enumerate(data_lista):
        try:
            print("Extracting data from header "+str(i+1)+" of "+str(len(data_lista)))

            #Lendo o cabeçalho da forma de onda diária do arclink
            st = client.get_waveforms(NETWORK_CODE,STATION_NAME,LOCATION,CHANNEL_CODE,j,j+86400)
            print(st)

            #----------------------------
            #Começo do dicionário

            time_dic['network'] = st[0].stats.network
            time_dic['station'] = st[0].stats.station
            time_dic['channel'] = st[0].stats.channel


            time_lst = []

            for t,trace in enumerate(st):
                starttime = trace.stats.starttime
                endtime = trace.stats.endtime
                time_lst.append(np.arange(starttime,endtime,60))

            flat_time_lst = [item for sublist in time_lst for item in sublist]

            time_dic[str(j)]['DATETIME'] = str(st[0].stats.starttime.year)+','+str(st[0].stats.starttime.month)+','+str(st[0].stats.starttime.day)

            #Contador da lista de horas
            time_flat_time_lst = [[]]*24
            for g,h in enumerate(np.arange(24)):
                lst_time = []
                for x,c in enumerate(flat_time_lst):
                    if c.hour == h:   
                        lst_time.append(c.hour)
                time_flat_time_lst[g] = lst_time

            for q,w in enumerate(np.arange(24)):
                time_dic[str(j)]['HOUR'][str(w)] = len(time_flat_time_lst[int(w)])


        except:
            print('No data to: '+str(j))

            time_dic[str(j)]['DATETIME'] = str(j.year)+','+str(j.month)+','+str(j.day) 
            for q,w in enumerate(np.arange(24)):            
                time_dic[str(j)]['HOUR'][str(w)] = 0
        #Fim do dicionário
        #----------------------------
    # ======================================================================
    # Alocando os valores do contador de horas na variável time_plot_lst 
    # ======================================================================

    time_plot_lst = [[]]*len(data_lista)
    for i,j in enumerate(data_lista):
        a = []
        for q,w in enumerate(np.arange(24)):
            a.append(time_dic[str(j)]['HOUR'][str(w)])
            time_plot_lst[i] = a


    # =========================================================================
    # Calculando as datas entre o período inicial e final para plotar no eixo x 
    # =========================================================================

    datatime_initial = datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day)

    datatime_final = datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day)

    datetime_lista = np.arange(datatime_initial, datatime_final, datetime.timedelta(days=1)).astype(datetime.datetime)


    xlim_initial = mdates.date2num(datatime_initial)
    xlim_final = mdates.date2num(datatime_final)
    
    #----------------------------
    #Checando se as datas do banco de dados batem com o período escolhido (INITIAL_DATE e FINAL_DATE)

    array_to_plot_by_xlim = []
    for x,c in enumerate(datetime_lista):
        lista_temp = []
        for i,j in enumerate(data_lista):
                if datetime.datetime(obspy.UTCDateTime(time_dic[str(j)]['DATETIME']).year,obspy.UTCDateTime(time_dic[str(j)]['DATETIME']).month,obspy.UTCDateTime(time_dic[str(j)]['DATETIME']).day) == c:
                        lista_temp.append(time_plot_lst[i])
        array_to_plot_by_xlim.append(lista_temp)

    data_x_axis = []
    for x,c in enumerate(array_to_plot_by_xlim):
        if c != []:
            data_x_axis.append(c[0][::-1])
        else:
            data_x_axis.append(np.zeros_like(np.arange(24)))

    data_x_axis = np.array(data_x_axis).T


    # ====================================
    # Function to plot DATA availability
    # ====================================

    #x axis parameters
    days1 = DayLocator(interval=1)   # every day
    days5 = DayLocator(interval=int(len(datetime_lista)*5/100))   # every day
    months = MonthLocator()  # every month
    yearsFmt = DateFormatter('%Y-%m-%d')

    #Matplotlib parameters
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,10))
    im = ax.imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.Greens,interpolation=None, vmin=0, vmax=60)
    ax.set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
    ax.yaxis.set_major_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(days5)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(days1)
    ax.tick_params(which='minor', length=4)
    ax.tick_params(which='major', length=10)
    ax.set_ylim(0,24)
    ax.set_ylabel('Hora do Dia')
    ax.set_xlabel('Tempo')
    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='k', linestyle='-')
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)
    ax.set_title(time_dic['network']+'.'+time_dic['station']+'.'+time_dic['channel'], fontsize=20   ,y=1.05)

    #criando a localização da barra de cores:
    axins = inset_axes(ax,
                       width="10%",  # width = 10% of parent_bbox width
                       height="5%",  # height : 50%
                       loc='upper left',
                       bbox_to_anchor=(0.85, 0.07, 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
    cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[0,30,60],label='Quantidade de arquivos por hora')
    cbar.ax.set_xticklabels(['0%','50%','100%'])
    os.makedirs(FIG_FOLDER_OUTPUT,exist_ok=True)
    fig.savefig(FIG_FOLDER_OUTPUT+time_dic['network']+'_'+time_dic['station']+'_'+time_dic['channel']+'_'+'COMPLETENESS.pdf',dpi=300)
    plt.show()
