'''
Script to estimate probabilistic power spectral densities for
one combination of network/station/location/channel/sampling_rate.
(https://docs.obspy.org/tutorial/code_snippets/probabilistic_power_spectral_density.html)

Calculations are based on the routine used by [McNamara2004]:
McNamara, D. E. and Buland, R. P. (2004),
Ambient Noise Levels in the Continental United States,
Bulletin of the Seismological Society of America, 94 (4), 1517-1527.
http://www.bssaonline.org/content/94/4/1517.abstract.


For information on New High/Low Noise Model see [Peterson1993]:
Peterson, J. (1993),
Observations and Modeling of Seismic Background Noise,
U.S. Geological Survey open-file report 93-322, Albuquerque, N.M.
http://ehp3-earthquake.wr.usgs.gov/regional/asl/pubs/files/ofr93-322.pdf
'''

from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
from obspy import UTCDateTime
import obspy
import glob
import json
import numpy as np
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib as mpl
import datetime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from obspy.core.util import AttribDict

from parameters_py.config import (
					OUTPUT_FIGURE_DIR,DIR_DATA,XML_FILE,OUTPUT_PSD_DIR,INITIAL_DATE,FINAL_DATE,
					TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR,
                    PERIOD_PSD,AMP_PSD_MIN,AMP_PSD_MAX,AMP_PSD_HYDROPHONE_MIN,AMP_PSD_HYDROPHONE_MAX
				   )

# ====================================
# Function to calculate PSD from file
# ====================================

def calc_PSD(file):
    st = obspy.read(file)
    l = st[0]

    sta_name = l.stats.station
    NETWORK_CODE = l.stats.network
    sta_channel = l.stats.channel

    time_data = l.stats.starttime
    time_data_year = '{:04}'.format(time_data.year)
    time_data_julday = '{:03}'.format(time_data.julday)
    time_data_hour = '{:02}'.format(time_data.hour)
    time_data_minute = '{:02}'.format(time_data.minute)

    inv = obspy.read_inventory(XML_FILE+'ON.'+sta_name+'.xml')

    if sta_channel == 'HHX':
        ppsd = PPSD(l.stats,inv,special_handling='hydrophone')
        ppsd.add(st)
        os.makedirs(OUTPUT_PSD_DIR+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/',exist_ok=True)
        ppsd.save_npz(OUTPUT_PSD_DIR+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+sta_channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

    else:
        ppsd = PPSD(l.stats, inv)
        ppsd.add(st)
        os.makedirs(OUTPUT_PSD_DIR+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/',exist_ok=True)
        ppsd.save_npz(OUTPUT_PSD_DIR+'/'+sta_name+'/'+sta_channel+'.PPSD'+'/'+NETWORK_CODE+'.'+sta_name+'..'+sta_channel+'.PPSD'+'.'+time_data_year+'.'+time_data_julday+'.npz')

    return 0


# ================================
# Ploting TOTAL PPSD DATA (SENSOR)
# ================================

def plot_PSD(directory):
	files = sorted(glob.glob(directory+'/*'))

	ppsd = PPSD.load_npz(files[0])

	[ppsd.add_npz(i) for i in files[1:]]

	ppsd.calculate_histogram(starttime=UTCDateTime(INITIAL_DATE),endtime=UTCDateTime(FINAL_DATE),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
	folder_output = OUTPUT_FIGURE_DIR+'WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd.station+'/'
	os.makedirs(folder_output,exist_ok=True)
	ppsd.plot(cmap=pqlx,show_coverage=False,filename=folder_output+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')

# =====================================
# Ploting TOTAL PPSD DATA (HYDROPHONE)
# =====================================

def plot_PSD_hydrophone(directory):
    files = sorted(glob.glob(directory+'/*'))

    ppsd = PPSD.load_npz(files[0])

    [ppsd.add_npz(i) for i in files[1:]]

    ppsd.calculate_histogram(starttime=UTCDateTime(INITIAL_DATE),endtime=UTCDateTime(FINAL_DATE),time_of_weekday=[(TIME_OF_WEEKDAY_DAY, TIME_OF_WEEKDAY_START_HOUR, TIME_OF_WEEKDAY_FINAL_HOUR)])
    folder_output = OUTPUT_FIGURE_DIR+'WINDOWED_'+str(int(TIME_OF_WEEKDAY_START_HOUR))+'_'+str(int(TIME_OF_WEEKDAY_FINAL_HOUR))+'/'+ppsd.station+'/'
    os.makedirs(folder_output,exist_ok=True)
    #ppsd.plot(cmap=pqlx,show_coverage=False,show_noise_models=False,filename=folder_output+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf')


    filename = folder_output+ppsd.network+'.'+ppsd.station+'.'+ppsd.channel+'.'+str(ppsd.times_processed[0].year)+'.pdf'

    percentiles=[0, 25, 50, 75, 100]
    period_lim=(0.01, 179)
    cumulative_number_of_colors=20

    fig = plt.figure()
    fig.ppsd = AttribDict()

    ax = fig.add_subplot(111)

    period_lim=(0.01, 179)

    fig.ppsd.cmap = pqlx
    fig.ppsd.label = "[%]"
    fig.ppsd.max_percentage = 30
    fig.ppsd.color_limits = (0, 30)

    title = "%s   %s -- %s  (%i/%i segments)"
    title = title % (ppsd.id,UTCDateTime(ns=ppsd._times_processed[0]).date,
                             UTCDateTime(ns=ppsd._times_processed[-1]).date,
                             ppsd.current_histogram_count,len(ppsd._times_processed))

    ax.set_title(title)
    ax.semilogx()
    ax.set_xlabel('Period [s]')
    ax.set_xlim(period_lim)
    ax.set_ylim(AMP_PSD_HYDROPHONE_MIN, AMP_PSD_HYDROPHONE_MAX)
    ax.set_ylabel('Amplitude [$m^2/s^4/Hz$] [dB]')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))

    data = (ppsd.current_histogram*100.0/(ppsd.current_histogram_count or 1))

    xedges = ppsd.period_xedges

    fig.ppsd.meshgrid = np.meshgrid(xedges, ppsd.db_bin_edges)
    ppsd = ax.pcolormesh(fig.ppsd.meshgrid[0], fig.ppsd.meshgrid[1], data.T,cmap=fig.ppsd.cmap, zorder=-1)
    fig.ppsd.quadmesh = ppsd

    cb = plt.colorbar(ppsd, ax=ax)
    #cb.set_clim(*fig.ppsd.color_limits)
    cb.set_label(fig.ppsd.label)
    fig.ppsd.colorbar = cb

    ppsd.set_clim(*fig.ppsd.color_limits)

    ax.grid(b=True, which="major")
    ax.grid(b=True, which="minor")
    plt.savefig(filename)

# ====================================
# Ploting PPSD DATA BY PERIOD (SENSOR)
# ====================================

def plot_PPSD_by_period_sensor(directory_data):

    data_lista = []

    print('Looking for data in the directory = '+directory_data)

    for root, dirs, files in os.walk(directory_data):
        for name in files:
            if name.endswith('.npz') and 'HHX' not in name:
                data_lista.append(os.path.join(root, name))


    data_lista = sorted(data_lista)

    dataframe_lista = []
    #create a empty dataframe with pandas
    print("Extracting data from PPSD header")

    for j in tqdm(data_lista):
        #Reading header from data
        ppsd = PPSD.load_npz(j,allow_pickle=True)

        #----------------------------
        #Dataframe starting

        network = ppsd.network
        station = ppsd.station
        channel = ppsd.channel

        flat_time_lst = ppsd.times_processed

        DATETIME = str(ppsd.times_data[0][0].year)+','+str(ppsd.times_data[0][0].month)+','+str(ppsd.times_data[0][0].day)

        #Contador da lista de horas
        time_flat_time_lst = [[]]*24
        for g,h in enumerate(np.arange(24)):
            lst_time = []
            for x,c in enumerate(flat_time_lst):
                if c.hour == h:
                    lst_time.append(ppsd.extract_psd_values(PERIOD_PSD)[0][x])
            time_flat_time_lst[g] = lst_time

        AMPLITUDE_HOUR = [[]]*24
        for q,w in enumerate(time_flat_time_lst):
            AMPLITUDE_HOUR[q] = np.mean(w)


        dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[AMPLITUDE_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','AMPLITUDE_HOUR']).T)

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

        def check_datetime_in_period(datetime_lst,df_DATETIME,df_AMPLITUDE_HOUR):

            array_to_plot_by_xlim = []
            for x,c in enumerate(datetime_lst):
                lista_temp = []
                for t,y in enumerate(df_DATETIME):
                    if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c:
                        lista_temp.append(df_AMPLITUDE_HOUR[df_DATETIME[df_DATETIME == y].index[0]])
                array_to_plot_by_xlim.append(lista_temp)

            data_x_axis = []
            for x,c in enumerate(array_to_plot_by_xlim):
                if c != []:
                    data_x_axis.append(c[0][::-1])
                else:
                    data_x_axis.append(np.full_like(np.arange(24),np.nan,dtype=np.double))

            data_x_axis = np.array(data_x_axis).T

            return data_x_axis

        # ====================================
        # Function to plot DATA availability
        # ====================================

        #x axis parameters

        days1 = DayLocator(interval=1)   # every day
        days5 = DayLocator(interval=5)   # every day
        months = MonthLocator()  # every month
        yearsFmt = DateFormatter('%Y-%m-%d')

        days1.MAXTICKS = 10000


        #Matplotlib parameters
        fig, ax = plt.subplots(nrows=len(channel_lista), ncols=1,sharex=True,sharey=True,figsize=(40,15))
        fig.suptitle(j,fontsize=25,y=0.9)
        for k,l in enumerate(channel_lista):

            df_ch = df_sta[df_sta['CHANNEL'] == l]

            data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['AMPLITUDE_HOUR'])

            im = ax[k].imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.viridis,interpolation=None,vmin=AMP_PSD_MIN,vmax=AMP_PSD_MAX)
            ax[k].set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
            ax[k].yaxis.set_major_locator(MultipleLocator(4))
            ax[k].yaxis.set_minor_locator(MultipleLocator(1))
            ax[k].xaxis.set_major_locator(days5)
            ax[k].xaxis.set_major_formatter(yearsFmt)
            ax[k].xaxis.set_minor_locator(days1)
            ax[k].tick_params(which='minor', length=4)
            ax[k].tick_params(which='major', length=10)
            ax[k].set_ylim(0,24)
            ax[k].set_yticklabels([' ',' ', '1h', '5h', '9h', '13h', '17h', '21h',' ']) #y axis according to Brazil UTC-3
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
        cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[AMP_PSD_MIN,np.mean([AMP_PSD_MIN,AMP_PSD_MAX]),AMP_PSD_MAX],label='Amplitude '+r'$[m^2/s^4/Hz][dB]$')

        os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
        fig.savefig(OUTPUT_FIGURE_DIR+j+'_'+'PSD_BY_PERIOD_'+str(PERIOD_PSD)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.pdf',dpi=500)

# ========================================
# Ploting PPSD DATA BY PERIOD (HYDROPHONE)
# ========================================

def plot_PPSD_by_period_hydrophone(directory_data):

    data_lista = []

    print('Looking for data in the directory = '+directory_data)

    for root, dirs, files in os.walk(directory_data):
        for name in files:
            if name.endswith('.npz') and 'HHX' in name:
                data_lista.append(os.path.join(root, name))


    data_lista = sorted(data_lista)

    dataframe_lista = []
    #create a empty dataframe with pandas
    print("Extracting data from PPSD header")

    for j in tqdm(data_lista):
        #Reading header from data
        ppsd = PPSD.load_npz(j,allow_pickle=True)

        #----------------------------
        #Dataframe starting

        network = ppsd.network
        station = ppsd.station
        channel = ppsd.channel

        flat_time_lst = ppsd.times_processed

        DATETIME = str(ppsd.times_data[0][0].year)+','+str(ppsd.times_data[0][0].month)+','+str(ppsd.times_data[0][0].day)

        #Contador da lista de horas
        time_flat_time_lst = [[]]*24
        for g,h in enumerate(np.arange(24)):
            lst_time = []
            for x,c in enumerate(flat_time_lst):
                if c.hour == h:
                    lst_time.append(ppsd.extract_psd_values(PERIOD_PSD)[0][x])
            time_flat_time_lst[g] = lst_time

        AMPLITUDE_HOUR = [[]]*24
        for q,w in enumerate(time_flat_time_lst):
            AMPLITUDE_HOUR[q] = np.mean(w)


        dataframe_lista.append(pd.DataFrame([[network],[station],[channel],[DATETIME],[AMPLITUDE_HOUR]], index=['NETWORK', 'STATION', 'CHANNEL', 'DATETIME','AMPLITUDE_HOUR']).T)

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

        def check_datetime_in_period(datetime_lst,df_DATETIME,df_AMPLITUDE_HOUR):

            array_to_plot_by_xlim = []
            for x,c in enumerate(datetime_lst):
                lista_temp = []
                for t,y in enumerate(df_DATETIME):
                    if datetime.datetime(obspy.UTCDateTime(y).year,obspy.UTCDateTime(y).month,obspy.UTCDateTime(y).day) == c:
                        lista_temp.append(df_AMPLITUDE_HOUR[df_DATETIME[df_DATETIME == y].index[0]])
                array_to_plot_by_xlim.append(lista_temp)

            data_x_axis = []
            for x,c in enumerate(array_to_plot_by_xlim):
                if c != []:
                    data_x_axis.append(c[0][::-1])
                else:
                    data_x_axis.append(np.full_like(np.arange(24),np.nan,dtype=np.double))

            data_x_axis = np.array(data_x_axis).T

            return data_x_axis

        # ====================================
        # Function to plot DATA availability
        # ====================================

        #x axis parameters

        days1 = DayLocator(interval=1)   # every day
        days5 = DayLocator(interval=5)   # every day
        months = MonthLocator()  # every month
        yearsFmt = DateFormatter('%Y-%m-%d')

        days1.MAXTICKS = 10000


        #Matplotlib parameters
        fig, ax = plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(40,15))
        fig.suptitle(j,fontsize=25,y=0.9)

        df_ch = df_sta[df_sta['CHANNEL'] == channel_lista[0]]

        data_x_axis = check_datetime_in_period(datetime_lista,df_ch['DATETIME'],df_ch['AMPLITUDE_HOUR'])

        im = ax.imshow(data_x_axis,extent = [xlim_initial,xlim_final,0,24],cmap=plt.cm.viridis,interpolation=None,vmin=AMP_PSD_HYDROPHONE_MIN,vmax=AMP_PSD_HYDROPHONE_MAX)
        ax.set_xlim(datetime.datetime(obspy.UTCDateTime(INITIAL_DATE).year,obspy.UTCDateTime(INITIAL_DATE).month,obspy.UTCDateTime(INITIAL_DATE).day),datetime.datetime(obspy.UTCDateTime(FINAL_DATE).year,obspy.UTCDateTime(FINAL_DATE).month,obspy.UTCDateTime(FINAL_DATE).day))
        ax.yaxis.set_major_locator(MultipleLocator(4))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_major_locator(days5)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(days1)
        ax.tick_params(which='minor', length=4)
        ax.tick_params(which='major', length=10)
        ax.set_ylim(0,24)
        ax.set_yticklabels([' ',' ', '1h', '5h', '9h', '13h', '17h', '21h',' ']) #y axis according to Brazil UTC-3
        ax.set_ylabel(channel_lista[0],fontsize=15)
        ax.grid(b=True, which='major', color='k', linestyle='-')
        ax.grid(b=True, which='minor', color='k', linestyle='-')


        plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=30)
        ax.set_xlabel('Time', fontsize=20)

        #criando a localização da barra de cores:
        axins = inset_axes(ax,
                            width="10%",  # width = 10% of parent_bbox width
                            height="5%",  # height : 50%
                            loc='upper left',
                            bbox_to_anchor=(0.85, 0.1, 1, 1),
                            bbox_transform=ax.transAxes,
                            borderpad=0,
                           )
        cbar = fig.colorbar(im, cax=axins, orientation="horizontal", ticklocation='top',ticks=[AMP_PSD_MIN,np.mean([AMP_PSD_MIN,AMP_PSD_MAX]),AMP_PSD_MAX],label='Amplitude '+r'$[m^2/s^4/Hz][dB]$')

        os.makedirs(OUTPUT_FIGURE_DIR,exist_ok=True)
        fig.savefig(OUTPUT_FIGURE_DIR+j+'_'+'PSD_BY_PERIOD_HYDROPHONE'+str(PERIOD_PSD)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).year)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).month)+'_'+str(obspy.UTCDateTime(INITIAL_DATE).day)+'_'+str(obspy.UTCDateTime(FINAL_DATE).year)+'_'+str(obspy.UTCDateTime(FINAL_DATE).month)+'_'+str(obspy.UTCDateTime(FINAL_DATE).day)+'.pdf',dpi=500)
