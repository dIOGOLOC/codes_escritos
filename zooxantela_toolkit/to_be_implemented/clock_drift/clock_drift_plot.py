 #!/usr/bin/env python
# coding: utf-8

import time
from tqdm import tqdm
from multiprocessing import Pool

import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
import matplotlib.dates as mdates
from matplotlib.patches import Ellipse
import matplotlib.cbook as cbook
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib.transforms import offset_copy

import obspy as op
from obspy import read,read_inventory, UTCDateTime, Stream, Trace
from obspy.io.xseed import Parser
from obspy.signal.filter import bandpass,lowpass
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.util import prev_pow_2
from obspy.signal.cross_correlation import correlate as obscorr
from obspy.signal.cross_correlation import xcorr_max
from obspy.core.util import AttribDict

import glob
import os
import numpy as np
from numpy.fft import rfft, irfft, fft, ifft, fftfreq
from itertools import combinations,product,compress
from numpy.lib.stride_tricks import as_strided
import pandas as pd
from scipy.signal import spectrogram, detrend, resample,savgol_filter,decimate
from scipy.linalg import norm
import pickle
import random
import collections
from copy import copy
import datetime

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import pyproj
from shapely.geometry import shape, Point, Polygon, MultiPolygon
import shapefile

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

from pyasdf import ASDFDataSet
import verde as vd
import geopandas as gpd

# ==================
# Configuration file
# ==================

# Folders input

MSEED_DIR_OBS = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/obs_data_MSEED/'

MSEED_DIR_STA = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/data/'

# -------------------------------

# Shapefile  boundary states input

BOUNDARY_STATES_SHP = '/media/diogoloc/Backup/dados_posdoc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/UFEBRASIL.shp'

OCEANO_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/area_batimetria/area_ON_projeto.shp'

# -------------------------------

# Stations and OBSs information

STATIONS_LST = ['ABR01','DUB01','MAN01','OBS20','OBS22','TER01','ALF01','GDU01','NAN01','TIJ01','CAJ01','GUA01','OBS17','PET01','TRI01','CAM01','JAC01','OBS18','RIB01','VAS01','CMC01','MAJ01','SLP01','PARB','CNLB','BSFB']
STATIONS_LST = sorted(STATIONS_LST)

STATIONXML_DIR = '/media/diogoloc/Backup/dados_posdoc/ON_MAR/XML_ON_OBS_CC/'

CHANNEL_LST = ['HHZ.D','HHN.D','HHE.D','HH1.D','HH2.D']

# -------------------------------

# Folders output

CLOCK_DRIFT_OUTPUT = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/FIGURAS/'

ASDF_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/ASDF_FILES/'

PICKLE_FILES = '/home/diogoloc/dados_posdoc/ON_MAR/CLOCK_DRIFT_OUTPUT/PICKLE_FILES/'

# -------------------------------

# Input parameters

FIRSTDAY = '2019-08-01'
LASTDAY = '2020-06-01'

#Each hour-long seismogram is amplitude clipped at twice its standard deviation of that hour-long time window.
CLIP_FACTOR = 2

MIN_WINDOWS = 30

WINDOW_LENGTH = 3600

#max time window (s) for cross-correlation
SHIFT_LEN = 1800

PERIOD_BANDS = [[2, 5], [7, 25], [20, 50], [50, 100]]
# (these bands focus on periods ~7, 15, 25 seconds)

FREQUENCY_BANDS = [[0.5, 1], [1, 2], [2, 3], [3, 4]]

# default parameters to define the signal and noise windows used to
# estimate the SNR:
# - the signal window is defined according to a min and a max velocity as:
#   dist/vmax < t < dist/vmin
# - the noise window has a fixed size and starts after a fixed trailing
#   time from the end of the signal window
SIGNAL_WINDOW_VMIN = 2.0
SIGNAL_WINDOW_VMAX = 4.0
SIGNAL2NOISE_TRAIL = 700.0
NOISE_WINDOW_SIZE = 700.0

#Returns pairs and spectral SNR array whose spectral SNRs are all >= minspectSNR
minspectSNR = 1

#RESAMPLING
NEW_SAMPLING_RATE = 2

# -------------------------------
# Grid tomography

LLCRNRLON_GRID_TOM = -52
URCRNRLON_GRID_TOM = -38
LLCRNRLAT_GRID_TOM = -30
URCRNRLAT_GRID_TOM = -12

GRID_TOM_STEP = 1

# -------------------------------

# Constants and parameters

ONESEC = datetime.timedelta(seconds=1)
ONEDAY = datetime.timedelta(days=1)

# -------------------------------

# MULTIPROCESSING

num_processes = 6

# =================
# Filtering by date
# =================

fday = UTCDateTime(FIRSTDAY)
lday = UTCDateTime(LASTDAY)
INTERVAL_PERIOD = [UTCDateTime(x.astype(str)) for x in np.arange(fday.datetime,lday.datetime+ONEDAY,ONEDAY)]
INTERVAL_PERIOD_DATE = [str(x.year)+'.'+"%03d" % x.julday for x in INTERVAL_PERIOD]

# =========
# Functions
# =========

# ------------------------------------------------------------------------------

def grid_maker(west, east, south, north, spacing,FILTER_GRID_BY_SHAPEFILE=False):
    '''
    Creating the grid points.
    '''

    region = (west, east, south, north)
    rows, cols = vd.grid_coordinates(region=region, spacing=spacing)

    grdx = rows.ravel()
    grdy = cols.ravel()

    if FILTER_GRID_BY_SHAPEFILE == True:
        polys = shapefile.Reader(OCEANO_SHAPEFILE)
        multi = shape(polys.shapeRecords()[0].shape.__geo_interface__)
        pontos = [Point(grdx[i],grdy[i]) for i,j in enumerate(grdy)]

        grdx = []
        grdy = []

        for i, j in enumerate(pontos):
            if j.within(multi):
                grdx.append(j.x)
                grdy.append(j.y)

    return grdx,grdy

# ------------------------------------------------------------------------------

def geodesic(coord1, coord2, npts):
    '''
    Returns a list of *npts* points along the geodesic between
    (and including) *coord1* and *coord2*, in an array of
    shape (*npts*, 2).
    @rtype: L{ndarray}
    '''

    # reference elipsoid to calculate distance
    wgs84 = pyproj.Geod(ellps='WGS84')

    if npts < 2:
        raise Exception('nb of points must be at least 2')

    path = wgs84.npts(lon1=coord1[0], lat1=coord1[1],
                      lon2=coord2[0], lat2=coord2[1],
                      npts=npts-2)
    return np.array([coord1] + path + [coord2])

# ------------------------------------------------------------------------------

def path_density(x_nodes, y_nodes, paths, window=(GRID_TOM_STEP, GRID_TOM_STEP)):
    '''
    Returns the path density, that is, on each node of the
    grid, the number of paths that cross the rectangular
    cell of size (window[0], window[1]) centered on the node.
    '''

    # initializing path density
    density = np.zeros((len(x_nodes),len(y_nodes)))

    # coordinates of grid nodes and associated windows
    lons_nodes = x_nodes
    lats_nodes = y_nodes

    df_grid = pd.DataFrame({'Latitude': y_nodes,'Longitude': x_nodes})
    gdf_grid = gpd.GeoDataFrame(df_grid, geometry=gdf.points_from_xy(df_grid.Longitude, df_grid.Latitude))

    for path in paths:

        #Figure
        fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(10,10))

        #Add shapefile coast
        reader_1_SHP = Reader(SHAPEFILE)
        shape_1_SHP = list(reader_1_SHP.geometries())
        plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))
        ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=1)
        ax.set_extent([LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT], crs=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))

        lons_path = np.array(path[:, 0])
        lats_path = np.array(path[:, 1])

        df_path = pd.DataFrame({'Latitude': y_nodes,'Longitude': x_nodes})
        gdf_path = gpd.GeoDataFrame(df, geometry=gdf_nodes.points_from_xy(df.Longitude, df.Latitude))


    '''
        # are points of paths in windows?
        # 1st dim = grid nodes; 2nd dim = points along path
        points_in_windows = (lons_path >= lons_min) & (lons_path <= lons_max) & \
                            (lats_path >= lats_min) & (lats_path <= lats_max)

        density += np.any(points_in_windows, axis=-1)

    return density
    '''
# ------------------------------------------------------------------------------

def plot_grid(grx,gry,density,LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT,SHAPEFILE,clock_drift_files_loc_sta1,clock_drift_files_loc_sta2,geopaths):
    #Figure
    fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(10,10))

    #Add shapefile coast
    reader_1_SHP = Reader(SHAPEFILE)
    shape_1_SHP = list(reader_1_SHP.geometries())
    plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))
    ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=1)

    ax.set_extent([LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT], crs=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))

    l3, = ax.plot(grx,gry, '.',markersize=2,markeredgecolor='k',markerfacecolor='k')

    # plotting path density
    #extent = (LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT)
    #d = density
    #m = ax.imshow(d,origin='lower',extent=extent,interpolation='bicubic',cmap='viridis',vmin=0,vmax=10)
    #c = plt.colorbar(m, ax=ax, orientation='horizontal', pad=0.1)
    #c.set_label('Path density')


    for i,j in enumerate(clock_drift_files_loc_sta1):
        ax.plot([clock_drift_files_loc_sta1[i][1],clock_drift_files_loc_sta2[i][1]],[clock_drift_files_loc_sta1[i][0],clock_drift_files_loc_sta2[i][0]],c='k',alpha=0.4)
        ax.scatter(clock_drift_files_loc_sta1[i][1],clock_drift_files_loc_sta1[i][0], marker='^',s=200,c='k',edgecolors='w')
        ax.scatter(clock_drift_files_loc_sta2[i][1],clock_drift_files_loc_sta2[i][0], marker='^',s=200,c='k',edgecolors='w')

    ax.set_xticks(np.arange(LLCRNRLON,URCRNRLON,4))
    ax.set_yticks(np.arange(LLCRNRLAT,URCRNRLAT,4))
    ax.tick_params(labeltop=True, labelbottom=True,labelright=True,labelleft=True,labelsize=15)
    ax.grid(color='grey', linestyle='--', linewidth=1)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    plt.show()

# ------------------------------------------------------------------------------

# ============
# Main program
# ============

print('\n')
print('=========================')
print('Clock Drift for each OBS:')
print('=========================')
print('\n')

OBS_LST = ['OBS17','OBS18','OBS20','OBS22']

clock_drift_files_lst = sorted(glob.glob(PICKLE_FILES+'/*'))
clock_drift_files_loc_sta1_all = []
clock_drift_files_loc_sta2_all = []
for l,k in enumerate(OBS_LST):
    clock_drift_files = [j for i,j in enumerate(clock_drift_files_lst) if k in j]

    clock_drift_files_name_sta1 = []
    clock_drift_files_name_sta2 = []
    clock_drift_files_loc_sta1 = []
    clock_drift_files_loc_sta2 = []

    clock_drift_files_date_to_plot_clock_True_static = []
    clock_drift_files_date_to_plot_clock_True_dynamic = []
    clock_drift_files_date_to_plot_clock_True_absolute = []

    clock_drift_files_data_to_plot_clock_True_static = []
    clock_drift_files_data_to_plot_clock_True_dynamic = []
    clock_drift_files_data_to_plot_clock_True_absolute = []

    for i in clock_drift_files:
        with open(i, 'rb') as f:
            dic_pickle = pickle.load(f)
        clock_drift_files_name_sta1.append(dic_pickle['name_sta1'])
        clock_drift_files_name_sta2.append(dic_pickle['name_sta2'])
        clock_drift_files_loc_sta1.append(dic_pickle['loc_sta1'])
        clock_drift_files_loc_sta2.append(dic_pickle['loc_sta2'])
        clock_drift_files_date_to_plot_clock_True_dynamic.append(dic_pickle['date_to_plot_clock_True_dynamic'])
        clock_drift_files_data_to_plot_clock_True_dynamic.append(dic_pickle['data_to_plot_clock_True_dynamic'])

    # ----------------------------------------------------------------------------------------------------

    clock_drift_files_date_to_plot_clock_True_dynamic = [item for sublist in clock_drift_files_date_to_plot_clock_True_dynamic for item in sublist]
    clock_drift_files_data_to_plot_clock_True_dynamic = [item for sublist in clock_drift_files_data_to_plot_clock_True_dynamic for item in sublist]

    # ----------------------------------------------------------------------------------------------------
    clock_drift_files_loc_sta1_all.append(clock_drift_files_loc_sta1)
    clock_drift_files_loc_sta2_all.append(clock_drift_files_loc_sta2)

'''
    #Creating the figure and plotting Clock-drift
    fig = plt.figure(figsize=(20, 10))
    fig.suptitle('Clock-drift total: '+OBS_LST[l],fontsize=20)

    gs = gridspec.GridSpec(1, 4,wspace=0.2, hspace=0.5)
    map_loc = fig.add_subplot(gs[0:2],projection=ccrs.PlateCarree())

    LLCRNRLON_LARGE = -52
    URCRNRLON_LARGE = -36
    LLCRNRLAT_LARGE = -30
    URCRNRLAT_LARGE = -10

    map_loc.set_extent([LLCRNRLON_LARGE,URCRNRLON_LARGE,LLCRNRLAT_LARGE,URCRNRLAT_LARGE])
    map_loc.yaxis.set_ticks_position('both')
    map_loc.xaxis.set_ticks_position('both')

    map_loc.set_xticks(np.arange(LLCRNRLON_LARGE,URCRNRLON_LARGE+2,2), crs=ccrs.PlateCarree())
    map_loc.set_yticks(np.arange(LLCRNRLAT_LARGE,URCRNRLAT_LARGE+2,2), crs=ccrs.PlateCarree())
    map_loc.tick_params(labelbottom=True, labeltop=True, labelleft=True, labelright=True, labelsize=12)
    map_loc.grid(True,which='major',color='k',linewidth=1,linestyle='-')

    reader_1_SHP = Reader(BOUNDARY_STATES_SHP)
    shape_1_SHP = list(reader_1_SHP.geometries())
    plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.PlateCarree())
    map_loc.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=0.5,zorder=-1)
    # Use the cartopy interface to create a matplotlib transform object
    # for the Geodetic coordinate system. We will use this along with
    # matplotlib's offset_copy function to define a coordinate system which
    # translates the text by 25 pixels to the left.
    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(map_loc)
    text_transform = offset_copy(geodetic_transform, units='dots', y=-5,x=80)

    for i,j in enumerate(clock_drift_files_name_sta1):
        map_loc.plot([clock_drift_files_loc_sta1[i][1],clock_drift_files_loc_sta2[i][1]],[clock_drift_files_loc_sta1[i][0],clock_drift_files_loc_sta2[i][0]],c='k', transform=ccrs.PlateCarree())
        map_loc.scatter(clock_drift_files_loc_sta1[i][1],clock_drift_files_loc_sta1[i][0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())
        map_loc.scatter(clock_drift_files_loc_sta2[i][1],clock_drift_files_loc_sta2[i][0], marker='^',s=200,c='k',edgecolors='w', transform=ccrs.PlateCarree())

    # ----------------------------------------------------------------------------------------------------

    days_major = DayLocator(interval=3)   # every day
    days_minor = DayLocator(interval=1)   # every day
    months = MonthLocator()  # every month
    yearsFmt = DateFormatter('%b-%Y')

    ax1 = fig.add_subplot(gs[2:])
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(yearsFmt)
    ax1.xaxis.set_minor_locator(days_minor)
    ax1.yaxis.set_major_locator(MultipleLocator(0.01))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.001))
    ax1.set_ylim(-0.025,0.025)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Dynamic drift (s)')

    sigma = 1 #70% of the data

    poly_reg = PolynomialFeatures(degree=4)
    X_poly = poly_reg.fit_transform(np.array(range(len(clock_drift_files_date_to_plot_clock_True_dynamic))).reshape(-1, 1))
    pol_reg = LinearRegression()
    pol_reg.fit(X_poly, clock_drift_files_data_to_plot_clock_True_dynamic)
    ax1.plot(clock_drift_files_date_to_plot_clock_True_dynamic, pol_reg.predict(poly_reg.fit_transform(np.array(range(len(clock_drift_files_data_to_plot_clock_True_dynamic))).reshape(-1, 1))), color='blue')

    for y,u in enumerate(clock_drift_files_data_to_plot_clock_True_dynamic):
        if np.mean(clock_drift_files_data_to_plot_clock_True_dynamic)-sigma*np.std(clock_drift_files_data_to_plot_clock_True_dynamic) <= u <= np.mean(clock_drift_files_data_to_plot_clock_True_dynamic)+sigma*np.std(clock_drift_files_data_to_plot_clock_True_dynamic):
            l1, = ax1.plot(clock_drift_files_date_to_plot_clock_True_dynamic[y],clock_drift_files_data_to_plot_clock_True_dynamic[y],'ok',ms=3)
        else:
            l2, = ax1.plot(clock_drift_files_date_to_plot_clock_True_dynamic[y],clock_drift_files_data_to_plot_clock_True_dynamic[y],'or',ms=3)
    ax1.legend((l1,l2),('%70 data','%30 data'),loc='upper right')

    fig.autofmt_xdate()

    # ----------------------------------------------------------------------------------------------------

    output_figure_CLOCK_DRIFT = CLOCK_DRIFT_OUTPUT+'CLOCK_DRIFT_TOTAL_FIGURES/'
    os.makedirs(output_figure_CLOCK_DRIFT,exist_ok=True)
    fig.savefig(output_figure_CLOCK_DRIFT+'CLOCK_DRIFT_SUM_TOTAL_'+OBS_LST[l]+'.png',dpi=300)
    plt.close()

print('\n')
print('=========================')
print('Clock Drift for each OBS:')
print('=========================')
print('\n')

'''

print('Creating GRID POINTS')
print('\n')

grdx, grdy = grid_maker(LLCRNRLON_GRID_TOM,URCRNRLON_GRID_TOM, LLCRNRLAT_GRID_TOM, URCRNRLAT_GRID_TOM,GRID_TOM_STEP)

paths = []
for i,j in enumerate(clock_drift_files_loc_sta1_all):
    geo_path = geodesic(clock_drift_files_loc_sta1_all[i], clock_drift_files_loc_sta2_all[i], 100)

    paths.append(geo_path)

density = path_density(grdx, grdy,paths, window=(GRID_TOM_STEP, GRID_TOM_STEP))

plot_grid(grdx,grdy,density,LLCRNRLON_GRID_TOM,URCRNRLON_GRID_TOM, LLCRNRLAT_GRID_TOM, URCRNRLAT_GRID_TOM,BOUNDARY_STATES_SHP,clock_drift_files_loc_sta1_all,clock_drift_files_loc_sta2_all,paths)
