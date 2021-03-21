#!/usr/bin/env python
# coding: utf-8

import time
from tqdm import tqdm
from multiprocessing import Pool
import numbers

import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.dates import YearLocator, MonthLocator, DayLocator, HourLocator, MinuteLocator, SecondLocator, DateFormatter
from matplotlib.patches import Circle,Rectangle

import json
import glob
import os
import numpy as np
import pandas as pd
import pyproj
import obspy as op

import datetime
import matplotlib.dates as mdates

import windrose
import cartopy.crs as ccrs
from cartopy.io import srtm
from cartopy.io import PostprocessedRasterSource, LocatedImage
from cartopy.io.srtm import SRTM3Source, SRTM1Source
from cartopy.io.shapereader import Reader
import cartopy.feature as cfeature
import shapefile
from shapely.geometry import shape, Point, Polygon, MultiPolygon


from osgeo import gdal, osr
from itertools import combinations

import verde as vd

# ==================
# Configuration file
# ==================

STATION_DISTANCE = 4.5 #degrees

EARTHQUAKES_MIN = 3.5

NUM_PROCESSES = 1

FILTER_GRID_BY_SHAPEFILE = True

LLCRNRLON = -52 #degrees
URCRNRLON = -28 #degrees
LLCRNRLAT = -28 #degrees
URCRNRLAT = -16 #degrees

STA_CONTINENT_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/rsbr_lat_lon_projeto.csv'
OBS_LOCATION_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/coord_obs_data.csv'
COAST_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/brasil_estados.shp'
OCEANO_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/area_batimetria/area_ON_projeto.shp'
EARTHQUAKES_SHAPEFILE =  '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/coord_terremotos_marinhos_brasil.csv'
SRTM_FILE =  '/home/diogoloc/SIG_dados/Projeto_ON_MAR/raster/projeto_dem/projeto_srtm.tif'

#--------------------------------------------------------------------------------------------------------------------

# =========
# Functions
# =========

def plot_grid_obs(plot_stats):
    #Figure
    fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(10,5))

    #Add srtm image
    gdal.UseExceptions()
    ds = gdal.Open(SRTM_FILE)
    data = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()

    inproj = osr.SpatialReference()
    inproj.ImportFromWkt(proj)

    # Add the shaded SRTM
    img = ax.imshow(data, extent=[LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT],cmap='terrain',origin='upper',alpha=0.5)

    #Add shapefile coast
    reader_1_SHP = Reader(COAST_SHAPEFILE)
    shape_1_SHP = list(reader_1_SHP.geometries())
    plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))
    ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=1)

    ax.set_extent([LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT], crs=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))

    l2, = ax.plot(terremoto_lon_lst,terremoto_lat_lst, '*',markersize=5,markeredgecolor='y',markerfacecolor='y')
    l1, = ax.plot(sta_terra_lon,sta_terra_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey')
    l3, = ax.plot(sta_OBS_lon,sta_OBS_lat, 's',markersize=10,markeredgecolor='k',markerfacecolor='w')

    for o,p in enumerate(sta_OBS_lon):
        circulo = Circle(radius=STATION_DISTANCE,xy=(sta_OBS_lon[o], sta_OBS_lat[o]),color='none', ec='k',linewidth=0.25,transform=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None),zorder=2)
        ax.add_patch(circulo)

    legend = ax.legend([l1,l2,l3],['Estações RSBR','Terremotos m$_{L}$>'+str(EARTHQUAKES_MIN),'Possíveis OBSs'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)
    ax.set_xticks(np.arange(LLCRNRLON,URCRNRLON,4))
    ax.set_yticks(np.arange(LLCRNRLAT,URCRNRLAT,4))
    ax.tick_params(labeltop=True, labelbottom=True,labelright=True,labelleft=True,labelsize=12)
    ax.grid(color='grey', linestyle='--', linewidth=1)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_title('Detecção de Terremotos de magnitude m$_{L}$>'+str(EARTHQUAKES_MIN),fontsize=15,y=1.05)

    if plot_stats == True:
        plt.show()
    else:
        pass

#------------------------------------------

def calc_angle(a,b,c):
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

#------------------------------------------

def calc_azimuthal_gap(input_lst):
    #input_list = [obs_xy[t],sta_obs_xy]

    #calculating the distance between the possible obs location and land stations with other obs locations:
    grx = input_lst[0][0]
    gry = input_lst[0][1]
    lst_sta_obs = input_lst[1]

    grid_sta_obs_lon = []
    grid_sta_obs_lat = []
    azimuth_lst = []
    azimuthal_gap_lst = []
    for i,sta_loc in enumerate(lst_sta_obs):
        dist,azimuth,_ = op.geodetics.base.gps2dist_azimuth(gry, grx, sta_loc[1], sta_loc[0], a=6378137.0, f=0.0033528106647474805)
        dist_degree = op.geodetics.base.kilometer2degrees(dist/1000)

        if dist_degree < STATION_DISTANCE:
            grid_sta_obs_lon.append(sta_loc[0])
            grid_sta_obs_lat.append(sta_loc[1])
            azimuth_lst.append(azimuth)

            azimuthal_gap_lst.append(1)

    return [[grx,gry],azimuthal_gap_lst,azimuth_lst]

# ============
# Main program
# ============

print(' ------------------------------ ')
print(' Importing variables and tables ')
print('\n')

#Coordenadas dos terremotos
terremoto_df = pd.read_csv(EARTHQUAKES_SHAPEFILE)
terremoto_lon = terremoto_df['LONG.']
terremoto_lat = terremoto_df['LAT.']
terremoto_mag = terremoto_df['Io']

terremoto_lon_lst = []
terremoto_lat_lst = []
for i,j in enumerate(terremoto_mag):
    try:
        tmp = float(j)
        if float(j) > EARTHQUAKES_MIN:
            terremoto_lon_lst.append(terremoto_lon[i])
            terremoto_lat_lst.append(terremoto_lat[i])
    except:
        pass

#Coordenadas das estações em terra
terra_STA_df = pd.read_csv(STA_CONTINENT_FILE)
sta_terra_lon = terra_STA_df['Longitude']
sta_terra_lat = terra_STA_df['Latitude']

#Coordenadas dos possíveis OBS
OBS_STA_df = pd.read_csv(OBS_LOCATION_FILE)
sta_OBS_lon = OBS_STA_df['LON']
sta_OBS_lat = OBS_STA_df['LAT']

#Plotting the map:
plot_grid_obs(True)

#Creating a list with input variables:
sta_xy = []
for t,u in enumerate(sta_terra_lon):
    sta_xy.append([sta_terra_lon[t],sta_terra_lat[t]])

obs_xy = []
for t,u in enumerate(sta_OBS_lon):
    obs_xy.append([sta_OBS_lon[t],sta_OBS_lat[t]])

input_var = []
for t,u in enumerate(obs_xy):
    temp = obs_xy
    temp.pop(t)
    sta_obs_xy = sta_xy+temp
    temp.insert(t, u)
    input_var.append([obs_xy[t],sta_obs_xy])

azimuthal_test_lst = []
with Pool(processes=NUM_PROCESSES) as p:
    max_ = len(input_var)
    with tqdm(total=max_, desc='Combinations loop') as pbar:
        for i, result in enumerate(p.imap_unordered(calc_azimuthal_gap, input_var)):
            azimuthal_test_lst.append(result)
            pbar.update()


#Figure final
fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(10,5))

#Add srtm image
gdal.UseExceptions()
ds = gdal.Open(SRTM_FILE)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()

inproj = osr.SpatialReference()
inproj.ImportFromWkt(proj)

# Add the shaded SRTM
img = ax.imshow(data, extent=[LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT],cmap='terrain',origin='upper',alpha=0.5)

#Add shapefile coast
reader_1_SHP = Reader(COAST_SHAPEFILE)
shape_1_SHP = list(reader_1_SHP.geometries())
plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))
ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=1)

ax.set_extent([LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT], crs=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))

ax.set_xticks(np.arange(LLCRNRLON,URCRNRLON,4))
ax.set_yticks(np.arange(LLCRNRLAT,URCRNRLAT,4))
ax.tick_params(labeltop=True, labelbottom=True,labelright=True,labelleft=True,labelsize=12)
ax.grid(color='grey', linestyle='--', linewidth=1)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_title('Resultado da cobertura azimutal da rede de OBSs para terremotos de m$_{L}$>'+str(EARTHQUAKES_MIN),fontsize=15,y=1.05)

l3, = ax.plot(sta_OBS_lon,sta_OBS_lat, 's',markersize=10,markeredgecolor='k',markerfacecolor='w',zorder=10)
l2, = ax.plot(terremoto_lon_lst,terremoto_lat_lst, '*',markersize=5,markeredgecolor='k',markerfacecolor='y',zorder=1)
l1, = ax.plot(sta_terra_lon,sta_terra_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',zorder=1)

legend = ax.legend([l1,l2,l3],['Estações RSBR','Terremotos m$_{L}$>'+str(EARTHQUAKES_MIN),'Possíveis OBSs'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)

for azi in azimuthal_test_lst:
    obs_gap_x = azi[0][0]
    obs_gap_y = azi[0][1]
    azi_gap = azi[1]
    azimuths = azi[2]
    # Inset axe it with a fixed size
    wrax_cham = inset_axes(ax,
        width=1,                             # size in inches
        height=1,                            # size in inches
        loc='center',                        # center bbox at given position
        bbox_to_anchor=(obs_gap_x, obs_gap_y), # position of the axe
        bbox_transform=ax.transData,    # use data coordinate (not axe coordinate)
        axes_class=windrose.WindroseAxes    # specify the class of the axe
        )

    wrax_cham.bar(azimuths,azi_gap, bins=np.arange(0, 20, 1), fc='k',ec='white',zorder=5,alpha=0.75)
    wrax_cham.axis('off')

plt.show()

#----------------------------------------------------------------------------
