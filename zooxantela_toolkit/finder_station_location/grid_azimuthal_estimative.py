#!/usr/bin/env python
# coding: utf-8

import time
from tqdm import tqdm
from multiprocessing import Pool
import numbers

import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm, colorbar, colors
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
from scipy.interpolate import griddata

# ==================
# Configuration file
# ==================

GRID_DISTANCE = 1.5 #degrees

GRID_SPACING = 0.25  #degrees

EARTHQUAKES_MIN = 2.5  #mag

NUM_PROCESSES = 8

FILTER_GRID_BY_SHAPEFILE = False

PLOT_STATS = False

LLCRNRLON = -52 #degrees
URCRNRLON = -28 #degrees
LLCRNRLAT = -28 #degrees
URCRNRLAT = -16 #degrees

STA_CONTINENT_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/rsbr_lat_lon_projeto.csv'
OBS_LOCATION_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/coord_obs_data_arr1.csv'
COAST_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/brasil_estados.shp'
OCEANO_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/area_batimetria/area_ON_projeto.shp'
EARTHQUAKES_SHAPEFILE =  '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/coord_terremotos_marinhos_brasil.csv'
SRTM_FILE =  '/home/diogoloc/SIG_dados/Projeto_ON_MAR/raster/projeto_dem/projeto_srtm.tif'

#--------------------------------------------------------------------------------------------------------------------

# =========
# Functions
# =========

#Function for creating GRID POINTS
def grid_maker(west, east, south, north):
    region = (west, east, south, north)
    rows, cols = vd.grid_coordinates(region=region, spacing=GRID_SPACING)

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

    l4, = ax.plot(grx,gry, '.',markersize=4,markeredgecolor='k',markerfacecolor='none')
    l2 = ax.scatter(np.array(terremoto_lon_lst),np.array(terremoto_lat_lst),s=np.array(terremoto_mag_lst)*10,c=np.array(terremoto_year_lst),vmin=2000,vmax=2020, cmap='RdGy_r')

    axins = inset_axes(ax,
                   width="20%",  # width = 10% of parent_bbox width
                   height="2.5%",  # height : 50%
                   loc='upper left',
                   bbox_to_anchor=(0.6, -0.95, 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
    plt.colorbar(l2, cax=axins, orientation="horizontal", ticklocation='top',label='Ano dos eventos')
    kw = dict(prop="sizes", num=4, color='k', fmt="{x:.1f}",func=lambda s:s/10)
    legend2 = ax.legend(*l2.legend_elements(**kw), loc="upper right", title="Magnitude dos eventos")
    ax.add_artist(legend2)

    l3, = ax.plot(sta_OBS_lon,sta_OBS_lat, 's',markersize=7,markeredgecolor='k',markerfacecolor='w')
    l1, = ax.plot(sta_terra_lon,sta_terra_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey')

    legend = ax.legend([l1,l3,l4],['Estações RSBR','Possíveis OBSs','Grade de eventos'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)
    ax.add_artist(legend)

    ax.set_xticks(np.arange(LLCRNRLON,URCRNRLON,4))
    ax.set_yticks(np.arange(LLCRNRLAT,URCRNRLAT,4))
    ax.tick_params(labeltop=True, labelbottom=True,labelright=True,labelleft=True,labelsize=12)
    ax.grid(color='grey', linestyle='--', linewidth=1)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_title('Grade para calcular o gap azimutal',fontsize=15,y=1.05)

    if PLOT_STATS == True:
        plt.show()
    else:
        pass

#------------------------------------------

def calc_azimuthal_gap(input_lst):
    #calculating the azimuth between the possible earthquakes and the stations (land and obs):
    gx = input_lst[0][0]
    gy = input_lst[0][1]
    lst_sta_obs = input_lst[1]

    azimuth_lst = []
    for i,sta_loc in enumerate(lst_sta_obs):
        dist,azim,bazim = op.geodetics.base.gps2dist_azimuth(gy,gx, sta_loc[1], sta_loc[0], a=6378137.0, f=0.0033528106647474805)
        dist_degree = op.geodetics.base.kilometer2degrees(dist/1000)
        if dist_degree <= GRID_DISTANCE:
            azimuth_lst.append(round(azim))

    return [[gx,gy],sorted(azimuth_lst)]

# ============
# Main program
# ============

print(' ------------------------------ ')
print(' Importing variables and tables ')

#Coordenadas dos terremotos
terremoto_df = pd.read_csv(EARTHQUAKES_SHAPEFILE)
terremoto_lon = terremoto_df['LONG.']
terremoto_lat = terremoto_df['LAT.']
terremoto_mag = terremoto_df['Io']
terremoto_year = terremoto_df['YEAR']

terremoto_lon_lst = []
terremoto_lat_lst = []
terremoto_mag_lst = []
terremoto_year_lst = []
for i,j in enumerate(terremoto_mag):
    try:
        tmp = float(j)
        if float(j) > 0.1:
            terremoto_lon_lst.append(terremoto_lon[i])
            terremoto_lat_lst.append(terremoto_lat[i])
            terremoto_mag_lst.append(tmp)
            terremoto_year_lst.append(terremoto_year[i])
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

#Creating the grid:
grx,gry = grid_maker(LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT)

#Plotting the map:
plot_grid_obs(True)

#Creating a list with input variables:
sta_xy = []
for t,u in enumerate(sta_terra_lon):
    sta_xy.append([sta_terra_lon[t],sta_terra_lat[t]])

obs_xy = []
for t,u in enumerate(sta_OBS_lon):
    obs_xy.append([sta_OBS_lon[t],sta_OBS_lat[t]])

sta_obs_xy = sta_xy+obs_xy

input_var = []
for t,u in enumerate(grx):
    input_var.append([[grx[t],gry[t]],sta_obs_xy])


#Calculating the azimuth gap between grid and stations (land/obs):

azimuthal_test_lst = []
with Pool(processes=NUM_PROCESSES) as p:
    max_ = len(input_var)
    with tqdm(total=max_, desc='Grid loop') as pbar:
        for i, result in enumerate(p.imap_unordered(calc_azimuthal_gap, input_var)):
            azimuthal_test_lst.append(result)
            pbar.update()

azimuthal_gap = []
for i,j in enumerate(azimuthal_test_lst):
    if len(j[1]) > 1:
        if j[1][0] == 0:
            pass
        else:
            j[1].insert(0, 0)

        if j[1][-1] == 360:
            pass
        else:
            j[1].append(360)

        azimuthal_gap.append(max(np.diff(j[1])))
    else:
        azimuthal_gap.append(360)

#Gridding data:
xi = np.linspace(min(grx), max(grx), 400)
yi = np.linspace(min(gry), max(gry), 400)
X, Y = np.meshgrid(xi, yi)

grid = griddata((np.array(grx), np.array(gry)), np.array(azimuthal_gap), (X, Y), method='nearest')

#-------------------------------------------------------------------

#Figure grid interpolation:
fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(20,10))

#Add srtm image
gdal.UseExceptions()
ds = gdal.Open(SRTM_FILE)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()

inproj = osr.SpatialReference()
inproj.ImportFromWkt(proj)

# Add the interpolated grid:
img2 = ax.imshow(grid,extent=[LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT],cmap='Spectral',origin='lower',vmin=0,vmax=360,alpha=0.75,zorder=-10,interpolation='kaiser')

axins = inset_axes(ax,
               width="20%",  # width = 10% of parent_bbox width
               height="2.5%",  # height : 50%
               loc='upper left',
               bbox_to_anchor=(0.6, -0.95, 1, 1),
               bbox_transform=ax.transAxes,
               borderpad=0,
               )
cb = plt.colorbar(img2, cax=axins, orientation="horizontal", ticklocation='top',label='Gap azimutal')
cb.ax.plot(180, 2, 'k|',ms=50)
cb.ax.annotate('180'+chr(176), xy=(180, 1), xytext=(145, 1))

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
ax.set_title('Cobertura azimutal da RSBR e da rede de OBSs para terremotos de m$_{L}$>'+str(EARTHQUAKES_MIN),fontsize=15,y=1.05)

l2, = ax.plot(sta_OBS_lon,sta_OBS_lat, 's',markersize=10,markeredgecolor='k',markerfacecolor='w',zorder=10)
l1, = ax.plot(sta_terra_lon,sta_terra_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey',zorder=1)
legend = ax.legend([l1,l2],['Estações RSBR','Possíveis OBSs'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)

plt.show()
