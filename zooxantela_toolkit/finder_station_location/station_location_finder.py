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

import datetime
import matplotlib.dates as mdates

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

STATION_SPACING = 2 #degrees
STATION_MAX_DISTANCE = 2 #degrees
NUMBER_STATIONS = 10

NUM_PROCESSES = 8

FILTER_GRID_BY_SHAPEFILE = True

AZIMUTHAL_GAP_VALUE = 180
SECONDARY_AZIMUTHAL_GAP_VALUE = 200

LLCRNRLON = -52 #degrees
URCRNRLON = -28 #degrees
LLCRNRLAT = -28 #degrees
URCRNRLAT = -16 #degrees

STA_CONTINENT_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/rsbr_lat_lon_projeto.csv'
OBS_LOCATION_FILE = '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/coord_potential_obs_location.csv'
COAST_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/brasil_estados/brasil_estados.shp'
OCEANO_SHAPEFILE = '/home/diogoloc/SIG_dados/Projeto_ON_MAR/shapefile/area_batimetria/area_ON_projeto.shp'
EARTHQUAKES_SHAPEFILE =  '/home/diogoloc/dados_posdoc/ON_MAR/sta_coord/coord_terremotos_marinhos_brasil.csv'
SRTM_FILE =  '/home/diogoloc/SIG_dados/Projeto_ON_MAR/raster/projeto_dem/projeto_srtm.tif'

#--------------------------------------------------------------------------------------------------------------------


# =========
# Functions
# =========

#Function for creating GRID POINTS
def grid_maker_circle(data_longitude, data_latitude, west, east, south, north):
    region = (west, east, south, north)

    # Generate the coordinates for a regular grid mask
    coordinates = vd.grid_coordinates(region, spacing=STATION_SPACING)

    # Generate a mask for points that are more than 2 grid spacings away from any data
    # point. The mask is True for points that are within the maximum distance. Distance
    # calculations in the mask are Cartesian only. We can provide a projection function to
    # convert the coordinates before distances are calculated (Mercator in this case). In
    # this case, the maximum distance is also Cartesian and must be converted from degrees
    # to meters.
    mask = vd.distance_mask(
        (data_longitude, data_latitude),
        maxdist=STATION_MAX_DISTANCE*111e3,
        coordinates=coordinates,
        projection=pyproj.Proj(proj="merc", lat_ts=np.mean(data_latitude)),
    )

    # Create a dummy grid with ones that we can mask to show the results.
    # Turn points that are too far into NaNs so they won't show up in our plot.
    dummy_data = np.ones_like(coordinates[0])
    dummy_data[~mask] = np.nan

    rows, cols = coordinates

    grdx = rows.ravel()
    grdy = cols.ravel()
    grdz = dummy_data.ravel()

    gridx = []
    gridy = []
    for x,c in enumerate(grdz):
        if np.isnan(c) == False:
            gridx.append(grdx[x])
            gridy.append(grdy[x])

    return gridx,gridy

#------------------------------------------

#Function for creating GRID POINTS
def grid_maker(west, east, south, north):
    region = (west, east, south, north)

    region = (west, east, south, north)
    rows, cols = vd.grid_coordinates(region=region, spacing=STATION_SPACING)

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

#------------------------------------------

def plot_grid_obs_circle(x_OBS,y_OBS,grx,gry):
    #Figure
    fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(10,10))

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

    l3, = ax.plot(grx,gry, 's',markersize=5,markeredgecolor='k',markerfacecolor='k')
    l2, = ax.plot(terremoto_lon_lst,terremoto_lat_lst, '*',markersize=5,markeredgecolor='y',markerfacecolor='y')
    l1, = ax.plot(sta_terra_lon,sta_terra_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey')

    circulo = Circle(radius=STATION_MAX_DISTANCE,xy=(sta_OBS_lon[i], sta_OBS_lat[i]),color='none', ec='k',linewidth=0.25,transform=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None),zorder=2)
    ax.add_patch(circulo)

    legend = ax.legend([l1,l2,l3,circulo],['Estações RSBR','Terremotos m$_{L}$>3','OBSs grid','Raio de '+str(round(STATION_MAX_DISTANCE))+' graus'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)
    ax.set_xticks(np.arange(LLCRNRLON,URCRNRLON,4))
    ax.set_yticks(np.arange(LLCRNRLAT,URCRNRLAT,4))
    ax.tick_params(labeltop=True, labelbottom=True,labelright=True,labelleft=True,labelsize=15)
    ax.grid(color='grey', linestyle='--', linewidth=1)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    plt.show()

#------------------------------------------

def plot_grid_obs(grx,gry):
    #Figure
    fig_grid, ax = plt.subplots(ncols=1, subplot_kw={'projection': ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None)},figsize=(10,10))

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
    reader_1_SHP = Reader(OCEANO_SHAPEFILE)
    shape_1_SHP = list(reader_1_SHP.geometries())
    plot_shape_1_SHP = cfeature.ShapelyFeature(shape_1_SHP, ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))
    ax.add_feature(plot_shape_1_SHP, facecolor='none', edgecolor='k',linewidth=1)

    ax.set_extent([LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT], crs=ccrs.Mercator(central_longitude=(LLCRNRLON-URCRNRLON)/2+LLCRNRLON, globe=None))

    l3, = ax.plot(grx,gry, 's',markersize=5,markeredgecolor='k',markerfacecolor='k')
    l2, = ax.plot(terremoto_lon_lst,terremoto_lat_lst, '*',markersize=5,markeredgecolor='y',markerfacecolor='y')
    l1, = ax.plot(sta_terra_lon,sta_terra_lat, '^',markersize=10,markeredgecolor='k',markerfacecolor='grey')

    legend = ax.legend([l1,l2,l3],['Estações RSBR','Terremotos m$_{L}$>3','OBSs grid'],scatterpoints=1, framealpha=1,labelspacing=1, loc='lower right',facecolor='w',fontsize='medium',markerscale=1.5)
    ax.set_xticks(np.arange(LLCRNRLON,URCRNRLON,4))
    ax.set_yticks(np.arange(LLCRNRLAT,URCRNRLAT,4))
    ax.tick_params(labeltop=True, labelbottom=True,labelright=True,labelleft=True,labelsize=15)
    ax.grid(color='grey', linestyle='--', linewidth=1)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    plt.show()

#------------------------------------------

def calc_angle(a,b,c):
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

#------------------------------------------

def calc_azimuthal_gap(d):

    #calculating azimuthal gap:
    comb1 = list(combinations(d, 3))
    azimuthal_gap_lst = []
    for k,l in enumerate(comb1):
        a = np.array([grx[l[0]],gry[l[0]]])
        b = np.array([grx[l[1]],gry[l[1]]])
        c = np.array([grx[l[2]],gry[l[2]]])
        azimuthal_gap_lst.append(calc_angle(a,b,c))
    azimuthal_gap = max(azimuthal_gap_lst)

    #calculating secondary azimuthal gap:
    secondary_azimuthal_gap_lst = []
    del_sta_lst = [del_sta for del_sta in range(len(d))]
    for del_index in del_sta_lst:
        new_d = list(d)
        del new_d[del_index]
        comb2 = list(combinations(new_d, 3))
        for k,l in enumerate(comb2):
            a2 = np.array([grx[l[0]],gry[l[0]]])
            b2 = np.array([grx[l[1]],gry[l[1]]])
            c2 = np.array([grx[l[2]],gry[l[2]]])
            secondary_azimuthal_gap_lst.append(calc_angle(a2,b2,c2))

        secondary_azimuthal_gap = max(secondary_azimuthal_gap_lst)

    if azimuthal_gap < AZIMUTHAL_GAP_VALUE and secondary_azimuthal_gap < SECONDARY_AZIMUTHAL_GAP_VALUE:
        return d

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
        if float(j) > 3:
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

print(' Calculating the network layout ')
print('\n')

#Creating the grid:
grx,gry = grid_maker(LLCRNRLON,URCRNRLON,LLCRNRLAT,URCRNRLAT)

#Plotting the map:
plot_grid_obs(grx,gry)

#Creating combinations between possible OBS locations:
lst_grid = [i for i in range(len(grx))]
comb = list(combinations(lst_grid, NUMBER_STATIONS))

print(' - Number of possible combinations: ',len(comb))

azimuthal_test_lst = []
with Pool(processes=NUM_PROCESSES) as p:
    max_ = len(comb)
    with tqdm(total=max_, desc='Combinations loop') as pbar:
        for i, result in enumerate(p.imap_unordered(calc_azimuthal_gap, comb)):
            azimuthal_test_lst.append(result)
            pbar.update()

for azi in tqdm(azimuthal_test_lst, desc='Approved combinations loop'):
        azi_gap_gx = [grx[i] for i in azi]
        azi_gap_gy = [gry[i] for i in azi]
        plot_grid_obs(azi_gap_gx,azi_gap_gy)
print('\n')

#for i,j in enumerate(sta_OBS_lon):
#    sta_OBS_lon[i], sta_OBS_lat[i]


#----------------------------------------------------------------------------
