#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 14:25:49 2025

@author: dlm356
"""

import geopandas as gpd
import os, sys
import matplotlib.pyplot as plt
import rasterio

# Read in NYC boroughs file with gpd

boros = gpd.read_file('/Volumes/NYC_geo/vectors/Borough Boundaries/geo_export_da133389-a6c6-45c3-a980-14295f0e4c2f.shp')

boros.plot("boro_name", cmap = "Set1", legend = True)

# Make plot with boro names centered
ax = boros.plot("boro_name", cmap = "Set1")
for idx, row in boros.iterrows():
    ax.text(row.geometry.centroid.x, row.geometry.centroid.y, row['boro_name'], fontsize=8, horizontalalignment = "center")

# calculates area, but in units of lat lon degrees so this doesn't actually make any sense
boros.area

boros_reproj = boros.to_crs(epsg = '26918') # Reproject to UTM Zone 18
boros_reproj.plot("boro_name", cmap = "Set1", legend = True)

boros_reproj.convex_hull.plot()

# load planet imagery
src = rasterio.open("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal/nyc_planet_composite_8band_20220710_nyccal20220710ref.tif")
array = src.read(8) # starts at 1
array.shape
