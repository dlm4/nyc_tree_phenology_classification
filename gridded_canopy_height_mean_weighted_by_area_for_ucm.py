#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 14:51:38 2025

@author: dlm356
"""

import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

with rasterio.open("/Users/dlm356/dlm356_files/nyc_trees/ucm/new_york_CONUS_LCZ_map.tif") as src:
    lcz = src.read(1) # this is just an array
    lcz_meta = src.meta
    lcz_crs = src.crs 
    res_x, res_y = src.res


plt.imshow(lcz)

# https://gis.stackexchange.com/questions/436022/finding-the-centroid-of-every-pixel-in-a-raster-python
path = '/Users/dlm356/dlm356_files/nyc_trees/ucm/new_york_CONUS_LCZ_map.tif'

src = rasterio.open(path)

with rasterio.open(path) as src:
    band1 = src.read(1)
    height = band1.shape[0]
    width = band1.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rasterio.transform.xy(src.transform, rows, cols)
    lons = np.array(xs)
    lats = np.array(ys)

    lon_lat_list = list(zip(lons.flatten(), lats.flatten()))
    
    lon_lat_point = [Point(x,y) for x,y in lon_lat_list] # had to modify this
    
    points = gpd.GeoSeries(lon_lat_point) # and this, .map(Point) doesn't work for me

    # use the feature loop in case shp is multipolygon
    geoms = points.values
    features = [i for i in range(len(geoms))]

    out = gpd.GeoDataFrame({'feature': features, 'geometry': geoms}, crs=src.crs)

out.plot("feature", cmap = "viridis")


# load tree GDB
tree_gdb_path = "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb"
tree_polys = gpd.read_file(tree_gdb_path, layer = 'treeobjects_2021_nyc')

# Use centroids, the SHAPE_Area column for each of this
tree_centroids = gpd.GeoDataFrame(tree_polys, geometry = tree_polys.centroid, crs = tree_polys.crs)

# reproject to lcz crs
tree_centroids_reproj = tree_centroids.to_crs(crs = src.crs)

# for each grid cell
# get x and y
# add res_x and res_y ranges to get min_x, max_x, min_y, max_y bounds
# get which tree points are within these bounds (logical)
# take weighted mean of Height weighted by SHAPE_Area
# Note that these are in feet, convert to m (/3.28)
# Note: will want to swap in our estimated height instead from lidar metrics