#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 09:20:37 2025

@author: dlm356
"""

import os, sys
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import colormaps
list(colormaps)

import pandas as pd
import xarray as xr

import numpy as np
from functools import partial

import geocube

from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_points_griddata, rasterize_points_radial
import rioxarray
# rasterizing points
# https://corteva.github.io/geocube/stable/examples/rasterize_point_data.html

# Summary variables at 100 m grid for Chenghao

# Land cover and trees from 2021 lidar
# https://zenodo.org/records/14053441
# doing this in R

# mean tree height from 2021 lidar
# Need to try CHM for all NYC in R

#%%

# Need to generate crown size weighted height at 100 m (LCZ map resolution)

with rasterio.open("/Users/dlm356/dlm356_files/nyc_trees/ucm/new_york_CONUS_LCZ_map.tif") as src:
    lcz = src.read(1) # this is just an array
    lcz_meta = src.meta
    lcz_crs = src.crs 
    res_x, res_y = src.res


plt.imshow(lcz)




#%%
# mean crown size from the tnc polygons
# Do this here

# load tif file basemap
lcz_file = rasterio.open("/Users/dlm356/dlm356_files/nyc_trees/ucm/new_york_CONUS_LCZ_map.tif")
lcz = lcz_file.read(1)
lcz_meta = lcz_file.meta
lcz.shape
lcz_file.crs

plt.imshow(lcz)

# load tree GDB
tree_gdb_path = "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb"
tree_polys = gpd.read_file(tree_gdb_path, layer = 'treeobjects_2021_nyc')

# Use centroids, the SHAPE_Area column for each of this
tree_centroids = gpd.GeoDataFrame(tree_polys, geometry = tree_polys.centroid, crs = tree_polys.crs)
tree_centroids.plot("SHAPE_Area", cmap = "viridis") # most trees are small, few big ones

# reproject to lcz crs
tree_centroids_reproj = tree_centroids.to_crs(crs = lcz_file.crs)
tree_centroids_reproj.plot("SHAPE_Area", cmap = "viridis") # most trees are small, few big ones

#res_x, res_y = lcz_file.res
lcz_rast_riox = rioxarray.open_rasterio("/Users/dlm356/dlm356_files/nyc_trees/ucm/new_york_CONUS_LCZ_map.tif")

out_grid = make_geocube(vector_data = tree_centroids_reproj,
                        measurements = ["SHAPE_Area"],
                        fill = 0,
                        like = lcz_rast_riox) # takes CRS, transform, resolution, extent from reference, needs rioxarray object for "like" call instead of regular rasterio


# Now no longer need to do this exactly...


#%% Junk that didn't work
# by default this does average

#averaged_raster = out_grid["SHAPE_Area"]

#out_grid.rio.to_raster("/Users/dlm356/dlm356_files/nyc_trees/ucm/tree_shape_area_mean.tif")

# this didn't take the mean, just took the last point
# redoing with a new zonal statistics workaround, it will be clunky

# lcz_vec = geocube.vector.vectorize(lcz_rast_riox)

# this is chatgpt to do the averaging outside of the raster as df and then turn back into raster
# # Convert the DataArray to a long-form DataFrame
# df = out_grid["SHAPE_Area"].to_dataframe().reset_index()
# df = df.drop(columns=["spatial_ref"])
# mean_df = df.groupby(["x", "y"], as_index=False).mean()

# # Get original raster shape and coordinates
# da = out_grid["SHAPE_Area"]
# shape = da.shape
# x_coords = da.coords["x"].values
# y_coords = da.coords["y"].values

# # Create an empty array filled with NaN
# mean_raster = np.full(shape, np.nan)

# # Build a quick lookup from coordinates to array indices
# coord_index = { (x, y): (i, j) 
#                 for i, y in enumerate(y_coords) 
#                 for j, x in enumerate(x_coords) }

# for _, row in mean_df.iterrows():
#     i, j = coord_index[(row["x"], row["y"])]
#     mean_raster[i, j] = row["SHAPE_Area"]
    
# mean_da = xr.DataArray(
#     mean_raster,
#     coords={"y": y_coords, "x": x_coords},
#     dims=("y", "x"),
#     name="SHAPE_Area"
# )

# mean_da = mean_da.rio.write_crs(out_grid.rio.crs)
# mean_da.rio.to_raster("/Users/dlm356/dlm356_files/nyc_trees/ucm/tree_shape_area_mean_3.tif")

# # checking
# 409.748
# 297.255

# This didn't work either

# # Chatgpt again
# # Load points and reference raster
# #points = gpd.read_file("tree_centroids_reproj.shp")
# points = tree_centroids_reproj
# ref_raster = rioxarray.open_rasterio("/Users/dlm356/dlm356_files/nyc_trees/ucm/new_york_CONUS_LCZ_map.tif")

# # Get raster info
# res_x, res_y = ref_raster.rio.resolution()
# x0, y0 = ref_raster.rio.bounds().left, ref_raster.rio.bounds().top # this breaks too
# ncols, nrows = ref_raster.rio.width, ref_raster.rio.height

# # Map each point to a raster row/col
# points['col'] = ((points.geometry.x - x0) // res_x).astype(int)
# points['row'] = ((y0 - points.geometry.y) // abs(res_y)).astype(int)

# # Compute mean per cell
# mean_df = points.groupby(['row', 'col'])['SHAPE_Area'].mean().reset_index()

# # Create empty array
# mean_raster = np.full((nrows, ncols), np.nan)

# # Fill raster
# for _, r in mean_df.iterrows():
#     mean_raster[r.row, r.col] = r.SHAPE_Area

# # Convert to DataArray
# mean_da = xr.DataArray(
#     mean_raster,
#     dims=('y', 'x'),
#     coords={
#         'y': ref_raster.y,
#         'x': ref_raster.x
#     },
#     name='SHAPE_Area'
# )

# # Assign CRS and export
# mean_da.rio.write_crs(ref_raster.rio.crs, inplace=True)
# mean_da.rio.to_raster("/Users/dlm356/dlm356_files/nyc_trees/ucm/tree_shape_area_mean_4.tif")

# Chatgpt is lost, no shortcuts here

#%%
# mean tree fraction (from whatever dataset is easiest)
# Doing this in R

# grass fraction
# Doing this in R