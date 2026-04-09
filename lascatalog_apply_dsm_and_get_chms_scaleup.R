

# Need to create digital surface model from lidar tiles
# Or directly created "urban canopy height model" which includes buildings
# Then use existing 2021 land cover map, aggregated at this scale 1.64 ft (0.5 m) to extract:
# Buildings for building surface model
# Trees for canopy height model

# Create DTM for all of NYC
# 1.64 ft spatial resolution (0.5 m, approx)

library(tidyverse)
library(sf)
library(lidR)
library(future)
library(terra) # for rasters
# Note: terra masks from lidR: area, crs, crs <-, is.empty, watershed


pts_classes <- c(LASGROUND, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASBUILDING,
                 LASWATER, LASRAIL, LASROADSURFACE, LASWIREGUARD, LASWIRECONDUCTOR, LASTRANSMISSIONTOWER,
                 LASBRIGDE) # add to readLAS filter instead of filterPOI, bridge typo is in package
filter_string <- paste("-drop_z_below 0 -keep_class", paste(as.character(sort(pts_classes)), collapse = " "))

grid_size <- 1.64 # ft approx 0.5 # half meter chm

# las catalog
setwd("/Volumes/DLM_backup/lidar_2021/NYC_2021/")
file_list <- list.files(pattern = glob2rx("*.las"))

test_ids <- file_list %>% strsplit('[.]') %>% sapply(`[`, 1)
test_files <- paste("/Volumes/DLM_backup/lidar_2021/NYC_2021/", test_ids, '.las', sep="")

ctg <- readLAScatalog(test_files, filter = filter_string)

plan(multisession, workers = 8) # don't use all the cores, easier on memory
temp_output_path <- "/Volumes/NYC_geo/processing_temporary/" # do this on the external disk

# set up new run for surface
opt_output_files(ctg) <- paste(temp_output_path, "surface_raster_1p64ft_{ORIGINALFILENAME}", sep = "")
opt_chunk_buffer(ctg) <- 328 # units in feet, approx 100 m

grid_size <- 1.64 # 0.5 # in feet instead of meters

rasterizeSurface <- function(chunk, grid_size) {
  
  thresh_ft <- c(0, 2, 5, 10, 15)*3.28
  max_edge_ft <- c(0, 1)*3.28
  
  # Load the chunk + buffer
  las <- readLAS(chunk)
  if (lidR::is.empty(las)) return(NULL)
  if (any(las$Classification == 2)){ # test for requirement that there be some bare ground (LASGROUND = 2)
    # do something
    output <- rasterize_canopy(chunk, grid_size, pitfree(thresholds = thresh_ft, max_edge = max_edge_ft))
    # this is called "canopy" but actually creates a digital surface model, DSM
    # will need to subtract terrain to create a urban canopy model (tree + building)
    
    # remove the buffer of the output
    output <- crop(output, st_bbox(chunk)) #crop for spatraster
    return(output) # do we want to temporarily write out the output instead and retain as a las catalog
  }
}

surface_raster_files <- catalog_apply(ctg, rasterizeSurface, grid_size)
# Note that not necessarily all of these files actually run, so need to do only do analysis for files that processed

gc()

#####
plan(sequential)
#surface_files <- surface_raster_files %>% sapply(`[`, 1) # might need to do a list.files() step instead
setwd("/Volumes/NYC_geo/processing_temporary/")
surface_files <- list.files(pattern = glob2rx("surface_raster*"))
#kept_ids <- surface_files %>% strsplit('[._]') %>% sapply(`[`, 6)
kept_ids <- surface_files %>% strsplit('[._]') %>% sapply(`[`, 4)
test_ids <- kept_ids
terrain_files <- paste0("/Volumes/NYC_geo/nyc_lidar_metrics/terrain_rasters/terrain_raster_1p64ft_", test_ids, ".tif")

# # Subtract terrain from surface to produce urban building + tree canopy height model
# surface_files <- paste0("/Volumes/NYC_geo/processing_temporary/surface_raster_1p64ft_", test_ids, ".tif")
# terrain_files <- paste0("/Volumes/NYC_geo/nyc_lidar_metrics/terrain_rasters/terrain_raster_1p64ft_", test_ids, ".tif")

for(i in 1:length(test_ids)){
  print(i)
  surface_rast <- rast(surface_files[i])
  terrain_rast <- rast(terrain_files[i])
  uchmbt_rast <- surface_rast - crop(terrain_rast, surface_rast) # uchmbt = urban canopy height model building terrain
  uchmbt_rast[uchmbt_rast < 0] <- 0
  writeRaster(uchmbt_rast,
              paste0("/Volumes/NYC_geo/processing_temporary/uchmbt_raster_1p64ft_", test_ids[i], ".tif"))
  rm(uchmbt_rast)
}

gc()

#####
# Take land cover data from 2021 6 inch
# Reproject (if needed)
# Subset to tile (with a 328 ft buffer)
# Coarsen to 1.64 ft as well (save this)
# take uchmbt file
# Mask uchmbt file to only keep:
# A. Buildings
# B. Trees
# Everything else becomes 0
# Save:
# A. Building height model
# B. Tree canopy height model
# Both will be in feet but can be converted to meters

lc_rast <- rast("/Volumes/NYC_geo/nyc_land_cover/landcover_nyc_2021_6in.tif")

surface_rast_reproj <- project(surface_rast, crs(lc_rast))

las_index_vect <- vect("/Volumes/DLM_backup/lidar_2021/NYC_2021/NYC2021_LAS_Index.shp")
las_index_vect_reproj <- project(las_index_vect, crs(lc_rast))

cls <- vect("/Volumes/NYC_geo/tree_classification/outputs_practical_v2_monthcomp/nyc_class_tree_genus_polygons_v2.gpkg")


# these have min of NA, actually want min of zero with NA outside of the study area
for (i in 1:length(test_ids)){
  #i <- 448
  print(i)
  las_index_vect_reproj_tile <- las_index_vect_reproj[las_index_vect_reproj$LAS_ID == test_ids[i]]
  las_index_vect_reproj_tile_buff <- buffer(las_index_vect_reproj_tile, width = 328) # 328 ft buffer, 100 m
  lc_rast_tile_buff <- crop(lc_rast, las_index_vect_reproj_tile_buff)
  lc_rast_tile_buff_reproj <- project(lc_rast_tile_buff, crs(surface_rast), method = "near")
  surface_rast <- rast(surface_files[i])
  lc_rast_tile_buff_reproj_resample <- resample(lc_rast_tile_buff_reproj, surface_rast, method = "modal")
  writeRaster(lc_rast_tile_buff_reproj_resample,
              paste0("/Volumes/NYC_geo/processing_temporary/landcover2021_raster_1p64ft_", test_ids[i], ".tif"))
  
  uchmbt_rast <- rast(paste0("/Volumes/NYC_geo/processing_temporary/uchmbt_raster_1p64ft_", test_ids[i], ".tif"))
  building_rast <- uchmbt_rast
  building_rast[lc_rast_tile_buff_reproj_resample != 5] <- NA # building only filter, need another filter to screen out areas beyond the land cover classification due to buffering
  building_rast[is.na(lc_rast_tile_buff_reproj_resample)] <- NA
  bmin <- minmax(building_rast)[1]
  bmax <- minmax(building_rast)[2]
  if (!is.na(bmin) & !is.na(bmax)){
    writeRaster(building_rast,
                paste0("/Volumes/NYC_geo/processing_temporary/buildingheight_raster_1p64ft_", test_ids[i], ".tif"), overwrite = TRUE)
  }
  tree_rast <- uchmbt_rast
  cls_tile_buff <- crop(cls, las_index_vect_reproj_tile_buff)
  cls_tile_buff_reproj <- project(cls_tile_buff, crs(uchmbt_rast))

  if (dim(cls_tile_buff_reproj)[1] > 0){
    cls_tile_height_rast <- rasterize(cls_tile_buff_reproj, uchmbt_rast, field = "height_max_ft")
    height_diff_rast <- cls_tile_height_rast - uchmbt_rast
    height_diff_rast[height_diff_rast < 0] <- NA
    tree_rast <- mask(uchmbt_rast, mask = height_diff_rast)

    writeRaster(tree_rast,
                paste0("/Volumes/NYC_geo/processing_temporary/treeheight_raster_1p64ft_", test_ids[i], ".tif"))

    tree_rast[tree_rast < 6.56] <- NA

    writeRaster(tree_rast,
                paste0("/Volumes/NYC_geo/processing_temporary/treeheightfilt_min6p56ft_raster_1p64ft_", test_ids[i], ".tif"))
  }
}

#####
# re do with min of zero
#for (i in 1:length(test_ids)){
for (i in 6:length(test_ids)){
  print(i)
  las_index_vect_reproj_tile <- las_index_vect_reproj[las_index_vect_reproj$LAS_ID == test_ids[i]]
  las_index_vect_reproj_tile_buff <- buffer(las_index_vect_reproj_tile, width = 328) # 328 ft buffer, 100 m
  lc_rast_tile_buff <- crop(lc_rast, las_index_vect_reproj_tile_buff)
  lc_rast_tile_buff_reproj <- project(lc_rast_tile_buff, crs(surface_rast), method = "near")
  surface_rast <- rast(surface_files[i])
  lc_rast_tile_buff_reproj_resample <- resample(lc_rast_tile_buff_reproj, surface_rast, method = "modal")
  # writeRaster(lc_rast_tile_buff_reproj_resample,
  #             paste0("/Volumes/NYC_geo/processing_temporary/landcover2021_raster_1p64ft_", test_ids[i], ".tif"))
  # don't need to write this out again
  
  uchmbt_rast <- rast(paste0("/Volumes/NYC_geo/processing_temporary/uchmbt_raster_1p64ft_", test_ids[i], ".tif"))
  building_rast <- uchmbt_rast
  building_rast[lc_rast_tile_buff_reproj_resample != 5] <- 0 # changed to zero for minimum value of zero
  building_rast[is.na(lc_rast_tile_buff_reproj_resample)] <- NA # keep as NA for areas beyond the edge of the map
  bmin <- minmax(building_rast)[1]
  bmax <- minmax(building_rast)[2]
  if (!is.na(bmin) & !is.na(bmax)){
    writeRaster(building_rast,
                paste0("/Volumes/NYC_geo/processing_temporary/buildingheight_raster_minzero_1p64ft_", test_ids[i], ".tif"), overwrite = TRUE)
  }
  tree_rast <- uchmbt_rast
  cls_tile_buff <- crop(cls, las_index_vect_reproj_tile_buff)
  cls_tile_buff_reproj <- project(cls_tile_buff, crs(uchmbt_rast))
  
  if (dim(cls_tile_buff_reproj)[1] > 0){
    cls_tile_height_rast <- rasterize(cls_tile_buff_reproj, uchmbt_rast, field = "height_max_ft")
    height_diff_rast <- cls_tile_height_rast - uchmbt_rast
    height_diff_rast[height_diff_rast < 0] <- NA
    tree_rast <- mask(uchmbt_rast, mask = height_diff_rast, updatevalue = 0) # change this to make the minimum value of zero
    tree_rast[is.na(lc_rast_tile_buff_reproj_resample)] <- NA # keep as NA for areas beyond the edge of the map
    
    writeRaster(tree_rast,
                paste0("/Volumes/NYC_geo/processing_temporary/treeheight_raster_minzero_1p64ft_", test_ids[i], ".tif"))
    
    tree_rast[tree_rast < 6.56] <- 0 # changed this to make minimum value of zero, check output to make sure this doesn't replace any NAs
    
    writeRaster(tree_rast,
                paste0("/Volumes/NYC_geo/processing_temporary/treeheightfilt_min6p56ft_minzero_raster_1p64ft_", test_ids[i], ".tif"))
  }
}

#####
# these are with min of NA, great for visualization in QGIS, but not good for aggregating
# this worked, now can create mosaics
# moving things over to: "/Volumes/NYC_geo/nyc_lidar_metrics/"

# land cover
setwd("/Volumes/NYC_geo/processing_temporary/")
file_prefix <- "landcover2021_raster_1p64ft"
file_list <- list.files(pattern = glob2rx(paste0(file_prefix, "*")))
rastcoll <- sprc(file_list)
merge(rastcoll, method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_mosaic.tif"))

# Building height
setwd("/Volumes/NYC_geo/processing_temporary/")
#file_prefix <- "buildingheight_raster_1p64ft"
#file_list <- list.files(pattern = glob2rx(paste0(file_prefix, "*")))
#rastcoll <- sprc(file_list)
# merge(rastcoll[1:50], method = "nearest",
#       filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_mosaic.tif"),
#       overwrite = TRUE)
# 
# merge(rastcoll[1:100], method = "nearest",
#       filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_mosaic.tif"),
#       scale = 0.01, datatype = "INT4S", overwrite = TRUE)
# mosaicking broke with a small empty tile that was all NaN, will need to filter this

# # Need to filter to remove any building tiles that are all NA
# rast_min <- rep(0, length(rastcoll))
# rast_max <- rep(0, length(rastcoll))
# for (i in 1:length(rastcoll)){
#   print(i)
#   mm <- minmax(rastcoll[i])
#   rast_min[i] <- mm[1]
#   rast_max[i] <- mm[2]
# }
# building_rast_minmax <- cbind.data.frame(rast_min, rast_max)
# write.csv(building_rast_minmax, "buildingheight_raster_minmax_vals.csv", row.names = FALSE)
file_prefix <- "buildingheight_raster_minzero_1p64ft"
file_list <- list.files(pattern = glob2rx(paste0(file_prefix, "*")))
rastcoll <- sprc(file_list)
#building_rast_minmax <- read.csv("buildingheight_raster_minmax_vals.csv")
#building_rast_withvals <- which(!is.na(building_rast_minmax$rast_min) & !is.na(building_rast_minmax$rast_max))
#rastcoll <- sprc(file_list[building_rast_withvals])
# mosaic in bigger chunks first, then do everything after. Too many tiles to do everything at once efficiently without risking crashing
merge(rastcoll[1:200], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1_200_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[201:400], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_201_400_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[401:600], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_401_600_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[601:800], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_601_800_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[801:1000], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_801_1000_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1001:1200], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1001_1200_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1201:1400], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1201_1400_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1401:1600], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1401_1600_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1601:1690], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1601_1690_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)

# mosaic everything together after this
setwd("/Volumes/NYC_geo/nyc_lidar_metrics/")
file_prefix <- "buildingheight_raster_minzero_1p64ft_scaled_"
file_list <- list.files(pattern = glob2rx(paste0(file_prefix, "*")))
rastcoll <- sprc(file_list)
merge(rastcoll, method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)


# Trees
# filt
setwd("/Volumes/NYC_geo/processing_temporary/")
#file_prefix <- "treeheightfilt_min6p56ft_minzero_raster_1p64ft"
file_prefix <- "treeheight_raster_minzero_1p64ft"
file_list <- list.files(pattern = glob2rx(paste0(file_prefix, "*")))
rastcoll <- sprc(file_list)
merge(rastcoll[1:200], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1_200_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[201:400], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_201_400_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[401:600], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_401_600_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[601:800], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_601_800_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[801:1000], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_801_1000_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1001:1200], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1001_1200_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1201:1400], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1201_1400_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1401:1600], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1401_1600_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
merge(rastcoll[1601:1674], method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "_scaled_1601_1674_mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)

# now merge all these together into one big mosaic
setwd("/Volumes/NYC_geo/nyc_lidar_metrics/")
#file_prefix <- "treeheightfilt_min6p56ft_minzero_raster_1p64ft_scaled_"
file_prefix <- "treeheight_raster_minzero_1p64ft_scaled_"
file_list <- list.files(pattern = glob2rx(paste0(file_prefix, "*")))
rastcoll <- sprc(file_list)
merge(rastcoll, method = "nearest",
      filename = paste0("/Volumes/NYC_geo/nyc_lidar_metrics/", file_prefix, "mosaic.tif"),
      scale = 0.01, datatype = "INT4S", overwrite = TRUE)
# changing name to move minzero after raster label

# test read in
chm <- rast('/Volumes/NYC_geo/nyc_lidar_metrics/treeheightfilt_min6p56ft_raster_minzero_1p64ft_scaled_mosaic.tif')
