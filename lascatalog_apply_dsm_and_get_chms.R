

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

#test_ids <-  "5240" # "2242"
# test_ids <- c("2242", "5242", "7242",
#               "2240", "5240", "7240",
#               "2237", "5237", "7237")
#test_files <- paste("/Volumes/DLM_backup/lidar_2021/NYC_2021/", test_ids, ".las", sep="") # will want to move this back over to SSD later

pts_classes <- c(LASGROUND, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASBUILDING,
                 LASWATER, LASRAIL, LASROADSURFACE, LASWIREGUARD, LASWIRECONDUCTOR, LASTRANSMISSIONTOWER,
                 LASBRIGDE) # add to readLAS filter instead of filterPOI, bridge typo is in package

#pts_classes <- c(LASGROUND, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASWATER, LASRAIL, LASROADSURFACE) # add to readLAS filter instead of filterPOI, bridge typo is in package

filter_string <- paste("-drop_z_below 0 -keep_class", paste(as.character(sort(pts_classes)), collapse = " "))

#las <- readLAS(test_files, filter = filter_string)

grid_size <- 1.64 # ft approx 0.5 # half meter chm

# make DTM
#dtm <- rasterize_terrain(las, grid_size, tin())

# las catalog
setwd("/Volumes/DLM_backup/lidar_2021/NYC_2021/")
#file_list <- list.files(pattern = glob2rx("*.las"))

test_ids <- c("2242", "5242", "7242",
              "2240", "5240", "7240",
              "2237", "5237", "7237")
#test_ids <- file_list
#test_files <- paste("/Volumes/DLM_backup/lidar_2021/NYC_2021/", test_ids[1:30], ".las", sep="") # will want to move this back over to SSD later

#test_files <- file_list
test_files <- paste0("/Volumes/DLM_backup/lidar_2021/NYC_2021/", test_ids, ".las")
ctg <- readLAScatalog(test_files, filter = filter_string)

plan(multisession, workers = 8) # don't use all the cores, easier on memory
temp_output_path <- "/Volumes/NYC_geo/processing_temporary/" # do this on the external disk

# set up new run for terrain, rpf is old naming convention
#ctg_rpf <- readLAScatalog(unlist(ctg))
#opt_output_files(ctg) <- paste(temp_output_path, "terrain_raster_1p64ft_{ORIGINALFILENAME}", sep = "")
opt_output_files(ctg) <- paste(temp_output_path, "surface_raster_1p64ft_{ORIGINALFILENAME}", sep = "")
opt_chunk_buffer(ctg) <- 328 # units in feet, approx 100 m

grid_size <- 1.64 # 0.5 # in feet instead of meters

# rasterizeTerrain <- function(chunk, grid_size) {
#   # Load the chunk + buffer
#   las <- readLAS(chunk)
#   if (lidR::is.empty(las)) return(NULL)
#   if (any(las$Classification == 2)){ # test for requirement that there be some bare ground (LASGROUND = 2)
#     # do something
#     output <- rasterize_terrain(chunk, grid_size, tin())
#     
#     # remove the buffer of the output
#     output <- crop(output, st_bbox(chunk)) #crop for spatraster
#     return(output) # do we want to temporarily write out the output instead and retain as a las catalog
#   }
# }
# 
# terrain_raster_files <- catalog_apply(ctg, rasterizeTerrain, grid_size)


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

terrain_raster_files <- catalog_apply(ctg, rasterizeSurface, grid_size)

#####
# Subtract terrain from surface to produce urban building + tree canopy height model
plan(sequential)
surface_files <- paste0("/Volumes/NYC_geo/processing_temporary/surface_raster_1p64ft_", test_ids, ".tif")
terrain_files <- paste0("/Volumes/NYC_geo/nyc_lidar_metrics/terrain_rasters/terrain_raster_1p64ft_", test_ids, ".tif")

for(i in 1:length(test_ids)){
  print(i)
  surface_rast <- rast(surface_files[i])
  terrain_rast <- rast(terrain_files[i])
  uchmbt_rast <- surface_rast - terrain_rast # uchmbt = urban canopy height model building terrain
  uchmbt_rast[uchmbt_rast < 0] <- 0
  writeRaster(uchmbt_rast,
              paste0("/Volumes/NYC_geo/processing_temporary/uchmbt_raster_1p64ft_", test_ids[i], ".tif"))
  rm(uchmbt_rast)
}

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

for (i in 2:9){
  #i <- 1
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
  building_rast[lc_rast_tile_buff_reproj_resample != 5] <- NA # building only filter
  writeRaster(building_rast,
              paste0("/Volumes/NYC_geo/processing_temporary/buildingheight_raster_1p64ft_", test_ids[i], ".tif"))
  
  tree_rast <- uchmbt_rast
  cls_tile_buff <- crop(cls, las_index_vect_reproj_tile_buff)
  cls_tile_buff_reproj <- project(cls_tile_buff, crs(uchmbt_rast))
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


