# Calculate structural variables from lidar for tree crowns

# Load in tree crown polygons
# Go polygon by polygon ?


# Set of structural variables, mostly all selected variables from Alonzo et al 2014

# Y: Max crown height, see what we get as a sanity check to the TNC data, good to have anyway
# Y: Median height of returns in crown
# Y: Crown width at median height of returns in crown
# Y: Ratio of crown height to width: median height
# Y: Average intensity below median height
# Crown surface intensity: 0.5* m spatial resolution (doing it a little coarser because we have lower points/m2)
# Surface heights (0.5* m)/surface heights (1 m)
# Y: Count of returns in 0.5 m vertical slice at 90th percentile divided by width at that height

#####

library(tidyverse)
library(sf)
library(lidR)
library(future)
library(terra) # for rasters
# Note: terra masks from lidR: area, crs, crs <-, is.empty, watershed


# Set this up to loop over sets of polygons
# Load in TNC polygons
# Load one lidar tile to do reprojection on TNC polygons
# Setup Loop, set ranges of 1000 polygons (or some similar amount, 500, 2000, something)
# Intersect with full las tile collection to know which las tiles we need
# Load in these few tiles as LAS (all together, ideally keep this to maybe 3 tiles max)
# Subset for this set of trees
# Clean out memory
# Do lidar metric processing steps
# Output data frame of tree metrics info df
# Repeat loop
# Note: figure out overhead to setup to run in parallel
# This also avoid edge effect wonkiness for lascatalog gridding with tree edges

# Possible issue: DTM tiling for subtraction
# Might want to generate the full DTM first
# Alternatively, might want to do it tile by tile to prevent big differences in DTM location
# update: created all DTM tiles, too big to do all in one file though

#plan(multisession, workers = 6) # don't use all the cores, easier on memory
#temp_output_path <- "/Volumes/NYC_geo/processing_temporary/" # do this on the external disk
# for catalog application

# Need to test it again without catalog application for center test_id and see what happens

test_ids <-  "5240" # "2242"
# test_ids <- c("2242", "5242", "7242",
#               "2240", "5240", "7240",
#               "2237", "5237", "7237")
test_files <- paste("/Volumes/DLM_backup/lidar_2021/NYC_2021/", test_ids, ".las", sep="") # will want to move this back over to SSD later

las <- readLAS(test_files, filter = "-drop_z_below 0")
#ctg <- readLAScatalog(test_files, filter = "-drop_z_below 0")
#plot(las)

pts_classes <- c(LASGROUND, LASWATER, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASBUILDING)
las2 <- filter_poi(las, Classification %in% c(pts_classes))
#plot(las2)
# 
# ctgFilterPoi <- function(chunk, pts_classes) {
#   # Load the chunk + buffer
#   las <- readLAS(chunk)
#   if (lidR::is.empty(las)) return(NULL)
#   
#   output <- filter_poi(las, Classification %in% c(pts_classes))
#   
#   # remove the buffer of the output
#   # output <- clip_roi(output, st_bbox(chunk))
#   # no need to remove buffer because buffer is 0 for this
#   
#   if (nrow(output@data) > 0){
#     return(output) # return the tiles if there are points, otherwise do nothing 
#   }
# }
# 
# opt_output_files(ctg) <- paste(temp_output_path, "filtpts_{ID}", sep = "") # temp output dir
# opt_chunk_buffer(ctg) <- 0 # no buffer for the filtering
# ctg_filt_lasfiles <- catalog_apply(ctg, ctgFilterPoi, pts_classes)


##### need to subtract off a dtm in order to get points as chm (height above ground only)
#make DSM
grid_size <- 1.64 # ft approx 0.5 # half meter chm
#dsm <- rasterize_canopy(las2, res = grid_size, pitfree(c(0, 2, 5, 10, 15))) # default, this might take longer than some other methods, overkill?

# make DTM
dtm <- rasterize_terrain(las2, grid_size, tin())


#####
# # Make DTM and DSM
# 
# # set up new run for terrain, rpf is old naming convention
# ctg_rpf <- readLAScatalog(unlist(ctg_filt_lasfiles))
# opt_output_files(ctg_rpf) <- paste(temp_output_path, "terrain_raster7_{ORIGINALFILENAME}", sep = "")
# opt_chunk_buffer(ctg_rpf) <- 328 # units in feet, approx 100 m
# 
# grid_size <- 1.64 # 0.5 # in feet instead of meters
# 
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
# #options <- list(automerge = TRUE) # this merges the output into a single raster # this option is probably the problem for memory
# # Do a SpatRaster merge() in terra instead!
# #terrain_raster <- catalog_apply(ctg_rpf, rasterizeTerrain, grid_size, .options = options) # this is a spatraster now with a pointer
# terrain_raster_files <- catalog_apply(ctg_rpf, rasterizeTerrain, grid_size)
# # Note: nearest neighbor interpolation warnings
# 
# plan(sequential) # turning this off kills the parallel process
# plan(multisession, workers = 6) # reduce number of workers to make it easier on memory
# 
# # combine raster files into a single file on the other side
# terrain_raster_sprc <- sprc(unlist(terrain_raster_files))
# terrain_raster <- merge(terrain_raster_sprc)

# to write out the 
#writeRaster(terrain_raster, paste("../sub1/bronx_terrain_test_yankeestadium_m.tif", sep=""))
# 
# # This is the DSM, I don't think we actually need this right now. Just need DTM to subtract from LAS files for canopy height estimate
# opt_output_files(ctg_rpf) <- paste(temp_output_path, "canopysurface3_{ORIGINALFILENAME}", sep="")
# 
# rasterizeCanopy <- function(chunk, grid_size) {
#   # Load the chunk + buffer
#   las <- readLAS(chunk)
#   if (lidR::is.empty(las)) return(NULL)
#   if (any(las$Classification == 2)) { # test for requirement that there be some bare ground (LASGROUND = 2)
#     # do something
#     output <- rasterize_canopy(las, res = grid_size, pitfree(c(0, 2, 5, 10, 15)))
#     
#     # remove the buffer of the output
#     #bbox <- sp::bbox(chunk)
#     #output <- lidR::remove_buffer(output, bbox)
#     output <- crop(output, st_bbox(chunk)) #crop for spatraster
#     return(output) # do we want to temporarily write out the output instead and retain as a las catalog 
#   }
# }
# 
# #options <- list(automerge = TRUE)
# #canopysurface_raster <- catalog_apply(ctg_rpf, rasterizeCanopy, grid_size, .options = options)
# canopysurface_raster_files <- catalog_apply(ctg_rpf, rasterizeCanopy, grid_size)
# 
# # reset memory to reclaim from multicore workers
# plan(sequential)
# plan(multisession, workers = 6) # reduce number of workers to make it easier on memory
# 
# canopysurface_raster_sprc <- sprc(unlist(canopysurface_raster_files))
# canopysurface_raster <- merge(canopysurface_raster_sprc) # this is a dsm




#####

# veg only las
pts_classes_veg <- c(LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION)
las3 <- filter_poi(las, Classification %in% c(pts_classes_veg))

# opt_output_files(ctg) <- paste(temp_output_path, "filtptsveg_{ID}", sep = "") # temp output dir
# opt_chunk_buffer(ctg) <- 0 # no buffer for the filtering
# ctg_filt_veg_lasfiles <- catalog_apply(ctg, ctgFilterPoi, pts_classes_veg)
# 
# ctg_filt_veg_raw <- readLAScatalog(unlist(ctg_filt_veg_lasfiles))
# 
las_veg_height <- las3 - dtm
# opt_output_files(ctg_filt_veg_raw) <- paste(temp_output_path, "filtptsvegnorm_{ID}", sep = "") # temp output dir
# ctg_filt_veg <- normalize_height(ctg_filt_veg_raw, terrain_raster) # requires a raster, not an sprc
# this is a canopy height version of the las file!

# how to create a lax file?
# this would speed up processing.

#####

# lidar tile index shapefile
tile_polys <- read_sf("/Volumes/DLM_backup/lidar_2021/NYC_2021/NYC2021_LAS_Index.shp")
tile_polys_sub <- tile_polys %>% filter(LAS_ID %in% test_ids)

# Tree polygon file
tree_poly_path_full <- "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb" # Note this is the FINAL TNC dataset
tnc_gdb_polys <- st_read(tree_poly_path_full, layer = "treeobjects_2021_nyc")
tnc_gdb_polys <- st_transform(tnc_gdb_polys, st_crs(tile_polys_sub)) # this is slow, but almost certainly faster than reprojecting the whole las dataset
tnc_gdb_polys_sub <- st_intersection(tnc_gdb_polys, tile_polys_sub)
tnc_gdb_polys_sub_sfc <- st_cast(tnc_gdb_polys_sub, "POLYGON")

# lidR user defined metrics
# https://r-lidar.github.io/lidRbook/metrics.html

# Pass on metrics from one function to the next in the implementation

# Average intensity above median height
fIntHigh <- function(z,i){
  med_height <- median(z)
  hi_inds <- which(z > med_height)
  return(mean(i[hi_inds]))
}

# Average intensity below median height
fIntLow <- function(z,i){
  med_height <- median(z)
  low_inds <- which(z < med_height)
  return(mean(i[low_inds]))
}

# Crown width at median height of returns in crown
# Implemented as get median, get all points +/- 1.64 ft (1 m bounds)
# then take max of x,y distances between these points, call this the width
fW1widthmax <- function(x,y,z){
  med_height <- median(z)
  ind_at_med <- which(z < med_height+1.64 & z > med_height-1.64)
  xy <- cbind(x[ind_at_med], y[ind_at_med])
  dist_xy <- dist(xy)
  return(max(dist_xy))
}

# Alternative version of the crown width
# create a convex hull of the points at the median height range
# calculate area, then calculate diameter of this area if it was a circle
fW1widthcircle <- function(x,y,z){
  med_height <- median(z)
  ind_at_med <- which(z < med_height+1.64 & z > med_height-1.64)
  xy <- cbind(x[ind_at_med], y[ind_at_med])
  xy_area <- st_area(st_convex_hull(st_multipoint(xy)))
  return(2*sqrt(xy_area/pi))
}

# Ratio of crown height to width at median height, max implementation
fHWrat2max <- function(x,y,z){
  med_height <- median(z)
  ind_at_med <- which(z < med_height+1.64 & z > med_height-1.64)
  xy <- cbind(x[ind_at_med], y[ind_at_med])
  dist_xy <- dist(xy)
  return(max(dist_xy))
}

# Ratio of crown height to width at median height, circle implementation
fHWrat2circle <- function(x,y,z){
  med_height <- median(z)
  ind_at_med <- which(z < med_height+1.64 & z > med_height-1.64)
  xy <- cbind(x[ind_at_med], y[ind_at_med])
  xy_area <- st_area(st_convex_hull(st_multipoint(xy)))
  med_width <- 2*sqrt(xy_area/pi)
  return(med_height/med_width)
}

# Count of returns in 0.5 m vertical slice at 90th percentile height divide by width at that height
fCP3circle <- function(x, y, z){
  q90_height <- unname(quantile(z, 0.9, na.rm = TRUE))
  if (is.na(q90_height)){
    return(NA)
  } else{
    ind_at_q90 <- which(z < q90_height+0.82 & z > q90_height-0.82) # half meter window
    slice_count <- length(ind_at_q90)
    
    xy <- cbind(x[ind_at_q90], y[ind_at_q90])
    xy_area <- st_area(st_convex_hull(st_multipoint(xy)))
    wid <- 2*sqrt(xy_area/pi)
    
    return(slice_count/wid)
  }
}

fCP3max <- function(x, y, z){
  q90_height <- unname(quantile(z, 0.9, na.rm = TRUE))
  if (is.na(q90_height)){
    return(NA)
  } else{
    ind_at_q90 <- which(z < q90_height+0.82 & z > q90_height-0.82) # half meter window
    slice_count <- length(ind_at_q90)
    
    xy <- cbind(x[ind_at_q90], y[ind_at_q90])
    dist_xy <- dist(xy)
    wid <- max(dist_xy)
    
    return(slice_count/wid)
  }
}

# Wrapper for all structural metric functions
calcMetrics <- function(x, y, z, i){
  list(height_max = max(z),
       height_median = median(z),
       intensity_mean = mean(i),
       intensity_mean_above_medianh = fIntHigh(z,i),
       intensity_mean_below_medianh = fIntLow(z,i),
       w_1_max = fW1widthmax(x,y,z),
       w_1_circle = fW1widthcircle(x,y,z),
       hw_rat_2max = fHWrat2max(x,y,z),
       hw_rat_2circle = fHWrat2circle(x,y,z),
       cp_3_max = fCP3max(x,y,z),
       cp_3_circle = fCP3circle(x,y,z))
}

calcMetrics2 <- function(las){
  x <- las$X
  y <- las$Y
  z <- las$Z
  i <- las$Intensity
  metrics <- list(height_max = max(z),
                  height_median = median(z),
                  intensity_mean = mean(i),
                  intensity_mean_above_medianh = fIntHigh(z,i),
                  intensity_mean_below_medianh = fIntLow(z,i),
                  w_1_max = fW1widthmax(x,y,z),
                  w_1_circle = fW1widthcircle(x,y,z),
                  hw_rat_2max = fHWrat2max(x,y,z),
                  hw_rat_2circle = fHWrat2circle(x,y,z),
                  cp_3_max = fCP3max(x,y,z),
                  cp_3_circle = fCP3circle(x,y,z))
  return(metrics)
}

# maybe flip this around so that the points are labeled for each tree polygon, then extract based on these unique labels??

las_veg_height_mask <- classify_poi(las_veg_height, as.integer(1), roi = tnc_gdb_polys_sub_sfc, inverse_roi = TRUE) # label anything 
las_veg_height_inpoly <- filter_poi(las_veg_height_mask, Classification %in% pts_classes_veg) # this removes all the points that are not within polygons
tnc_gdb_polys_sub_sfc$Poly_ID <- as.integer(rownames(tnc_gdb_polys_sub_sfc)) # append the unique ID for each Poly, make sure this the correct form of the indexing!
las_veg_height_inpoly_labeled <- merge_spatial(las_veg_height_inpoly, tnc_gdb_polys_sub_sfc, attribute = "Poly_ID")
# now we have each point retained labeled as part of a polygon
# Will need to buffer or figure out borders somehow...

#tnc_gdb_polys_sub_metrics <- plot_metrics(las_veg_height, calcMetrics(X, Y, Z, Intensity), tnc_gdb_polys_sub_sfc)
# this works fine (with acceptable warnings), so some issue with setting it up for the catalog processing?

tnc_metrics_test <- crown_metrics(las_veg_height_inpoly_labeled, ~calcMetrics(X, Y, Z, Intensity), attribute = "Poly_ID")
# this is much faster too

# if can do the merge_spatial for the entire data set (or for parts of it, chunked perhaps), could likely run this as a catalog process for entire study area 



# check function list and simplify to reduce repeated calculations
# setup to work with chunked buffers, test to see how it's extracting
# if it duplicates, also extract n points and keep the very of the poly extract with most points
ctgCalcMetrics <- function(las, polygons) {
  las <- readLAS(las)
  if (lidR::is.empty(las)) return(NULL)
  output <- plot_metrics(las, calcMetrics(X, Y, Z, Intensity), polygons)
  return(output)
}

# Setup new catalog
setwd("/Volumes/NYC_geo/processing_temporary/")
normlasfiles <- list.files(pattern = glob2rx("*filtptsvegnorm*.las"))
ctg_filt_veg2 <- readLAScatalog(normlasfiles, filter = "-drop_z_below 0")

# appends onto existing sfc_POLYGONS file
las_ex <- readLAS(normlasfiles[9])
#tnc_gdb_polys_sub_metrics <- plot_metrics(las_veg_height, calcMetrics(X, Y, Z, Intensity), tnc_gdb_polys_sub_sfc)
test_metrics_las_ex <- plot_metrics(las_ex, calcMetrics(X, Y, Z, Intensity), tnc_gdb_polys_sub_sfc)
# it may not be an issue with catalog_apply() step, might be something in the plot_metrics behavior

opt_output_files(ctg_filt_veg) <- paste(temp_output_path, "treemetrics_{ID}", sep = "") # temp output dir
opt_chunk_buffer(ctg_filt_veg) <- 328 # units in feet, approx 100 m
tnc_gdb_polys_sub_metrics <- catalog_apply(ctg_filt_veg, ctgCalcMetrics, tnc_gdb_polys_sub_sfc)

test_output <- polygon_metrics(las = ctg_filt_veg, func = ~calcMetrics2(.x), geometry = tnc_gdb_polys_sub_sfc)
# wrap this as a catalog_apply() function, crashes when run directly

ggplot(tnc_gdb_polys_sub_metrics) +
  geom_point(aes(x = median, y = max)) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dotted") +
  coord_equal()

ggplot(tnc_gdb_polys_sub_metrics) +
  geom_histogram(aes(x = log10(Height - max)))
                         