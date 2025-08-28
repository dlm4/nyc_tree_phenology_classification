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

test_ids <- "5240"
# test_ids <- c("2242", "5242", "7242",
#               "2240", "5240", "7240",
#               "2237", "5237", "7237")
test_files <- paste("/Volumes/DLM_backup/lidar_2021/NYC_2021/", test_ids, ".las", sep="") # will want to move this back over to SSD later

las <- readLAS(test_files, filter = "-drop_z_below 0")
plot(las)

pts_classes <- c(LASGROUND, LASWATER, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASBUILDING)
las2 <- filter_poi(las, Classification %in% c(pts_classes))
plot(las2)

##### need to subtract off a dtm in order to get points as chm (height above ground only)
# make DSM
grid_size <- 1.64 # ft approx 0.5 # half meter chm
dsm <- rasterize_canopy(las2, res = grid_size, pitfree(c(0, 2, 5, 10, 15))) # default, this might take longer than some other methods, overkill?

# make DTM
dtm <- rasterize_terrain(las2, grid_size, tin())


# 
# # Prep stuff to make cleaner canopy height model (CHM)
# 
# # Mask out buildings
# mask_base <- dsm*0
# pts_classes <- c(LASBUILDING)
# las_buildings <- filter_poi(las_reproj_z2m, Classification %in% c(pts_classes)) # still needed las_reproj_z2m...
# building_density <- rasterize_density(las_buildings, grid_size) %>% extend(mask_base, fill = 0) %>% crop(mask_base) # extend and then crop to force exact match with the extent of the ground
# mask_building <- mask(mask_base, building_density, inverse = TRUE, maskvalues = 0, updatevalue = 1)
# mask_building_smooth <- focal(x = mask_building, w = 7, fun = "modal") # use 'focal' fn in terra for smoothing, 11 for 0.25, use 7 for 0.5
# mask_building_smooth[is.na(mask_building_smooth)] <- 0 # remove NAs from the mask
# 
# dsm_masked <- mask(dsm, mask_building_smooth, maskvalues = 1, updatevalue = 0)
# dtm_masked <- mask(dtm, mask_building_smooth, inverse = TRUE, maskvalues = 1, updatevalue = 0)
# 
# #new_dsm <- dsm_masked + dtm_masked
# 
# # Inverted density mask for tree CHM to remove everything else
# mask_base <- dsm*0
# pts_classes <- c(LASGROUND, LASWATER, LASHIGHVEGETATION, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASBUILDING)
# las_fullcls <- filter_poi(las_reproj_z2m, Classification %in% c(pts_classes)) # still needed las_reproj_z2m...
# fullcls_density <- rasterize_density(las_fullcls, grid_size) %>% extend(mask_base, fill = 0) %>% crop(mask_base) # extend and then crop to force exact match with the extent of the ground
# mask_fullcls <- mask(mask_base, fullcls_density, inverse = TRUE, maskvalues = 0, updatevalue = 1)
# mask_fullcls_smooth <- focal(mask_fullcls, w = 7, fun = "modal", na.rm = T) # need to remove NAs otherwise get square holes
# mask_fullcls_smooth[is.na(mask_fullcls_smooth)] <- 0 # remove NAs from the mask
# dsm_masked <- mask(dsm_masked, mask_fullcls_smooth, maskvalues = 0, updatevalue = 0) # invert mask
# dtm_masked2 <- mask(dtm, mask_fullcls_smooth, maskvalues = 1, updatevalue = 0) # this is the fullcls version of dtm, mask again with buildings
# dtm_masked3 <- mask(dtm_masked2, mask_building_smooth, maskvalues = 1, updatevalue = 0)
# dtm_masked_sum <- dtm_masked + dtm_masked3
# new_dsm <- dsm_masked + dtm_masked_sum
# 
# # Make CHM
# chm <- new_dsm - dtm # try with new one; 0.25 m might be too small for this point density, might need to do 0.5 m do be safe
# chm[chm < 0] <- 0 # remove negative values

#####

# veg only las
pts_classes <- c(LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION)
las3 <- filter_poi(las, Classification %in% c(pts_classes))

las_veg_height <- las3 - dtm
# this is a canopy height version of the las file

# lidar tile index shapefile
tile_polys <- read_sf("/Volumes/DLM_backup/lidar_2021/NYC_2021/NYC2021_LAS_Index.shp")
tile_polys_sub <- tile_polys %>% filter(LAS_ID %in% test_ids)

# Tree polygon file
tree_poly_path_full <- "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb" # Note this is the FINAL TNC dataset
tnc_gdb_polys <- st_read(tree_poly_path_full, layer = "treeobjects_2021_nyc")
tnc_gdb_polys <- st_transform(tnc_gdb_polys, st_crs(tile_polys_sub)) # this is slow, so maybe should do it the other way, hmmmm
tnc_gdb_polys_sub <- st_intersection(tnc_gdb_polys, tile_polys_sub)
tnc_gdb_polys_sub_sfc <- st_cast(tnc_gdb_polys_sub, "POLYGON")

# lidR user defined metrics
# https://r-lidar.github.io/lidRbook/metrics.html

# Pass on metrics from one function to the next in the implementation

# Average intensity below median height
fInt2 <- function(z,i){
  med_height <- median(z)
  low_inds <- which(z < med_height)
  return(mean(z[low_inds]))
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

# Ratio of crown to width at median height
fHWrat2 <- function(x,y,z){
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

# Need to figure out the surface metrics next... not sure how to do those yet, voxels??
# Crown surface intensity

# Get first return vegetation points
# Voxelize_points at resolution
# get the average intensity

# Wrapper for all structural metric functions
calcMetrics = function(x, y, z, i){
  list(height_max = max(z),
       height_median = median(z),
       intensity_mean = mean(i),
       int_2 = fInt2(z,i),
       w_1_max = fW1widthmax(x,y,z),
       w_1_circle = fW1widthcircle(x,y,z),
       hw_rat_2 = fHWrat2(x,y,z),
       cp_3_circle = fCP3circle(x,y,z))
}

calcMetricsSet2 = function(x, y, z, i){
  list(crown_surf_intensity = mean(i),
       crown_surf_height05 = mean(z))
}

calcMetricsSet3 = function(z){
  list(crown_surf_height1 = mean(z))
}

# appends onto existing sfc_POLYGONS file
tnc_gdb_polys_sub_metrics <- plot_metrics(las_veg_height, calcMetrics(X, Y, Z, Intensity), tnc_gdb_polys_sub_sfc)
las_veg_height_surf_vox05 <- voxelize_points(las_veg_height, 1.64) # don't need first return filter because we're doing max anyway
las_veg_height_surf_vox1 <- voxelize_points(las_veg_height, 3.28/4)

#tnc_gdb_polys_sub_metrics <- plot_metrics(las_veg_height_surf_vox05, calcMetricsSet2(X, Y, Z, Intensity), tnc_gdb_polys_sub_metrics)
#tnc_gdb_polys_sub_metrics <- plot_metrics(las_veg_height_surf_vox1, calcMetricsSet3(X, Y, Z, Intensity), tnc_gdb_polys_sub_metrics)

tnc_gdb_polys_sub_metrics <- plot_metrics(las_veg_height_surf_vox05, calcMetricsSet2(X, Y, Z, Intensity), tnc_gdb_polys_sub_metrics)

# this issue with the 1 m resolution is not an empty las issue, this is handled fine at 0.5 m resolution. Really not sure why this is failing then
m05 <- clip_roi(las_veg_height_surf_vox05, tnc_gdb_polys_sub_sfc)
m1 <- clip_roi(las_veg_height_surf_vox1, tnc_gdb_polys_sub_sfc)
#good_poly_inds <- intersect(which(unlist(sapply(m1, nrow)) > 0), which(sapply(m1, class) == "LAS")) # get joint grouping of both
#tnc_sub <- tnc_gdb_polys_sub_metrics[good_poly_inds,]
test1 <- plot_metrics(las_veg_height_surf_vox1, calcMetricsSet3(Z), tnc_sub) # this still does not work even with filtering, not sure why voxelization is not working at this resolution.

ggplot(tnc_gdb_polys_sub_metrics) +
  geom_point(aes(x = median, y = max)) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dotted") +
  coord_equal()

ggplot(tnc_gdb_polys_sub_metrics) +
  geom_histogram(aes(x = log10(Height - max)))
                         