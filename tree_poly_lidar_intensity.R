# Extract lidar intensity values for vegetation points within TNC tree crown polygons

library(terra)
library(sf)
library(lidR)
library(tidyverse)
library(purrr)

# Load lidar data (might need to move it back over to SSD)
# will need to apply across catalog
# likely will need to chunk this up by tree polygon id

# Load tnc polygon data gdb

# Load las file index poly

# Fix projections
# Intersect to get list of las files needed for a given subset area

# for each tree crown polygon, loop or catalog or alt parallel process
# Subset for just vegetation classified lidar points (low, medium, high)
# intensity
# To start: get mean of all vegetation (over 2 m height?)
# Alt: get mean all vegetation points, 95th percentile, SD
# retain these values for each tree id
#

# Export dataframe with these intensity values

las_tile_id <- "10225"

# test single las tile
# lidar
setwd("/Users/dlm356/dlm356_files/nyc_trees/nyc_lidar_2021/") # a couple tester files are here, others on big external disk
las_filename <- paste0(las_tile_id, ".las")
las <- readLAS(las_filename, select = "*", filter = "-drop_z_below 0")
pts_classes <- c(LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION) # LASGROUND, LASWATER
las_veg <- filter_poi(las, Classification %in% c(pts_classes))
las_lowveg <- filter_poi(las, Classification %in% c(LASLOWVEGETATION))
las_mediumveg <- filter_poi(las, Classification %in% c(LASMEDIUMVEGETATION))
las_highveg <- filter_poi(las, Classification %in% c(LASHIGHVEGETATION))

# load las index poly
las_tiles <- st_read("NYC2021_LAS_Index.shp")
las_tile <- las_tiles %>% filter(LAS_ID == las_tile_id)

# Load tree polys
tree_poly_path_full <- "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb"
tnc_gdb_polys <- st_read(tree_poly_path_full, layer = "treeobjects_2021_nyc")
tnc_gdb_polys$Object_ID <- seq(1,nrow(tnc_gdb_polys))

# subset tree polys to just the tile area for now
las_tile_reproj <- st_transform(las_tile, st_crs(tnc_gdb_polys))
tnc_gdb_polys_sub <- st_intersection(tnc_gdb_polys, las_tile_reproj)

# extract intensity
tnc_gdb_polys_sub_reproj <- st_transform(tnc_gdb_polys_sub, st_crs(las_veg))
las_in_polys <- lidR::clip_roi(las_veg, geometry = tnc_gdb_polys_sub_reproj)

#df_means <- map_dbl(las_in_polys, ~mean(.x[[1]]$Intensity, na.rm = TRUE)) # this is likely faster than loop but doesn't work yet

# this works to get mean of intensity for all vegetation
tnc_gdb_polys_sub$veg_intensity_mean <- NA
for (i in 1:length(las_in_polys)){
  tnc_gdb_polys_sub$veg_intensity_mean[i] <- mean(las_in_polys[[i]]$Intensity)
}

plot(tnc_gdb_polys_sub["veg_intensity_mean"])

ggplot(tnc_gdb_polys_sub) +
  geom_sf(aes(fill = veg_intensity_mean)) +  # Use the column you want
  scale_fill_viridis_c() +     # Optional: Improve color scale
  theme_minimal()

quantile(tnc_gdb_polys_sub$veg_intensity_mean, probs = seq(0,1,0.1), na.rm = TRUE)
