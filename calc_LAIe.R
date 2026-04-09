library(tidyverse)
library(lidR)
library(terra)
library(sf)

calcLAIe <- function(beta, n_ground, n_total){
  lai_e <- -beta * log(n_ground/n_total)
  return(lai_e)
}

pts_range <- 0:1000
#plot(pts_range/max(pts_range)*100, calcLAIe(2, pts_range, max(pts_range)))

ggplot() +
  geom_point(aes(x = pts_range/max(pts_range)*100, y = calcLAIe(2, pts_range, max(pts_range)))) +
  labs(x = "Ground Return %", y = "LAIe") +
  theme_bw()

#####

# apply LAIe calculation to lidar tiles

# Ground class
# Vegetation classes > 1 m in height
# 987202.las

#tile_ids <- 987202
tile_ids <- 990190

# Load las tile vector
tile_polys <- read_sf("/Volumes/DLM_backup/lidar_2021/NYC_2021/NYC2021_LAS_Index.shp")
tile_polys <- tile_polys[which(tile_polys$LAS_ID %in% tile_ids),] # Do ID intersect to only work with tiles where a DTM raster could be created

dtm <- rast(paste0("/Volumes/NYC_geo/nyc_lidar_metrics/terrain_rasters/terrain_raster_1p64ft_", tile_ids,".tif"))

pts_classes <- c(LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION) # load only veg classes!
filter_string <- paste("-drop_z_below 0 -keep_class", paste(as.character(sort(pts_classes)), collapse = " "))
las_veg <- readLAS(paste0("/Volumes/DLM_backup/lidar_2021/NYC_2021/", tile_ids, ".las"), filter_string)

las_veg_height <- las_veg - dtm
rm(las_veg)

#pts_classes <- c(LASGROUND) # load only ground class (might need to include water at some point later too)
#filter_string <- paste("-drop_z_below 0 -keep_class LASGROUND")
las_ground <- readLAS(paste0("/Volumes/DLM_backup/lidar_2021/NYC_2021/", tile_ids, ".las"), "-drop_z_below 0 -keep_class LASGROUND")


# Load TNC polygons
tree_poly_path_full <- "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb" # Note this is the FINAL TNC dataset
tnc_gdb_polys <- st_read(tree_poly_path_full, layer = "treeobjects_2021_nyc")
tnc_gdb_polys$Poly_ID <- 1:nrow(tnc_gdb_polys)
# reproject
tnc_gdb_polys_reproj <- st_transform(tnc_gdb_polys, st_crs(dtm))
tnc_gdb_polys_reproj_centroid <- st_centroid(tnc_gdb_polys_reproj)
# Get polygon subset from relevant centroids
tnc_gdb_polys_reproj_centroid_sub <- st_intersection(tnc_gdb_polys_reproj_centroid, tile_polys)
tnc_gdb_polys_reproj_sub <- tnc_gdb_polys_reproj[which(tnc_gdb_polys_reproj$Poly_ID %in% tnc_gdb_polys_reproj_centroid_sub$Poly_ID),]


calcNpts <- function(z){
  list(npts = length(z))
}

# Prep LAS files
las_veg_height_mask <- classify_poi(las_veg_height, as.integer(1), roi = tnc_gdb_polys_reproj_sub, inverse_roi = TRUE)
las_veg_height_inpoly <- filter_poi(las_veg_height_mask, Classification %in% pts_classes, Z > (1/0.3048)) # this removes all the points that are not within polygons



# if this is empty, then no tree points within the remaining las, and this won't work
#if (!lidR::is.empty(las_veg_height_inpoly)){
las_veg_height_inpoly_labeled <- merge_spatial(las_veg_height_inpoly, tnc_gdb_polys_reproj_sub, attribute = "Poly_ID")
  
# Calculate crown metrics and strip geometry
tnc_metrics_sub <- crown_metrics(las_veg_height_inpoly_labeled, ~calcNpts(Z), attribute = "Poly_ID") %>% st_drop_geometry()
colnames(tnc_metrics_sub)[2] <- "npts_veg"

# Prep LAS files
las_ground_mask <- classify_poi(las_ground, as.integer(1), roi = tnc_gdb_polys_reproj_sub, inverse_roi = TRUE)
las_ground_inpoly <- filter_poi(las_ground_mask, Classification == 2) # this removes all the points that are not within polygons
las_ground_inpoly_labeled <- merge_spatial(las_ground_inpoly, tnc_gdb_polys_reproj_sub, attribute = "Poly_ID")

# Calculate crown metrics and strip geometry
tnc_metrics_sub_ground <- crown_metrics(las_ground_inpoly_labeled, ~calcNpts(Z), attribute = "Poly_ID") %>% st_drop_geometry()
colnames(tnc_metrics_sub_ground)[2] <- "npts_ground"

tnc_merged <- merge(tnc_metrics_sub, tnc_metrics_sub_ground, by = "Poly_ID")

tnc_merged <- tnc_merged %>% mutate(lai_e <- calcLAIe(2, npts_ground, (npts_ground + npts_veg)))
colnames(tnc_merged)[4] <- "LAIe"

tnc_merged_sub_poly <- merge(tnc_gdb_polys_reproj_sub, tnc_merged, by = "Poly_ID")

tnc_merged_sub_poly$LAIe <- round(tnc_merged_sub_poly$LAIe, 2)

plot(tnc_merged_sub_poly["LAIe"])

plot(tnc_merged_sub_poly$LAIe, tnc_merged_sub_poly$Height)

setwd("/Volumes/NYC_geo/lai/laie_test")
st_write(tnc_merged_sub_poly, paste0("laie_test_", tile_ids, ".geojson"))
