

library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(lubridate)
`%notin%` <- Negate(`%in%`)
# Ref targets

#

# Prep input stuff
num_bands <- 8 # 4 or 8
path_planet_mosaics <- "/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly"

# to write out
output_prefix <- "all_nyc_8b_"
output_target_refl_dir <- "/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal" # no trailing "/"
#...

# Apply empirical line corrections
# setup

setwd(output_target_refl_dir)
output_empline_f2 <- read.csv(paste0(output_prefix, "_nyc_daily_stack_highsunonly_calibrated_to_nyc_ground_ref_targets_revised_empline.csv")) # this will change between runs
output_empline_f2$date <- ymd(output_empline_f2$date)

# Loop through output_empline_f2 for each image
# apply slope and intercept for each band
# output each calibrated image (bronxcal20220710ref)

img_list_f2 <- unique(output_empline_f2$img_file)
cal_name <- "nyccal20220710ref"
setwd(path_planet_mosaics)
#band_list <- c("blue", "green", "red", "nir")

# subset image list to start new run (if needed)
setwd(output_target_refl_dir)
outputted_imgs <- list.files(pattern = glob2rx("*.tif"))
img_list_prefix <- sapply(strsplit(img_list_f2, "\\."), `[`, 1)
imgs_to_output <- paste0(img_list_prefix, "_", cal_name, ".tif")
img_list_f2_sub <- img_list_f2[which(imgs_to_output %notin% outputted_imgs)]
img_list_f2 <- img_list_f2_sub
#output_empline_f2_2 <- output_empline_f2 %>% filter(img_file %in% img_list_f2)

setwd(path_planet_mosaics)
for (i in 1:length(img_list_f2)) { # test range
#for (i in 1:10) { # test range
  print(i)
  img_prefix <- unlist(strsplit(img_list_f2[i], "\\."))[1]
  
  print(img_prefix)
  img <- rast(img_list_f2[i])
  output_empline_f2_sub <- output_empline_f2[which(output_empline_f2$img_file == img_list_f2[i]),]
  
  # calcCor <- function(img, output_empline_f2_sub, b){
  #   b_ind <- which(output_empline_f2_sub$band == b)
  #   slope <- output_empline_f2_sub$slope[b_ind]
  #   intercept <- output_empline_f2_sub$intercept[b_ind]
  #   new_img_layer <- round(img[[b]]*slope + intercept, 0)
  #   return(new_img_layer)
  # }
  
  calcCor <- function(img, output_empline_f2_sub, b){
    b_ind <- which(output_empline_f2_sub$band == b)
    slope <- output_empline_f2_sub$slope[b_ind]
    intercept <- output_empline_f2_sub$intercept[b_ind]
    new_img_layer <- round(img[[b]]*slope + intercept, 0)
    return(new_img_layer)
  }
  
  # 4 band
  img2 <- img
  img2[["blue"]] <- calcCor(img, output_empline_f2_sub, "blue")
  img2[["green"]] <- calcCor(img, output_empline_f2_sub, "green")
  img2[["red"]] <- calcCor(img, output_empline_f2_sub, "red")
  img2[["nir"]] <- calcCor(img, output_empline_f2_sub, "nir")
  
  # 8 band includes 4 band plus 4 more
  if (num_bands == 8){
    img2[["coastal_blue"]] <- calcCor(img, output_empline_f2_sub, "coastal_blue")
    img2[["green_i"]] <- calcCor(img, output_empline_f2_sub, "green_i")
    img2[["yellow"]] <- calcCor(img, output_empline_f2_sub, "yellow")
    img2[["rededge"]] <- calcCor(img, output_empline_f2_sub, "rededge")
  }
  
  #output_dir <- "/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_highsunonly_cal/"
  #output_dir <- "/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_v2_test/" # new output dir for tests
  output_dir <- output_target_refl_dir
  output_filename <- paste0(output_dir, "/", img_prefix, "_", cal_name, ".tif")
  
  img2[img2 > 10000] <- 10000 # make sure reflectance is capped at 10000 (1 x 10000) after calibrating
  img2[img2 < 0] <- 0 # make sure no reflectance is below zero after calibrating
  
  writeRaster(img2, filename = output_filename, datatype = "INT2S")
  
  print(paste0("Written out raster: ", output_filename))
  
  # Cleanup
  rm(img)
  rm(img2)
  gc()
  tmpFiles(remove=TRUE, orphan = TRUE) # remove temporary files
}
