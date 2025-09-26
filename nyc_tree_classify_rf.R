# this is actually compositing script!

library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(ranger)
library(future)
library(future.apply)

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
}

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

#####
# Prep stuff

#####
# Get imagery for all 8 bands
# Can make this more efficient to only search one Objset for a single id version - could check range before running instead of a loop?
getTreePlanetSpectra <- function(ids_all){
  for (i in 1:7){
    print(paste0("Objset ", i)) # looping over each of the Objset groupings
    objset_text <- paste0("*objset", i, "_*.csv")
    file_list_sub <- list.files(pattern = glob2rx(objset_text))
    names(file_list_sub) <- strsplit(file_list_sub, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 27) %>% ymd() # name the list
    tree_df_sub <- purrr::map_df(file_list_sub, fread, .id = 'date') 
    tree_df_sub2 <- tree_df_sub %>% filter(Object_ID %in% ids_all)
    if (i == 1){
      tree_df_agg <- tree_df_sub2
    } else {
      tree_df_agg <- bind_rows(tree_df_agg, tree_df_sub2)
    }
    rm(tree_df_sub)
    rm(tree_df_sub2)
  }
  gc()
  
  return(tree_df_agg)
}


getTreePlanetSpectra2 <- function(ids_all){
  
  objset_inds <- c()
  if (any(ids_all %in% 1:300000)){objset_inds <- c(1)}
  if (any(ids_all %in% 300001:600000)){objset_inds <- c(objset_inds, 2)}
  if (any(ids_all %in% 600001:900000)){objset_inds <- c(objset_inds, 3)}
  if (any(ids_all %in% 900001:1200000)){objset_inds <- c(objset_inds, 4)}
  if (any(ids_all %in% 1200001:1500000)){objset_inds <- c(objset_inds, 5)}
  if (any(ids_all %in% 1500001:1800000)){objset_inds <- c(objset_inds, 6)}
  if (any(ids_all > 1800000)){objset_inds <- c(objset_inds, 7)}
  
  check_inds <- 1
  
  for (i in objset_inds){
    print(paste0("Objset ", i)) # looping over each of the Objset groupings
    objset_text <- paste0("*objset", i, "_*.csv")
    file_list_sub <- list.files(pattern = glob2rx(objset_text))
    names(file_list_sub) <- strsplit(file_list_sub, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 27) %>% ymd() # name the list
    tree_df_sub <- purrr::map_df(file_list_sub, fread, .id = 'date') 
    tree_df_sub2 <- tree_df_sub %>% filter(Object_ID %in% ids_all)
    if (check_inds == 1){
      tree_df_agg <- tree_df_sub2
    } else {
      tree_df_agg <- bind_rows(tree_df_agg, tree_df_sub2)
    }
    rm(tree_df_sub)
    rm(tree_df_sub2)
    check_inds <- 2
  }
  gc()
  
  return(tree_df_agg)
}


# get full list of tree IDs
tree_ids <- unique(tree_pheno$Poly_ID)

stepsize <- 5000
top <- floor(length(tree_ids)/stepsize)*5000+1
step_ranges <- seq(1, top, stepsize)

# thiw would be parallel setup, but doing as old fashioned loop instead to prevent memory crash for now
#calcMonthMeans8b <- function(idset, tree_ids, stepsize, step_ranges){
for (idset in 1:length(step_ranges)){

  start_time <- Sys.time()
  range_min <- step_ranges[idset]
  idset_range <- seq(range_min, range_min + stepsize - 1)
  tree_ids_sub <- tree_ids[idset_range]
  print(idset)
  # Might need to set this up to do smaller groups of 1000 trees or something because it crashes with too many
  
  setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
  tree_id_spectra <- getTreePlanetSpectra2(tree_ids_sub)
  
  # take mean for each month (across all years)
  tree_id_spectra$date <- ymd(tree_id_spectra$date)
  
  #####
  # OMNBR step
  
  tree_id_spectra <- as.data.frame(tree_id_spectra)
  
  band_list <- c("coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge", "nir")
  colname_list <- c()
  
  for (i in 1:(length(band_list)-1)){
    b1 <- band_list[i]
    for (j in (i+1):length(band_list)){
      b2 <- band_list[j]
      
      # Do longer wavelength first since we might expect longer to be greater in NDVI form
      band_name <- paste("nd", b2, b1, sep = "_")
      #print(band_name)
      
      colname_list <- c(colname_list, band_name)
      
      #nd_index <- calcNormDif(tree_df[,b2], tree_df[,b1]) # note this only works on data frames, not tibbles
      nd_index <- calcNormDif(tree_id_spectra[,b2], tree_id_spectra[,b1]) # note this only works on data frames, not tibbles
      
      if (i == 1 & j == 2){
        omnbr <- nd_index
      } else {
        omnbr <- cbind.data.frame(omnbr, nd_index)
        colnames(omnbr) <- colname_list
      }
      
    }
  }
  
  #df_spectra_omnbr <- bind_cols(tree_df, omnbr)
  df_spectra_omnbr <- bind_cols(tree_id_spectra, omnbr)
  
  #####
  # Need to do z-scaling here
  
  tree_id_list <- unique(df_spectra_omnbr$Object_ID)
  
  df_spectra_omnbr_zscaled <- as.data.frame(df_spectra_omnbr)
  colnames(df_spectra_omnbr_zscaled)[3:ncol(df_spectra_omnbr)] <- paste(colnames(df_spectra_omnbr)[3:ncol(df_spectra_omnbr)], "_zscaled", sep = "")
  
  for (i in 1:length(tree_id_list)){
    #print(i) 
    tree_inds <- which(df_spectra_omnbr$Object_ID == tree_id_list[i])
    nd_ts_zscaled <- scale(df_spectra_omnbr[tree_inds, 3:ncol(df_spectra_omnbr)])
    df_spectra_omnbr_zscaled[tree_inds, 3:ncol(df_spectra_omnbr_zscaled)] <- nd_ts_zscaled # can't do this with too many rows, doesn't work
  }
  # Monthly averaging of z-scaled here
  # Only include 2020-2024 for consistency (only have a couple months of 2025)
  
  yr <- year(df_spectra_omnbr_zscaled$date)
  df_spectra_omnbr_zscaled <- df_spectra_omnbr_zscaled[which(yr %in% 2020:2024),] # filtering to keep only 2020-2024
  df_spectra_omnbr <- df_spectra_omnbr[which(yr %in% 2020:2024),] # filtering to keep only 2020-2024
  
  # aggregate mean by month
  mo <- month(df_spectra_omnbr_zscaled$date)
  df_spectra_omnbr_zscaled_monthmean <- aggregate(df_spectra_omnbr_zscaled[,3:ncol(df_spectra_omnbr_zscaled)], by = list(mo, df_spectra_omnbr_zscaled$Object_ID), FUN = "mean", na.rm = TRUE)
  colnames(df_spectra_omnbr_zscaled_monthmean)[1:2] <- c("Month", "Object_ID")
  
  mo <- month(df_spectra_omnbr$date)
  df_spectra_omnbr_monthmean <- aggregate(df_spectra_omnbr[,3:ncol(df_spectra_omnbr)], by = list(mo, df_spectra_omnbr$Object_ID), FUN = "mean", na.rm = TRUE)
  colnames(df_spectra_omnbr_monthmean)[1:2] <- c("Month", "Object_ID")
  
  # make long, append month, make wide again
  df_spectra_omnbr_zscaled_monthmean_long <- df_spectra_omnbr_zscaled_monthmean %>% pivot_longer(cols = colnames(df_spectra_omnbr_zscaled_monthmean)[3:ncol(df_spectra_omnbr_zscaled_monthmean)])
  df_spectra_omnbr_zscaled_monthmean_long <- df_spectra_omnbr_zscaled_monthmean_long %>% mutate(name_month = paste0(name, "_mean_", Month))
  df_spectra_omnbr_zscaled_monthmean_wide <- df_spectra_omnbr_zscaled_monthmean_long %>% pivot_wider(id_cols = c(Object_ID), names_from = name_month, values_from = value)
  
  df_spectra_omnbr_monthmean_long <- df_spectra_omnbr_monthmean %>% pivot_longer(cols = colnames(df_spectra_omnbr_monthmean)[3:ncol(df_spectra_omnbr_monthmean)])
  df_spectra_omnbr_monthmean_long <- df_spectra_omnbr_monthmean_long %>% mutate(name_month = paste0(name, "_mean_", Month))
  df_spectra_omnbr_monthmean_wide <- df_spectra_omnbr_monthmean_long %>% pivot_wider(id_cols = c(Object_ID), names_from = name_month, values_from = value)
  
  outpath <- "/Volumes/NYC_geo/tree_classification/extracted_8band_traintest"
  write.csv(df_spectra_omnbr_zscaled_monthmean_wide,
            paste0(outpath, "/zscaled/tree_points_8b_spectra_zscaled_monthmean_idset", idset, ".csv"),
            row.names = FALSE)
  write.csv(df_spectra_omnbr_monthmean_wide,
            paste0(outpath, "/raw/tree_points_8b_spectra_raw_monthmean_idset", idset, ".csv"),
            row.names = FALSE)
  print(Sys.time() - start_time)
}

# test
#idset <- 4
#calcMonthMeans8b(idset, tree_ids, stepsize, step_ranges)

# maybe don't parallelize this because I don't want to crash the laptop again
# 
# plan(multisession, workers = 6) # start parallel again each loop
# idset_range <- 1:length(step_ranges)
# future_lapply(idset_range, FUN = calcMonthMeans8b, tree_ids, stepsize, step_ranges)
# plan(sequential)
