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


tree_ids <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures/tree_id_list_95pct_2week_intersect.csv")
#date_ranges_4b <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures/tree_id_date_ranges_95pct_2week_4b.csv")
date_ranges_8b <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures/tree_id_date_ranges_95pct_2week_8b.csv")

date_ranges <- date_ranges_8b

#tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

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
#tree_ids <- unique(tree_pheno$Poly_ID)

tree_ids <- tree_ids$x

stepsize <- 2000
top <- floor(length(tree_ids)/stepsize)*stepsize+1
step_ranges <- seq(1, top, stepsize)

#for (idset in 1:length(step_ranges)){
for (idset in 11:length(step_ranges)){
  start_time <- Sys.time()
  range_min <- step_ranges[idset]
  idset_range <- seq(range_min, range_min + stepsize - 1)
  tree_ids_sub <- tree_ids[idset_range]
  print(idset)
  
  setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point") # 8 b setup
  tree_id_spectra <- getTreePlanetSpectra2(tree_ids_sub)
  
  # take mean for each month (across all years)
  tree_id_spectra$date <- ymd(tree_id_spectra$date)
  
  # drop tree ids dates not in the date range, and define which row in the date_ranges 
  tree_id_spectra$in_date_range <- NA
  for (i in 1:nrow(date_ranges)){
    # get rows within each date range, and assign a number based on the row from the date_ranges file df (as a reference)
    tree_id_spectra$in_date_range[which(ymd(tree_id_spectra$date) %within% interval(date_ranges$start_date[i], date_ranges$end_date[i]))] <- i
  }
  
  tree_id_spectra <- na.omit(tree_id_spectra) # drop the non-included dates
  tree_id_spectra <- as.data.frame(tree_id_spectra)
  # Do 2 week range averaging, then do band combinations, then do Z-scaling
  band_list <- c("coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge", "nir")
  tree_id_spectra_date_range_mean <- aggregate(tree_id_spectra[,band_list], by = list(tree_id_spectra$in_date_range, tree_id_spectra$Object_ID), FUN = "mean", na.rm = TRUE)
  colnames(tree_id_spectra_date_range_mean)[1:2] <- c("Date_Range_Group", "Object_ID")
  # would need to add the date back into the date_range_group
  
  #####
  # OMNBR step (all band combinations, not true omnbr but would be a setup step)
  
  #tree_id_spectra <- as.data.frame(tree_id_spectra)
  
  #band_list <- c("coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge", "nir")
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
      #nd_index <- calcNormDif(tree_id_spectra[,b2], tree_id_spectra[,b1]) # note this only works on data frames, not tibbles
      nd_index <- calcNormDif(tree_id_spectra_date_range_mean[,b2], tree_id_spectra_date_range_mean[,b1]) # note this only works on data frames, not tibbles
      
      if (i == 1 & j == 2){
        omnbr <- nd_index
      } else {
        omnbr <- cbind.data.frame(omnbr, nd_index)
        colnames(omnbr) <- colname_list
      }
      
    }
  }
  
  #df_spectra_omnbr <- bind_cols(tree_df, omnbr)
  df_spectra_omnbr <- bind_cols(tree_id_spectra_date_range_mean, omnbr)
  
  subDateChar <- function(i){
    return(gsub("-", "", as.character(date_ranges$start_date[i])))
  }
  # replace numbers with actual dates (as character)
  df_spectra_omnbr$Date_Range_Group <- sapply(df_spectra_omnbr$Date_Range_Group, FUN = subDateChar)
  
  
  # df_spectra_omnbr_zscaled %>% filter(Object_ID == 547) %>%
  #   ggplot() +
  #   geom_histogram(aes(nd_nir_red_zscaled))
  
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
  
  # make long, append date_range_group and band type, make wide again
  df_spectra_omnbr_zscaled_long <- df_spectra_omnbr_zscaled %>% pivot_longer(cols = colnames(df_spectra_omnbr_zscaled)[3:ncol(df_spectra_omnbr_zscaled)])
  df_spectra_omnbr_zscaled_long <- df_spectra_omnbr_zscaled_long %>% mutate(name_drg = paste0(name, "_mean_", Date_Range_Group))
  df_spectra_omnbr_zscaled_wide <- df_spectra_omnbr_zscaled_long %>% pivot_wider(id_cols = c(Object_ID), names_from = name_drg, values_from = value)
  
  df_spectra_omnbr_long <- df_spectra_omnbr %>% pivot_longer(cols = colnames(df_spectra_omnbr)[3:ncol(df_spectra_omnbr)])
  df_spectra_omnbr_long <- df_spectra_omnbr_long %>% mutate(name_month = paste0(name, "_mean_", Date_Range_Group))
  df_spectra_omnbr_wide <- df_spectra_omnbr_long %>% pivot_wider(id_cols = c(Object_ID), names_from = name_month, values_from = value)
  
  # UPDATE THESE OUTPUTS
  outpath <- "/Volumes/NYC_geo/tree_classification/extracted_8band_2week"
  write.csv(df_spectra_omnbr_zscaled_wide,
            paste0(outpath, "/zscaled/tree_points_elegant_8b_spectra_zscaled_2week_idset", idset, ".csv"),
            row.names = FALSE)
  write.csv(df_spectra_omnbr_wide,
            paste0(outpath, "/raw/tree_points_elegant_8b_spectra_raw_2week_idset", idset, ".csv"),
            row.names = FALSE)
  print(Sys.time() - start_time)
}
