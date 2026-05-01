# tree point variable set from:
# tree_counts_in_planetscope3_updated_mset.R

library(tidyverse)
library(lubridate)
library(purrr)
library(data.table)
library(arrow)

agg_fun <- "median" # MEDIAN OR MEAN

mset <- readRDS("/Volumes/NYC_geo/tree_classification/mset_v2_wlidarvars.rds")

mset$genus %>% table()


# 2018 - 2024
# 4 band setup, but make possible to do 8 band

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as range-year combinations, maybe do by DOY

# 1 day
# will need to sort out these available images

getRanges4b <- function(week_var, n_shifted_day){
  x <- c(seq(ymd('2018-01-01'), ymd('2018-12-31'), by = week_var), # will need to drop the last day of the year
        seq(ymd('2018-01-01') + days(n_shifted_day), ymd('2018-12-31'), by = week_var), # shifted by +3 days
        
        seq(ymd('2019-01-01'), ymd('2019-12-31'), by = week_var),
        seq(ymd('2019-01-01') + days(n_shifted_day), ymd('2019-12-31'), by = week_var),
        
        seq(ymd('2020-01-01'), ymd('2020-12-31'), by = week_var),
        seq(ymd('2020-01-01') + days(n_shifted_day), ymd('2020-12-31'), by = week_var),
        
        seq(ymd('2021-01-01'), ymd('2021-12-31'), by = week_var),
        seq(ymd('2021-01-01') + days(n_shifted_day), ymd('2021-12-31'), by = week_var),
        
        seq(ymd('2022-01-01'), ymd('2022-12-31'), by = week_var),
        seq(ymd('2022-01-01') + days(n_shifted_day), ymd('2022-12-31'), by = week_var),
        
        seq(ymd('2023-01-01'), ymd('2023-12-31'), by = week_var),
        seq(ymd('2023-01-01') + days(n_shifted_day), ymd('2023-12-31'), by = week_var),
        
        seq(ymd('2024-01-01'), ymd('2024-12-31'), by = week_var),
        seq(ymd('2024-01-01') + days(n_shifted_day), ymd('2024-12-31'), by = week_var))
  x <- sort(x)
  x2 <- x[which(yday(x) < 364)] # remove dates at the end of the year, not needed and overlapping
  return(x2)
}

# 1 week
wk1_start_days <- getRanges4b("1 week", 3)

# 2 week
wk2_start_days <- getRanges4b("2 week", 7)

# 4 week
wk4_start_days <- getRanges4b("4 week", 14)

# 8 week
wk8_start_days <- getRanges4b("8 week", 28)

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# Modifying from: 4ms_tree_spectra_agg_monthcomp_4b.R

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
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


# get full list of tree IDs (this was originally on the entire polygon file, but only need it for our sample here)
tree_ids <- mset$Poly_ID

stepsize <- 5000
top <- floor(length(tree_ids)/stepsize)*stepsize+1
step_ranges <- seq(1, top, stepsize)

# this is all setup for 1 week
for (idset in 1:length(step_ranges)){
  #idset <- 30
  start_time <- Sys.time()
  range_min <- step_ranges[idset]
  idset_range <- seq(range_min, range_min + stepsize - 1)
  tree_ids_sub <- tree_ids[idset_range]
  print(idset)
  
  setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point") # 4 b setup
  tree_id_spectra <- getTreePlanetSpectra2(tree_ids_sub)
  
  # take mean for each date range
  # take mean for each month (across all years)
  tree_id_spectra$date <- ymd(tree_id_spectra$date)
  
  for (j in c(1, 2, 4, 8)){
    n_weeks_ts <- j
    if (j == 1){
      unq_names <- wk1_start_days
    } else if (j == 2){
      unq_names <- wk2_start_days
    } else if (j == 4){
      unq_names <- wk4_start_days
    } else if (j == 8){
      unq_names <- wk8_start_days
    } else{
      print("Invalid date range")
      break
    }
    
    init <- 0 # initializer for appending
    for (i in 1:length(unq_names)){
      #for (i in 1:12){
      print(i)
      #i <- 50
      #i <- 1
      time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i]) + weeks(n_weeks_ts) - days(1)) # a one week time span, inclusive
      tree_id_spectra_sub <- tree_id_spectra %>% filter(date %within% time_span) %>% as.data.frame()
      # if statement for if there is any data in the range
      if (nrow(tree_id_spectra_sub) > 0){
        print(paste0(time_span, " Tree spectra found"))
        # need to aggregate here
        band_list <- c("blue", "green", "red", "nir")
        tree_id_spectra_date_range_mean <- aggregate(tree_id_spectra_sub[,band_list], by = list(tree_id_spectra_sub$Object_ID), FUN = agg_fun, na.rm = TRUE) 
        # THIS IS NOW MEDIAN! Slower than mean.
        # Median chosen based on McMahon et al. 2024 on handling cloudy bad data over ranges
        # Mean would be more likely to emphasize outlier values, both good and bad, so would need to think about it
        
        colnames(tree_id_spectra_date_range_mean)[1] <- c("Object_ID")
        tree_id_spectra_date_range_mean$Date_Range_Group <- paste0(n_weeks_ts, "wk_", unq_names[i]) # organized by X week count, start date (since week time ranges overlap)
        # for multiyear collapse, need to decompose this variable to get DOY, then collapse mean on that. Mean of means, ie each year is equally weighted.
        
        tree_id_spectra_date_range_mean$n_img_days <- table(tree_id_spectra_sub$Object_ID)
        
        
        # OMNBR step (all band combinations, not true omnbr but would be a setup step)
        # do this after loading extracted data in next script, would be quick and that would save space on output here (1/6 size as RDS)
        # colname_list <- c()
        # 
        # for (i in 1:(length(band_list)-1)){
        #   b1 <- band_list[i]
        #   for (j in (i+1):length(band_list)){
        #     b2 <- band_list[j]
        #     
        #     # Do longer wavelength first since we might expect longer to be greater in NDVI form
        #     band_name <- paste("nd", b2, b1, sep = "_")
        #     #print(band_name)
        #     
        #     colname_list <- c(colname_list, band_name)
        #     
        #     nd_index <- calcNormDif(tree_id_spectra_date_range_mean[,b2], tree_id_spectra_date_range_mean[,b1]) # note this only works on data frames, not tibbles
        #     
        #     if (i == 1 & j == 2){
        #       omnbr <- nd_index
        #     } else {
        #       omnbr <- cbind.data.frame(omnbr, nd_index)
        #       colnames(omnbr) <- colname_list
        #     }
        #     
        #   }
        # }
        
        #df_spectra_omnbr <- bind_cols(tree_id_spectra_date_range_mean, omnbr)
        #df_spectra_omnbr <- tree_id_spectra_date_range_mean # doesn't have omnbr, just a tester
        
        if (init == 0){
          #df_spectra_all <- df_spectra_omnbr
          df_spectra_all <- tree_id_spectra_date_range_mean
          init <- 1
        } else {
          #df_spectra_all <- bind_rows(df_spectra_all, df_spectra_omnbr)
          df_spectra_all <- bind_rows(df_spectra_all, tree_id_spectra_date_range_mean)
        }
      }
    }
    
    df_spectra_all <- df_spectra_all %>% relocate(Date_Range_Group, .after = Object_ID) %>% relocate(n_img_days, .after = Date_Range_Group)
    
    # do the pivot_long to pivot_wide flip in another step, will generate many NAs in data frame otherwise and need to filter later
    #df_spectra_omnbr_long <- df_spectra_omnbr %>% pivot_longer(cols = colnames(df_spectra_omnbr)[3:ncol(df_spectra_omnbr)])
    #df_spectra_omnbr_long <- df_spectra_omnbr_long %>% mutate(name_drg = paste0(name, "_mean_", Date_Range_Group))
    #df_spectra_omnbr_wide <- df_spectra_omnbr_long %>% pivot_wider(id_cols = c(Object_ID), names_from = name_drg, values_from = value)
    
    outpath <- paste0("/Volumes/NYC_geo/tree_classification/extracted_timeranges/", agg_fun, "/single_years") # MEDIAN OR MEAN
    saveRDS(df_spectra_all,
            paste0(outpath, "/", n_weeks_ts, "week/tree_points_4b_spectra_raw_", n_weeks_ts, "weekcomplong_idset", idset, ".rds"))
    # can make output files smaller by multiplier for integer and/or saving to RDS (instead of csv), but I don't know the correct value ranges that will be important for the band indices...
    # rds is much smaller, check parquet
    # write_parquet(df_spectra_all,
    #         paste0(outpath, "/", n_weeks_ts, "week/tree_points_4b_spectra_raw_", n_weeks_ts, "weekcomplong_idset", idset, ".parquet")) # parquet is bigger than RDS
  }
}

#####
# Collapse all into collapsed year value, for example 1wk_DOY. Use median for median and use mean for mean for cross year collapse

#n_weeks_ts <- 1

for (n_weeks_ts in c(1, 2, 4, 8)){
  
  setwd(paste0("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median/single_years/", n_weeks_ts, "week/"))
  rds_files <- list.files()
  all_single_year <- purrr::map(rds_files, readRDS) %>% list_rbind()
  
  #n_weeks_ts <- 1
  
  all_single_year$doy <- str_split(all_single_year$Date_Range_Group, "_") %>% lapply('[[', 2) %>% ymd() %>% yday()
  timerange_label <- str_split(all_single_year$Date_Range_Group[1], "_")[[1]][1]
  
  all_single_year$n_img_days_num <- as.numeric(all_single_year$n_img_days)
  
  # need to figure out how to do this efficiently if this is the way to do this
  unique_doy <- all_single_year$doy %>% unique() %>% sort()
  band_list <- c("blue", "green", "red", "nir")
  init <- 0
  for (i in 1:length(unique_doy)){
    print(i)
    #i <- 1
    
    all_collapsed_year_sub <- aggregate(all_single_year[unique_doy_inds,band_list], by = list(all_single_year$Object_ID[unique_doy_inds], all_single_year$doy[unique_doy_inds]), FUN = agg_fun, na.rm = TRUE) 
    all_collapsed_year_img_sums_sub <- aggregate(all_single_year$n_img_days_num[unique_doy_inds], by = list(all_single_year$Object_ID[unique_doy_inds], all_single_year$doy[unique_doy_inds]), FUN = "sum", na.rm = TRUE) 
    all_collapsed_year_sub$n_img_counts_sum <- all_collapsed_year_img_sums_sub$x
    
    if (init == 0){
      all_collapsed_year <- all_collapsed_year_sub
      init <- 1
    } else {
      all_collapsed_year <- rbind.data.frame(all_collapsed_year, all_collapsed_year_sub)
    }
  }
  
  colnames(all_collapsed_year)[1:2] <- c("Object_ID", "Agg_StartDOY")
  all_collapsed_year$Agg_StartDOY <- paste0(timerange_label, "_doy", sprintf("%03d", all_collapsed_year$Agg_StartDOY))
  
  outpath <- paste0("/Volumes/NYC_geo/tree_classification/extracted_timeranges/", agg_fun, "/collapsed_years") # MEDIAN OR MEAN
  saveRDS(all_collapsed_year,
          paste0(outpath, "/tree_points_4b_spectra_raw_", n_weeks_ts, "weekcomplong_multiyearcollapsed.rds"))
  
}
