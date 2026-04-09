library(data.table)
library(tidyverse)
library(lubridate)
library(phenofit)
library(purrr)
library(future)
library(future.apply)
library(weights) # for weighted correlation

# v6 - adding in R2 fits and NDVI amplitude for spring greenup curves

# v5 - v4 but revised to include columns for number of observations near SOS and EOS dates
# n_SOS_3day : +/- 3 days (within week)
# n_SOS_7day : +/- 7 days  (within 2 weeks)
# n_SOS_14day :+/- 14 days (within ~month)\
# n_SOS_20_80 : Also do within the SOS_20 and SOS_80 dates, and EOS_80 and EOS_20 (will need to retain these values too)
# as well as EOS flavors of all of the above

# try for tree_id_list
# tree_id_list <- unique(tree_df_all$Object_ID[1:2000])

# would need to loop over for mean, ndvi75, point
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
agg_method <- "point" # "mean", "ndvi75", "point"

#output_dir <- "/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/raw_pheno_extract" # with no trailing "/"
output_dir <- "/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/raw_pheno_extract_springfit" # with no trailing "/"

calcNDVI <- function(NIR, red){
  return((NIR - red)/(NIR + red))
}

# USE THIS
setupPhenoElmore <- function(tree_ts_df_complete, year_val){
  # make sure the labeling of the "date" column is formatted correctly
  tree_ts_df_complete_oneyear <- tree_ts_df_complete[which(year(tree_ts_df_complete$date) ==  year_val),] # get which values are in this year
  t <- tree_ts_df_complete_oneyear$date
  vi <- tree_ts_df_complete_oneyear$NDVI
  vi_checked <- check_input(t, vi)
  tout <- seq(1, 366, 1) # need to be the same length, so assume a year of 366 days
  c_out <- curvefit(vi_checked$y, yday(vi_checked$t), tout, w = vi_checked$w, methods = c("Elmore")) # needed the t_out variable, this fixed it. Not sure why
  rm(tree_ts_df_complete_oneyear, t, vi, vi_checked, tout)
  #gc()
  return(c_out)
}

# subs in the correct year instead of 2000, not sure why that is the default
# this is revised to correctly account for feb 29
getCorYear <- function(l_pheno_df, colname){
  new_date <- as.Date(yday(l_pheno_df[[colname]]) - 1, origin = paste0(l_pheno_df$flag, "-01-01"))
  return(new_date)
}

getNumDays <- function(df_tree, l_pheno_sub, colname1, colname2, n, yr_ind){
  return(length(which(df_tree$date %within% interval(l_pheno_sub[[colname1]][yr_ind] - days(n), l_pheno_sub[[colname2]][yr_ind] + days(n)))))
}

#df_tree <- tree_df_sub1tree
runPhenoElmore <- function(df_tree){
  # run setup for each year desired
  
  fits <- list('2017' = setupPhenoElmore(df_tree, 2017),
               '2018' = setupPhenoElmore(df_tree, 2018),
               '2019' = setupPhenoElmore(df_tree, 2019),
               '2020' = setupPhenoElmore(df_tree, 2020),
               '2021' = setupPhenoElmore(df_tree, 2021),
               '2022' = setupPhenoElmore(df_tree, 2022),
               '2023' = setupPhenoElmore(df_tree, 2023),
               '2024' = setupPhenoElmore(df_tree, 2024))
  l_pheno   <- get_pheno(fits, method = "Elmore", TRS = c(0.2, 0.5, 0.8), asymmetric = TRUE)
  
  l_pheno_sub <- l_pheno$date$Elmore[,c("flag", "origin", "TRS5.sos", "TRS5.eos",  "TRS2.sos", "TRS2.eos",  "TRS8.sos", "TRS8.eos")]
  
  # get SOS and EOS value ranges
  l_pheno_sub$TRS5.sos <- getCorYear(l_pheno_sub, "TRS5.sos")
  l_pheno_sub$TRS5.eos <- getCorYear(l_pheno_sub, "TRS5.eos")
  
  l_pheno_sub$TRS2.sos <- getCorYear(l_pheno_sub, "TRS2.sos")
  l_pheno_sub$TRS2.eos <- getCorYear(l_pheno_sub, "TRS2.eos")
  
  l_pheno_sub$TRS8.sos <- getCorYear(l_pheno_sub, "TRS8.sos")
  l_pheno_sub$TRS8.eos <- getCorYear(l_pheno_sub, "TRS8.eos")
  
  l_pheno_sub$SOS_50 <-  yday(l_pheno_sub$TRS5.sos)
  l_pheno_sub$EOS_50 <- yday(l_pheno_sub$TRS5.eos)
  
  l_pheno_sub$SOS_20 <-  yday(l_pheno_sub$TRS2.sos)
  l_pheno_sub$EOS_20 <- yday(l_pheno_sub$TRS2.eos)
  
  l_pheno_sub$SOS_80 <-  yday(l_pheno_sub$TRS8.sos)
  l_pheno_sub$EOS_80 <- yday(l_pheno_sub$TRS8.eos)
  
  # get number of observations within time ranges
  
  l_pheno_sub$n_SOS_50_3day <- NA
  l_pheno_sub$n_EOS_50_3day <- NA
  l_pheno_sub$n_SOS_50_7day <- NA
  l_pheno_sub$n_EOS_50_7day <- NA
  l_pheno_sub$n_SOS_50_14day <- NA
  l_pheno_sub$n_EOS_50_14day <- NA
  l_pheno_sub$n_SOS_20_80 <- NA
  l_pheno_sub$n_EOS_20_80 <- NA
  
  # get results of the fits for each year
  l_pheno_sub$NDVI_fit_SOS_min <- NA
  l_pheno_sub$NDVI_fit_SOS_max <- NA
  l_pheno_sub$NDVI_fit_SOS_20 <- NA
  l_pheno_sub$NDVI_fit_SOS_50 <- NA
  l_pheno_sub$NDVI_fit_SOS_80 <- NA
  l_pheno_sub$SOS_R2_20_80 <- NA
  l_pheno_sub$SOS_R2_wtd_20_80 <- NA
  #l_pheno_sub$RMSE_20_80 <- NA
  
  
  for (yr_ind in 1:nrow(l_pheno_sub)) {
    #print(yr_ind)
    l_pheno_sub$n_SOS_50_3day[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS5.sos", "TRS5.sos", 3, yr_ind)
    l_pheno_sub$n_EOS_50_3day[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS5.eos", "TRS5.eos", 3, yr_ind)
    
    l_pheno_sub$n_SOS_50_7day[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS5.sos", "TRS5.sos", 7, yr_ind)
    l_pheno_sub$n_EOS_50_7day[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS5.eos", "TRS5.eos", 7, yr_ind)
    
    l_pheno_sub$n_SOS_50_14day[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS5.sos", "TRS5.sos", 14, yr_ind)
    l_pheno_sub$n_EOS_50_14day[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS5.eos", "TRS5.eos", 14, yr_ind)
    
    l_pheno_sub$n_SOS_20_80[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS2.sos", "TRS8.sos", 0, yr_ind)
    l_pheno_sub$n_EOS_20_80[yr_ind] <- getNumDays(df_tree, l_pheno_sub, "TRS8.eos", "TRS2.eos", 0, yr_ind) # flipped for EOS
    
    # Adding spring fit characteristics 
    model_ndvi <- fits[[yr_ind]]$model$Elmore$zs$iter2
    l_pheno_sub$NDVI_fit_SOS_min[yr_ind] <- model_ndvi[1] %>% round(6) # pick the first element
    l_pheno_sub$NDVI_fit_SOS_max[yr_ind] <- max(model_ndvi) %>% round(6) 
    
    sos20 <- l_pheno_sub$SOS_20[yr_ind]
    l_pheno_sub$NDVI_fit_SOS_20[yr_ind] <- model_ndvi[sos20] %>% round(6) 
    
    sos50 <- l_pheno_sub$SOS_50[yr_ind]
    l_pheno_sub$NDVI_fit_SOS_50[yr_ind] <- model_ndvi[sos50] %>% round(6) 
    
    sos80 <- l_pheno_sub$SOS_80[yr_ind]
    l_pheno_sub$NDVI_fit_SOS_80[yr_ind] <- model_ndvi[sos80] %>% round(6) 
    
    if (!is.na(sos20) && !is.na(sos50) && !is.na(sos80)){
      ps_20_80 <- fits[[yr_ind]]$data %>% filter(t >= sos20, t <= sos80)
      w_inds <- which(fits[[yr_ind]]$data$t %in% ps_20_80$t)
      if(nrow(ps_20_80) > 1){ # need at least two elements to do R2
        l_pheno_sub$SOS_R2_20_80[yr_ind] <- cor(ps_20_80$y, model_ndvi[ps_20_80$t])^2 %>% round(6) 
        weights <- fits[[yr_ind]]$model$Elmore$ws$iter2[w_inds]
        l_pheno_sub$SOS_R2_wtd_20_80[yr_ind] <- wtd.cor(ps_20_80$y, model_ndvi[ps_20_80$t], weight = weights)[1]^2 %>% round(6) 
        # could also add similar values too 
        #l_pheno_sub$RMSE_20_80[yr_ind] <- NA
      }
    }
  }
  
  
  gof <- get_GOF(fits)
  pheno_output<- merge(l_pheno_sub[, c("flag", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80", 
                                       "n_SOS_50_3day", "n_EOS_50_3day", "n_SOS_50_7day", "n_EOS_50_7day",
                                       "n_SOS_50_14day", "n_EOS_50_14day", "n_SOS_20_80", "n_EOS_20_80",
                                       "NDVI_fit_SOS_min", "NDVI_fit_SOS_max", "NDVI_fit_SOS_20", "NDVI_fit_SOS_50", "NDVI_fit_SOS_80",
                                       "SOS_R2_20_80", "SOS_R2_wtd_20_80")], 
                       gof, 
                       by = "flag")
  pheno_output$Object_ID <- df_tree$Object_ID[1] # will only be one tree ID
  
  # get number of observations (images) within ranges of SOS_50 and EOS_50
  
  rm(fits, l_pheno, l_pheno_sub, gof)
  #gc()
  
  return(pheno_output)
}

#n_year <- 8
calcTreePheno <- function(i, idrange, set_num, agg_method, tree_df_all, n_year, output_dir){
  
  row_start <- i
  row_end <- i + idrange - 1
  
  tree_df_sub <- tree_df_all[which(tree_df_all$Object_ID %in% seq(row_start, row_end)),]
  
  if (agg_method == "mean"){
    tree_df_sub$mean.blue[tree_df_sub$mean.blue < 0] <- 0
    tree_df_sub$mean.green[tree_df_sub$mean.green < 0] <- 0
    tree_df_sub$mean.red[tree_df_sub$mean.red < 0] <- 0
    tree_df_sub$mean.nir[tree_df_sub$mean.nir < 0] <- 0
    tree_df_sub$NDVI <- calcNDVI(tree_df_sub$mean.nir, tree_df_sub$mean.red)
  } else if (agg_method == "ndvi75"){
    tree_df_sub$NDVI <- tree_df_sub$q75/10000
  } else if (agg_method == "point"){ # only using point method now, but could do other methods
    tree_df_sub$blue[tree_df_sub$blue < 0] <- 0
    tree_df_sub$green[tree_df_sub$green < 0] <- 0
    tree_df_sub$red[tree_df_sub$red < 0] <- 0
    tree_df_sub$nir[tree_df_sub$nir < 0] <- 0
    tree_df_sub$NDVI <- calcNDVI(tree_df_sub$nir, tree_df_sub$red)
  } else {
    print("Not a valid method")
  }
  
  #tree_df_sub$date <- ymd(tree_df_sub$date) # doing this outside the function now, redundant
  tree_id_list <- unique(tree_df_sub$Object_ID)
  row_ind <- 1
  
  # preallocate a big data frame to fill with all the polygon information
  # will have to update the ncol, this is hard coded
  #all_pheno_output <- as.data.frame(matrix(data = NA, nrow = n_year*length(tree_id_list), ncol = 23)) # update the number of columns here, ncol = 11 -> 23 with 20 and 80 thresholds and n obs ranges
  all_pheno_output <- as.data.frame(matrix(data = NA, nrow = n_year*length(tree_id_list), ncol = 30)) # now 30 columns
  
  for (j in 1:length(tree_id_list)){
    print(j)
    tree_df_sub1tree <- tree_df_sub[which(tree_df_sub$Object_ID == tree_id_list[j])]
    #pheno_output <- runPhenoElmore(tree_df_sub1tree) # update this function here
    pheno_output <- tryCatch({runPhenoElmore(tree_df_sub1tree)}, 
                             error = function(e) {
                               # if there's an error with runPhenoElmore, skip that tree and leave it as NA
                               print("Error in runPhenoElmore, filling with NA")
                               pheno_output_base <- as.data.frame(matrix(data = NA, nrow = 8, ncol = 30))
                               colnames(pheno_output_base) <-  c("flag", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80", "n_SOS_50_3day", "n_EOS_50_3day", "n_SOS_50_7day",
                                                                 "n_EOS_50_7day", "n_SOS_50_14day", "n_EOS_50_14day", "n_SOS_20_80", "n_EOS_20_80", "NDVI_fit_SOS_min", "NDVI_fit_SOS_max", "NDVI_fit_SOS_20",  "NDVI_fit_SOS_50", "NDVI_fit_SOS_80", 
                                                                 "SOS_R2_20_80", "SOS_R2_wtd_20_80", "meth", "R2", "NSE", "R", "RMSE", "pvalue", "n_sim", "Object_ID")
                               pheno_output_base$flag <- as.character(2017:2024)
                               return(pheno_output_base)
                             }) # update this function here
    pheno_output$Object_ID <- unique(tree_df_sub1tree$Object_ID)
    if (j == 1){
      colnames(all_pheno_output) <- colnames(pheno_output)
    }
    all_pheno_output[(row_ind):(row_ind + n_year - 1),] <- pheno_output # shift by year range
    row_ind <- row_ind + n_year # shift by year range
    rm(pheno_output)
  }
  
  all_pheno_output$SOS_50_date <- as.Date(all_pheno_output$SOS_50-1, origin = paste0(all_pheno_output$flag, "-01-01"))
  all_pheno_output$EOS_50_date <- as.Date(all_pheno_output$EOS_50-1, origin = paste0(all_pheno_output$flag, "-01-01"))
  
  all_pheno_output$SOS_20_date <- as.Date(all_pheno_output$SOS_20-1, origin = paste0(all_pheno_output$flag, "-01-01"))
  all_pheno_output$EOS_20_date <- as.Date(all_pheno_output$EOS_20-1, origin = paste0(all_pheno_output$flag, "-01-01"))
  
  all_pheno_output$SOS_80_date <- as.Date(all_pheno_output$SOS_80-1, origin = paste0(all_pheno_output$flag, "-01-01"))
  all_pheno_output$EOS_80_date <- as.Date(all_pheno_output$EOS_80-1, origin = paste0(all_pheno_output$flag, "-01-01"))
  
  colnames(all_pheno_output)[1] <- "Year"
  
  all_pheno_output$agg_method <- agg_method
  
  write.csv(all_pheno_output, 
            paste0(output_dir, "/trees_pheno_output_springfit_objset", #"/trees_pheno_output_objset", # output path dir with no trailing "/" # adding "springfit" here
                   set_num,"_", as.character(as.integer(row_start)), "_", as.character(as.integer(row_end)), "_", agg_method, ".csv"), 
            row.names = FALSE) 
  
  # thinks that it runs out of application memory so clean up at the end of each run
  rm(tree_df_sub, tree_df_sub1tree, all_pheno_output, row_ind, tree_id_list)
  gc()
}

#plan(multisession, workers = 10) # trying 10, usually at 8 to prevent memory blow up, start parallel inside loop instead
#set_num in 5:7
#set_num <- 1

#for (set_num in 1:7){
for (set_num in 7:7){
  #set_num <- 2
  
  #plan(multisession, workers = availableCores(omit = 2)) # start parallel again each loop
  plan(multisession, workers = 6) # start parallel again each loop
  
  tree_file_list <- list.files(pattern = glob2rx(paste0("nyc_trees_objset", as.character(set_num), "*.csv"))) # stored as different sets by object id
  names(tree_file_list) <- strsplit(tree_file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 27) %>% ymd() # giving a date name for each
  # potentially big memory pull
  tree_df_all <- purrr::map_df(tree_file_list, fread, .id = 'date') # works quite fast to get all files, cannot be loaded for point, these need to be split up to objset too
  # this might be causing R to abort randomly in certain runs, perhaps issue with data.table library
  #print(paste("Loaded set num ", set_num))
  
  tree_df_all$date <- ymd(tree_df_all$date)
  
  idrange <- 1000
  
  tree_ids_to_apply <- unique(tree_df_all$Object_ID)
  tree_id_ranges <- seq(min(tree_ids_to_apply), max(tree_ids_to_apply), idrange)
  
  future_lapply(tree_id_ranges, FUN = calcTreePheno, idrange, set_num, "point", tree_df_all, 8, output_dir, future.seed = TRUE)
  
  rm(tree_df_all, tree_file_list, tree_ids_to_apply, tree_id_ranges)
  gc()
  plan(sequential) # remove parallelization with each loop to get fresh memory
}