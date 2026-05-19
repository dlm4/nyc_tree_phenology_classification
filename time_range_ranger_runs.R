# Setup random forest runs
# use ranger

library(tidyverse)
library(ranger)
library(caret)

# Do ranger RF runs for each individual aggregated set
# All single year
# All collapsed year
# All together

# Likely need to do multiple RF runs to check stability of importance variables on output

# tree information, mset
mset <- readRDS("/Volumes/NYC_geo/tree_classification/mset_v2_wlidarvars.rds")

#for (input_file in 1:8){
for (input_file in 9:11){
  # Aggragated (collapsed)
  # 1 week
  if (input_file == 1){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_rf_tests.csv"
  }
  
  # 2 week
  if (input_file == 2){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_rf_tests.csv"
  }
  
  # 4 week
  if (input_file == 3){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_rf_tests.csv"
  }
  
  # 8 week
  if (input_file == 4){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_rf_tests.csv"
  }
  
  # Single years
  # 1 week
  if (input_file == 5){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_rf_tests.csv"
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
  }
  
  # 2 week
  if (input_file == 6){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_singleyear_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_singleyear_rf_tests.csv"
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
  }
  
  # 4 week
  if (input_file == 7){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_singleyear_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_singleyear_rf_tests.csv"
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
  }
  
  # 8 week
  if (input_file == 8){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_singleyear_subset_wide.rds")
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_singleyear_rf_tests.csv"
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
  }
  
  # All collapsed
  if (input_file == 9){
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    rds_files <- list.files(pattern = glob2rx("*weekcomplong_multiyearcollapsed_subset_wide.rds"))
    df_spectra <- purrr::map(rds_files, readRDS) %>% list_rbind()
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_all_multiyearcollapsed_rf_tests.csv"
  }
  
  # All single year
  # This will take a very very long time to run (multiple hours, for each iteration)
  if (input_file == 10){
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    rds_files <- list.files(pattern = glob2rx("*weekcomplong_singleyear_subset_wide.rds"))
    df_spectra <- purrr::map(rds_files, readRDS) %>% list_rbind()
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_all_singleyear_rf_tests.csv"
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
  }
  
  # Everything!
  # haven't started this yet
  if (input_file == 11){
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    rds_files <- list.files(pattern = glob2rx("*weekcomplong_*_subset_wide.rds"))
    df_spectra <- purrr::map(rds_files, readRDS) %>% list_rbind()
    output_file_path <- "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_everything_rf_tests.csv"
  }
  
  # This was setup stuff
  # Aggregated spectra data
  # df1c <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_subset_wide.rds")
  # df2c <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_subset_wide.rds")
  # df4c <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_subset_wide.rds")
  # df8c <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds") # this is the biggest aggregation, simplest
  # 
  # df1s <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_subset_wide.rds")
  # df2s <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_singleyear_subset_wide.rds")
  # df4s <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_singleyear_subset_wide.rds")
  # df8s <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_singleyear_subset_wide.rds")
  
  #names(df2s) <- gsub("-", "", names(df2s)) # drop hyphens, ranger can't handle these
  #names(df8s) <- gsub("-", "", names(df8s)) # drop hyphens, ranger can't handle these
  
  #df_spectra <- df8c
  # df_spectra <- df2s
  # df_spectra <- df8s
  
  mset_sub <- mset %>% filter(Poly_ID %in% df_spectra$Object_ID)
  
  gen_counts <- table(mset_sub$genus)
  gen_flist <- names(gen_counts[which(gen_counts >= 1000)]) # minimum 1000 trees in a genus to include it in the RF model
  
  df_spectra <- df_spectra %>% select(-contains("img"))
  genus_id <- mset_sub %>% select(Poly_ID, genus)
  genus_id <- genus_id %>% subset(genus %in% gen_flist)
  
  df_spectra <- merge(df_spectra, genus_id, by.x = "Object_ID", by.y = "Poly_ID")
  
  df_spectra <- df_spectra %>% relocate(genus, .after = Object_ID)
  
  set.seed(321) # keep training and validation the same for every run
  # sample size fraction for training
  samp_frac <- 0.8
  train_samp_size <- round(nrow(df_spectra)*samp_frac, 0)
  all_row_inds <- 1:nrow(df_spectra)
  train_inds <- sample(all_row_inds, train_samp_size)
  val_inds <- all_row_inds[-train_inds]
  df_spectra$genus <- as.factor(df_spectra$genus)
  train_spectra <- df_spectra[train_inds,] %>% select(-contains("Object_ID"))
  val_spectra <- df_spectra[val_inds,] %>% select(-contains("Object_ID"))
  
  # set.seed(i) # here for ranger seed, not just sampling seed
  #rf_trained <- ranger(genus ~ ., train_spectra, importance = 'permutation', local.importance = TRUE, num.threads = 10, num.trees = 500) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
  
  # test with just two, setup for longer version later, more loops
  i_range <- 1:2
  out_df <- data.frame(i_range, 
                       rep(NA, length(i_range)),
                       rep(NA, length(i_range)),
                       rep(NA, length(i_range)))
  colnames(out_df) <- c("Seed", "OOB_error", "OA", "Kappa")
  
  # Can setup these outputs to also include users and producers accuracies, counts for each genus
  # Or save out all of confusion matrices in list object
  
  for (i in i_range){
    print(i)
    start_time <- Sys.time()
    set.seed(i) # here for ranger seed, not just sampling seed
    rf_trained <- ranger(genus ~ ., train_spectra, num.threads = 10, num.trees = 500, seed = i) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
    val_predicted <- predict(rf_trained, data = val_spectra)
    conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)
    
    out_df$Seed[i] <- i
    out_df$OOB_error[i] <- rf_trained$prediction.error
    out_df$OA[i] <- conf_mat$overall[1]
    out_df$Kappa[i] <- conf_mat$overall[2]
    
    print(Sys.time() - start_time)
  }
  
  write.csv(out_df, output_file_path, row.names = FALSE)
  rm(df_spectra)
  gc()
}

# 


# permutation importance can take much longer to calculate than the ranger object itself, becomes hours with thousands of variables
# Should likely do some recursive variable removal

# https://www.r-bloggers.com/2018/06/be-aware-of-bias-in-rf-variable-importance-metrics/
# May be a good idea to do varselrf() to reduce number of input features, so that remaining importance scores are more interepretable
# Can also use tuneRanger (with mlr) to estimate more optimal hyperparameters, if these likely should be varied with different runs

# Save the RF object
# saveRDS(rf_trained, "/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/rf_test_s321_nt500_2week_single_year.rds")

# check importance (from impurity, for speed) with multiple seeds and if the importance changes around, then need to increase the number of rf trees
# conditional permutation importance - how to handle multiple collinear variables...

# need to save the rf model each time, and will need to evaluate if 500 trees is robust or should try something else?

# overall importance
var_imp <- rf_trained$variable.importance %>% sort(decreasing = TRUE)
var_imp_df <- cbind.data.frame(names(var_imp), var_imp)
colnames(var_imp_df) <- c("var", "imp")

ggplot(var_imp_df[1:25,]) + 
  geom_point(aes(x = imp, y = reorder(var, imp))) + 
  theme_bw()

# local importance
# https://stackoverflow.com/questions/61340327/how-to-obtain-feature-importance-by-class-using-ranger
library(data.table)    
genus_var_imp <- as.data.table(rf_trained$variable.importance.local)[,genus := train_spectra$genus][,lapply(.SD,mean),by=genus]

max_cols <- max.col(genus_var_imp[,2:ncol(genus_var_imp)]) + 1
colnames(genus_var_imp)[max_cols]
max_vals <- as.matrix(genus_var_imp[,..max_cols]) %>% diag()

genus_max_imp_df <- data.frame(genus_var_imp[,1], colnames(genus_var_imp)[max_cols], max_vals)


# validation prediction
val_predicted <- predict(rf_trained, data = val_spectra)
conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)

# To save from each run:
# rf_trained object
# train_inds or val_inds index object
# confusion matrix object
# plot for overall importance
# individual 
