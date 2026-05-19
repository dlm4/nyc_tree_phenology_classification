# Setup random forest runs
# use ranger

# v2 using variable feature selection, backwards recusion

library(tidyverse)
library(ranger)
library(caret) # for confusion matrix
#library(varSelRF) # tried this to see what it comes up with, was taking too long

# tree information, mset
mset <- readRDS("/Volumes/NYC_geo/tree_classification/mset_v2_wlidarvars.rds")

#for (input_file in 1:8){
#for (input_file in c(1,2,5,6,7)){
for (input_file in c(11)){
  #input_file <- 4
  
  # Aggragated (collapsed)
  # 1 week
  if (input_file == 1){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_subset_wide.rds")
    feature_set_name <- "1week_collapsed"
  }
  
  # 2 week
  if (input_file == 2){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_subset_wide.rds")
    feature_set_name <- "2week_collapsed"
  }
  
  # 4 week
  if (input_file == 3){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_subset_wide.rds")
    feature_set_name <- "4week_collapsed"
  }
  
  # 8 week
  if (input_file == 4){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds")
    feature_set_name <- "8week_collapsed"
  }
  
  # Single years
  # 1 week
  if (input_file == 5){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_subset_wide.rds")
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
    feature_set_name <- "1week_singleyears"
  }
  
  # 2 week
  if (input_file == 6){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_singleyear_subset_wide.rds")
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
    feature_set_name <- "2week_singleyears"
  }
  
  # 4 week
  if (input_file == 7){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_singleyear_subset_wide.rds")
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
    feature_set_name <- "4week_singleyears"
  }
  
  # 8 week
  if (input_file == 8){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_singleyear_subset_wide.rds")
    names(df_spectra) <- gsub("-", "", names(df_spectra)) # drop hyphens, ranger can't handle these
    feature_set_name <- "8week_singleyears"
  }
  
  getSelVars <- function(file_path_rds){
    x <- readRDS(file_path_rds)
    sel_vars <- x$vars_rounds$var[which(x$vars_rounds$max_round >= x$iter_round_selected)]
    return(sel_vars)
  }
  
  # All collapsed
  if (input_file == 9){
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/1week_collapsed_selected_features_info.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/2week_collapsed_selected_features_info.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/4week_collapsed_selected_features_info.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/8week_collapsed_selected_features_info.rds")
    
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    df1wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v1wk)))
    df2wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v2wk)))
    df4wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v4wk)))
    df8wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v8wk)))
    
    df_spectra <- merge(df1wk, df2wk) %>% merge(df4wk) %>% merge(df8wk)
    
    feature_set_name <- "allweek_collapsed"
  }
  
  # All single year
  if (input_file == 10){
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/1week_singleyears_selected_features_info.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/2week_singleyears_selected_features_info.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/4week_singleyears_selected_features_info.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/8week_singleyears_selected_features_info.rds")
    
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    df1wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_subset_wide.rds")
    names(df1wk) <- gsub("-", "", names(df1wk))
    df2wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_singleyear_subset_wide.rds")
    names(df2wk) <- gsub("-", "", names(df2wk))
    df4wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_singleyear_subset_wide.rds")
    names(df4wk) <- gsub("-", "", names(df4wk))
    df8wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_singleyear_subset_wide.rds")
    names(df8wk) <- gsub("-", "", names(df8wk))
    
    df1wk <- df1wk %>% select(all_of(c("Object_ID", v1wk)))
    df2wk <- df2wk %>% select(all_of(c("Object_ID", v2wk)))
    df4wk <- df4wk %>% select(all_of(c("Object_ID", v4wk)))
    df8wk <- df8wk %>% select(all_of(c("Object_ID", v8wk)))
    
    df_spectra <- merge(df1wk, df2wk) %>% merge(df4wk) %>% merge(df8wk)
    
    feature_set_name <- "allweek_singleyear"
  }
  
  # Everything combined!
  if (input_file == 11){
    # Collapsed
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/1week_collapsed_selected_features_info.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/2week_collapsed_selected_features_info.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/4week_collapsed_selected_features_info.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/8week_collapsed_selected_features_info.rds")
    
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    df1wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v1wk)))
    df2wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v2wk)))
    df4wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v4wk)))
    df8wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v8wk)))
    
    df_spectra_collapsed <- merge(df1wk, df2wk) %>% merge(df4wk) %>% merge(df8wk)
    
    # Single year
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/1week_singleyears_selected_features_info.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/2week_singleyears_selected_features_info.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/4week_singleyears_selected_features_info.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/8week_singleyears_selected_features_info.rds")
    
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    df1wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_subset_wide.rds")
    names(df1wk) <- gsub("-", "", names(df1wk))
    df2wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_singleyear_subset_wide.rds")
    names(df2wk) <- gsub("-", "", names(df2wk))
    df4wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_singleyear_subset_wide.rds")
    names(df4wk) <- gsub("-", "", names(df4wk))
    df8wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_singleyear_subset_wide.rds")
    names(df8wk) <- gsub("-", "", names(df8wk))
    
    df1wk <- df1wk %>% select(all_of(c("Object_ID", v1wk)))
    df2wk <- df2wk %>% select(all_of(c("Object_ID", v2wk)))
    df4wk <- df4wk %>% select(all_of(c("Object_ID", v4wk)))
    df8wk <- df8wk %>% select(all_of(c("Object_ID", v8wk)))
    
    df_spectra_singleyear <- merge(df1wk, df2wk) %>% merge(df4wk) %>% merge(df8wk)
    
    # Select and merge both from previously reduced for each version
    v_collapsed <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/allweek_collapsed_selected_features_info.rds")
    v_singleyear <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/allweek_singleyear_selected_features_info.rds")
    
    df_spectra_collapsed <- df_spectra_collapsed %>% select(all_of(c("Object_ID", v_collapsed)))
    df_spectra_singleyear <- df_spectra_singleyear %>% select(all_of(c("Object_ID", v_singleyear)))
    
    df_spectra <- merge(df_spectra_collapsed, df_spectra_singleyear)
                                                              
    feature_set_name <- "everything"
  }
  
  
  print(feature_set_name)
  
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
  
  # This is very slow with no updates 
  #vsrf <- varSelRF(train_spectra[,2:ncol(train_spectra)], Class = train_spectra$genus, verbose = TRUE)
  
  oob_out <- c()
  oa_out <- c()
  kappa_out <- c()
  cilo_out <- c()
  cihi_out <- c()
  nvar_out <- c()
  n_vars <- ncol(train_spectra) - 1
  
  vars <- colnames(train_spectra)[2:ncol(train_spectra)]
  vars_rounds <- data.frame(vars, rep(1, length(vars)))
  colnames(vars_rounds) <- c("var", "max_round")
  
  conf_mat_list <- list()
  
  # variable selection
  round_iter <- 2
  round_out <- c(1)
  while(n_vars > 5){
    
    print(paste0("Round iteration: ", round_iter - 1, ",  n_vars: ", n_vars))
    i <- 1
    set.seed(i) # here for ranger seed, not just sampling seed
    rf_trained <- ranger(genus ~ ., train_spectra, importance = 'permutation', num.threads = 10, num.trees = 500, seed = 1) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
    
    oob_out <- c(oob_out, rf_trained$prediction.error)
    
    val_predicted <- predict(rf_trained, data = val_spectra)
    conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)
    oa_out <- c(oa_out, conf_mat$overall[1])
    kappa_out <- c(kappa_out, conf_mat$overall[2])
    cilo_out <- c(cilo_out, conf_mat$overall[3])
    cihi_out <- c(cihi_out, conf_mat$overall[4])
    nvar_out <- c(nvar_out, n_vars)
    round_out <- c(round_out, round_iter)
    
    n_el_retain <- round(length(rf_trained$variable.importance)*0.8, 0)
    
    sorted_vars <- sort(rf_trained$variable.importance, decreasing = TRUE)
    
    kept_vars <- sorted_vars[1:n_el_retain]
    
    conf_mat_list[[round_iter - 1]] <- conf_mat
    
    vars_rounds$max_round[which(vars_rounds$var %in% names(kept_vars))] <- round_iter
    round_iter <- round_iter + 1
    
    train_spectra <- train_spectra[, c("genus", names(kept_vars))]
    n_vars <- ncol(train_spectra) - 1
    
    # could also implement a way for this to break earlier instead of running all the way down to only 5...
  }
  
  # 95% CI is ~ double the SD for a normal distribution...
  # For our purposes, would be reasonable to choose which has 95% CI overlapping with 
  
  data_out <- data.frame(round_out[1:(length(round_out)-1)], nvar_out, oob_out, oa_out, kappa_out, cilo_out, cihi_out)
  colnames(data_out)[1] <- "Round" 
  
  # which is the highest round (ie lowest number of variables) st the 95% overlaps with the maximum OA
  iter_round_selected <- max(data_out$Round[which(data_out$cihi_out > max(data_out$oa_out))])
  #data_out[iter_round_selected,]
  varsel_f <- vars_rounds$var[which(vars_rounds$max_round >= iter_round_selected)] # these are the selected variables.
  max_acc <- data_out$oa_out[which.max(data_out$oa_out)]
  
  ggplot(data_out) +
    geom_vline(xintercept = length(varsel_f), color = "navy", linetype = "dotted") +
    geom_hline(yintercept = max_acc, color = "red", linetype = "dotted") +
    geom_linerange(aes(nvar_out, ymin = cilo_out, ymax = cihi_out)) +
    geom_point(aes(nvar_out, oa_out)) +
    labs(x = "Number of features", y = "Overall Accuracy", title = feature_set_name) +
    theme_bw()
  
  # To save on output
  
  setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/")
  
  # plot
  ggsave(paste0(feature_set_name, "_feature_selection.png"), width = 6, height = 5, units = "in")
  
  # confusion matrix list object, needed if we want to plot based on accuracy for individual genera later
  saveRDS(conf_mat_list, paste0(feature_set_name, "_conf_mat_list.rds"))
  
  # data_out
  write.csv(data_out, paste0(feature_set_name, "_feature_selection_data.csv"), row.names = FALSE)
  # iter_round_selected
  # var_rounds df
  # selected vars (just for ease of use)
  out_list <- list(iter_round_selected = iter_round_selected, 
                   vars_rounds = vars_rounds, 
                   varsel_f = varsel_f)
  saveRDS(out_list, paste0(feature_set_name, "_selected_features_info.rds"))
  
}

#####
# Need to setup alternate version to do 10 fold cross validation for each round
# For each iteration, setup 10 folds (could do 5 to start for checking)
#   Sample dataset without replacement to make index column
# For each cross val run
# Pull from index column each time to setup training and validation data
# Run training through ranger rf
# retain all values as earlier version
# Then, take mean of importance scores for all variables (for 10 cross val)
# Remove least important, bottom 20% of the variables
# Retain OOB, OA, Kappa, CI_LO and CI_HI for each run
# Take mean of each of these for decision making on "best" model run
# Repeat loop

# At end, select which



#####

  # 
  # # I don't think ranger outputs the OOB SD, I was calculating this incorrectly here
  # oob_se <- c() # not sure if this is the correct format
  # n_vars <- c()
  # for (i in 1:length(oob_out)){
  #   n_var <- length(which(vars_rounds$max_round >= i))
  #   oob_se <- c(oob_se, sqrt(oob_out[i]*(1-oob_out[i]) * 1/n_var))
  #   n_vars <- c(n_vars, n_var)
  # }
  # 
  # oob_se <- sqrt(oob_out[i]*(1-oob_out[i]) * 1/n_var)
  # 
  # min(oob_out) + oob_se - oob_out
  
  # Alternative, multiple runs at a given level, multiple seeds
  # Take SD of importance across multiple runs
  # Remove 20% with lowest SD of importance, because high SD = more possibility of importance
  

  # # test with just two, setup for longer version later, more loops
  # i_range <- 1:2
  # out_df <- data.frame(i_range, 
  #                      rep(NA, length(i_range)),
  #                      rep(NA, length(i_range)),
  #                      rep(NA, length(i_range)))
  # colnames(out_df) <- c("Seed", "OOB_error", "OA", "Kappa")
  # 
  # # Can setup these outputs to also include users and producers accuracies, counts for each genus
  # # Or save out all of confusion matrices in list object
  # 
  # for (i in i_range){
  #   print(i)
  #   start_time <- Sys.time()
  #   set.seed(i) # here for ranger seed, not just sampling seed
  #   rf_trained <- ranger(genus ~ ., train_spectra, num.threads = 10, num.trees = 500, seed = i) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
  #   val_predicted <- predict(rf_trained, data = val_spectra)
  #   conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)
  #   
  #   out_df$Seed[i] <- i
  #   out_df$OOB_error[i] <- rf_trained$prediction.error
  #   out_df$OA[i] <- conf_mat$overall[1]
  #   out_df$Kappa[i] <- conf_mat$overall[2]
  #   
  #   print(Sys.time() - start_time)
  # }
  
  #write.csv(out_df, output_file_path, row.names = FALSE)
  #rm(df_spectra)
  #gc()
#}

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
# 
# # overall importance
# var_imp <- rf_trained$variable.importance %>% sort(decreasing = TRUE)
# var_imp_df <- cbind.data.frame(names(var_imp), var_imp)
# colnames(var_imp_df) <- c("var", "imp")
# 
# ggplot(var_imp_df[1:25,]) + 
#   geom_point(aes(x = imp, y = reorder(var, imp))) + 
#   theme_bw()
# 
# # local importance
# # https://stackoverflow.com/questions/61340327/how-to-obtain-feature-importance-by-class-using-ranger
# library(data.table)    
# genus_var_imp <- as.data.table(rf_trained$variable.importance.local)[,genus := train_spectra$genus][,lapply(.SD,mean),by=genus]
# 
# max_cols <- max.col(genus_var_imp[,2:ncol(genus_var_imp)]) + 1
# colnames(genus_var_imp)[max_cols]
# max_vals <- as.matrix(genus_var_imp[,..max_cols]) %>% diag()
# 
# genus_max_imp_df <- data.frame(genus_var_imp[,1], colnames(genus_var_imp)[max_cols], max_vals)
# 
# 
# # validation prediction
# val_predicted <- predict(rf_trained, data = val_spectra)
# conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)

# To save from each run:
# rf_trained object
# train_inds or val_inds index object
# confusion matrix object
# plot for overall importance
# individual 
