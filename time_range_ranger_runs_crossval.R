# Setting up 10 fold cv for ranger runs

# currently using 95% CI error bar overlaps for selection
# but could also use OOB error directly, and error on this, if this is the more appropriate thing to do

library(tidyverse)
library(ranger)
library(caret) # for confusion matrix and 10 fold cross val
`%notin%` <- Negate(`%in%`)

for (input_file in c(3, 2, 1, 5, 6, 7, 8)){
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
  
  if (input_file == 4){
    df_spectra <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds")
    #feature_set_name <- "8week_collapsed_test"
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
  print(feature_set_name)
  
  # combined versions
  
  
  
  # tree information, mset
  mset <- readRDS("/Volumes/NYC_geo/tree_classification/mset_v2_wlidarvars.rds")
  mset_sub <- mset %>% filter(Poly_ID %in% df_spectra$Object_ID)
  
  gen_counts <- table(mset_sub$genus)
  gen_flist <- names(gen_counts[which(gen_counts >= 1000)]) # minimum 1000 trees in a genus to include it in the RF model
  
  df_spectra <- df_spectra %>% select(-contains("img"))
  genus_id <- mset_sub %>% select(Poly_ID, genus)
  genus_id <- genus_id %>% subset(genus %in% gen_flist)
  
  df_spectra <- merge(df_spectra, genus_id, by.x = "Object_ID", by.y = "Poly_ID")
  
  df_spectra <- df_spectra %>% relocate(genus, .after = Object_ID)
  
  # make it smaller for testing - remove and do for real!
  # df_spectra <- df_spectra[1:5000,]
  
  set.seed(321) # keep training and validation the same for every run
  
  # this is where we setup the cross val
  # ctrl <- trainControl(method = "cv", number = 10)
  # 
  # rf_trained_cv <- train(genus ~ ., 
  #       data = df_spectra,
  #       method = "ranger",
  #       importance = 'permutation', num.threads = 10, num.trees = 500, seed = 1)
  # Not sure what is going to come out of this?
  # No updates, will need to check this later
  
  # sample size fraction for training
  # Repeated for 10-fold cross validation, set it up myself because of what we need for outputs
  samp_frac <- 0.1
  val_samp_size <- round(nrow(df_spectra)*samp_frac, 0)
  all_row_inds <- 1:nrow(df_spectra)
  val_samps <- list()
  for (i in 1:9){
    val_inds <- sample(all_row_inds, val_samp_size)
    val_samps[[i]] <- val_inds
    all_row_inds <- all_row_inds[which(all_row_inds %notin% val_inds)]
    print(length(all_row_inds))
  }
  val_samps[[10]] <- all_row_inds # add whatever is left
  
  # loop over all of these iterations
  # take average of accuracies and average of CIs, but save everything else
  # put this inside the while loop
  
  df_spectra <- df_spectra %>% select(-contains("Object_ID"))
  
  # Setup while loop here?
  # while()
  
  #conf_mat_list <- list()
  oob_out <- c()
  oa_out <- c()
  kappa_out <- c()
  cilo_out <- c()
  cihi_out <- c()
  nvar_out <- c()
  n_vars <- ncol(df_spectra) - 1 # remove genus part
  
  all_row_inds <- 1:nrow(df_spectra) # for training resample
  
  # variable selection
  round_iter <- 2
  round_out <- c()
  
  conf_mat_biglist <- list()
  
  var_imp_list <- list()
  
  vars <- colnames(df_spectra)[2:ncol(df_spectra)] # these will be the same within a cross val set
  vars_rounds <- data.frame(vars, rep(1, length(vars)))
  colnames(vars_rounds) <- c("var", "max_round")
  
  #for(round_iter in 2:5){ # just a tester
  while(n_vars > 5){
    
    print(paste0("Round iteration: ", round_iter - 1, ",  n_vars: ", n_vars))
    
    conf_mat_list <- list()
    
    #var_imp <- matrix(NA, nrow = n_vars, ncol = 11) # 10 cols, last one is the mean
    var_imp <- as.data.frame(colnames(df_spectra[2:(n_vars + 1)]))
    colnames(var_imp)[1] <- "var_names"
    
    # Then cross validation inner loop
    for (cv_fold in 1:10){
      
      print(paste0("CV fold: ", cv_fold))
      #cv_fold <- 1
      
      #train_inds <- sample(all_row_inds, train_samp_size)
      val_inds <- val_samps[[cv_fold]] #all_row_inds[-train_inds]
      train_inds <- all_row_inds[-val_inds]
      
      df_spectra$genus <- as.factor(df_spectra$genus)
      train_spectra <- df_spectra[train_inds,] %>% select(-contains("Object_ID"))
      val_spectra <- df_spectra[val_inds,] %>% select(-contains("Object_ID"))
      
      set.seed(1) # here for ranger seed, not just sampling seed
      rf_trained <- ranger(genus ~ ., train_spectra, importance = 'permutation', num.threads = 10, num.trees = 500, seed = 1) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
      
      oob_out <- c(oob_out, rf_trained$prediction.error)
      
      val_predicted <- predict(rf_trained, data = val_spectra)
      conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)
      oa_out <- c(oa_out, conf_mat$overall[1])
      kappa_out <- c(kappa_out, conf_mat$overall[2])
      cilo_out <- c(cilo_out, conf_mat$overall[3])
      cihi_out <- c(cihi_out, conf_mat$overall[4])
      nvar_out <- c(nvar_out, n_vars)
      round_out <- c(round_out, round_iter-1) # this is looping incorrectly now
      # these are now all organized within a single cross val set
      
      # conf_mat
      conf_mat_list[[cv_fold]] <- conf_mat
      
      # variable importance
      rf_imp <- cbind.data.frame(names(rf_trained$variable.importance), rf_trained$variable.importance)
      colnames(rf_imp) <- c("var_names", paste0("imp_", cv_fold))
      var_imp <- merge(var_imp, rf_imp)
      
      rm(rf_trained)
    }
    
    gc()
    
    conf_mat_biglist[[round_iter - 1]] <- conf_mat_list
    # Need to keep average of accuracy variables for each iteration too
    
    # Then:
    # this is the dimension reduction part
    
    n_el_retain <- round(n_vars*0.8, 0) # keep 80% of the vars, rounded
    
    var_imp$imp_mean <- rowMeans(var_imp[,2:11])
    
    var_imp_list[[round_iter - 1]] <- var_imp
    
    imp_mean <- var_imp$imp_mean
    names(imp_mean) <- var_imp$var_names
    
    sorted_vars <- sort(imp_mean, decreasing = TRUE)
    
    kept_vars <- sorted_vars[1:n_el_retain]
    
    vars_rounds$max_round[which(vars_rounds$var %in% names(kept_vars))] <- round_iter
    round_iter <- round_iter + 1
    
    # will need to do with this with the original df_spectra because it's getting repeatedly resampled with the cross validation
    df_spectra <- df_spectra[, c("genus", names(kept_vars))]
    n_vars <- ncol(df_spectra) - 1
    
  }
  
  # combine data out
  data_out <- data.frame(round_out, nvar_out, oob_out, oa_out, kappa_out, cilo_out, cihi_out)
  colnames(data_out)[1] <- "Round" 
  
  # aggregate mean values for output information afterwards
  data_out_m <- aggregate(data_out, by = list(data_out$Round), FUN = "mean")
  data_out_m <- data_out_m[,2:ncol(data_out_m)]
  
  # which is the highest round (ie lowest number of variables) st the 95% overlaps with the maximum OA
  iter_round_selected <- max(data_out_m$Round[which(data_out_m$cihi_out > max(data_out_m$oa_out))])
  #data_out[iter_round_selected,]
  varsel_f <- vars_rounds$var[which(vars_rounds$max_round >= iter_round_selected)] # these are the selected variables.
  max_acc <- data_out_m$oa_out[which.max(data_out_m$oa_out)]
  
  ggplot(data_out_m) +
    geom_vline(xintercept = length(varsel_f), color = "navy", linetype = "dotted") +
    geom_hline(yintercept = max_acc, color = "red", linetype = "dotted") +
    geom_linerange(aes(nvar_out, ymin = cilo_out, ymax = cihi_out)) +
    geom_point(aes(nvar_out, oa_out)) +
    labs(x = "Number of features", y = "Overall Accuracy", title = feature_set_name) +
    theme_bw()
  
  # To save on output
  
  setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/")
  
  # plot
  ggsave(paste0(feature_set_name, "_feature_selection_cv.png"), width = 6, height = 5, units = "in")
  
  # confusion matrix list object, needed if we want to plot based on accuracy for individual genera later
  saveRDS(conf_mat_biglist, paste0(feature_set_name, "_conf_mat_biglist.rds"))
  
  # variable importance
  saveRDS(var_imp_list, paste0(feature_set_name, "_var_imp_list.rds"))
  
  # data_out and data_out_m
  write.csv(data_out, paste0(feature_set_name, "_feature_selection_data_full.csv"), row.names = FALSE)
  write.csv(data_out_m, paste0(feature_set_name, "_feature_selection_data_mean.csv"), row.names = FALSE)
  
  # iter_round_selected
  # var_rounds df
  # selected vars (just for ease of use)
  out_list <- list(iter_round_selected = iter_round_selected, 
                   vars_rounds = vars_rounds, 
                   varsel_f = varsel_f)
  saveRDS(out_list, paste0(feature_set_name, "_selected_features_info_cv.rds"))
  
  
}


# to get list again
#vars_rounds$var[which(vars_rounds$max_round >= iter_round_selected)]

#cfx <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/8week_collapsed_conf_mat_biglist.rds")
#sfx <- readRDS('/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/8week_collapsed_selected_features_info_cv.rds')

# get list
#sfx$vars_rounds$var[which(sfx$vars_rounds$max_round >= sfx$iter_round_selected)]

