
# Re-running selected main iterations of random forest ranger in order to get importance by genus

# Run particular "variable selection" random forest runs to get variable importance for individual genera

# Setting up 10 fold cv for ranger runs

# currently using 95% CI error bar overlaps for selection
# but could also use OOB error directly, and error on this, if this is the more appropriate thing to do

library(tidyverse)
library(ranger)
library(caret) # for confusion matrix and 10 fold cross val
`%notin%` <- Negate(`%in%`)
library(data.table)
library(tidytext)

for (input_file in 2:11){  
  #input_file <- 11
  
  # Aggregated (collapsed)
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
  
  # combined versions
  
  getSelVars <- function(file_path_rds){
    x <- readRDS(file_path_rds)
    sel_vars <- x$vars_rounds$var[which(x$vars_rounds$max_round >= x$iter_round_selected)]
    return(sel_vars)
  }
  
  # All collapsed
  if (input_file == 9){
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/1week_collapsed_selected_features_info_cv.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/2week_collapsed_selected_features_info_cv.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/4week_collapsed_selected_features_info_cv.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/8week_collapsed_selected_features_info_cv.rds")
    
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
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/1week_singleyears_selected_features_info_cv.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/2week_singleyears_selected_features_info_cv.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/4week_singleyears_selected_features_info_cv.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/8week_singleyears_selected_features_info_cv.rds")
    
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
  
  # just doing this one for now.
  
  # Everything combined!
  if (input_file == 11){
    # Collapsed
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/1week_collapsed_selected_features_info_cv.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/2week_collapsed_selected_features_info_cv.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/4week_collapsed_selected_features_info_cv.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/8week_collapsed_selected_features_info_cv.rds")
    
    setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
    df1wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v1wk)))
    df2wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_2weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v2wk)))
    df4wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_4weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v4wk)))
    df8wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_8weekcomplong_multiyearcollapsed_subset_wide.rds") %>% select(all_of(c("Object_ID", v8wk)))
    
    df_spectra_collapsed <- merge(df1wk, df2wk) %>% merge(df4wk) %>% merge(df8wk)
    
    # Single year
    v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/1week_singleyears_selected_features_info_cv.rds")
    v2wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/2week_singleyears_selected_features_info_cv.rds")
    v4wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/4week_singleyears_selected_features_info_cv.rds")
    v8wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/8week_singleyears_selected_features_info_cv.rds")
    
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
    v_collapsed <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/allweek_collapsed_selected_features_info_cv.rds")
    v_singleyear <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/allweek_singleyear_selected_features_info_cv.rds")
    
    df_spectra_collapsed <- df_spectra_collapsed %>% select(all_of(c("Object_ID", v_collapsed)))
    df_spectra_singleyear <- df_spectra_singleyear %>% select(all_of(c("Object_ID", v_singleyear)))
    
    df_spectra <- merge(df_spectra_collapsed, df_spectra_singleyear)
    
    feature_set_name <- "everything"
  }
  
  print(feature_set_name)
  
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
  
  #conf_mat_list <- list()
  oob_out <- c()
  oa_out <- c()
  kappa_out <- c()
  cilo_out <- c()
  cihi_out <- c()
  nvar_out <- c()
  #n_vars <- ncol(df_spectra) - 1 # remove genus part
  
  all_row_inds <- 1:nrow(df_spectra) # for training resample
  
  # New data analysis
  # load outputted information
  setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv")
  sel_feat_filename <- list.files(pattern = glob2rx(paste0(feature_set_name, "*_selected_features_info_cv.rds")))
  sel_feat <- readRDS(sel_feat_filename)
  sel_feat_itervec <- sel_feat$vars_rounds$var[which(sel_feat$vars_rounds$max_round >= sel_feat$iter_round_selected)]
  
  df_spectra <- df_spectra[, c("genus", sel_feat_itervec)] # subset df_spectra to these variables
  
  n_vars <- length(sel_feat_itervec)
  
  # variable selection
  # load the round_iter from the output file, this should be reduced by 1 because of some artifact of how this script was setup
  round_iter <- sel_feat$iter_round_selected + 1 #2 
  round_out <- c()
  
  conf_mat_biglist <- list()
  
  var_imp_list <- list()
  
  vars <- colnames(df_spectra)[2:ncol(df_spectra)] # these will be the same within a cross val set
  vars_rounds <- data.frame(vars, rep(1, length(vars)))
  colnames(vars_rounds) <- c("var", "max_round")
  
  #for(round_iter in 2:5){ # just a tester
  #while(n_vars > 5){
    
    print(paste0("Round iteration: ", round_iter - 1, ",  n_vars: ", n_vars))
    
    conf_mat_list <- list()
    
    local_varimp_list <- list()
    
    #var_imp <- matrix(NA, nrow = n_vars, ncol = 11) # 10 cols, last one is the mean
    var_imp <- as.data.frame(colnames(df_spectra[2:(n_vars + 1)]))
    colnames(var_imp)[1] <- "var_names"
    
    # Then cross validation inner loop
    for (cv_fold in 1:10){
      
      #cv_fold <- 1 # tester
      
      print(paste0("CV fold: ", cv_fold))
      
      #train_inds <- sample(all_row_inds, train_samp_size)
      val_inds <- val_samps[[cv_fold]] #all_row_inds[-train_inds]
      train_inds <- all_row_inds[-val_inds]
      
      df_spectra$genus <- as.factor(df_spectra$genus)
      train_spectra <- df_spectra[train_inds,] %>% select(-contains("Object_ID"))
      val_spectra <- df_spectra[val_inds,] %>% select(-contains("Object_ID"))
      
      set.seed(1) # here for ranger seed, not just sampling seed
      rf_trained <- ranger(genus ~ ., train_spectra, importance = 'permutation', local.importance = TRUE, num.threads = 10, num.trees = 500, seed = 1) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
      # added local.importance, get outputs of this
      
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
      
      # add local variable importance for each cv fold like for conf_mat, as list
      local_varimp <- cbind.data.frame(train_spectra$genus, rf_trained$variable.importance.local)
      colnames(local_varimp)[1] <- "genus"
      local_varimp_list[[cv_fold]] <- local_varimp
      
      # variable importance
      rf_imp <- cbind.data.frame(names(rf_trained$variable.importance), rf_trained$variable.importance)
      colnames(rf_imp) <- c("var_names", paste0("imp_", cv_fold))
      var_imp <- merge(var_imp, rf_imp)
      
      rm(rf_trained)
    }
    
    gc()
    
    #conf_mat_biglist[[1]] <- conf_mat_list
    # Need to keep average of accuracy variables for each iteration too
    
    # Then:
    # this is the dimension reduction part
    
    #n_el_retain <- round(n_vars*0.8, 0) # keep 80% of the vars, rounded
    
    var_imp$imp_mean <- rowMeans(var_imp[,2:11])
    
    var_imp_list[[1]] <- var_imp
    
    #imp_mean <- var_imp$imp_mean
    #names(imp_mean) <- var_imp$var_names
    
    #sorted_vars <- sort(imp_mean, decreasing = TRUE)
    
    #kept_vars <- sorted_vars[1:n_el_retain]
    
    #vars_rounds$max_round[which(vars_rounds$var %in% names(kept_vars))] <- round_iter
    #round_iter <- round_iter + 1
    
    # will need to do with this with the original df_spectra because it's getting repeatedly resampled with the cross validation
    #df_spectra <- df_spectra[, c("genus", names(kept_vars))]
    n_vars <- ncol(df_spectra) - 1
    
  #}
  
  # combine data out
  data_out <- data.frame(round_out, nvar_out, oob_out, oa_out, kappa_out, cilo_out, cihi_out)
  colnames(data_out)[1] <- "Round" 
  
  # aggregate mean values for output information afterwards
  data_out_m <- aggregate(data_out, by = list(data_out$Round), FUN = "mean")
  data_out_m <- data_out_m[,2:ncol(data_out_m)]
  
  # which is the highest round (ie lowest number of variables) st the 95% overlaps with the maximum OA
  iter_round_selected <- max(data_out_m$Round[which(data_out_m$cihi_out > max(data_out_m$oa_out))])
  #data_out[iter_round_selected,]
  #varsel_f <- vars_rounds$var[which(vars_rounds$max_round >= iter_round_selected)] # these are the selected variables.
  #max_acc <- data_out_m$oa_out[which.max(data_out_m$oa_out)]
  
  # won't need this
  # ggplot(data_out_m) +
  #   geom_vline(xintercept = length(varsel_f), color = "navy", linetype = "dotted") +
  #   geom_hline(yintercept = max_acc, color = "red", linetype = "dotted") +
  #   geom_linerange(aes(nvar_out, ymin = cilo_out, ymax = cihi_out)) +
  #   geom_point(aes(nvar_out, oa_out)) +
  #   labs(x = "Number of features", y = "Overall Accuracy", title = feature_set_name) +
  #   theme_bw()
  
  # To save on output
  
  setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/genus_importance/")
  
  # Will need to be a new suite of files coming out of here
  # plot
  #ggsave(paste0(feature_set_name, "_feature_selection_cv.png"), width = 6, height = 5, units = "in")
  
  # confusion matrix list object, needed if we want to plot based on accuracy for individual genera later
  saveRDS(conf_mat_list, paste0(feature_set_name, "_conf_mat_list_mainiter.rds"))
  # this will just be a single conf_mat to check that it matches the original runs (it should, based on seed)
  
  
  # variable importance
  saveRDS(var_imp_list, paste0(feature_set_name, "_var_imp_list_mainiter.rds"))
  # keep this as well, but will just be for a single iterations
  
  # local variable importance
  saveRDS(local_varimp_list, paste0(feature_set_name, "_local_var_imp_list_mainiter.rds"))
  
  # data_out and data_out_m
  write.csv(data_out, paste0(feature_set_name, "_feature_selection_data_full_mainiter.csv"), row.names = FALSE)
  write.csv(data_out_m, paste0(feature_set_name, "_feature_selection_data_mean_mainiter.csv"), row.names = FALSE)
  # the mean is for collapsing each of the 10-fold cv rounds
  
  # iter_round_selected
  # var_rounds df
  # selected vars (just for ease of use)
  # out_list <- list(iter_round_selected = iter_round_selected, 
  #                  vars_rounds = vars_rounds, 
  #                  varsel_f = varsel_f)
  # saveRDS(out_list, paste0(feature_set_name, "_selected_features_info_cv.rds"))
  
  # Local importance summary figure
  
  local_varimp_giantdf <- rbind.data.frame(local_varimp_list[[1]], local_varimp_list[[2]], local_varimp_list[[3]], local_varimp_list[[4]], local_varimp_list[[5]],
                                           local_varimp_list[[6]], local_varimp_list[[7]], local_varimp_list[[8]], local_varimp_list[[9]], local_varimp_list[[10]])
  local_varimp_giantdf_summary <- as.data.table(local_varimp_giantdf[,2:ncol(local_varimp_giantdf)])[,genus := local_varimp_giantdf$genus][,lapply(.SD,mean),by=genus]
  
  # get means of each test, so we can get SDs of these
  for (i in 1:10){
    local_varimp_test <- local_varimp_list[[i]]
    local_varimp_summary <- as.data.table(local_varimp_test[,2:ncol(local_varimp_test)])[,genus := local_varimp_test$genus][,lapply(.SD,mean),by=genus]
    if (i == 1){
      local_varimp_summary_mean_each <- local_varimp_summary
    } else {
      local_varimp_summary_mean_each <- rbind.data.frame(local_varimp_summary_mean_each, local_varimp_summary)
    }
  }
  
  local_varimp_giantdf_summary_SD <- local_varimp_summary_mean_each %>% group_by(genus) %>% summarize_all(sd)
  
  local_varimp_giantdf_summary_long <- local_varimp_giantdf_summary %>% as_tibble() %>% pivot_longer(cols = !genus, names_to = "time_range_feature", values_to = "value")
  
  local_varimp_giantdf_summary_SD_long <- local_varimp_giantdf_summary_SD %>% as_tibble() %>% pivot_longer(cols = !genus, names_to = "time_range_feature", values_to = "value_sd")
  
  
  topx <- local_varimp_giantdf_summary_long %>%
    group_by(genus) %>%
    top_n(10, value)
  
  topx_sd <- inner_join(topx, local_varimp_giantdf_summary_SD_long)
  
  p_imp <- ggplot(topx_sd) + 
    geom_point(aes(x = value, y = reorder_within(time_range_feature, value, genus)), size = 1) +
    geom_linerange(aes(xmin = value - value_sd, xmax = value + value_sd, y = reorder_within(time_range_feature, value, genus)), linewidth = 0.3) +
    facet_wrap(~genus, scales = "free", ncol = 3) +
    scale_y_reordered() +
    labs(x = "Mean Permutation Importance", y = "Most Important Features by Genus") +
    theme_bw() +
    theme(panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
          panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
          axis.title = element_text(size = 7, color = "black"),
          axis.text = element_text(size = 5, color = "black"),
          strip.text = element_text(size = 7, color = "black", face = "italic"),
          strip.background = element_rect(fill = "white", color = NA))
  ggsave(paste0(feature_set_name, "_local_imp.jpg"), p_imp, width = 190, height = 150, units = "mm")
  
}