library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(caret)
library(xgboost)
library(ggrepel)


#####

out_path <- "/Volumes/NYC_geo/tree_classification/year_test_outputs/"
name_prefix <- "test1_" # if any name modifications are desired appended to beginning of each output, make to sure include trailing "_"

tree_df_for_cls_clean <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file.rds") # not too bad for loading in, < 30 seconds on SSD

# looping for different years

col_list <- colnames(tree_df_for_cls_clean)

for (yr in 2018:2024){
  #yr <- 2024
  col_list_sub <- col_list[grep(paste0("_", yr), col_list)]
  
  kept_cols <- c("Poly_ID", "genus", "npts", "height_max_ft", "height_med_ft", "int_mean", "int_mean_above_medh", "int_mean_below_medh", "wid_medh_max_ft", 
                 "wid_medh_circlerad_ft", "hw_rat_2_max", "hw_rat_2_circlerad", "cp_3_max", "cp_3_circlerad")
  
  cols_to_use <- c(kept_cols, col_list_sub)
  
  tree_df_for_cls_clean_sub <- tree_df_for_cls_clean[,..cols_to_use] #.. for tibble
  
  #
  
  # Do an xgboost run
  set.seed(14)
  
  colnames(tree_df_for_cls_clean_sub)[1] <- "Object_ID"
  
  df_cls_agg_mean_cc <- tree_df_for_cls_clean_sub[sample(1:nrow(tree_df_for_cls_clean_sub), 2000, 0),] # try a sample of 2000 trees for a test run
  #df_cls_agg_mean_cc <- tree_df_for_cls_clean_sub
  
  genus_list <- sort(unique(df_cls_agg_mean_cc$genus))
  
  df_cls_agg_mean_cc$genus <- as.numeric(as.factor(as.character(df_cls_agg_mean_cc$genus))) - 1 # this scrambles the order of the genera unless df_cls_agg_mean_cc is redefined
  train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
  train_data <- df_cls_agg_mean_cc[train_index, ]
  test_data <- df_cls_agg_mean_cc[-train_index, ]
  
  train_data_noid <- train_data %>% select(!Object_ID)
  test_data_noid <- test_data %>% select(!Object_ID)
  
  # https://rpubs.com/mharris/multiclass_xgboost
  # might need to normalize classes
  
  train_data_input <- train_data_noid %>% select(!genus) %>% as.matrix() 
  train_matrix <- xgb.DMatrix(data = train_data_input, label = train_data_noid$genus)
  
  test_data_input <- test_data_noid %>% select(!genus) %>% as.matrix() 
  test_matrix <- xgb.DMatrix(data = test_data_input, label = test_data_noid$genus)
  
  #
  numberOfClasses <- length(unique(train_data_noid$genus))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "nthread" = 8)
  nround    <- 999999 # number of XGBoost rounds, was originally 50 but setting to giant number to see if it does in fact stop
  esr <- 10 # early stopping rounds, need to figure out how to use this
  cv.nfold  <- 5
  
  # this is for cross validation folding (not needed right now)
  # Fit cv.nfold * cv.nround XGB models and save OOF predictions
  # try with verbose?
  # start_time <- Sys.time()
  # cv_model <- xgb.cv(params = xgb_params,
  #                    data = train_matrix,
  #                    nrounds = nround,
  #                    nfold = cv.nfold,
  #                    verbose = TRUE,
  #                    prediction = TRUE)
  # t_dif <- Sys.time() - start_time
  # print(t_dif) # 8.4 hours to run...
  # 
  # OOF_prediction <- data.frame(cv_model$pred) %>%
  #   mutate(max_prob = max.col(., ties.method = "last"),
  #          label = train_data_noid$genus + 1)
  # head(OOF_prediction)
  # # confusion matrix
  # confusionMatrix(factor(OOF_prediction$max_prob),
  #                 factor(OOF_prediction$label),
  #                 mode = "everything")
  # 
  # xgb_impcv <- xgb.importance(model = cv_model)
  # 
  # ggplot(xgb_impcv[1:25,]) +
  #   geom_col(aes(x = Gain, y = fct_reorder(Feature, Gain))) +
  #   labs(x = "Gain", y = "Variable", title = "Top 25 most important variables")
  
  
  wl <- list(train = train_matrix, eval = test_matrix) # sometimes crashes and R aborts with a watchlist, didn't happen before
  start_time <- Sys.time()
  bst_model <- xgb.train(params = xgb_params, data = train_matrix, nrounds = nround, watchlist = wl, verbose = TRUE, early_stopping_rounds = esr) # added verbose setting and early_stopping_rounds = 10
  t_dif <- Sys.time() - start_time
  print(t_dif)
  
  # output xgb model
  setwd(out_path)
  xgb_outname <- paste0(name_prefix, "xgb_", yr, "_only.model")
  xgb.save(bst_model, xgb_outname)
  #test_load <- xgb.load(xgb_outname)
  
  # Predict hold-out test set
  test_pred <- predict(bst_model, newdata = test_matrix)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses) %>%
    t() %>%
    data.frame() %>%
    mutate(label = test_data_noid$genus + 1,
           max_prob = max.col(., "last"))
  # confusion matrix of test set
  # confusionMatrix(factor(test_prediction$max_prob),
  #                 factor(test_prediction$label),
  #                 mode = "everything")
  # so xgboost is not a panacea, but does work better, especially for some small classes. 
  # 50 rounds takes 8x longer to run than ranger with 500 rf trees though
  
  xgb_imp <- xgb.importance(model = bst_model)
  xgb_imp_df <- paste0(name_prefix, "xgb_", yr, "_importance.csv")
  write.csv(xgb_imp, xgb_imp_df, row.names = FALSE)
  
  ggplot(xgb_imp[1:25,]) +
    geom_col(aes(x = Gain, y = fct_reorder(Feature, Gain))) +
    labs(x = "Gain", y = "Variable", title = "Top 25 most important variables")
  xgb_imp_fig <- paste0(name_prefix, "xgb_", yr, "_importance.png")
  ggsave(xgb_imp_fig, height = 5, width = 8, units = "in")
  
  #
  xgb_cm <- confusionMatrix(factor(test_prediction$max_prob),
                            factor(test_prediction$label),
                            mode = "everything")
  xgb_cm_out <- paste0(name_prefix, "xgb_", yr, "_confmat.rds")
  saveRDS(xgb_cm, xgb_cm_out)
  
  # save test prediction for referencing
  test_prediction_path <- paste0(name_prefix, "xgb_", yr, "_predictions.csv")
  write.csv(cbind.data.frame(test_data[,1:2], test_prediction), test_prediction_path, row.names = FALSE)
  
  table(test_data$genus, test_prediction$max_prob)
  
  prod_acc <- unname(xgb_cm$byClass[,1])
  user_acc <- unname(xgb_cm$byClass[,3])
  class_num_xgb <- sort(unique(test_data$genus))
  classes <- genus_list[class_num_xgb + 1] # add 1 because xgb numbering for classes starts at 0
  cm_results <- cbind.data.frame(classes, class_num_xgb, prod_acc, user_acc)
  
  test_prediction_accuracies <- paste0(name_prefix, "xgb_", yr, "_accuracies.csv")
  write.csv(cm_results, test_prediction_accuracies, row.names = FALSE)
  
  cm_results$classes <- as.factor(cm_results$classes)
  
  ggplot(cm_results) +
    geom_point(aes(x = prod_acc, y = user_acc)) +
    geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
    coord_equal() +
    ylim(0,1) + xlim(0,1) +
    labs(x = "Producer's Accuracy", y = "User's Accuracy") +
    theme_bw()
  xgb_acc_fig <- paste0(name_prefix, "xgb_", yr, "_acc.png")
  ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")
}

