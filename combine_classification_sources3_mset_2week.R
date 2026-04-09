library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(caret)
#library(ranger)
#library(missRanger)
#library(ranger)
#library(tuneRanger)
library(xgboost)
#library(ggrepel)

# Tree genera labels and SOS values
tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

# good tree ids have already been pre-filtered in previous scripts

# Aggregated spectra
# 4 band
setwd("/Volumes/NYC_geo/tree_classification/extracted_4band_2week/zscaled/")
file_list <- list.files()
tree_spectra_4band <- purrr::map_df(file_list, fread)
colnames(tree_spectra_4band)[2:ncol(tree_spectra_4band)] <- paste0("4b_", colnames(tree_spectra_4band)[2:ncol(tree_spectra_4band)])

tree_spectra_4band %>% na.omit() %>% nrow()

# # figuring out which are na, for example
# which(is.na(tree_spectra_4band$`4b_blue_zscaled_mean_20180312`)) %>% head()
# #239 240 660 661 665 669
# tree_spectra_4band$Object_ID[239] # 10343

# 8 band
setwd("/Volumes/NYC_geo/tree_classification/extracted_8band_2week/zscaled/")
file_list <- list.files()
tree_spectra_8band <- purrr::map_df(file_list, fread)
colnames(tree_spectra_8band)[2:ncol(tree_spectra_8band)] <- paste0("8b_", colnames(tree_spectra_8band)[2:ncol(tree_spectra_8band)])

tree_spectra_8band %>% na.omit() %>% nrow()

# Combine and remove individual versions
tree_spectra <- merge(tree_spectra_4band, tree_spectra_8band, by = c("Object_ID"))
rm(tree_spectra_4band)
rm(tree_spectra_8band)

# add genus labels to tree spectra
tree_labels <- tree_pheno %>% filter(Poly_ID %in% tree_spectra$Object_ID) %>% select("Poly_ID", "genus")
tree_labels <- tree_labels[match(tree_spectra$Object_ID, tree_labels$Poly_ID),]
tree_spectra_wlabels <- merge(tree_labels, tree_spectra, by.x = "Poly_ID", by.y = "Object_ID")
rm(tree_spectra)

tree_spectra_wlabels %>% na.omit() %>% nrow()

# Add in SOS variables
pheno_vars <- tree_pheno %>% filter(Poly_ID %in% tree_spectra_wlabels$Poly_ID) %>% select("Poly_ID", "Year", "SOS_50", "SOS_20", "SOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)

tree_spectra_wlabels_sos <- merge(tree_spectra_wlabels, pheno_vars_wide)
rm(pheno_vars)
rm(pheno_vars_long)
rm(tree_spectra_wlabels)
gc()

#tree_spectra_wlabels_sos %>% na.omit() %>% nrow()

# add in variables from our lidar crown metrics
setwd("/Volumes/NYC_geo/nyc_lidar_metrics/crown_metrics")
crown_metrics_file_list <- list.files()
tree_crown_metrics <- purrr::map_df(crown_metrics_file_list, fread)
tree_spectra_wlabels_sos_lcm <- merge(tree_spectra_wlabels_sos, tree_crown_metrics) # lcm = lidar crown metrics

#tree_spectra_wlabels_sos_lcm %>% na.omit %>% nrow()
tree_df_for_cls <- na.omit(tree_spectra_wlabels_sos_lcm) # losing 3 trees that have random bad data values in the zscaled band combinations, will remove from analysis.

tree_cls_mat <- as.matrix(tree_df_for_cls[,3:ncol(tree_df_for_cls)])
summed_rows <- rowSums(tree_cls_mat) # much faster on matrix
tree_df_for_cls_clean <- tree_df_for_cls[is.finite(summed_rows), ] # removing any more trees with missing data. 16 total removed from dataset because of Infinite or NaN band combination values
# the Object_ID on this should be the final set of trees we use, 225027 trees
#table(tree_df_for_cls_clean$genus)

# remove unneeded variables to converve memory
rm(tree_crown_metrics)
rm(tree_cls_mat)
rm(tree_labels)
rm(tree_spectra_wlabels_sos)
rm(tree_spectra_wlabels_sos_lcm)
rm(tree_pheno)
rm(pheno_vars_wide)
rm(tree_df_for_cls)
gc()

saveRDS(tree_df_for_cls_clean, "/Volumes/NYC_geo/tree_classification/cls_prepped_file.rds")

# #####
# # Test RF
# 
# # Skipping RF run for now, just doing XGBoost
# # RF Run
# 
# # remove NA
# #df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno_filt) # do this cleaning for prep, could fold this into previous dplyr exclusion fn
# 
# set.seed(14)
# train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
# train_data <- df_cls_agg_mean_cc[train_index, ]
# test_data <- df_cls_agg_mean_cc[-train_index, ]
# 
# train_data_noid <- train_data %>% select(!Object_ID) # %>% as.data.frame()
# test_data_noid <- test_data %>% select(!Object_ID) # %>% as.data.frame()
# 
# # might want to do tuneRanger or similar to adjust the hyperparameters used for the ranger RF
# # hit a fatal error in R when tried to do this last time
# # treecls_task <- makeClassifTask(data = train_data_noid, target = "genus")
# # estimateTimeTuneRanger(treecls_task)
# # # Tuning
# # res = tuneRanger(iris.task, measure = list(multiclass.brier), num.trees = 1000, 
# #                  num.threads = 2, iters = 70, save.file.path = NULL)
# 
# # Also need to consider rebalancing class levels using SMOTE or some other rebalancing algorithm
# # Should also measure runtime for ranger model run
# 
# # <20 minutes with 500 trees
# rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity", num.trees = 500, num.threads = 10)
# #rf_model <- ranger(genus ~ ., train_data_noid, importance = "permutation", local.importance = TRUE) # Doing local variable importance, trying to see what is driving each tree class
# p_rf <- predict(rf_model, test_data_noid[,2:ncol(test_data_noid)])
# confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)
# 
# var_imp <- as.data.frame(rf_model$variable.importance)
# var_imp$name <- names(rf_model$variable.importance)
# colnames(var_imp)[1] <- "importance"
# 
# var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]
# 
# ggplot(var_imp_ordered[1:25,]) +
#   geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
#   labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")
# 
# # This works, but very different for order of importance likely because this is not really a representative sample
# 
# # Make plot from confusion matrix for Sensitivity (Producer's Acc) vs. Pos Pred Value (User's Acc)
# rf_cm <- confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)
# 
# prod_acc <- unname(rf_cm$byClass[,1])
# user_acc <- unname(rf_cm$byClass[,3])
# classes <- rownames(rf_cm$table)
# cm_results <- cbind.data.frame(classes, prod_acc, user_acc)
# cm_results$classes <- as.factor(cm_results$classes)
# 
# ggplot(cm_results) +
#   geom_point(aes(x = prod_acc, y = user_acc)) +
#   geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
#   coord_equal() +
#   theme_bw()

#####

tree_df_for_cls_clean <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file.rds") # not too bad for loading in, < 30 seconds on SSD
# Do an xgboost run
set.seed(14)

colnames(tree_df_for_cls_clean)[1] <- "Object_ID"

df_cls_agg_mean_cc <- tree_df_for_cls_clean[sample(1:nrow(tree_df_for_cls_clean), 2000, 0),] # try a sample of 2000 trees for a test run

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
#bst_model <- xgb.train(params = xgb_params, data = train_matrix, nrounds = nround)
bst_model <- xgb.train(params = xgb_params, data = train_matrix, nrounds = nround, watchlist = wl, verbose = TRUE, early_stopping_rounds = esr) # added verbose setting and early_stopping_rounds = 10
t_dif <- Sys.time() - start_time
print(t_dif)

# Predict hold-out test set
test_pred <- predict(bst_model, newdata = test_matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test_data_noid$genus + 1,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")
# so xgboost is not a panacea, but does work better, especially for some small classes. 
# 50 rounds takes 8x longer to run than ranger with 500 rf trees though

xgb_imp <- xgb.importance(model = bst_model)

ggplot(xgb_imp[1:25,]) +
  geom_col(aes(x = Gain, y = fct_reorder(Feature, Gain))) +
  labs(x = "Gain", y = "Variable", title = "Top 25 most important variables")

#
xgb_cm <- confusionMatrix(factor(test_prediction$max_prob),
                          factor(test_prediction$label),
                          mode = "everything")

prod_acc <- unname(xgb_cm$byClass[,1])
user_acc <- unname(xgb_cm$byClass[,3])
classes <- rownames(rf_cm$table)
cm_results <- cbind.data.frame(classes, prod_acc, user_acc)
cm_results$classes <- as.factor(cm_results$classes)

ggplot(cm_results) +
  geom_point(aes(x = prod_acc, y = user_acc)) +
  geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  coord_equal() +
  ylim(0,1) + xlim(0,1) +
  theme_bw()


col_list <- colnames(tree_df_for_cls_clean)

yr <- 2024
col_list_sub <- col_list[grep(paste0("_", yr), col_list)]

kept_cols <- c("Poly_ID", "genus", "npts", "height_max_ft", "height_med_ft", "int_mean", "int_mean_above_medh", "int_mean_below_medh", "wid_medh_max_ft", 
  "wid_medh_circlerad_ft", "hw_rat_2_max", "hw_rat_2_circlerad", "cp_3_max", "cp_3_circlerad")

cols_to_use <- c(kept_cols, col_list_sub)

tree_df_for_cls_clean_sub <- tree_df_for_cls_clean[,..cols_to_use] #.. for tibble

# feed this into model