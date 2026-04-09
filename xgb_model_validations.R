library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(caret)
library(xgboost)
library(sf)
library(ggrepel)
library(wconf)

# Validation datasets

# Check test dataset

#out_path <- "E:/nyc/classification/outputs_practical_v2_monthcomp"
out_path <- "/Volumes/NYC_geo/tree_classification/outputs_practical_v2_monthcomp"
name_prefix <- "practical_v2_monthcomp_"

#tree_df_for_cls_clean <- readRDS("E:/nyc/classification/cls_prepped_file_v2_monthcomp.rds")
tree_df_for_cls_clean <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file_v2_monthcomp.rds")

#####
# Evaluate the xgboost run
set.seed(14)

colnames(tree_df_for_cls_clean)[1] <- "Object_ID"

df_cls_agg_mean_cc <- tree_df_for_cls_clean # full version

genus_list <- sort(unique(df_cls_agg_mean_cc$genus))

rm(tree_df_for_cls_clean)

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

# #
numberOfClasses <- length(unique(train_data_noid$genus))
# xgb_params <- list("objective" = "multi:softprob",
#                    "eval_metric" = "mlogloss",
#                    "num_class" = numberOfClasses,
#                    "nthread" = 8)
# nround    <- 999999 # number of XGBoost rounds, was originally 50 but setting to giant number to see if it does in fact stop
# esr <- 10 # early stopping rounds, need to figure out how to use this
# cv.nfold  <- 5

#wl <- list(train = train_matrix, eval = test_matrix) # sometimes crashes and R aborts with a watchlist, didn't happen before
#start_time <- Sys.time()
#bst_model <- xgb.train(params = xgb_params, data = train_matrix, nrounds = nround, watchlist = wl, verbose = TRUE, early_stopping_rounds = esr) # added verbose setting and early_stopping_rounds = 10
#t_dif <- Sys.time() - start_time
#print(t_dif)

# Save with xgb.save()
# Reload with xgb.load()

bst_model <- xgb.load("/Volumes/NYC_geo/tree_classification/outputs_practical_v2_monthcomp/practical_v2_monthcomp_xgb_full.model")

# output xgb model
#setwd(out_path)
#xgb_outname <- paste0(name_prefix, "xgb_full.model")
#xgb.save(bst_model, xgb_outname)
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
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")
# this is the same as the outputted validation dataset from earlier


# xgb_imp <- xgb.importance(model = bst_model)
# xgb_imp_df <- paste0(name_prefix, "xgb_full_importance.csv")
# write.csv(xgb_imp, xgb_imp_df, row.names = FALSE)

# # Additional code for if there are less than 25 variables
# nvar <- nrow(xgb_imp)
# if(nvar > 25){
#   var_plotted <- 25
# } else {
#   var_plotted <- nvar
# }
# 
# ggplot(xgb_imp[1:var_plotted,]) +
#   geom_col(aes(x = Gain, y = fct_reorder(Feature, Gain))) +
#   labs(x = "Gain", y = "Variable", title = paste0("Top ", var_plotted, " most important variables"))
# xgb_imp_fig <- paste0(name_prefix, "xgb_full.jpg")
# ggsave(xgb_imp_fig, height = 5, width = 8, units = "in")

setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v2_monthcomp/validations")

eval_type <- "hq"

xgb_cm <- confusionMatrix(factor(test_prediction$max_prob),
                          factor(test_prediction$label),
                          mode = "everything")
xgb_cm_out <- paste0(name_prefix, "xgb_full_confmat_", eval_type, ".rds")
saveRDS(xgb_cm, xgb_cm_out)

# save test prediction for referencing
test_prediction_path <- paste0(name_prefix, "xgb_full_predictions_", eval_type, ".csv")
write.csv(cbind.data.frame(test_data[,1:2], test_prediction), test_prediction_path, row.names = FALSE)

prod_acc <- unname(xgb_cm$byClass[,1])
user_acc <- unname(xgb_cm$byClass[,3])
class_num_xgb <- sort(unique(test_data$genus))
classes <- genus_list[class_num_xgb + 1] # add 1 because xgb numbering for classes starts at 0
cm_results <- cbind.data.frame(classes, class_num_xgb, prod_acc, user_acc)

test_prediction_accuracies <- paste0(name_prefix, "xgb_full_accuracies_", eval_type, ".csv")
write.csv(cm_results, test_prediction_accuracies, row.names = FALSE)

cm_results$classes <- as.factor(cm_results$classes)

ggplot(cm_results) +
  geom_point(aes(x = prod_acc, y = user_acc)) +
  geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  coord_equal() +
  ylim(0,1) + xlim(0,1) +
  labs(x = "Producer's Accuracy", y = "User's Accuracy") +
  theme_bw()
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_", eval_type, "_nstems_wtext.eps")
ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_", eval_type, "_nstems_wtext.jpg")
ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")

ggplot(cm_results) +
  geom_point(aes(x = prod_acc, y = user_acc), shape = 1, size = 1.5, color = "gray10") +
  #geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  #geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  coord_equal() +
  #ylim(0,1) + xlim(0,1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0.1, 0.9, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0.1, 0.9, 0.2), limits = c(0, 1)) +
  labs(x = "Producer's Accuracy", y = "User's Accuracy") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))
#xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_", eval_type,"_nstems_notext.eps")
#xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_", eval_type,"_nstems_notext.jpg")
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_", eval_type,"_nstems_notext.svg")
#xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_", eval_type,"_nstems_wtext.jpg")
ggsave(xgb_acc_fig, height = 4, width = 4, units = "in")
#ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")

# might need to update the settings for text and output sizing

# test_data$Object_ID is the validation data that should be included in the balanced dataset

# Next, evaluate this same set but by canopy area
cls <- st_read("/Volumes/NYC_geo/tree_classification/outputs_practical_v2_monthcomp/nyc_class_tree_genus_polygons_v2.gpkg")

cls_test <- cls[1:5,]
cls_test$TNC_SHAPE_Area
st_area(cls_test)

# confusion matrix weighted by area
# xgb_cm <- confusionMatrix(factor(test_prediction$max_prob),
#                           factor(test_prediction$label),
#                           mode = "everything")
area_weights <- cls %>% filter(Poly_ID %in% test_data$Object_ID) %>% select(Poly_ID, TNC_SHAPE_Area, Genus_Ref, Genus_Predicted)
#all(area_weights$Poly_ID, test_data$Object_ID) # TRUE, same order
xgb_wcm <- wconfusionmatrix(xgb_cm, weight.type = "custom", custom.weights = area_weights$TNC_SHAPE_Area)

# Area weighted accuracies got worse... classifier struggling in big, mixed up crowns? Will need to think about this more

# Pretty sure what I did last time was weight the final accuracies by the mean or sum of the areas, not the accuracy within class per tree
# Let's see what that accuracy would be
area_by_genus_ref <- area_weights %>% 
  st_drop_geometry() %>%
  group_by(Genus_Ref) %>% 
  summarize(area_sum = sum(TNC_SHAPE_Area)) %>%
  mutate(area_frac = area_sum/sum(area_sum))
# these are the fractions but need to think about how to include them to upweight or downweight this

# Just going to stick to stems for now, then circle back to the area weighting

#####

# Filter based on minimum summer NDVI threshold for all years
setwd('/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract')
mean_ndvi_all <- fread("mean_summer_ndvi.csv")
ndvi_min <- 0.3 # Thapa

# drop 2017, only need them to be green 2018-2024 for classification, 2017 is weird? Or did we use 2017 for the classfication?
mean_ndvi_all_live <- mean_ndvi_all %>% select(!ndvi_2017) %>% filter_at(vars(starts_with("ndvi")), all_vars(. > ndvi_min))

cls_sub <- cls %>% filter(height_max_ft > 10 & Poly_ID %in% mean_ndvi_all_live$Object_ID &
                            !is.na(Genus_Ref) & !is.na(Genus_Predicted))
xgb_cm2 <- confusionMatrix(factor(cls_sub$Genus_Predicted),
                          factor(cls_sub$Genus_Ref),
                          mode = "everything")
# OA = 0.8255, accuracy is higher maybe because lots of training data is in here?

cls_sub2 <- cls %>% filter(!is.na(Genus_Ref) & !is.na(Genus_Predicted))
xgb_cm3 <- confusionMatrix(factor(cls_sub2$Genus_Predicted),
                           factor(cls_sub2$Genus_Ref),
                           mode = "everything")
# OA = 0.7293, this includes little trees and trees with low NDVI values

#
cls_sub3 <- cls %>% filter(height_max_ft > 10 & Poly_ID %in% mean_ndvi_all_live$Object_ID & Poly_ID %notin% train_data$Object_ID &
                            !is.na(Genus_Ref) & !is.na(Genus_Predicted))
xgb_cm4 <- confusionMatrix(factor(cls_sub3$Genus_Predicted),
                           factor(cls_sub3$Genus_Ref),
                           mode = "everything")
xgb_cm4
# OA = 0.5501, but this distribution is not correct for some reason because the train_data counts should be in here

cls_sub_test <- cls %>% filter(height_max_ft > 10 & Poly_ID %in% mean_ndvi_all_live$Object_ID & Poly_ID %in% test_data$Object_ID &
                             !is.na(Genus_Ref) & !is.na(Genus_Predicted))

cls_sub_train <- cls %>% filter(height_max_ft > 10 & Poly_ID %in% mean_ndvi_all_live$Object_ID & Poly_ID %in% train_data$Object_ID &
                                              !is.na(Genus_Ref) & !is.na(Genus_Predicted))

cls_sub_notin_test_or_train <- cls %>% filter(height_max_ft > 10 & Poly_ID %in% mean_ndvi_all_live$Object_ID & Poly_ID %notin% test_data$Object_ID & Poly_ID %notin% train_data$Object_ID &
                                 !is.na(Genus_Ref) & !is.na(Genus_Predicted))

# Test should be repeated to be 5x copies of itself, then appended to notin_test_or_train
# This accounts for the missing length of training data, without the training data itself being included

cls_sub_new_comb <- bind_rows(cls_sub_notin_test_or_train, cls_sub_test, cls_sub_test, cls_sub_test, cls_sub_test, cls_sub_test)
# very close to length of cls_sub, so this is basically accounting for not included training data
xgb_cm_new_comb <- confusionMatrix(factor(cls_sub_new_comb$Genus_Predicted),
                           factor(cls_sub_new_comb$Genus_Ref),
                           mode = "everything")
xgb_cm_new_comb
# OA = 0.7211, and this still includes genera that the model does not know about and cannot possibly get correct
# I feel like this is reasonable

prod_acc_new_comb <- unname(xgb_cm_new_comb$byClass[,1])
user_acc_new_comb <- unname(xgb_cm_new_comb$byClass[,3])
classes_xgb_new_comb <- sort(unique(cls_sub_new_comb$Genus_Ref))
cm_results_new_comb <- cbind.data.frame(classes_xgb_new_comb, prod_acc_new_comb, user_acc_new_comb)

test_prediction_accuracies_new_comb <- paste0(name_prefix, "xgb_full_accuracies_new_comb.csv")
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v2_monthcomp/validations")
write.csv(cm_results_new_comb, test_prediction_accuracies_new_comb, row.names = FALSE)

xgb_cm_new_comb$table[2,2]

for(gen in genus_list){
  print(paste0(gen, ": ", xgb_cm_new_comb$table[gen, gen]))
}

#

cm_results_new_comb2 <- cm_results_new_comb %>% filter(classes_xgb_new_comb %in% genus_list)
cm_results_new_comb2$classes <- as.factor(cm_results_new_comb2$classes_xgb_new_comb)

ggplot(cm_results_new_comb2) +
  geom_point(aes(x = prod_acc_new_comb, y = user_acc_new_comb)) +
  geom_text_repel(aes(x = prod_acc_new_comb, y = user_acc_new_comb, label = classes)) +
  coord_equal() +
  ylim(0,1) + xlim(0,1) +
  labs(x = "Producer's Accuracy", y = "User's Accuracy") +
  theme_bw()
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_newval_nstems_wtext.eps")
ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_newval_nstems_wtext.jpg")
ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")

ggplot(cm_results_new_comb2) +
  geom_point(aes(x = prod_acc_new_comb, y = user_acc_new_comb), shape = 1, size = 1.5, color = "gray10") +
  #geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  coord_equal() +
  #ylim(0,1) + xlim(0,1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0.1, 0.9, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0.1, 0.9, 0.2), limits = c(0, 1)) +
  labs(x = "Producer's Accuracy", y = "User's Accuracy") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_newval_nstems_notext.svg")
ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")

ggplot(cm_results_new_comb2) +
  geom_point(aes(x = prod_acc_new_comb, y = user_acc_new_comb), shape = 1, size = 1.5, color = "gray10") +
  #geom_text_repel(aes(x = prod_acc_new_comb, y = user_acc_new_comb, label = classes)) +
  coord_equal() +
  #ylim(0,1) + xlim(0,1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0.1, 0.9, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0.1, 0.9, 0.2), limits = c(0, 1)) +
  labs(x = "Producer's Accuracy", y = "User's Accuracy") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_newval_nstems_notext.svg")
ggsave(xgb_acc_fig, height = 4, width = 4, units = "in")
# xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc_newval_nstems_wtext.jpg")
# ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")
