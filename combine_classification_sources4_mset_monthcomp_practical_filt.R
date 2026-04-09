library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(caret)
library(xgboost)
library(sf)

# good_trees_prep <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file.rds")
# good_trees_prep_ids <- good_trees_prep %>% select(Poly_ID, genus)
# rm(good_trees_prep)
# gc()
# 
# # Aggregated spectra
# # 4 band
# setwd("/Volumes/NYC_geo/tree_classification/extracted_4band_monthcomp/zscaled/")
# file_list <- list.files()
# tree_spectra_4band <- purrr::map_df(file_list, fread)
# colnames(tree_spectra_4band)[2:ncol(tree_spectra_4band)] <- paste0("4b_", colnames(tree_spectra_4band)[2:ncol(tree_spectra_4band)])
# tree_spectra_4band_sub <- tree_spectra_4band %>% filter(Object_ID %in% good_trees_prep_ids$Poly_ID)
# 
# tree_spectra_4band_sub %>% na.omit() %>% nrow()
# 
# rm(tree_spectra_4band)
# gc()
# 
# # 8 band
# setwd("/Volumes/NYC_geo/tree_classification/extracted_8band_monthcomp/zscaled/")
# file_list <- list.files()
# tree_spectra_8band <- purrr::map_df(file_list, fread)
# colnames(tree_spectra_8band)[2:ncol(tree_spectra_8band)] <- paste0("8b_", colnames(tree_spectra_8band)[2:ncol(tree_spectra_8band)])
# tree_spectra_8band_sub <- tree_spectra_8band %>% filter(Object_ID %in% good_trees_prep_ids$Poly_ID)
# 
# tree_spectra_8band_sub %>% na.omit() %>% nrow()
# 
# rm(tree_spectra_8band)
# gc()
# 
# # Combine and remove individual versions
# tree_spectra <- merge(tree_spectra_4band_sub, tree_spectra_8band_sub, by = c("Object_ID"))
# rm(tree_spectra_4band_sub)
# rm(tree_spectra_8band_sub)
# 
# # add genus labels to tree spectra
# #tree_labels <- tree_pheno %>% filter(Poly_ID %in% tree_spectra$Object_ID) %>% select("Poly_ID", "genus")
# #tree_labels <- tree_labels[match(tree_spectra$Object_ID, good_trees_prep_ids$Poly_ID),]
# tree_spectra_wlabels <- merge(good_trees_prep_ids, tree_spectra, by.x = "Poly_ID", by.y = "Object_ID")
# rm(tree_spectra)
# 
# tree_spectra_wlabels %>% na.omit() %>% nrow()
# 
# # # Add in SOS variables
# # pheno_vars <- tree_pheno %>% filter(Poly_ID %in% tree_spectra_wlabels$Poly_ID) %>% select("Poly_ID", "Year", "SOS_50", "SOS_20", "SOS_80")
# # pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
# # pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
# # pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
# # pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)
# # 
# # tree_spectra_wlabels_sos <- merge(tree_spectra_wlabels, pheno_vars_wide)
# # rm(pheno_vars)
# # rm(pheno_vars_long)
# # rm(tree_spectra_wlabels)
# # gc()
# 
# #tree_spectra_wlabels_sos %>% na.omit() %>% nrow()
# 
# # add in variables from our lidar crown metrics
# setwd("/Volumes/NYC_geo/nyc_lidar_metrics/crown_metrics")
# crown_metrics_file_list <- list.files()
# tree_crown_metrics <- purrr::map_df(crown_metrics_file_list, fread)
# #tree_spectra_wlabels_sos_lcm <- merge(tree_spectra_wlabels_sos, tree_crown_metrics) # lcm = lidar crown metrics
# tree_spectra_wlabels_lcm <- merge(tree_spectra_wlabels, tree_crown_metrics) # lcm = lidar crown metrics
# 
# 
# #tree_spectra_wlabels_sos_lcm %>% na.omit %>% nrow()
# #tree_df_for_cls <- na.omit(tree_spectra_wlabels_sos_lcm) # losing 3 trees that have random bad data values in the zscaled band combinations, will remove from analysis.
# tree_df_for_cls <- na.omit(tree_spectra_wlabels_lcm)
# 
# tree_cls_mat <- as.matrix(tree_df_for_cls[,3:ncol(tree_df_for_cls)])
# summed_rows <- rowSums(tree_cls_mat) # much faster on matrix
# tree_df_for_cls_clean <- tree_df_for_cls[is.finite(summed_rows), ] # removing any more trees with missing data.
# 
# # remove unneeded variables to converve memory
# rm(tree_crown_metrics)
# rm(tree_cls_mat)
# #rm(tree_labels)
# #rm(tree_spectra_wlabels_sos)
# #rm(tree_spectra_wlabels_sos_lcm)
# #rm(tree_pheno)
# #rm(pheno_vars_wide)
# rm(tree_df_for_cls)
# gc()
# 
# saveRDS(tree_df_for_cls_clean, "/Volumes/NYC_geo/tree_classification/cls_prepped_file_practical_v1_monthcomp.rds")

# Filter for not dead and reasonable crowns only
tree_df_for_cls_clean <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file_practical_v1_monthcomp.rds")
ref_trees <- fread('/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv')

ref_trees_sub <- ref_trees %>% filter(tpconditio %in% c("Excellent", "Good", "Fair"),
                     tpstructur %in% c("Full"))

tree_df_for_cls_clean_filt <- tree_df_for_cls_clean %>% filter(Poly_ID %in% unique(ref_trees_sub$Poly_ID))
# from 225014 trees to 217609, not much but it's not nothing
saveRDS(tree_df_for_cls_clean_filt, "/Volumes/NYC_geo/tree_classification/cls_prepped_file_practical_v1_monthcomp_filt.rds")

#####

### THIS SECTION WAS RUN ON BIG COMPUTER, needs editing

out_path <- "E:/nyc/classification/outputs_practical_v1_monthcomp"
#name_prefix <- "no_platanus_"
name_prefix <- "practical_v1_monthcomp_"
#name_prefix <- "everything_no_sos_" # if any name modifications are desired appended to beginning of each output, make to sure include trailing "_"
#name_prefix <- "lidar_only_"
#name_prefix <- "sos_lidar_only_"

#tree_df_for_cls_clean <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file.rds") # not too bad for loading in, < 30 seconds on SSD
tree_df_for_cls_clean <- readRDS("E:/nyc/classification/cls_prepped_file_practical_v1_monthcomp.rds")

# Drop SOS from 2017, not using 2017 for this analysis
# col_list <- colnames(tree_df_for_cls_clean)
# col_list_sub <- col_list[!grepl("_2017", col_list)]
# tree_df_for_cls_clean <- tree_df_for_cls_clean[,..col_list_sub]

# Section for more filtering
# col_list <- colnames(tree_df_for_cls_clean)
# # lidar cols only (if needed)
# kept_cols <- c("Poly_ID", "genus", "npts", "height_max_ft", "height_med_ft", "int_mean", "int_mean_above_medh", "int_mean_below_medh", "wid_medh_max_ft", 
#                "wid_medh_circlerad_ft", "hw_rat_2_max", "hw_rat_2_circlerad", "cp_3_max", "cp_3_circlerad")
# # cols_to_use <- kept_cols # for lidar only run
# 
# # # Filtering to remove SOS (if desired, otherwise comment out) or keep
# # col_list_sub_4b <- col_list[grep("4b_", col_list)]
# # col_list_sub_8b <- col_list[grep("8b_", col_list)]
# # cols_to_use <- c(kept_cols, col_list_sub_4b, col_list_sub_8b)
# 
# # Filter to only use SOS (+lidar)
# col_list_sos <- col_list[grep("SOS_", col_list)]
# cols_to_use <- c(kept_cols, col_list_sos) # sos + lidar only run
# 
# tree_df_for_cls_clean <- tree_df_for_cls_clean[,..cols_to_use] #.. for tibble

#####
# Do an xgboost run
set.seed(14)
#set.seed(26) # CHANGED SEED FOR THIS RUN, otherwise was 14

colnames(tree_df_for_cls_clean)[1] <- "Object_ID"

#df_cls_agg_mean_cc <- tree_df_for_cls_clean[sample(1:nrow(tree_df_for_cls_clean), 2000, 0),] # try test run length
df_cls_agg_mean_cc <- tree_df_for_cls_clean # full version
#df_cls_agg_mean_cc <- tree_df_for_cls_clean %>% filter(genus != "Platanus") # try a no platanus run

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

# Save with xgb.save()
# Reload with xgb.load()

# output xgb model
setwd(out_path)
xgb_outname <- paste0(name_prefix, "xgb_full.model")
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
confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")
# so xgboost is not a panacea, but does work better, especially for some small classes. 
# 50 rounds takes 8x longer to run than ranger with 500 rf trees though

xgb_imp <- xgb.importance(model = bst_model)
xgb_imp_df <- paste0(name_prefix, "xgb_full_importance.csv")
write.csv(xgb_imp, xgb_imp_df, row.names = FALSE)

# Additional code for if there are less than 25 variables
nvar <- nrow(xgb_imp)
if(nvar > 25){
  var_plotted <- 25
} else {
  var_plotted <- nvar
}

ggplot(xgb_imp[1:var_plotted,]) +
  geom_col(aes(x = Gain, y = fct_reorder(Feature, Gain))) +
  labs(x = "Gain", y = "Variable", title = paste0("Top ", var_plotted, " most important variables"))
xgb_imp_fig <- paste0(name_prefix, "xgb_full.jpg")
ggsave(xgb_imp_fig, height = 5, width = 8, units = "in")

#
xgb_cm <- confusionMatrix(factor(test_prediction$max_prob),
                          factor(test_prediction$label),
                          mode = "everything")
xgb_cm_out <- paste0(name_prefix, "xgb_full_confmat.rds")
saveRDS(xgb_cm, xgb_cm_out)

# save test prediction for referencing
test_prediction_path <- paste0(name_prefix, "xgb_full_predictions.csv")
write.csv(cbind.data.frame(test_data[,1:2], test_prediction), test_prediction_path, row.names = FALSE)

prod_acc <- unname(xgb_cm$byClass[,1])
user_acc <- unname(xgb_cm$byClass[,3])
class_num_xgb <- sort(unique(test_data$genus))
classes <- genus_list[class_num_xgb + 1] # add 1 because xgb numbering for classes starts at 0
cm_results <- cbind.data.frame(classes, class_num_xgb, prod_acc, user_acc)

test_prediction_accuracies <- paste0(name_prefix, "xgb_full_accuracies.csv")
write.csv(cm_results, test_prediction_accuracies, row.names = FALSE)

cm_results$classes <- as.factor(cm_results$classes)

ggplot(cm_results) +
  geom_point(aes(x = prod_acc, y = user_acc)) +
  geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  coord_equal() +
  ylim(0,1) + xlim(0,1) +
  labs(x = "Producer's Accuracy", y = "User's Accuracy") +
  theme_bw()
xgb_acc_fig <- paste0(name_prefix, "xgb_full_acc.jpg")
ggsave(xgb_acc_fig, height = 6, width = 6, units = "in")


#####

# Updating from here

# Predict back onto extract tree data

# Load xgb model
bst_model <- xgb.load("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_1_monthcomp_filt/practical_v1_1_monthcomp_filt_xgb_full.model")

# Load crown metrics
setwd("/Volumes/NYC_geo/nyc_lidar_metrics/crown_metrics")
crown_metrics_file_list <- list.files()
tree_crown_metrics <- purrr::map_df(crown_metrics_file_list, fread)

# Loop (if needed)

# Load each ID set for monthcomp:
# Load 4 b
# Load 8 b

# Aggregated spectra
# 4 band
setwd("/Volumes/NYC_geo/tree_classification/extracted_4band_monthcomp/zscaled/")
file_list <- list.files()
tree_spectra_4band <- purrr::map_df(file_list, fread)
colnames(tree_spectra_4band)[2:ncol(tree_spectra_4band)] <- paste0("4b_", colnames(tree_spectra_4band)[2:ncol(tree_spectra_4band)])
tree_spectra_4band_sub <- tree_spectra_4band %>% filter(Object_ID %in% tree_crown_metrics$Poly_ID)
rm(tree_spectra_4band)
gc()

# 8 band
setwd("/Volumes/NYC_geo/tree_classification/extracted_8band_monthcomp/zscaled/")
file_list <- list.files()
tree_spectra_8band <- purrr::map_df(file_list, fread)
colnames(tree_spectra_8band)[2:ncol(tree_spectra_8band)] <- paste0("8b_", colnames(tree_spectra_8band)[2:ncol(tree_spectra_8band)])
tree_spectra_8band_sub <- tree_spectra_8band %>% filter(Object_ID %in% tree_crown_metrics$Poly_ID)
rm(tree_spectra_8band)
gc()

# merge crown metrics, 4 b, and 8 b into single file
tree_spectra <- merge(tree_spectra_4band_sub, tree_spectra_8band_sub, by = c("Object_ID"))
rm(tree_spectra_4band_sub)
rm(tree_spectra_8band_sub)

tree_spectra_lcm <- merge(tree_spectra, tree_crown_metrics, by.x = "Object_ID", by.y = "Poly_ID")
rm(tree_spectra)
rm(tree_crown_metrics)
gc()

# na omit
tree_spectra_lcm <- na.omit(tree_spectra_lcm) 
tree_spectra_lcm_mat <- as.matrix(tree_spectra_lcm)
summed_rows <- rowSums(tree_spectra_lcm_mat) # much faster on matrix
tree_spectra_lcm_clean <- tree_spectra_lcm[is.finite(summed_rows), ] # removing any more trees with missing or infinite data.
# losing 4180 trees because some data is missing, 1810991 to 1806811
# would be able to fill in later versions

# Setup loop

tree_ids <- tree_spectra_lcm_clean$Object_ID

#i <- 1
stepsize <- 50000
top <- floor(length(tree_ids)/stepsize)*stepsize+1
step_ranges <- seq(1, top, stepsize)

# for genus names
tree_df_for_cls_clean <- readRDS("/Volumes/NYC_geo/tree_classification/cls_prepped_file_practical_v1_1_monthcomp_filt.rds")
genus_list <- tree_df_for_cls_clean$genus %>% unique() %>% sort()
rm(tree_df_for_cls_clean)
gc()

for (idset in 1:length(step_ranges)){
  
  print(idset)
  
  range_min <- step_ranges[idset]
  idset_range <- seq(range_min, range_min + stepsize - 1)
  tree_ids_sub <- tree_ids[idset_range]
  
  tree_spectra_lcm_clean_sub <- tree_spectra_lcm_clean %>% filter(Object_ID %in% tree_ids_sub)
  
  # setup as xgb.DMatrix
  tree_spectra_lcm_noid <- tree_spectra_lcm_clean_sub %>% select(!Object_ID)
  tree_spectra_lcm_input <- tree_spectra_lcm_noid %>% as.matrix() # memory limit hit here if doing the whole thing, loop subset
  test_matrix <- xgb.DMatrix(data = tree_spectra_lcm_input) 
  
  numberOfClasses <- 18 # hard coded, 18 classes
  
  # predict on the model
  test_pred <- predict(bst_model, newdata = test_matrix)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses) %>%
    t() %>%
    data.frame() 
  
  colnames(test_prediction) <- genus_list
  
  test_prediction <- test_prediction %>% mutate(max_prob = max.col(., "last"))
  
  test_prediction <- cbind.data.frame(tree_spectra_lcm_clean_sub$Object_ID, genus_list[test_prediction$max_prob], test_prediction[,1:numberOfClasses])
  colnames(test_prediction)[1:2] <- c("Poly_ID", "Genus_Predicted")
  
  output_name <- paste0("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_1_monthcomp_filt/predictions/", "practical_xgb_v1_1_predictions_set", idset, ".csv")
  write.csv(test_prediction, output_name, row.names = FALSE)
}

rm(tree_spectra_lcm_mat)
gc()

#####
# merge into tnc polygons
# include known polygons
# Genus_Merged
# Genus_Predicted
# Genus_Known
# Species
# Cultivar

# Modified from: merge_street_tree_pheno2.R

# Load TNC polygons
tree_poly_path_full <- "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb" # Note this is the FINAL TNC dataset
tnc_gdb_polys <- st_read(tree_poly_path_full, layer = "treeobjects_2021_nyc")
tnc_gdb_polys$Poly_ID <- 1:nrow(tnc_gdb_polys)

# Load tree points
tree_pts <- st_read('/Users/dlm356/dlm356_files/nyc_trees/Forestry Tree Points_20240228/geo_export_3a66edb6-f0e4-44ec-9a82-e8b36349c231.shp')

# Remove points with bad tpconditio and bad tpstructur
# this won't be the same strict set as training input, but captures all decent stuff
tree_pts <- tree_pts %>% filter(tpconditio %in% c("Excellent", "Good", "Fair", "Poor", "Critical", "Unknown"), # only remove dead
                                      tpstructur %in% c("Full")) # only want Full trees, not stumps or only shafts or retired points

tree_pts_reproj <- st_transform(tree_pts, st_crs(tnc_gdb_polys))

colnames(tree_pts_reproj)[which(colnames(tree_pts_reproj) == "objectid")] <- "Point_ID"

tree_intersect <- st_intersection(tnc_gdb_polys, tree_pts_reproj) # intersection step

# which polygons only appear once after intersection, so only 1 stem
poly_ids_ints <- table(tree_intersect$Poly_ID) %>% as.data.frame()
poly_ids_1unique <- poly_ids_ints$Var1[which(poly_ids_ints$Freq == 1)] %>% as.character() %>% as.numeric()
# these are the trees we actually have unique labels for!

# String splitting genus + species + cultivar information from source tree points data set
tree_intersect_1unique_df <- tree_intersect %>% filter(Poly_ID %in% poly_ids_1unique) %>% st_drop_geometry()
front <- sapply(strsplit(tree_intersect_1unique_df$genusspeci, " - "), '[', 1)
species <- sapply(strsplit(front, " '"), '[', 1)
species2 <- sapply(strsplit(species, " var. inermis"), '[', 1) # gleditsia split
genus <- sapply(strsplit(species, " "), '[', 1)
species2[which(species2 == "Platanus x acerfolia")] <- "Platanus x acerifolia" # fix typo in original dataset

# merge altogether, these are filtered!
tree_intersect_1unique_df <- cbind.data.frame(genus, species2, tree_intersect_1unique_df)

# Load predicted genus info
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_1_monthcomp_filt/predictions")
file_list <- list.files()
predicted_trees <- purrr::map_df(file_list, fread)

tree_intersect_1unique_df <- tree_intersect_1unique_df %>% select(-Height, -Radius, -SHAPE_Length, -SHAPE_Area) # drop these columns from TNC, will add this in for every row

# do two full joins to merge everything
predictions_w_tree_points <- full_join(predicted_trees, tree_intersect_1unique_df, by = "Poly_ID")
predictions_w_tree_points_and_polys <- full_join(predictions_w_tree_points, tnc_gdb_polys, by = "Poly_ID")


known_genus_inds <- which(!is.na(predictions_w_tree_points_and_polys$genus))
predictions_w_tree_points_and_polys$Genus_Merged <- predictions_w_tree_points_and_polys$Genus_Predicted
predictions_w_tree_points_and_polys$Genus_Merged[known_genus_inds] <- predictions_w_tree_points_and_polys$genus[known_genus_inds]

# Move columns
out_geodf <- predictions_w_tree_points_and_polys %>% relocate(Genus_Merged, .after = Poly_ID)
out_geodf <- out_geodf %>% relocate(genus, .after = Genus_Predicted)
out_geodf <- out_geodf %>% relocate(species2, .after = genus)
out_geodf <- out_geodf %>% relocate(genusspeci, .after = species2)
colnames(out_geodf)[4:6] <- c("Genus_Ref", "Species_Ref", "FullCultivar_Ref")

out_geodf <- out_geodf %>% mutate(across(all_of(genus_list), round, 4)) %>% 
  mutate(across(all_of(genus_list), format, scientific = FALSE)) # reformat accuracy estimates, 4 digit, no scientific notation

# Write this whole thing out as a geopackage
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_1_monthcomp_filt/")
st_write(out_geodf, "cls_poly_practical_v1_1_monthcomp_filt.gpkg")

# load confusion matrix
confmat <- readRDS("practical_v1_1_monthcomp_filt_xgb_full_confmat.rds")

#####
# Merge lidar data into new geopackage with existing classification
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_1_monthcomp_filt/")
out_geodf <- st_read("cls_poly_practical_v1_1_monthcomp_filt.gpkg")

setwd("/Volumes/NYC_geo/nyc_lidar_metrics/crown_metrics")
crown_metrics_file_list <- list.files()
tree_crown_metrics <- purrr::map_df(crown_metrics_file_list, fread)
tree_crown_metrics <- tree_crown_metrics[order(tree_crown_metrics$Poly_ID),]

out_geodf <- out_geodf %>% rename("TNC_Height" = "Height")
out_geodf <- out_geodf %>% rename("TNC_Radius" = "Radius")
out_geodf <- out_geodf %>% rename("TNC_SHAPE_Length" = "SHAPE_Length")
out_geodf <- out_geodf %>% rename("TNC_SHAPE_Area" = "SHAPE_Area")

tree_crown_metrics[,2:ncol(tree_crown_metrics)] <- round(tree_crown_metrics[,2:ncol(tree_crown_metrics)], 4)

#tree_spectra_wlabels_sos_lcm <- merge(tree_spectra_wlabels_sos, tree_crown_metrics) # lcm = lidar crown metrics
out_geodf_2 <- left_join(out_geodf, tree_crown_metrics, by = "Poly_ID") # lcm = lidar crown metrics
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_1_monthcomp_filt/")
st_write(out_geodf_2, "cls_poly_practical_v1_1_monthcomp_filt_w_crown_metrics.gpkg")


table(out_geodf_2$Genus_Predicted)


#####
# Merge SOS 50 dates into new geopackages with existing classification
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_monthcomp/")
out_geodf <- st_read("cls_poly_practical_v1_monthcomp_w_crown_metrics.gpkg")

file_list <- list.files("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/raw_pheno_extract/")
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/raw_pheno_extract/")
big_sos <- purrr::map_df(file_list, fread)

sos_2018 <- big_sos %>% filter(Year == 2018) %>% select(Object_ID, SOS_50)
sos_2018 <- sos_2018 %>% rename(SOS_2018 = SOS_50)

sos_2019 <- big_sos %>% filter(Year == 2019) %>% select(Object_ID, SOS_50)
sos_2019 <- sos_2019 %>% rename(SOS_2019 = SOS_50)

sos_2020 <- big_sos %>% filter(Year == 2020) %>% select(Object_ID, SOS_50)
sos_2020 <- sos_2020 %>% rename(SOS_2020 = SOS_50)

sos_2021 <- big_sos %>% filter(Year == 2021) %>% select(Object_ID, SOS_50)
sos_2021 <- sos_2021 %>% rename(SOS_2021 = SOS_50)

sos_2022 <- big_sos %>% filter(Year == 2022) %>% select(Object_ID, SOS_50)
sos_2022 <- sos_2022 %>% rename(SOS_2022 = SOS_50)

sos_2023 <- big_sos %>% filter(Year == 2023) %>% select(Object_ID, SOS_50)
sos_2023 <- sos_2023 %>% rename(SOS_2023 = SOS_50)

sos_2024 <- big_sos %>% filter(Year == 2024) %>% select(Object_ID, SOS_50)
sos_2024 <- sos_2024 %>% rename(SOS_2024 = SOS_50)

sos_joined <- sos_2018 %>% 
  full_join(sos_2019) %>%
  full_join(sos_2020) %>%
  full_join(sos_2021) %>%
  full_join(sos_2022) %>%
  full_join(sos_2023) %>%
  full_join(sos_2024)

sos_joined <- sos_joined %>% rename(Poly_ID = Object_ID)
out_geodf_2 <- out_geodf %>% left_join(sos_joined)

setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_monthcomp/")
st_write(out_geodf_2, "cls_poly_practical_v1_monthcomp_w_crown_metrics_sos.gpkg")


######
# Quercus palustris only
setwd("/Volumes/NYC_geo/tree_classification/outputs_practical_v1_monthcomp/")
out_geodf <- st_read("cls_poly_practical_v1_monthcomp_w_crown_metrics_sos_qupa.gpkg")

out_geodf_qupa <- out_geodf %>% filter(Species_Ref == "Quercus palustris")
st_write(out_geodf_qupa, "cls_poly_practical_v1_monthcomp_w_crown_metrics_sos_qupa.gpkg")

#####
# Quercus palustris high quality 2019

sos_qupa <- big_sos %>% 
  filter(Object_ID %in% out_geodf_qupa$Poly_ID,
                   Year == 2019,
                   R2 > 0.8)

out_geodf_qupa_filt <- out_geodf_qupa %>% filter(Poly_ID %in% sos_qupa$Object_ID)
st_write(out_geodf_qupa_filt, "cls_poly_practical_v1_monthcomp_w_crown_metrics_sos_qupa_R2_08.gpkg")
