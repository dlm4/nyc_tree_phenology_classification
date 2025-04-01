# Make a script that loads in a subset of the tree database for a neighborhood with likely a lot of Platanus unmapped
# Use existing extracted point phenology data + lidar polygon variables (and temperature?? perhaps not needed) to train random forest
# Ranger (random forest) training for binary classifier of Platanus / not Platanus (very simple)
# Apply back onto unmapped trees in neighborhood for binary classifier demonstration

library(tidyverse)
library(data.table)
library(lubridate)
library(reshape2)
library(ranger)
library(missRanger)
library(terra)
library(sf)
`%notin%` <- Negate(`%in%`)
library(purrr)

pheno_output <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output.csv")

# Filter based on R2? 0.8 or higher might be needed
# Also likely filter based on tpconditio: c('Excellent', 'Good', 'Fair')
# And tpstructur: c('Full') only!

# setup column for genus and then reassigned to 'Platanus' and 'OtherGenus'

# Variables

# To predict:
# TreeLabel: c('Platanus', 'OtherGenus')

# Inputs:
# PlanetScope derived:
# SOS_50_2018
# SOS_50_2019
# SOS_50_2020
# SOS_50_2021
# SOS_50_2022
# SOS_50_2023
# EOS_50_2018
# EOS_50_2019
# EOS_50_2020
# EOS_50_2021
# EOS_50_2022
# EOS_50_2023

# TNC Polygons (from lidar)
# Height
# Radius
# SHAPE_Length
# SHAPE_Area

# MAYBE, Location
# Lon
# Lat


# # SOS
# sos_vals <- pheno_output[, c("Poly_ID", "Year", "SOS_50")] %>% pivot_wider(names_from = Year, values_from = SOS_50, names_prefix = "SOS_50_", names_sort = TRUE, values_fill = NA)
# 
# # EOS
# eos_vals <- pheno_output[, c("Poly_ID", "Year", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = EOS_50, names_prefix = "EOS_50_", names_sort = TRUE, values_fill = NA)
# 
# # Both at once
# pheno_output[, c("Poly_ID", "Year", "SOS_50", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = c(SOS_50, EOS_50), names_sort = TRUE, values_fill = NA)
# 
# pheno_output[, c("Poly_ID", "genus", "Height", "Radius", "SHAPE_Length", "SHAPE_Area", "Year", "SOS_50", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = c(SOS_50, EOS_50), names_sort = TRUE, values_fill = NA)

pheno_output$TreeLabel <- pheno_output$genus
pheno_output$TreeLabel[which(pheno_output$TreeLabel != "Platanus")] <- "Other_Genus"

#pheno_output[, c("TreeLabel", "Height", "Radius", "SHAPE_Length", "SHAPE_Area", "Year", "SOS_50", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = c(SOS_50, EOS_50), names_sort = TRUE, values_fill = NA)

pheno_output_sub <- pheno_output %>% filter(tpconditio %in% c('Excellent', 'Good', 'Fair') & tpstructur == "Full" & R2 > 0.8)
pheno_output_ranger_df <- pheno_output_sub[, c("TreeLabel", "Height", "Radius", "SHAPE_Length", "SHAPE_Area", "Year", "SOS_50", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = c(SOS_50, EOS_50), names_sort = TRUE, values_fill = NA)
cc_inds <- complete.cases(pheno_output_ranger_df)
pheno_output_ranger_df_noNA <- pheno_output_ranger_df[cc_inds,]

# to retain random samples
set.seed(14)
rf_input_full <- pheno_output_ranger_df_noNA # could revise this
rf_input_full$TreeLabel <- as.factor(rf_input_full$TreeLabel)
train_fraction <- 0.8
n_elements <- nrow(rf_input_full)
ntrain <- floor(n_elements*train_fraction)
nval <- n_elements - ntrain
element_ids <- 1:n_elements
train_sample <- sample(element_ids, ntrain)
val_sample <- element_ids[element_ids %notin% train_sample]
rf_input_train <- rf_input_full[train_sample,]
rf_input_val <- rf_input_full[val_sample,]

true_cc <- which(cc_inds == TRUE)
pheno_output_ranger_df_nottrain <- pheno_output_ranger_df[-true_cc[train_sample],] # exclude the values that were included in training

# to retain rf setup
set.seed(14)
rf_trees <- ranger(TreeLabel ~ ., data = rf_input_train, importance = "permutation") # only works when there isn't missing data...

# prediction
trees_predicted <- predict(rf_trees, data = rf_input_val)
rf_input_val$TreeLabel_o <- paste0(as.character(rf_input_val$TreeLabel), "_o")
pred_table <- table(rf_input_val$TreeLabel_o, trees_predicted$predictions)
# original is rows, predictions are columns

#             Other_Genus   Platanus
#Other_Genus       12915      389
#Platanus            906     4335

val_correct <- pred_table[1,1] + pred_table[2,2]
val_total <- sum(pred_table)
val_acc <- val_correct/val_total # 93% accuracy for Platanus on separate validation data

as.data.frame(rf_trees$variable.importance)

# To fill missing values, would need to run missRanger to predict on all the values that don't exist first...
# Or have to accept that some trees are going to be missing and cannot be labeled

# # This is much worse for Platanus if there is no R2 fitting filter first
# pheno_output_sub2 <- pheno_output %>% filter(tpconditio %in% c('Excellent', 'Good', 'Fair') & tpstructur == "Full")
# pheno_output_ranger_df2 <- pheno_output_sub2[, c("TreeLabel", "Height", "Radius", "SHAPE_Length", "SHAPE_Area", "Year", "SOS_50", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = c(SOS_50, EOS_50), names_sort = TRUE, values_fill = NA)
# pheno_output_ranger_df2_noNA <- pheno_output_ranger_df2[complete.cases(pheno_output_ranger_df2),]
# 
# # to retain random samples
# set.seed(14)
# rf_input_full <- pheno_output_ranger_df2_noNA # could revise this
# rf_input_full$TreeLabel <- as.factor(rf_input_full$TreeLabel)
# train_fraction <- 0.8
# n_elements <- nrow(rf_input_full)
# ntrain <- floor(n_elements*train_fraction)
# nval <- n_elements - ntrain
# element_ids <- 1:n_elements
# train_sample <- sample(element_ids, ntrain)
# val_sample <- element_ids[element_ids %notin% train_sample]
# rf_input_train <- rf_input_full[train_sample,]
# rf_input_val <- rf_input_full[val_sample,]
# 
# # to retain rf setup
# set.seed(14)
# rf_trees <- ranger(TreeLabel ~ ., data = rf_input_train, importance = "permutation") # only works when there isn't missing data...
# trees_predicted <- predict(rf_trees, data = rf_input_val)
# pred_table <- table(rf_input_val$TreeLabel, trees_predicted$predictions)
# 
# val_correct <- pred_table[1,1] + pred_table[2,2]
# val_total <- sum(pred_table)
# val_acc <- val_correct/val_total # 93% accuracy for Platanus on separate validation data

#####


# Need to impute SOS and EOS for poor fits or otherwise missing data
# Fallback is polygon information only, no pheno
pheno_output_ranger_df_nottrain <- pheno_output_ranger_df_nottrain %>% mutate(across(6:17, as.numeric))
mr_pheno <- missRanger(pheno_output_ranger_df_nottrain) # this can take awhile
# if this works OK, would need to impute all missing values first, then run classifier on it
# will need to setup in a way to retain accuracy of the random forests used to fill gaps
# estimates SOS and EOS dates as decimals, so might need to round for mapping (and maybe for classifying too)
mr_pheno_trees_predicted <- predict(rf_trees, data = mr_pheno)

table(mr_pheno$TreeLabel, mr_pheno_trees_predicted$predictions)
#               Other_Genus Platanus
#Other_Genus      139482     2133
#Platanus           7974    22617

# Overall accuracy is still 94%

# Classification does a much better job on trees with good fits rather than those with imputed SOS and EOS values, but at least it works and is able to make a prediction
# Alternative would be just using the lidar derived variables


#####
# Set it up to classify all unlabeled trees in the Bronx as Platanus or not Platanus

# Read in polygons
# get poly object ids
tree_poly_path_full <- "/Volumes/NYC_geo/tree_polygons/tnc_2021/Trees_Centroids_Crown_Objects_2021.gdb"
tnc_gdb_polys <- st_read(tree_poly_path_full, layer = "treeobjects_2021_nyc")
tnc_gdb_polys$Object_ID <- seq(1,nrow(tnc_gdb_polys))

# load in bronx borough
# limit polygons to just bronx borough
boros <- st_read("/Volumes/NYC_geo/vectors/Borough Boundaries/geo_export_da133389-a6c6-45c3-a980-14295f0e4c2f.shp")
bronx <- boros %>% filter(boro_name == "Bronx")

bronx_reproj <- st_transform(bronx, st_crs(tnc_gdb_polys))

tnc_gdb_polys_bronx <- st_intersection(tnc_gdb_polys, bronx_reproj)

# Next time, retain all the IDs for the Bronx and then use these as reference back to the full polygon map so that polygons don't get cut by borough boundaries

# Read in extracted and labeled SOS and EOS data
pheno_output <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output.csv")

# Read in extracted and unlabeled polygon SOS and EOS data
# Read in phenology data
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_highsunonly_cal_pheno/")
pheno_file_list <- list.files(pattern = glob2rx("trees_pheno_output_objset*point.csv")) # stored as different sets by object id
pheno_file_all <- purrr::map_df(pheno_file_list, fread, .id = 'object_id_group') 


# Anything that is already labeled as Platanus or another genus retains that label
pheno_output_sublabels <- pheno_output[which(pheno_output$Year == 2018), c("Poly_ID", "genus")]

tnc_gdb_polys_bronx_labeled <- merge(tnc_gdb_polys_bronx, pheno_output_sublabels, by.x = "Object_ID", by.y = "Poly_ID", all.x = TRUE)
tnc_gdb_polys_bronx_labeled$TreeLabel <- NA
tnc_gdb_polys_bronx_labeled$TreeLabel[!is.na(tnc_gdb_polys_bronx_labeled$genus)] <- "Other_Genus"
tnc_gdb_polys_bronx_labeled$TreeLabel[tnc_gdb_polys_bronx_labeled$genus == "Platanus"] <- "Platanus"

# Anything that is NOT already labeled as these we apply the classifier, using the rf_trees we made before
unlabeled_ids <- tnc_gdb_polys_bronx_labeled$Object_ID[is.na(tnc_gdb_polys_bronx_labeled$TreeLabel)] # get IDs

pheno_file_all_unlabeled <- pheno_file_all %>% filter(Object_ID %in% unlabeled_ids)
tnc_poly_info_unlabeled <- tnc_gdb_polys_bronx_labeled %>% filter(Object_ID %in% unlabeled_ids) %>% st_drop_geometry()

merged_poly_info <- merge(pheno_file_all_unlabeled, tnc_poly_info_unlabeled, by = "Object_ID") # this is setup to be labeled
# convert to wide on only columns that are needed, needed Object_ID to keep unique rows and for casting back to polygon set
merged_poly_info_ranger_df <- merged_poly_info[, c("Object_ID", "TreeLabel", "Height", "Radius", "SHAPE_Length", "SHAPE_Area", "Year", "SOS_50", "EOS_50")] %>% pivot_wider(names_from = Year, values_from = c(SOS_50, EOS_50), names_sort = TRUE, values_fill = NA)

# BUT, in order to use it, need to use missRanger first to fill in the few SOS and EOS data that totally failed so the classifier can be applied (it can't work with missing values)
miss_ranger_inputs <- merged_poly_info_ranger_df[, !names(merged_poly_info_ranger_df) %in% c("Object_ID", "TreeLabel")]
miss_ranger_inputs <- miss_ranger_inputs %>% mutate(across(5:16, as.numeric))
miss_ranger_output <- missRanger(miss_ranger_inputs) # this is very slow with this many data points, but it's what we've got. Could reduce the overall size to just fill with a few to make this faster, or reduce number of iterations

miss_ranger_output$TreeLabel <- merged_poly_info_ranger_df$TreeLabel
tnc_predictions <- predict(rf_trees, data = miss_ranger_output) # this is the random forest output

# fill into polygons
merged_poly_info_ranger_df$TreeLabel <- tnc_predictions$predictions
tnc_gdb_polys_bronx_labeled$TreeLabel[which(tnc_gdb_polys_bronx_labeled$Object_ID %in% merged_poly_info_ranger_df$Object_ID)] <- as.character(merged_poly_info_ranger_df$TreeLabel)

# Write out to shapefile because geodatabase is not supported by gdal in my install right now
setwd("/Volumes/NYC_geo/tree_polygons/classifications")
st_write(tnc_gdb_polys_bronx_labeled[,c("Object_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area", "TreeLabel")], dsn = "platanus_bronx_class_test.shp")

