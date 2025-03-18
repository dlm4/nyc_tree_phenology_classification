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
pheno_output_ranger_df_noNA <- pheno_output_ranger_df[complete.cases(pheno_output_ranger_df),]

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

# to retain rf setup
set.seed(14)
rf_trees <- ranger(TreeLabel ~ ., data = rf_input_train, importance = "permutation") # only works when there isn't missing data...
trees_predicted <- predict(rf_trees, data = rf_input_val)
pred_table <- table(rf_input_val$TreeLabel, trees_predicted$predictions)

val_correct <- pred_table[1,1] + pred_table[2,2]
val_total <- sum(pred_table)
val_acc <- val_correct/val_total # 93% accuracy for Platanus on separate validation data

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
