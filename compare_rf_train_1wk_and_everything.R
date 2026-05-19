# Compare outputs of ranger dimension reduction runs

library(tidyverse)
library(ranger)
library(caret) # for confusion matrix
library(data.table)
# 1 week vs everything

getSelVars <- function(file_path_rds){
  x <- readRDS(file_path_rds)
  sel_vars <- x$vars_rounds$var[which(x$vars_rounds$max_round >= x$iter_round_selected)]
  return(sel_vars)
}

# tree information, mset
mset <- readRDS("/Volumes/NYC_geo/tree_classification/mset_v2_wlidarvars.rds")

v1wk <- getSelVars("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs/1week_singleyears_selected_features_info.rds")
setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/")
df1wk <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/tree_points_4b_spectra_raw_1weekcomplong_singleyear_subset_wide.rds")
names(df1wk) <- gsub("-", "", names(df1wk))
df1wk <- df1wk %>% select(all_of(c("Object_ID", v1wk)))

df_spectra <- df1wk


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

i <- 1
set.seed(i) # here for ranger seed, not just sampling seed
rf_trained <- ranger(genus ~ ., train_spectra, importance = 'permutation', num.threads = 10, num.trees = 500, seed = 1, local.importance = TRUE) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
val_predicted <- predict(rf_trained, data = val_spectra)
conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)

rf_trained_1wk <- rf_trained
val_predicted_1wk <- val_predicted
conf_mat_1wk <- conf_mat
df_spectra_1wk <- df_spectra
train_spectra_1wk <- train_spectra
genus_var_imp_1wk <- as.data.table(rf_trained$variable.importance.local)[,genus := train_spectra$genus][,lapply(.SD,mean),by=genus]

#####
# Repeat for everything, not just 1 week single year
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

i <- 1
set.seed(i) # here for ranger seed, not just sampling seed
rf_trained <- ranger(genus ~ ., train_spectra, importance = 'permutation', num.threads = 10, num.trees = 500, seed = 1, local.importance = TRUE) # default settings, ranger, permutation importance to allow local.importance, unscaled for now (should this be scaled though?)
val_predicted <- predict(rf_trained, data = val_spectra)
conf_mat <- confusionMatrix(val_predicted$predictions, val_spectra$genus)

rf_trained_everything <- rf_trained
val_predicted_everything <- val_predicted
conf_mat_everything <- conf_mat
df_spectra_everything <- df_spectra
train_spectra_everything <- train_spectra
genus_var_imp_everything <- as.data.table(rf_trained$variable.importance.local)[,genus := train_spectra$genus][,lapply(.SD,mean),by=genus]


# Accuracies are almost exactly the same...
conf_mat_1wk$byClass[,1]
conf_mat_everything$byClass[,1]

conf_mat_1wk$byClass[,3]
conf_mat_everything$byClass[,3]

# Variable importance
sort(rf_trained_1wk$variable.importance, decreasing = TRUE)[1:15]
sort(rf_trained_everything$variable.importance, decreasing = TRUE)[1:15]
# these are different lists, the everything version is not all 1 week

# Importance comparison for individual genera
#library(data.table)
#genus_var_imp <- as.data.table(rf_trained$variable.importance.local)[,genus := train_spectra$genus][,lapply(.SD,mean),by=genus]

# This is the maximum
max_cols_1wk <- max.col(genus_var_imp_1wk[,2:ncol(genus_var_imp_1wk)]) + 1
colnames(genus_var_imp_1wk)[max_cols_1wk]
max_vals <- as.matrix(genus_var_imp[,..max_cols]) %>% diag() # these are the importance values

max_cols_everything <- max.col(genus_var_imp_everything[,2:ncol(genus_var_imp_everything)]) + 1
colnames(genus_var_imp_everything)[max_cols_everything]

genus_var_imp_everything_df <- data.frame(genus_var_imp_everything$genus, colnames(genus_var_imp_everything)[max_cols_everything])
colnames(genus_var_imp_everything_df) <- c("genus", "MaxImp_everything")
genus_var_imp_1wk_df <- data.frame(genus_var_imp_1wk$genus, colnames(genus_var_imp_1wk)[max_cols_1wk])
colnames(genus_var_imp_1wk_df) <- c("genus", "MaxImp_1wk")
genus_var_imp_max <- merge(genus_var_imp_1wk_df, genus_var_imp_everything_df)
