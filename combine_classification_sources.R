library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(caret)
library(ranger)
library(missRanger)
library(ranger)
library(tuneRanger)
library(xgboost)
library(ggrepel)

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

# will need to load this with purrr
#tree_spectra <- fread("/Volumes/NYC_geo/tree_classification/extracted_8band_traintest/zscaled/tree_points_8b_spectra_zscaled_monthmean_idset1.csv")
setwd("/Volumes/NYC_geo/tree_classification/extracted_8band_traintest/zscaled/")
file_list <- list.files()
tree_spectra <- purrr::map_df(file_list, fread)

gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")

# Filter based on minimum summer NDVI threshold for all years
setwd('/Volumes/NYC_geo/tree_mortality')
mean_ndvi_all <- fread("mean_summer_ndvi.csv")
ndvi_min <- 0.3 # just picked this number for now, Thapa
#ndvi_min <- 0.5 # just picked this number for now, Alonzo
#mean_ndvi_all_live <- mean_ndvi_all %>% filter(ndvi_2024 > ndvi_min) # accidentally tested this first, will do it the other way next
mean_ndvi_all_live <- mean_ndvi_all %>% filter_at(vars(starts_with("ndvi")), all_vars(. > ndvi_min))

# combine everything to setup an RF run like previous organization
#tree_pheno_sub <- tree_pheno %>% filter(dbh > 4 & tpstructur == "Full" & genus %in% gen_list & Poly_ID %in% unique(tree_spectra$Object_ID))
tree_pheno_sub <- tree_pheno %>% filter(dbh > 4 & tpstructur == "Full" & genus %in% gen_list & Poly_ID %in% unique(tree_spectra$Object_ID) & Poly_ID %in% unique(mean_ndvi_all_live$Object_ID)) 

#add in lidar variables
lidar_vars <- tree_pheno_sub[match(tree_spectra$Object_ID, tree_pheno_sub$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
df_cls_agg_mean_lidar <- merge(tree_spectra, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

# Add in SOS and EOS variables
pheno_vars <- tree_pheno_sub %>% filter(Poly_ID %in% tree_spectra$Object_ID) %>% select("Poly_ID", "Year", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)

# # Just a few missing values due to infeasible pheno fits
# # Maybe these should not be included in training, but should be included in validation??
# mr_pheno_vars_wide <- missRanger(pheno_vars_wide[, 2:ncol(pheno_vars_wide)]) # don't include Poly_ID in the missRanger filling inputs
# pheno_vars_wide_filled <- round(mr_pheno_vars_wide, 0)
# pheno_vars_wide_filled <- cbind.data.frame(pheno_vars_wide$Poly_ID, pheno_vars_wide_filled)
# colnames(pheno_vars_wide_filled)[1] <- "Poly_ID"

#length(which(complete.cases(pheno_vars_wide))) # 323479
# 323479/325083, 99.5% are complete, let's just use these for now.

# Merge all together
#df_cls_agg_mean_lidar_pheno <- merge(df_cls_agg_mean_lidar, pheno_vars_wide_filled, by.x = "Object_ID", by.y = "Poly_ID")

pheno_vars_wide_sub <- na.omit(pheno_vars_wide)
df_cls_agg_mean_lidar_pheno <- merge(df_cls_agg_mean_lidar, pheno_vars_wide_sub, by.x = "Object_ID", by.y = "Poly_ID") #subsetting for complete is built-in on the merge


# add genus
genus_vars <- tree_pheno_sub[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno_sub$Poly_ID),] %>% select("Poly_ID", "genus")
colnames(genus_vars)[1] <- "Object_ID"
df_cls_agg_mean_lidar_pheno <- merge(genus_vars, df_cls_agg_mean_lidar_pheno)


# # Swap in species for Acer
# sp_acer <- c('Acer platanoides', 'Acer rubrum', 'Acer saccharinum', 'Acer',
#              'Acer campestre', 'Acer pseudoplatanus', 'Acer saccharum', 'Acer ginnala',
#              'Acer palmatum', 'Acer griseum', 'Acer negundo', 'Acer x freemanii')
# 
# 
# tree_pheno_sub2 <- tree_pheno_sub[match(df_cls_agg_mean_lidar_pheno$Object_ID[df_cls_agg_mean_lidar_pheno$genus == 'Acer'], tree_pheno_sub$Poly_ID)]
# tree_pheno_sub2 <- tree_pheno_sub2 %>% filter(species %in% sp_acer)
# sp_acer_vars <- tree_pheno_sub2[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno_sub2$Poly_ID),] %>% select("Poly_ID", "species") %>% na.omit()
# sp_acer_vars$sub_ind <- match(sp_acer_vars$Poly_ID, df_cls_agg_mean_lidar_pheno$Object_ID)
# df_cls_agg_mean_lidar_pheno$genus <- replace(df_cls_agg_mean_lidar_pheno$genus, sp_acer_vars$sub_ind, sp_acer_vars$species)

# Change to be as.factor()
df_cls_agg_mean_lidar_pheno$genus <- as.factor(df_cls_agg_mean_lidar_pheno$genus)



#####
# Test RF
# RF Run
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno) # do this cleaning just for prep

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID) # %>% as.data.frame()
test_data_noid <- test_data %>% select(!Object_ID) # %>% as.data.frame()

# might want to do tuneRanger or similar to adjust the hyperparameters used for the ranger RF
# hit a fatal error in R when tried to do this last time
# treecls_task <- makeClassifTask(data = train_data_noid, target = "genus")
# estimateTimeTuneRanger(treecls_task)
# # Tuning
# res = tuneRanger(iris.task, measure = list(multiclass.brier), num.trees = 1000, 
#                  num.threads = 2, iters = 70, save.file.path = NULL)

# Also need to consider rebalancing class levels using SMOTE or some other rebalancing algorithm
# Should also measure runtime for ranger model run

# <20 minutes with 500 trees
rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity", num.trees = 500, num.threads = 10)
#rf_model <- ranger(genus ~ ., train_data_noid, importance = "permutation", local.importance = TRUE) # Doing local variable importance, trying to see what is driving each tree class
p_rf <- predict(rf_model, test_data_noid[,2:ncol(test_data_noid)])
confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

# This works, but very different for order of importance likely because this is not really a representative sample

# Make plot from confusion matrix for Sensitivity (Producer's Acc) vs. Pos Pred Value (User's Acc)
rf_cm <- confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)

prod_acc <- unname(rf_cm$byClass[,1])
user_acc <- unname(rf_cm$byClass[,3])
classes <- rownames(rf_cm$table)
cm_results <- cbind.data.frame(classes, prod_acc, user_acc)
cm_results$classes <- as.factor(cm_results$classes)

ggplot(cm_results) +
  geom_point(aes(x = prod_acc, y = user_acc)) +
  geom_text_repel(aes(x = prod_acc, y = user_acc, label = classes)) +
  coord_equal() +
  theme_bw()

#####
# Try an xgboost run
set.seed(14)

# test on a much smaller sample, will need to try to scale this up if we want to use this.
#df_cls_agg_mean_cc <- df_cls_agg_mean_cc[sample(1:nrow(df_cls_agg_mean_cc)),]
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno) # no difference here, old code setup
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
nround    <- 50 # number of XGBoost rounds
esr <- 10 # early stopping rounds, need to figure out how to use this
cv.nfold  <- 5

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

#####
# Comparing RF and XGB
classes <- rownames(rf_cm$table)
rf_prod_acc <- unname(rf_cm$byClass[,1])
rf_user_acc <- unname(rf_cm$byClass[,3])
xgb_prod_acc <- unname(xgb_cm$byClass[,1])
xgb_user_acc <- unname(xgb_cm$byClass[,3])

xgb_n_oa <- c()
rf_n_oa <- c()
for (i in 1:nrow(xgb_cm$table)){
  xgb_n_oa <- c(xgb_n_oa, xgb_cm$table[i,i])
  rf_n_oa <- c(rf_n_oa, rf_cm$table[i,i])
}

rf_xgb_compare <- cbind.data.frame(classes, rf_prod_acc, xgb_prod_acc, rf_user_acc, xgb_user_acc, rf_n_oa, xgb_n_oa)


#####

# Get area-weighted OA
area_fraction <- test_data_noid$SHAPE_Area/sum(test_data_noid$SHAPE_Area)
area_weighted_oa <- sum(area_fraction[which(test_prediction$label == test_prediction$max_prob)])

# Confidence in correct calls, generate a box plot from this
prob_correct_all <- as.data.frame(matrix(data = NA, nrow = length(area_fraction[which(test_prediction$label == test_prediction$max_prob)]), ncol = 2))
colnames(prob_correct_all) <- c("ID_Cor", "Prob")
row_track <- 1
for (i in 1:numberOfClasses){
  print(i)
  result_sub <- test_prediction %>% filter(label == i & max_prob == i)
  prob_correct <- result_sub[,i]
  len_newcor <- length(prob_correct)
  prob_correct_all$ID_Cor[row_track:(row_track+len_newcor-1)] <- i
  prob_correct_all$Prob[row_track:(row_track+len_newcor-1)] <- prob_correct
  row_track <- row_track + len_newcor
}

prob_correct_all$classes <- factor(prob_correct_all$ID_Cor, labels = classes)

ggplot(prob_correct_all) +
  geom_boxplot(aes(x = classes, y = Prob), outlier.size = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Genus", y = "Probability when correct")


#####
# Testing to do, RF and XGB
# Remove small classes and test accuracy
# Tune input hyperparameters
# Split genus into species, start with Acer
# likely need to try some of these setup with smaller subset of total to improve efficiency in testing


