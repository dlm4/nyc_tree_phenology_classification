
# Make time series plots across all bands + combinations for a known ginkgo, likely 402050 since this one is big and should be clear signal
# Do the same thing for a Liquidambar? Or other tree?


library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(rpart)
library(rpart.plot)
library(caret)
library(ranger)
library(paletteer)
library(missRanger)

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
}

# Can make this more efficient to only search one Objset for a single id version
getTreePlanetSpectra <- function(ids_all){
  for (i in 1:7){
    print(paste0("Objset ", i))
    objset_text <- paste0("*objset", i, "_*.csv")
    file_list_sub <- list.files(pattern = glob2rx(objset_text))
    names(file_list_sub) <- strsplit(file_list_sub, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 27) %>% ymd() # name the list
    tree_df_sub <- purrr::map_df(file_list_sub, fread, .id = 'date') 
    tree_df_sub2 <- tree_df_sub %>% filter(Object_ID %in% ids_all)
    if (i == 1){
      tree_df_agg <- tree_df_sub2
    } else {
      tree_df_agg <- bind_rows(tree_df_agg, tree_df_sub2)
    }
    rm(tree_df_sub)
    rm(tree_df_sub2)
  }
  gc()
  
  return(tree_df_agg)
}

#

tree_id <- c(402050, 691424)
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
tree_df <- getTreePlanetSpectra(tree_id)
tree_df <- as.data.frame(tree_df)

band_list <- c("coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge", "nir")
colname_list <- c()

for (i in 1:(length(band_list)-1)){
  b1 <- band_list[i]
  for (j in (i+1):length(band_list)){
    b2 <- band_list[j]
    
    # Do longer wavelength first since we might expect longer to be greater in NDVI form
    band_name <- paste("nd", b2, b1, sep = "_")
    print(band_name)
    
    colname_list <- c(colname_list, band_name)
    
    nd_index <- calcNormDif(tree_df[,b2], tree_df[,b1]) # note this only works on data frames, not tibbles
    
    if (i == 1 & j == 2){
      omnbr <- nd_index
    } else {
      omnbr <- cbind.data.frame(omnbr, nd_index)
      colnames(omnbr) <- colname_list
    }
    
  }
}

df_spectra_omnbr <- bind_cols(tree_df, omnbr)

#

df_spectra_omnbr$date <- as.Date(df_spectra_omnbr$date)
cols_set <- c("nd_green_blue", "nd_nir_green", "nd_yellow_blue", "nd_nir_yellow")

df_spectra_omnbr[,c("date", "Object_ID", cols_set)] %>%
  pivot_longer(cols = all_of(cols_set)) %>%
  filter(date %within% interval(ymd("2023-01-01"), ymd("2023-12-31"))) %>%
ggplot(aes(x = date, y = value, color = as.factor(Object_ID))) +
  geom_vline(xintercept = c(ymd("2023-09-14"), ymd("2023-11-23")), color = "gray80") +
  geom_point() +
  geom_line(alpha = 0.3) +
  facet_wrap(~name) +
  theme_bw()


#####
# Try doing this style of plotting for many different trees of varying quality
# Check ginkgos that the classifier got right vs ones that it got wrong
# See how close the color peak is to the EOS date and whether this can be used as a guide
# Check out possible time windows for anomalies in color and consistency
# baseline during summer vs consistent change in the shoulder seasons

#####

# Load in tree spectra, get gingkos from list for 2023, see how color peak compares to EOS?
tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")
tree_pheno_sub <- tree_pheno %>% filter(genus != "Ginkgo" & genus != "Unknown") # Ginkgo

other_tree_ids <- unique(tree_pheno_sub$Poly_ID)
other_tree_ids_sub <- sample(other_tree_ids, 13841) # made this sample the same length as the ginkgo sample

# ggplot(tree_pheno_sub) +
#   geom_point(aes(x = EOS_80, y = EOS_20, color = as.factor(Year), alpha = R2, size = R2)) +
#   geom_abline(slope = 1, intercept = 0)

# Try getting spectra info for all ginkgos
tree_id <- other_tree_ids_sub
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
tree_df <- getTreePlanetSpectra(tree_id)
tree_df <- as.data.frame(tree_df)

band_list <- c("coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge", "nir")
colname_list <- c()

for (i in 1:(length(band_list)-1)){
  b1 <- band_list[i]
  for (j in (i+1):length(band_list)){
    b2 <- band_list[j]
    
    # Do longer wavelength first since we might expect longer to be greater in NDVI form
    band_name <- paste("nd", b2, b1, sep = "_")
    print(band_name)
    
    colname_list <- c(colname_list, band_name)
    
    nd_index <- calcNormDif(tree_df[,b2], tree_df[,b1]) # note this only works on data frames, not tibbles
    
    if (i == 1 & j == 2){
      omnbr <- nd_index
    } else {
      omnbr <- cbind.data.frame(omnbr, nd_index)
      colnames(omnbr) <- colname_list
    }
    
  }
}

df_spectra_omnbr <- bind_cols(tree_df, omnbr)

#
df_spectra_omnbr$date <- as.Date(df_spectra_omnbr$date)

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/")
#write.csv(df_spectra_omnbr, "ginkgo_8b_indices.csv", row.names = FALSE)
#df_spectra_omnbr$Object_ID %>% unique() %>% length() # 13841
write.csv(df_spectra_omnbr, "other_trees_sample_8b_indices.csv", row.names = FALSE)

cols_set <- c("nd_green_blue", "nd_nir_green", "nd_yellow_blue", "nd_nir_yellow")

df_input_polys <- tree_pheno_sub %>% filter(dbh > 18)

df_spectra_omnbr[,c("date", "Object_ID", cols_set)] %>%
  #filter(Object_ID %in% df_input_polys$Poly_ID) %>%
  pivot_longer(cols = all_of(cols_set)) %>%
  filter(date %within% interval(ymd("2023-09-01"), ymd("2023-12-31"))) %>%
ggplot(aes(x = date, y = value, group = as.factor(Object_ID))) +
  geom_vline(xintercept = c(ymd("2023-09-14"), ymd("2023-11-23")), color = "gray80") +
  geom_point(alpha = 0.3) +
  #geom_line(alpha = 0.3) +
  facet_wrap(~name) +
  theme_bw() + theme(legend.position = "none")

#####

# Normalize time series for individual trees
setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/")
df_spectra_omnbr <- fread("ginkgo_8b_indices.csv")
df_spectra_omnbr$date <- as.Date(df_spectra_omnbr$date)

tree_inds <- which(df_spectra_omnbr$Object_ID == 402050) #df_spectra_omnbr$Object_ID[1]
nd_ts_scaled <- scale(df_spectra_omnbr[tree_inds, 3:ncol(df_spectra_omnbr)])


nd_ts_scaled <- scale(nd_ts)
hist(nd_ts, breaks = 30)
hist(nd_ts_scaled, breaks = 30)

plot(df_spectra_omnbr$date[tree_inds], nd_ts_scaled)


df_spectra_omnbr_scaled <- as.data.frame(df_spectra_omnbr)
colnames(df_spectra_omnbr_scaled)[3:ncol(df_spectra_omnbr)] <- paste(colnames(df_spectra_omnbr)[3:ncol(df_spectra_omnbr)], "_zscaled", sep = "")

tree_id_list <- unique(df_spectra_omnbr$Object_ID)
for (i in 1:length(tree_id_list)){
  print(i) # this will take 10+ minutes for 13841 trees
  tree_inds <- which(df_spectra_omnbr$Object_ID == tree_id_list[i])
  nd_ts_scaled <- scale(df_spectra_omnbr[tree_inds, 3:ncol(df_spectra_omnbr)])
  df_spectra_omnbr_scaled[tree_inds, 3:ncol(df_spectra_omnbr_scaled)] <- nd_ts_scaled # doesn't work with tibble
}
# this does all the indices and bands, which is overkill but at least we have it

write.csv(df_spectra_omnbr_scaled, "ginkgo_8b_indices_zscaled.csv", row.names = FALSE)

# then evaluate the peaks for key indices that are likely to increase
# this seems promising
id <- 402050 # 836013 # 630759 #108001 #402050
df_spectra_omnbr_scaled %>% filter(Object_ID == id) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_yellow_blue_zscaled)) +
  geom_hline(yintercept = median(df_spectra_omnbr_scaled$nd_yellow_blue_zscaled[which(month(df_spectra_omnbr_scaled$date) == 11 & df_spectra_omnbr_scaled$Object_ID == id)]), linetype = "dashed") +
  labs(title = id)

df_spectra_omnbr_scaled %>% filter(Object_ID == id) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_green_blue_zscaled)) +
  geom_hline(yintercept = median(df_spectra_omnbr_scaled$nd_green_blue_zscaled[which(month(df_spectra_omnbr_scaled$date) == 11 & df_spectra_omnbr_scaled$Object_ID == id)]), linetype = "dashed") +
  labs(title = id)


# Do normalization again for non-ginkgo trees
# Normalize time series for individual trees
setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/")
df_spectra_omnbr <- fread("other_trees_sample_8b_indices.csv")
df_spectra_omnbr$date <- as.Date(df_spectra_omnbr$date)
df_spectra_omnbr_scaled <- as.data.frame(df_spectra_omnbr)
colnames(df_spectra_omnbr_scaled)[3:ncol(df_spectra_omnbr)] <- paste(colnames(df_spectra_omnbr)[3:ncol(df_spectra_omnbr)], "_zscaled", sep = "")

tree_id_list <- unique(df_spectra_omnbr$Object_ID)
for (i in 1:length(tree_id_list)){
  print(i) # this will take 10+ minutes for 13841 trees
  tree_inds <- which(df_spectra_omnbr$Object_ID == tree_id_list[i])
  nd_ts_scaled <- scale(df_spectra_omnbr[tree_inds, 3:ncol(df_spectra_omnbr)])
  df_spectra_omnbr_scaled[tree_inds, 3:ncol(df_spectra_omnbr_scaled)] <- nd_ts_scaled # doesn't work with tibble
}
# this does all the indices and bands, which is overkill but at least we have it

write.csv(df_spectra_omnbr_scaled, "other_trees_sample_8b_indices_zscaled.csv", row.names = FALSE)

#
id <- 93788
df_spectra_omnbr_scaled %>% filter(Object_ID == id) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_yellow_blue_zscaled)) +
  geom_hline(yintercept = median(df_spectra_omnbr_scaled$nd_yellow_blue_zscaled[which(month(df_spectra_omnbr_scaled$date) == 11 & df_spectra_omnbr_scaled$Object_ID == id)]), linetype = "dashed") +
  labs(title = id)

df_spectra_omnbr_scaled %>% filter(Object_ID == id) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_green_blue_zscaled)) +
  geom_hline(yintercept = median(df_spectra_omnbr_scaled$nd_green_blue_zscaled[which(month(df_spectra_omnbr_scaled$date) == 11 & df_spectra_omnbr_scaled$Object_ID == id)]), linetype = "dashed") +
  labs(title = id)

# other trees also can have relatively high values in Nov, but perhaps not as high and not as consistently

#####
# Take Ginkgo and Other trees
# get monthly median of nd_yellow_blue_scaled and nd_green_blue_scaled
# Differences *should* occur in month 11 (and 10)
# include all monthly medians (probably) in an RF run
# Could include other normalized scaled indices too but these likely tell most of the story that we're looking for??

ginkgo <- fread("ginkgo_8b_indices_zscaled.csv") %>% as.data.frame()
other <- fread("other_trees_sample_8b_indices_zscaled.csv") %>% as.data.frame()


ginkgo_ids <- unique(ginkgo$Object_ID)
ginkgo_ids_samp <- sample(ginkgo_ids, 30)

ginkgo$date <- as.Date(ginkgo$date)

ginkgo %>% filter(Object_ID %in% ginkgo_ids_samp) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_yellow_blue_zscaled)) +
  facet_wrap(~as.factor(Object_ID))

ginkgo %>% filter(Object_ID %in% ginkgo_ids_samp) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_green_blue_zscaled)) +
  facet_wrap(~as.factor(Object_ID))


other_ids <- unique(other$Object_ID)
other_ids_samp <- sample(other_ids, 30)

other$date <- as.Date(other$date)

other %>% filter(Object_ID %in% other_ids_samp) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_yellow_blue_zscaled)) +
  facet_wrap(~as.factor(Object_ID))

other %>% filter(Object_ID %in% other_ids_samp) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(month(date)), y = nd_green_blue_zscaled)) +
  facet_wrap(~as.factor(Object_ID))

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv") %>% as.data.frame()

# Try an RF run

set.seed(14)

ginkgo$genus <- "Ginkgo"
other$genus <- "Other"

df_cls <- rbind.data.frame(ginkgo, other)
#df_cls <- df_cls[which(df_cls$Object_ID %in% tree_pheno$Poly_ID[which(tree_pheno$dbh > 18)]),] # can test subsetting here
#df_cls <- df_cls[which(df_cls$Object_ID %in% tree_pheno$Poly_ID[which(tree_pheno$R2 > 0.7)]),] # can test subsetting here
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus
col_sub_medians <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "median")
colnames(col_sub_medians)[1:3] <- c("month", "Object_ID", "genus")
col_sub_medians_long <- col_sub_medians %>% pivot_longer(cols = colnames(col_sub_medians)[4:ncol(col_sub_medians)])
col_sub_medians_long <- col_sub_medians_long %>% mutate(name_month = paste0(name, "_", month))
df_cls_agg <- col_sub_medians_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
  
df_cls_agg$genus <- as.factor(df_cls_agg$genus)
df_cls_agg <- df_cls_agg %>% select(!Object_ID)

df_cls_agg_cc <- df_cls_agg[complete.cases(df_cls_agg),]

train_index <- sample(1:nrow(df_cls_agg_cc), round(0.8*nrow(df_cls_agg_cc), 0))
train_data <- df_cls_agg_cc[train_index, ]
test_data <- df_cls_agg_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

#
var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

library(forcats)

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")


#####

# Try doing a similar setup with other aggregation functions
# Originally doing median, also do: mean, 90% quantile, and max, and can combine and use all of these too
set.seed(14)

ginkgo$genus <- "Ginkgo"
other$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo, other)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Median
col_sub_medians <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "median", na.rm = TRUE)
colnames(col_sub_medians)[1:3] <- c("month", "Object_ID", "genus")
col_sub_medians_long <- col_sub_medians %>% pivot_longer(cols = colnames(col_sub_medians)[4:ncol(col_sub_medians)])
col_sub_medians_long <- col_sub_medians_long %>% mutate(name_month = paste0(name, "_median_", month))
df_cls_agg_median <- col_sub_medians_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

# Max
col_sub_max <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "max", na.rm = TRUE)
colnames(col_sub_max)[1:3] <- c("month", "Object_ID", "genus")
col_sub_max_long <- col_sub_max %>% pivot_longer(cols = colnames(col_sub_max)[4:ncol(col_sub_max)])
col_sub_max_long <- col_sub_max_long %>% mutate(name_month = paste0(name, "_max_", month))
df_cls_agg_max <- col_sub_max_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

# Quantile 90%
col_sub_q90 <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "quantile", probs = seq(0.9), na.rm = TRUE) # this takes much much longer than other functions
colnames(col_sub_q90)[1:3] <- c("month", "Object_ID", "genus")
col_sub_q90_long <- col_sub_q90 %>% pivot_longer(cols = colnames(col_sub_q90)[4:ncol(col_sub_q90)])
col_sub_q90_long <- col_sub_q90_long %>% mutate(name_month = paste0(name, "_q90_", month))
df_cls_agg_q90 <- col_sub_q90_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_merge <- merge(merge(merge(df_cls_agg_median, df_cls_agg_mean), df_cls_agg_max), df_cls_agg_q90)

df_cls_merge$genus <- as.factor(df_cls_merge$genus)
df_cls_merge <- df_cls_merge %>% select(!Object_ID)
df_cls_merge_cc <- df_cls_merge[complete.cases(df_cls_merge),]

train_index <- sample(1:nrow(df_cls_merge_cc), round(0.8*nrow(df_cls_merge_cc), 0))
train_data <- df_cls_merge_cc[train_index, ]
test_data <- df_cls_merge_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

#
var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#####
# Now need to set it up to organize by individual year 
# try with mean values because that worked the best

set.seed(14)

ginkgo$genus <- "Ginkgo"
other$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo, other)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
year_vals <- year(df_cls$date)
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals, year_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:4] <- c("month", "Object_ID", "genus", "year")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[5:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month, "_", year))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)
df_cls_agg_mean <- df_cls_agg_mean %>% select(!Object_ID)

#df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),]
# No complete cases... values are always missing for some trees somewhere at this time scale
# Need to use missRanger to fill in the gaps if we want to do this
start <- Sys.time()
#mrf <- missRanger(df_cls_agg_mean, seed = 14) # would take hours to run, pausing it
stop <- Sys.time()
stop - start

# trying with just a few months for efficiency
df_cls_agg_mean <- col_sub_means_long %>% filter(month %in% c(4, 10, 11)) %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)
df_cls_agg_mean <- df_cls_agg_mean %>% select(!Object_ID)

start <- Sys.time()
mrf <- missRanger(df_cls_agg_mean, seed = 14) # this took 3 h 10 m to run, not reasonable to do this regularly
stop <- Sys.time()
stop - start

#df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),] # this works OK

train_index <- sample(1:nrow(mrf), round(0.8*nrow(mrf), 0))
train_data <- mrf[train_index, ]
test_data <- mrf[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

###

test_df <- col_sub_means_long %>% filter(name == "blue_zscaled")
num_obs <- as.data.frame(matrix(NA, nrow = 6*12, ncol = 2))
colnames(num_obs) <- c("Date", "n")
i <- 1
for (yr in 2020:2025){
  for (mo in 1:12){
    num_obs$n[i] <- test_df %>% filter(month == mo & year == yr) %>% nrow()
    num_obs$Date[i] <- ymd(paste(yr, mo, "01", sep = "-"))
    i <- i + 1
  }
}

ggplot(num_obs) +
  geom_col(aes(x = Date, y = n)) +
  scale_x_date(breaks = "month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Can do this by year
# choosing all months that exceed 27000 samples, could do something procedural based on quantiles

num_obs_sub <- num_obs %>% filter(n > 27000)
#paste0(substr(as.character(df_cls$date),1,7), "-01") %in% num_obs_sub$Date
#df_cls_sub <- df_cls %>% filter(month(date) %in% month(num_obs_sub$Date) & year(date) %in% year(num_obs_sub$Date))

set.seed(14)

ginkgo$genus <- "Ginkgo"
other$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo, other)
df_cls <- df_cls %>% filter(month(date) %in% month(num_obs_sub$Date) & year(date) %in% year(num_obs_sub$Date))
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
year_vals <- year(df_cls$date)
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals, year_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:4] <- c("month", "Object_ID", "genus", "year")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[5:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month, "_", year))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)
df_cls_agg_mean <- df_cls_agg_mean %>% select(!Object_ID)

df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),]

train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")


#####

# Try combining zscaled and raw values

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/")

ginkgo_raw <- fread("ginkgo_8b_indices.csv")
ginkgo_zs <- fread("ginkgo_8b_indices_zscaled.csv")
ginkgo_merge <- merge(ginkgo_raw, ginkgo_zs)
ginkgo_merge$genus <- "Ginkgo"

other_raw <- fread("other_trees_sample_8b_indices.csv")
other_zs <- fread("other_trees_sample_8b_indices_zscaled.csv")
other_merge <- merge(other_raw, other_zs)
other_merge$genus <- "Other"

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_merge, other_merge)
df_cls$date <- as.Date(df_cls$date)
df_cls <- df_cls %>% filter(month(date) %in% month(num_obs_sub$Date) & year(date) %in% year(num_obs_sub$Date))
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
year_vals <- year(df_cls$date)
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals, year_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:4] <- c("month", "Object_ID", "genus", "year")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[5:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month, "_", year))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)
df_cls_agg_mean <- df_cls_agg_mean %>% select(!Object_ID)

df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),]

train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#####
# Do monthly means for raw (not z-score normalized) because we haven't tried that yet

set.seed(14)

ginkgo_raw$genus <- "Ginkgo"
other_raw$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo_raw, other_raw)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# # Median
# col_sub_medians <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "median", na.rm = TRUE)
# colnames(col_sub_medians)[1:3] <- c("month", "Object_ID", "genus")
# col_sub_medians_long <- col_sub_medians %>% pivot_longer(cols = colnames(col_sub_medians)[4:ncol(col_sub_medians)])
# col_sub_medians_long <- col_sub_medians_long %>% mutate(name_month = paste0(name, "_median_", month))
# df_cls_agg_median <- col_sub_medians_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

# # Mean
# col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
# colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
# col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
# col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
# df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

# # Max
# col_sub_max <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "max", na.rm = TRUE)
# colnames(col_sub_max)[1:3] <- c("month", "Object_ID", "genus")
# col_sub_max_long <- col_sub_max %>% pivot_longer(cols = colnames(col_sub_max)[4:ncol(col_sub_max)])
# col_sub_max_long <- col_sub_max_long %>% mutate(name_month = paste0(name, "_max_", month))
# df_cls_agg_max <- col_sub_max_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
# 
# Quantile 90%
col_sub_q90 <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "quantile", probs = seq(0.9), na.rm = TRUE) # this takes much much longer than other functions
colnames(col_sub_q90)[1:3] <- c("month", "Object_ID", "genus")
col_sub_q90_long <- col_sub_q90 %>% pivot_longer(cols = colnames(col_sub_q90)[4:ncol(col_sub_q90)])
col_sub_q90_long <- col_sub_q90_long %>% mutate(name_month = paste0(name, "_q90_", month))
df_cls_agg_q90 <- col_sub_q90_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

#df_cls_merge <- merge(merge(merge(df_cls_agg_median, df_cls_agg_mean), df_cls_agg_max), df_cls_agg_q90)

df_cls_merge <- df_cls_agg_q90

df_cls_merge$genus <- as.factor(df_cls_merge$genus)
df_cls_merge <- df_cls_merge %>% select(!Object_ID)
df_cls_merge_cc <- df_cls_merge[complete.cases(df_cls_merge),]

train_index <- sample(1:nrow(df_cls_merge_cc), round(0.8*nrow(df_cls_merge_cc), 0))
train_data <- df_cls_merge_cc[train_index, ]
test_data <- df_cls_merge_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

#
var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#####
# Try other temporal ranges
# Options: 2 week, 2 week sliding (week overlap), 1 week
# Could also try 1 month sliding (or 30 day) with half overlap, but not sure how much this would help

set.seed(14)

ginkgo_raw$genus <- "Ginkgo"
other_raw$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo_raw, other_raw)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)

# setup temporal values as two week ranges, with last two week range a day (or two) longer to get the end of the year
two_week_ranges <- seq(1, 365, 14)
two_week_ranges[length(two_week_ranges)] <- 367 # make this longer than possible to set this interval
temporal_vals <- rep(NA, nrow(df_cls)) # fill this up
doys <- yday(df_cls$date)
for (i in 1:(length(two_week_ranges)-1)){
  doy_range <- seq(two_week_ranges[i], two_week_ranges[i+1]-1)
  temporal_vals[which(doys %in% doy_range)] <- i
}
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(temporal_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("two_week_interval", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_interval = paste0(name, "_mean_", two_week_interval))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_interval, values_from = value)

df_cls_merge <- df_cls_agg_mean

df_cls_merge$genus <- as.factor(df_cls_merge$genus)
df_cls_merge <- df_cls_merge %>% select(!Object_ID)
df_cls_merge_cc <- df_cls_merge[complete.cases(df_cls_merge),]

train_index <- sample(1:nrow(df_cls_merge_cc), round(0.8*nrow(df_cls_merge_cc), 0))
train_data <- df_cls_merge_cc[train_index, ]
test_data <- df_cls_merge_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

#
var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#####
# Try one week, then switch to multiple years

set.seed(14)

ginkgo_raw$genus <- "Ginkgo"
other_raw$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo_raw, other_raw)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)

# setup temporal values as two week ranges, with last two week range a day (or two) longer to get the end of the year
one_week_ranges <- seq(1, 365, 7)
one_week_ranges[length(one_week_ranges)] <- 367 # make this longer than possible to set this interval
temporal_vals <- rep(NA, nrow(df_cls)) # fill this up
doys <- yday(df_cls$date)
for (i in 1:(length(one_week_ranges)-1)){
  doy_range <- seq(one_week_ranges[i], one_week_ranges[i+1]-1)
  temporal_vals[which(doys %in% doy_range)] <- i
}
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(temporal_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("two_week_interval", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_interval = paste0(name, "_mean_", two_week_interval))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_interval, values_from = value)

df_cls_merge <- df_cls_agg_mean

df_cls_merge$genus <- as.factor(df_cls_merge$genus)
df_cls_merge <- df_cls_merge %>% select(!Object_ID)
df_cls_merge_cc <- df_cls_merge[complete.cases(df_cls_merge),]

train_index <- sample(1:nrow(df_cls_merge_cc), round(0.8*nrow(df_cls_merge_cc), 0))
train_data <- df_cls_merge_cc[train_index, ]
test_data <- df_cls_merge_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

#
var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")


#####
# try two week at month-yearly resolution just to check

set.seed(14)

ginkgo_raw$genus <- "Ginkgo"
other_raw$genus <- "Other"
df_cls <- rbind.data.frame(ginkgo_raw, other_raw)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)

# setup temporal values as two week ranges, with last two week range a day (or two) longer to get the end of the year
two_week_ranges <- seq(1, 365, 14)
two_week_ranges[length(two_week_ranges)] <- 367 # make this longer than possible to set this interval
temporal_vals <- rep(NA, nrow(df_cls)) # fill this up
doys <- yday(df_cls$date)
for (i in 1:(length(two_week_ranges)-1)){
  doy_range <- seq(two_week_ranges[i], two_week_ranges[i+1]-1)
  temporal_vals[which(doys %in% doy_range)] <- i
}
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus
year_vals <- year(df_cls$date)

# Mean
col_sub_means <- aggregate(col_sub, by = list(temporal_vals, id_vals, genus_vals, year_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:4] <- c("two_week_interval", "Object_ID", "genus", "year")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[5:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_interval = paste0(name, "_mean_", two_week_interval, "_", year))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_interval, values_from = value)

range_len <- length(two_week_ranges) - 1
test_df <- col_sub_means_long %>% filter(name == "blue")
num_obs <- as.data.frame(matrix(NA, nrow = range_len*12, ncol = 4))
colnames(num_obs) <- c("Date", "n", "year", "two_week_range")

i <- 1
for (yr in 2020:2025){
  for (m in 1:range_len){
    num_obs$n[i] <- test_df %>% filter(two_week_interval == m & year == yr) %>% nrow()
    num_obs$Date[i] <- as.Date(two_week_ranges[m] - 1, origin = paste0(yr, "-01-01"))
    num_obs$year[i] <- yr
    num_obs$two_week_range[i] <- m
    i <- i + 1
  }
}

num_obs$Date <- as.Date(num_obs$Date) # I don't know why this doesn't take within the loop

ggplot(num_obs) +
  geom_col(aes(x = Date, y = n)) +
  scale_x_date(breaks = "month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

quantile(num_obs$n, probs = seq(0, 1, 0.1), na.rm = TRUE)
# try 27000 as limit again for test
num_obs_sub <- num_obs %>% filter(n > 27000)

#df_cls <- df_cls %>% filter(month(date) %in% month(num_obs_sub$Date) & year(date) %in% year(num_obs_sub$Date))
# filtering step
df_cls$date <- as.Date(df_cls$date)
df_cls$year_ind <- paste0(year(df_cls$date), "_", temporal_vals)
num_obs_sub$year_ind <- paste0(num_obs_sub$year, "_", num_obs_sub$two_week_range)
#df_cls <- df_cls[which(year(df_cls$date) %in% num_obs_sub$year & temporal_vals %in% num_obs_sub$two_week_range),] # this is not filtering like it's supposed to, that's the issue
df_cls <- df_cls[which(df_cls$year_ind %in% num_obs_sub$year_ind),]
#col_sub <- df_cls %>% select(!date & !Object_ID & !genus) # this is now a subsetted df_cls
col_sub <- df_cls %>% select(!date & !Object_ID & !genus & !year_ind) # this is now a subsetted df_cls

# setup temporal values as two week ranges, with last two week range a day (or two) longer to get the end of the year
two_week_ranges <- seq(1, 365, 14)
two_week_ranges[length(two_week_ranges)] <- 367 # make this longer than possible to set this interval
temporal_vals <- rep(NA, nrow(df_cls)) # fill this up
doys <- yday(df_cls$date)
for (i in 1:(length(two_week_ranges)-1)){
  doy_range <- seq(two_week_ranges[i], two_week_ranges[i+1]-1)
  temporal_vals[which(doys %in% doy_range)] <- i
}
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus
year_vals <- year(df_cls$date)

# Mean
col_sub_means <- aggregate(col_sub, by = list(temporal_vals, id_vals, genus_vals, year_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:4] <- c("two_week_interval", "Object_ID", "genus", "year")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[5:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_interval = paste0(name, "_mean_", two_week_interval, "_", year))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_interval, values_from = value)

df_cls_merge <- df_cls_agg_mean

df_cls_merge$genus <- as.factor(df_cls_merge$genus)
df_cls_merge <- df_cls_merge %>% select(!Object_ID)
df_cls_merge_cc <- df_cls_merge[complete.cases(df_cls_merge),] # this works now

train_index <- sample(1:nrow(df_cls_merge_cc), round(0.8*nrow(df_cls_merge_cc), 0))
train_data <- df_cls_merge_cc[train_index, ]
test_data <- df_cls_merge_cc[-train_index, ]

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#####

# Setup monthly mean classification, and then be able to reference back to Object_ID (Poly_ID) in the original data frame and see where we're getting predictions right and wrong

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv") %>% as.data.frame()

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/")

ginkgo_raw <- fread("ginkgo_8b_indices.csv")
#ginkgo_zs <- fread("ginkgo_8b_indices_zscaled.csv")
#ginkgo_merge <- merge(ginkgo_raw, ginkgo_zs)
ginkgo_raw$genus <- "Ginkgo"

other_raw <- fread("other_trees_sample_8b_indices.csv")
#other_zs <- fread("other_trees_sample_8b_indices_zscaled.csv")
#other_merge <- merge(other_raw, other_zs)
other_raw$genus <- "Other"

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_raw, other_raw)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
#year_vals <- year(df_cls$date)
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)
df_cls_agg_mean <- df_cls_agg_mean

df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),]

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
p_rf <- predict(rf_model, test_data_noid[,2:ncol(test_data_noid)])
confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#
# Characteristics of the object IDs that were predicted correctly or incorrectly.

cor_ids <- test_data$Object_ID[which(p_rf$predictions == test_data$genus)]
inc_ids <- test_data$Object_ID[which(p_rf$predictions != test_data$genus)]

#
tree_pheno_inc_sub <- tree_pheno %>% filter(Poly_ID %in% inc_ids)

unq_first_inds_inc <- match(unique(tree_pheno_inc_sub$Poly_ID), tree_pheno_inc_sub$Poly_ID)

tree_pheno_inc_sub[unq_first_inds_inc, ] %>% 
  ggplot() +
  geom_histogram(aes(dbh))

tree_pheno_inc_sub$dbh[unq_first_inds_inc] %>% mean(na.rm = T) # 8.4 in
tree_pheno_inc_sub$dbh[unq_first_inds_inc] %>% sd(na.rm = T) # 6.0 in

tree_pheno_inc_sub$Height[unq_first_inds_inc] %>% mean(na.rm = T) # 33.1 ft
tree_pheno_inc_sub$Height[unq_first_inds_inc] %>% sd(na.rm = T) # 19.0 ft

tree_pheno_inc_sub$Radius[unq_first_inds_inc] %>% mean(na.rm = T) # 12.4 ft
tree_pheno_inc_sub$Radius[unq_first_inds_inc] %>% sd(na.rm = T) # 5.9 ft

#
tree_pheno_cor_sub <- tree_pheno %>% filter(Poly_ID %in% cor_ids)

unq_first_inds_cor <- match(unique(tree_pheno_cor_sub$Poly_ID), tree_pheno_cor_sub$Poly_ID)

tree_pheno_cor_sub[unq_first_inds_cor, ] %>% 
  ggplot() +
  geom_histogram(aes(dbh))

tree_pheno_cor_sub$dbh[unq_first_inds_cor] %>% mean(na.rm = T) # 12.0 in
tree_pheno_cor_sub$dbh[unq_first_inds_cor] %>% sd(na.rm = T) # 7.8 in
# correct trees are much larger on average
tree_pheno_cor_sub$Height[unq_first_inds_cor] %>% mean(na.rm = T) # 37.7 ft
tree_pheno_cor_sub$Height[unq_first_inds_cor] %>% sd(na.rm = T) # 18 ft
# correct trees are taller

tree_pheno_cor_sub$Radius[unq_first_inds_cor] %>% mean(na.rm = T) # 14.0 ft
tree_pheno_cor_sub$Radius[unq_first_inds_cor] %>% sd(na.rm = T) # 6.2 ft

library(sf)
boro_sf <- read_sf('/Volumes/NYC_geo/vectors/Borough Boundaries/geo_export_da133389-a6c6-45c3-a980-14295f0e4c2f.shp')

# map them
ggplot() +
  geom_sf(data = boro_sf) +
  geom_point(data = tree_pheno_cor_sub[unq_first_inds_cor,], 
             aes(x = Lon, y = Lat), col = "blue") +
  geom_point(data = tree_pheno_inc_sub[unq_first_inds_inc,], 
             aes(x = Lon, y = Lat), col = "red") +
  theme_bw()
# Doesn't seem to have any geographic pattern to the correct/incorrect of the validation points

# Is there a difference in which species are more often confused with Ginkgo?

# Not normal, stacked up at lower values, so doing rank-sum test comparison
wilcox.test(tree_pheno_cor_sub$dbh[unq_first_inds_cor], tree_pheno_inc_sub$dbh[unq_first_inds_inc]) # sig different, p < 0.001
wilcox.test(tree_pheno_cor_sub$Height[unq_first_inds_cor], tree_pheno_inc_sub$Height[unq_first_inds_inc]) # sig different, p < 0.001
wilcox.test(tree_pheno_cor_sub$Radius[unq_first_inds_cor], tree_pheno_inc_sub$Radius[unq_first_inds_inc]) # sig different, p < 0.001
#

ggplot() +
  geom_point(data = tree_pheno_cor_sub[unq_first_inds_cor,], 
             aes(x = dbh, y = Height), col = "blue") +
  geom_point(data = tree_pheno_inc_sub[unq_first_inds_inc,], 
             aes(x = dbh, y = Height), col = "red") +
  scale_x_continuous(limits = c(0, 60)) +
  scale_y_continuous(limits = c(0, 160)) +
  theme_bw()

ggplot() +
  geom_point(data = tree_pheno_cor_sub[unq_first_inds_cor,], 
             aes(x = dbh, y = Height), col = "blue") +
  scale_x_continuous(limits = c(0, 60)) +
  scale_y_continuous(limits = c(0, 160)) +
  theme_bw()

ggplot() +
  geom_point(data = tree_pheno_inc_sub[unq_first_inds_inc,], 
             aes(x = dbh, y = Height), col = "red") +
  scale_x_continuous(limits = c(0, 60)) +
  scale_y_continuous(limits = c(0, 160)) +
  theme_bw()

#####

# Which trees are getting mixed up the most with Ginkgo?

cor_trees <- tree_pheno_cor_sub[unq_first_inds_cor,]
inc_trees <- tree_pheno_inc_sub[unq_first_inds_inc,]


table(cor_trees$species)
table(inc_trees$species)


#####

# Try separating additional genera (maybe species, but more likely genus)
# Add lidar derived variables
# Add pheno variables, SOS 50 and EOS 50 for each year, and maybe SOS 20, 80, and EOS 20, 80

# Genera, State of the Urban Forest of NYC
# Street Trees: Platanus, Gleditsia, Pyrus, Quercus, Acer, Tilia, Prunus, Zelkova, Ginkgo, Styphnolobium (10)
# Parkland Trees: + Ulmus, Liquidambar, Malus (+3)
# Natural Areas: + Sassafras, Robinia, Liriodendron, Betula (+4)
# 17 total classes if we get everything, could have one additional class for evergreen conifer, but likely very small.
# From Urban forest of New York City (Nowak):
# + Thuja, Ailanthus, Myrica (no Myrica bayberry in street tree database, so can't do it)
# maybe consider Thuja all evergreen conifers, because that's likely what's going to happen...

# based on what is in the street tree db, won't likely get:
# Thuja 133, Sassafras 155, maybe Ailanthus 796

# Should consider adding:
# Fraxinus 12552

# So the core list should be:
# Platanus, Gleditsia, Pyrus, Quercus, Acer, Tilia, Prunus, Zelkova, Ginkgo, Styphnolobium, Ulmus, Liquidambar, Malus, Robinia, Liriodendron, Betula (16)
# With additional: Ailanthus, Fraxinus (+2)
# And other: Evergreen Conifer, Other Deciduous? (+2?) Could also add deciduous conifer if we want Metasequoia and Taxodium (and Larix) and be complete about the trees we have
# 20 or 21 classes total at maximum, we won't get all these but we'll see what we get
# So, start with the genus list, then add in more classes as desired

# Platanus, Gleditsia, Pyrus, Quercus, Acer, Tilia, Prunus, Zelkova, Ginkgo, Styphnolobium, Ulmus, Liquidambar, Malus, Robinia, Liriodendron, Betula, Ailanthus, Fraxinus


table(tree_pheno$genus[match(unique(other_raw$Object_ID), tree_pheno$Poly_ID)]) %>% sort(decreasing = TRUE)
# can do an initial test with the other_raw extracted data that I have already, will likely need to sample Ginkgo down to a more reasonable number.
# will need to revise other_raw genus names


# Setup
# ginkgo_raw # genus is Ginkgo
gen_lookup <- tree_pheno[match(unique(other_raw$Object_ID), tree_pheno$Poly_ID),] %>% select("Poly_ID", "genus")
other_raw_wgenus <- merge(select(other_raw, -genus), gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")
other_raw_wgenus_sub <- other_raw_wgenus %>% filter(genus %in% gen_list)

set.seed(14)
table(other_raw_wgenus_sub$genus[match(unique(other_raw_wgenus_sub$Object_ID), other_raw_wgenus_sub$Object_ID)]) # number of individual trees is messy, but 500 is probably a good starting place for n Ginkgo
ginkgo_samp_ids <- sample(unique(ginkgo_raw$Object_ID), 500)
ginkgo_raw_sub <- ginkgo_raw %>% filter(Object_ID %in% ginkgo_samp_ids)

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_raw_sub, other_raw_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
#year_vals <- year(df_cls$date)
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)
df_cls_agg_mean <- df_cls_agg_mean

df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),]

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
p_rf <- predict(rf_model, test_data_noid[,2:ncol(test_data_noid)])
confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")


#####
# See if IAV matters for this now...

test_df <- col_sub_means_long %>% filter(name == "blue")
num_obs <- as.data.frame(matrix(NA, nrow = 6*12, ncol = 2))
colnames(num_obs) <- c("Date", "n")
i <- 1
for (yr in 2020:2025){
  for (mo in 1:12){
    num_obs$n[i] <- test_df %>% filter(month == mo & year == yr) %>% nrow()
    num_obs$Date[i] <- ymd(paste(yr, mo, "01", sep = "-"))
    i <- i + 1
  }
}

num_obs$Date <- as.Date(num_obs$Date)

ggplot(num_obs) +
  geom_col(aes(x = Date, y = n)) +
  scale_x_date(breaks = "month") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

num_obs_sub <- num_obs %>% filter(n > 12500)
#paste0(substr(as.character(df_cls$date),1,7), "-01") %in% num_obs_sub$Date
#df_cls_sub <- df_cls %>% filter(month(date) %in% month(num_obs_sub$Date) & year(date) %in% year(num_obs_sub$Date))

set.seed(14)
df_cls <- rbind.data.frame(ginkgo_raw_sub, other_raw_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
df_cls <- df_cls %>% filter(month(date) %in% month(num_obs_sub$Date) & year(date) %in% year(num_obs_sub$Date))
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
year_vals <- year(df_cls$date)
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals, year_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:4] <- c("month", "Object_ID", "genus", "year")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[5:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month, "_", year))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)

df_cls_agg_mean_cc <- df_cls_agg_mean[complete.cases(df_cls_agg_mean),]
df_cls_agg_mean_cc$genus <- as.factor(df_cls_agg_mean_cc$genus)

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
p_rf <- predict(rf_model, test_data_noid[,2:ncol(test_data_noid)])
confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")


#####

gen_lookup <- tree_pheno[match(unique(other_raw$Object_ID), tree_pheno$Poly_ID),] %>% select("Poly_ID", "genus")
other_raw_wgenus <- merge(select(other_raw, -genus), gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")
other_raw_wgenus_sub <- other_raw_wgenus %>% filter(genus %in% gen_list)

set.seed(14)
table(other_raw_wgenus_sub$genus[match(unique(other_raw_wgenus_sub$Object_ID), other_raw_wgenus_sub$Object_ID)]) # number of individual trees is messy, but 500 is probably a good starting place for n Ginkgo
ginkgo_samp_ids <- sample(unique(ginkgo_raw$Object_ID), 500)
ginkgo_raw_sub <- ginkgo_raw %>% filter(Object_ID %in% ginkgo_samp_ids)

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_raw_sub, other_raw_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Do monthly mean
# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)

#add in lidar variables
lidar_vars <- tree_pheno[match(df_cls_agg_mean$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
df_cls_agg_mean_lidar <- merge(df_cls_agg_mean, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

# RF Run
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar)

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
p_rf <- predict(rf_model, test_data_noid[,2:ncol(test_data_noid)])
confusionMatrix(p_rf$predictions, reference = test_data_noid$genus)

var_imp <- as.data.frame(rf_model$variable.importance)
var_imp$name <- names(rf_model$variable.importance)
colnames(var_imp)[1] <- "importance"

var_imp_ordered <- var_imp[order(var_imp$importance, decreasing = TRUE),]

ggplot(var_imp_ordered[1:25,]) +
  geom_col(aes(x = importance, y = fct_reorder(name, importance))) +
  labs(x = "Gini importance", y = "Variable (name + month)", title = "Top 25 most important variables")

#####

# Add in SOS and EOS variables for each year


gen_lookup <- tree_pheno[match(unique(other_raw$Object_ID), tree_pheno$Poly_ID),] %>% select("Poly_ID", "genus")
other_raw_wgenus <- merge(select(other_raw, -genus), gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")
other_raw_wgenus_sub <- other_raw_wgenus %>% filter(genus %in% gen_list)

set.seed(14)
table(other_raw_wgenus_sub$genus[match(unique(other_raw_wgenus_sub$Object_ID), other_raw_wgenus_sub$Object_ID)]) # number of individual trees is messy, but 500 is probably a good starting place for n Ginkgo
ginkgo_samp_ids <- sample(unique(ginkgo_raw$Object_ID), 500)
ginkgo_raw_sub <- ginkgo_raw %>% filter(Object_ID %in% ginkgo_samp_ids)

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_raw_sub, other_raw_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Do monthly mean
# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)

#add in lidar variables
lidar_vars <- tree_pheno[match(df_cls_agg_mean$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
df_cls_agg_mean_lidar <- merge(df_cls_agg_mean, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

# Add in SOS and EOS variables
pheno_vars <- tree_pheno %>% filter(Poly_ID %in% df_cls_agg_mean$Object_ID) %>% select("Poly_ID", "Year", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)

# Just a few missing values due to infeasible pheno fits, missRanger fills these reasonably quickly (few short minutes)
mr_pheno_vars_wide <- missRanger(pheno_vars_wide)
pheno_vars_wide_filled <- round(mr_pheno_vars_wide, 0)

#
df_cls_agg_mean_lidar_pheno <- merge(df_cls_agg_mean_lidar, pheno_vars_wide_filled, by.x = "Object_ID", by.y = "Poly_ID")

# # Filter trees based on size
# # Dan does >= 10 cm (~4 in) DBH, and >= 5 m (16.4 ft) tall
# size_vars <- tree_pheno[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "dbh")
# size_vars_sub <- size_vars %>% filter(dbh > 4 & Height > 16.4)
# df_cls_agg_mean_lidar_pheno <- df_cls_agg_mean_lidar_pheno %>% filter(Object_ID %in% size_vars_sub$Poly_ID)

# RF Run
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno)

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
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

#####
# 
local_imp <- rf_model$variable.importance.local
#https://stackoverflow.com/questions/61340327/how-to-obtain-feature-importance-by-class-using-ranger
#as.data.table(rf.iris$variable.importance.local)[,Species := iris$Species][,lapply(.SD,mean),by=Species]
rf_local_imp_summary <- as.data.table(rf_model$variable.importance.local)[,Genus := train_data_noid$genus][,lapply(.SD,mean),by=Genus]

rf_local_imp_summary_long <- rf_local_imp_summary %>% pivot_longer(cols = 2:ncol(rf_local_imp_summary))

genus_list <- unique(rf_local_imp_summary_long$Genus)
for (i in 1:length(genus_list)){
}

i <- 1
genus_inds <- which(rf_local_imp_summary_long$Genus == genus_list[i])
frank(rf_local_imp_summary_long$value[genus_inds])

# would need to filter these for top 25
ggplot(rf_local_imp_summary_long) +
  geom_point(aes(x = name, y = value)) + 
  facet_wrap(~Genus)
# looks like the importance is mostly tuned towards Platanus and less so for other variables

#####

# Swap in normalized data

ginkgo_zs <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/ginkgo_8b_indices_zscaled.csv")
other_zs <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/other_trees_sample_8b_indices_zscaled.csv")

ginkgo_zs$genus <- "Ginkgo"

gen_lookup <- tree_pheno[match(unique(other_zs$Object_ID), tree_pheno$Poly_ID),] %>% select("Poly_ID", "genus")
other_zs_wgenus <- merge(other_zs, gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")
other_zs_wgenus_sub <- other_zs_wgenus %>% filter(genus %in% gen_list)

set.seed(14)
table(other_zs_wgenus_sub$genus[match(unique(other_zs_wgenus_sub$Object_ID), other_zs_wgenus_sub$Object_ID)]) # number of individual trees is messy, but 500 is probably a good starting place for n Ginkgo
ginkgo_samp_ids <- sample(unique(ginkgo_zs$Object_ID), 500)
ginkgo_zs_sub <- ginkgo_zs %>% filter(Object_ID %in% ginkgo_samp_ids)

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_zs_sub, other_zs_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Do monthly mean
# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)

#add in lidar variables
lidar_vars <- tree_pheno[match(df_cls_agg_mean$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
df_cls_agg_mean_lidar <- merge(df_cls_agg_mean, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

# Add in SOS and EOS variables
pheno_vars <- tree_pheno %>% filter(Poly_ID %in% df_cls_agg_mean$Object_ID) %>% select("Poly_ID", "Year", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)

# Just a few missing values due to infeasible pheno fits, missRanger fills these reasonably quickly (few short minutes)
mr_pheno_vars_wide <- missRanger(pheno_vars_wide)
pheno_vars_wide_filled <- round(mr_pheno_vars_wide, 0)

#
df_cls_agg_mean_lidar_pheno <- merge(df_cls_agg_mean_lidar, pheno_vars_wide_filled, by.x = "Object_ID", by.y = "Poly_ID")

# # Filter trees based on size
# # Dan does >= 10 cm (~4 in) DBH, and >= 5 m (16.4 ft) tall
# size_vars <- tree_pheno[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "dbh")
# size_vars_sub <- size_vars %>% filter(dbh > 4 & Height > 16.4)
# df_cls_agg_mean_lidar_pheno <- df_cls_agg_mean_lidar_pheno %>% filter(Object_ID %in% size_vars_sub$Poly_ID)

# RF Run
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno)

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
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

######
# Testing for filtering of source

ginkgo_zs <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/ginkgo_8b_indices_zscaled.csv")
other_zs <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/other_trees_sample_8b_indices_zscaled.csv")

ginkgo_zs$genus <- "Ginkgo"

gen_lookup <- tree_pheno[match(unique(other_zs$Object_ID), tree_pheno$Poly_ID),] %>% select("Poly_ID", "genus")
other_zs_wgenus <- merge(other_zs, gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")
other_zs_wgenus_sub <- other_zs_wgenus %>% filter(genus %in% gen_list)

set.seed(14)
table(other_zs_wgenus_sub$genus[match(unique(other_zs_wgenus_sub$Object_ID), other_zs_wgenus_sub$Object_ID)]) # number of individual trees is messy, but 500 is probably a good starting place for n Ginkgo
ginkgo_samp_ids <- sample(unique(ginkgo_zs$Object_ID), 500)
ginkgo_zs_sub <- ginkgo_zs %>% filter(Object_ID %in% ginkgo_samp_ids)

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_zs_sub, other_zs_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Do monthly mean
# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)

#add in lidar variables
lidar_vars <- tree_pheno[match(df_cls_agg_mean$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
df_cls_agg_mean_lidar <- merge(df_cls_agg_mean, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

# Add in SOS and EOS variables
pheno_vars <- tree_pheno %>% filter(Poly_ID %in% df_cls_agg_mean$Object_ID) %>% select("Poly_ID", "Year", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)

# Just a few missing values due to infeasible pheno fits, missRanger fills these reasonably quickly (few short minutes)
mr_pheno_vars_wide <- missRanger(pheno_vars_wide)
pheno_vars_wide_filled <- round(mr_pheno_vars_wide, 0)

#
df_cls_agg_mean_lidar_pheno <- merge(df_cls_agg_mean_lidar, pheno_vars_wide_filled, by.x = "Object_ID", by.y = "Poly_ID")

# Filter trees based on size
# Dan does >= 10 cm (~4 in) DBH, and >= 5 m (16.4 ft) tall
#size_vars <- tree_pheno[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "dbh")
#size_vars_sub <- size_vars %>% filter(dbh > 4 & Height > 16.4)
#df_cls_agg_mean_lidar_pheno <- df_cls_agg_mean_lidar_pheno %>% filter(Object_ID %in% size_vars_sub$Poly_ID)

size_vars <- tree_pheno[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "dbh", "tpstructur", "tpconditio")
size_vars_sub <- size_vars %>% filter(dbh > 4 & tpstructur %in% c("Full")) # & tpconditio %in% c("Excellent", "Good")) #dbh > 4) #) # tpstructur %in% c("Full")) # & tpconditio %in% c("Excellent", "Good", "Fair", "Poor")) & dbh > 4 & Height > 16.4 
df_cls_agg_mean_lidar_pheno <- df_cls_agg_mean_lidar_pheno %>% filter(Object_ID %in% size_vars_sub$Poly_ID)

# RF Run
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno)

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
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

#####

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv") %>% as.data.frame()

table(tree_pheno$genus)

#####


ginkgo_zs <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/ginkgo_8b_indices_zscaled.csv")
other_zs <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/other_trees_sample_8b_indices_zscaled.csv")

ginkgo_raw <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/ginkgo_8b_indices.csv")
other_raw <- fread("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/other_trees_sample_8b_indices.csv")

ginkgo_zs$genus <- "Ginkgo"
ginkgo_raw$genus <- "Ginkgo"

gen_lookup <- tree_pheno[match(unique(other_zs$Object_ID), tree_pheno$Poly_ID),] %>% select("Poly_ID", "genus")
other_zs_wgenus <- merge(other_zs, gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")
other_zs_wgenus_sub <- other_zs_wgenus %>% filter(genus %in% gen_list)

other_raw_wgenus <- merge(other_raw, gen_lookup, by.x = "Object_ID", by.y = "Poly_ID")
other_raw_wgenus_sub <- other_raw_wgenus %>% filter(genus %in% gen_list)

set.seed(14)
table(other_zs_wgenus_sub$genus[match(unique(other_zs_wgenus_sub$Object_ID), other_zs_wgenus_sub$Object_ID)]) # number of individual trees is messy, but 500 is probably a good starting place for n Ginkgo
ginkgo_samp_ids <- sample(unique(ginkgo_zs$Object_ID), 500)
ginkgo_zs_sub <- ginkgo_zs %>% filter(Object_ID %in% ginkgo_samp_ids)
ginkgo_raw_sub <- ginkgo_raw %>% filter(Object_ID %in% ginkgo_samp_ids)

ginkgo_merge_sub <- cbind.data.frame(ginkgo_zs_sub, ginkgo_raw_sub[,3:(ncol(ginkgo_raw_sub)-1)])
other_merge_wgenus_sub <- cbind.data.frame(other_zs_wgenus_sub, other_raw_wgenus_sub[,3:(ncol(other_raw_wgenus_sub)-1)])

# try for years, need num_obs_sub from before for same month-year samples
df_cls <- rbind.data.frame(ginkgo_merge_sub, other_merge_wgenus_sub)
df_cls$date <- as.Date(df_cls$date)
col_sub <- df_cls %>% select(!date & !Object_ID & !genus)
month_vals <- month(df_cls$date) # can revise this part of the aggregation step
id_vals <- df_cls$Object_ID
genus_vals <- df_cls$genus

# Do monthly mean
# Mean
col_sub_means <- aggregate(col_sub, by = list(month_vals, id_vals, genus_vals), FUN = "mean", na.rm = TRUE)
colnames(col_sub_means)[1:3] <- c("month", "Object_ID", "genus")
col_sub_means_long <- col_sub_means %>% pivot_longer(cols = colnames(col_sub_means)[4:ncol(col_sub_means)])
col_sub_means_long <- col_sub_means_long %>% mutate(name_month = paste0(name, "_mean_", month))
df_cls_agg_mean <- col_sub_means_long %>% pivot_wider(id_cols = c(Object_ID, genus), names_from = name_month, values_from = value)
df_cls_agg_mean$genus <- as.factor(df_cls_agg_mean$genus)

#add in lidar variables
lidar_vars <- tree_pheno[match(df_cls_agg_mean$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
df_cls_agg_mean_lidar <- merge(df_cls_agg_mean, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

# Add in SOS and EOS variables
pheno_vars <- tree_pheno %>% filter(Poly_ID %in% df_cls_agg_mean$Object_ID) %>% select("Poly_ID", "Year", "SOS_50", "EOS_50", "SOS_20", "EOS_20", "SOS_80", "EOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)

# Just a few missing values due to infeasible pheno fits, missRanger fills these reasonably quickly (few short minutes)
mr_pheno_vars_wide <- missRanger(pheno_vars_wide)
pheno_vars_wide_filled <- round(mr_pheno_vars_wide, 0)

#
df_cls_agg_mean_lidar_pheno <- merge(df_cls_agg_mean_lidar, pheno_vars_wide_filled, by.x = "Object_ID", by.y = "Poly_ID")

# Filter trees based on size and quality
size_vars <- tree_pheno[match(df_cls_agg_mean_lidar_pheno$Object_ID, tree_pheno$Poly_ID),] %>% select("Poly_ID", "Height", "dbh", "tpstructur", "tpconditio")
size_vars_sub <- size_vars %>% filter(dbh > 4 & tpstructur %in% c("Full"))
df_cls_agg_mean_lidar_pheno <- df_cls_agg_mean_lidar_pheno %>% filter(Object_ID %in% size_vars_sub$Poly_ID)

# RF Run
df_cls_agg_mean_cc <- na.omit(df_cls_agg_mean_lidar_pheno)

set.seed(14)
train_index <- sample(1:nrow(df_cls_agg_mean_cc), round(0.8*nrow(df_cls_agg_mean_cc), 0))
train_data <- df_cls_agg_mean_cc[train_index, ]
test_data <- df_cls_agg_mean_cc[-train_index, ]

train_data_noid <- train_data %>% select(!Object_ID)
test_data_noid <- test_data %>% select(!Object_ID)

rf_model <- ranger(genus ~ ., train_data_noid, importance = "impurity")
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