# Number of trees by date

#####
# check extracted number of trees

library(tidyverse)
library(purrr)
library(data.table)
'%notin%' <- Negate('%in%')

# 4 band
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
#names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 27) %>% ymd() # name the list
# names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 17, end = 26)  # name the list
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

unq_names <- unique(names(file_list))

tree_counts <- cbind.data.frame(unq_names, rep(NA, length(unq_names)))
colnames(tree_counts) <- c("Date", "num_trees")

for (i in 1:length(unq_names)){
  print(i)
  #i <- 1
  file_list_sub <- file_list[which(names(file_list) == unq_names[i])]
  tree_spectra <- purrr::map_df(file_list_sub, fread, .id = 'date')
  tree_counts$num_trees[i] <- nrow(tree_spectra)
}

tree_counts_4b <- tree_counts
tree_counts_4b$n_bands <- "4band"

#
# 8 band
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

unq_names <- unique(names(file_list))

tree_counts <- cbind.data.frame(unq_names, rep(NA, length(unq_names)))
colnames(tree_counts) <- c("Date", "num_trees")

for (i in 1:length(unq_names)){
  print(i)
  #i <- 1
  file_list_sub <- file_list[which(names(file_list) == unq_names[i])]
  tree_spectra <- purrr::map_df(file_list_sub, fread, .id = 'date')
  tree_counts$num_trees[i] <- nrow(tree_spectra)
}

tree_counts_8b <- tree_counts
tree_counts_8b$n_bands <- "8band"

#####

colnames(tree_counts_4b)[2] <- "n_trees_4b"
colnames(tree_counts_8b)[2] <- "n_trees_8b"

tree_counts_merged <- full_join(tree_counts_4b[,1:2], tree_counts_8b[,1:2])
tree_counts_merged_ordered <- tree_counts_merged[order(tree_counts_merged$Date),]

tree_counts_merged_ordered$Date <- ymd(tree_counts_merged_ordered$Date)

tree_counts_merged_ordered_long <- tree_counts_merged_ordered %>% pivot_longer(cols = c("n_trees_4b", "n_trees_8b"), names_to = "n_trees_source")
tree_counts_merged_ordered_long$Year <- year(tree_counts_merged_ordered_long$Date)

tree_counts_merged_ordered_long$Year <- as.factor(tree_counts_merged_ordered_long$Year)

ggplot(tree_counts_merged_ordered_long) +
  geom_hline(yintercept = 2127625) +
  geom_point(aes(x = Date, y = value, alpha = value, color = n_trees_source)) +
  scale_x_date(breaks = "1 month", date_labels = "%m") +
  facet_wrap(~Year, scales = "free_x", nrow = 5) +
  theme_bw()
setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
ggsave("tree_counts_fullset_daily.png", width = 12, height = 10, units = "in")


#####

# Repeat, but filter for the trees in what we can call the Maximum Set of Elegant Trees (MSET)
library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

# will need to load this with purrr
#tree_spectra <- fread("/Volumes/NYC_geo/tree_classification/extracted_8band_traintest/zscaled/tree_points_8b_spectra_zscaled_monthmean_idset1.csv")
# setwd("/Volumes/NYC_geo/tree_classification/extracted_8band_traintest/zscaled/")
# file_list <- list.files()
# tree_spectra <- purrr::map_df(file_list, fread)

gen_list <- c("Platanus", "Gleditsia", "Pyrus", "Quercus", "Acer", "Tilia", 
              "Prunus", "Zelkova", "Ginkgo", "Styphnolobium", "Ulmus", 
              "Liquidambar", "Malus", "Robinia", "Liriodendron", "Betula", "Ailanthus", "Fraxinus")

# Filter based on minimum summer NDVI threshold for all years
setwd('/Volumes/NYC_geo/tree_mortality')
mean_ndvi_all <- fread("mean_summer_ndvi.csv")
ndvi_min <- 0.3 # Thapa
#ndvi_min <- 0.5 #  Alonzo
#mean_ndvi_all_live <- mean_ndvi_all %>% filter(ndvi_2024 > ndvi_min) # accidentally tested this first, will do it the other way next

# drop 2017
mean_ndvi_all_live <- mean_ndvi_all %>% select(!ndvi_2017) %>% filter_at(vars(starts_with("ndvi")), all_vars(. > ndvi_min))

# combine everything to setup an RF run like previous organization
#tree_pheno_sub <- tree_pheno %>% filter(dbh > 4 & tpstructur == "Full" & genus %in% gen_list & Poly_ID %in% unique(tree_spectra$Object_ID))
tree_pheno_sub <- tree_pheno %>% filter(dbh > 4 & tpstructur == "Full" & genus %in% gen_list & Poly_ID %in% unique(mean_ndvi_all_live$Object_ID)) 

#add in lidar variables - not including these, using our own crown metrics
#lidar_vars <- tree_pheno_sub[match(tree_spectra$Object_ID, tree_pheno_sub$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
#lidar_vars <- tree_pheno_sub[match(unique(tree_pheno_sub$Poly_ID), tree_pheno_sub$Poly_ID),] %>% select("Poly_ID", "Height", "Radius", "SHAPE_Length", "SHAPE_Area")
#df_cls_agg_mean_lidar <- merge(tree_spectra, lidar_vars, by.x = "Object_ID", by.y = "Poly_ID")

#pheno_vars <- tree_pheno_sub %>% filter(Poly_ID %in% tree_spectra$Object_ID) %>% select("Poly_ID", "Year", "SOS_50", "SOS_20", "SOS_80")
pheno_vars <- tree_pheno_sub %>% select("Poly_ID", "Year", "SOS_50", "SOS_20", "SOS_80")
pheno_vars_long <- pheno_vars %>% pivot_longer(cols = colnames(pheno_vars)[3:ncol(pheno_vars)])
pheno_vars_long <- pheno_vars_long %>% mutate(name_year = paste0(name, "_", Year))
pheno_vars_long$value <- as.numeric(pheno_vars_long$value)
pheno_vars_wide <- pheno_vars_long %>% pivot_wider(id_cols = c(Poly_ID), names_from = name_year, values_from = value)
#pheno_vars_wide_sub <- na.omit(pheno_vars_wide)

# add in variables from our lidar crown metrics
# later version, should drop variables from tree polygon variables (probably)
setwd("/Volumes/NYC_geo/nyc_lidar_metrics/crown_metrics")
crown_metrics_file_list <- list.files()
tree_crown_metrics <- purrr::map_df(crown_metrics_file_list, fread)

pheno_vars_wmetrics <- merge(pheno_vars_wide, tree_crown_metrics, by = "Poly_ID")

mset_vars <- na.omit(pheno_vars_wmetrics)

mset_genus <- tree_pheno[match(mset_vars$Poly_ID, tree_pheno$Poly_ID), c("Poly_ID", "genus")]
mset <- merge(mset_genus, mset_vars, by = "Poly_ID")
unique(mset$genus) %>% sort()
mset %>% na.omit() %>% nrow() # 243675, same length
table(mset$genus) # counts for each genus

#####


#####

# 4 band
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

unq_names <- unique(names(file_list))

tree_counts <- cbind.data.frame(unq_names, rep(NA, length(unq_names)))
colnames(tree_counts) <- c("Date", "num_trees")

for (i in 1:length(unq_names)){
  print(i)
  #i <- 1
  file_list_sub <- file_list[which(names(file_list) == unq_names[i])]
  tree_spectra <- purrr::map_df(file_list_sub, fread, .id = 'date')
  tree_spectra_sub <- tree_spectra %>% filter(Object_ID %in% mset$Poly_ID)# filtering for mset
  tree_counts$num_trees[i] <- nrow(tree_spectra_sub)
}

tree_counts_4b <- tree_counts
tree_counts_4b$n_bands <- "4band"

#
# 8 band
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

unq_names <- unique(names(file_list))

tree_counts <- cbind.data.frame(unq_names, rep(NA, length(unq_names)))
colnames(tree_counts) <- c("Date", "num_trees")

for (i in 1:length(unq_names)){
  print(i)
  #i <- 1
  file_list_sub <- file_list[which(names(file_list) == unq_names[i])]
  tree_spectra <- purrr::map_df(file_list_sub, fread, .id = 'date')
  tree_spectra_sub <- tree_spectra %>% filter(Object_ID %in% mset$Poly_ID)# filtering for mset
  tree_counts$num_trees[i] <- nrow(tree_spectra_sub)
}

tree_counts_8b <- tree_counts
tree_counts_8b$n_bands <- "8band"

# Merging
colnames(tree_counts_4b)[2] <- "n_trees_4b"
colnames(tree_counts_8b)[2] <- "n_trees_8b"

tree_counts_merged <- full_join(tree_counts_4b[,1:2], tree_counts_8b[,1:2])
tree_counts_merged_ordered <- tree_counts_merged[order(tree_counts_merged$Date),]

tree_counts_merged_ordered$Date <- ymd(tree_counts_merged_ordered$Date)

tree_counts_merged_ordered_long <- tree_counts_merged_ordered %>% pivot_longer(cols = c("n_trees_4b", "n_trees_8b"), names_to = "n_trees_source")
tree_counts_merged_ordered_long$Year <- year(tree_counts_merged_ordered_long$Date)

tree_counts_merged_ordered_long$Year <- as.factor(tree_counts_merged_ordered_long$Year)

ggplot(tree_counts_merged_ordered_long) +
  geom_hline(yintercept = 243675) +
  geom_point(aes(x = Date, y = value, alpha = value, color = n_trees_source)) +
  scale_x_date(breaks = "1 month", date_labels = "%m") +
  facet_wrap(~Year, scales = "free_x", nrow = 5) +
  theme_bw()
setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
ggsave("tree_counts_mset_daily.png", width = 12, height = 10, units = "in")


#####



#####

# Now do summarizing by date ranges (4 band)

# Monthly summaries
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as month-year combinations
unq_names <- seq(ymd('2018-01-01'), ymd('2024-12-31'), by = "1 month")

mset_month <- cbind.data.frame(mset_genus, as.data.frame(matrix(NA, nrow = nrow(mset_genus), ncol = length(unq_names))))
colnames(mset_month) <- c(colnames(mset_genus), as.character(unq_names))

for (i in 1:length(unq_names)){
  #i <- 1
  print(unq_names[i])
  time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i]) + months(1) - days(1)) # a one month time span
  file_inds <- which(ymd(names(file_list)) %within% time_span)
  if (length(file_inds) > 0){
    tree_spectra <- purrr::map_df(file_list[file_inds], fread, .id = 'date') %>% filter(Object_ID %in% mset_genus$Poly_ID)
    
    tree_obs_counts <- as.data.frame(table(tree_spectra$Object_ID))
    colnames(tree_obs_counts)[1] <- "Poly_ID"
    tree_obs_counts$Poly_ID <- as.numeric(as.character(tree_obs_counts$Poly_ID))
    mset_joined_span <- full_join(mset_genus, tree_obs_counts) # check if factor indexing doesn't get changed 
    
    mset_month[, i+2] <- mset_joined_span$Freq
  } else {
    print("No extracted data in time range")
  }
}

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
write.csv(mset_month, "mset_counts_month.csv", row.names = FALSE)

# Note: no extracted trees in Feb 2021, poor images throughout that month...

# Can do all month related outputs from this resulting data frame

# Composite across years
mset_month_colinds <- c(NA, NA, month(ymd(colnames(mset_month)[3:ncol(mset_month)])))
mset_month_composite <- cbind.data.frame(mset_month[,1:2], as.data.frame(matrix(NA, nrow = nrow(mset_month), ncol = 12)))
colnames(mset_month_composite)[3:14] <- paste0("month_", 1:12)
for (i in 1:12){
  mset_month_composite[,i+2] <- rowSums(mset_month[,which(mset_month_colinds == i)], na.rm = TRUE)
}

mset_month_composite %>% na.omit %>% nrow()
table(na.omit(mset_month_composite)$genus)
# no drop off in poly ids when we do monthly composites across years

# Drop February for all years
mset_month_nofeb <- cbind.data.frame(mset_month[,1:2], mset_month[,which(mset_month_colinds != 2)]) # NAs are omitted from which()!
mset_month_nofeb %>% na.omit %>% nrow() # 169657
table(na.omit(mset_month_nofeb)$genus)

# Drop January and February for all years
mset_month_nojf <- cbind.data.frame(mset_month[,1:2], mset_month[,which(mset_month_colinds %notin% 1:2)]) # NAs are omitted from which()!
mset_month_nojf %>% na.omit %>% nrow()
table(na.omit(mset_month_nojf)$genus)

# Drop December, January, and February for all years
mset_month_nodjf <- cbind.data.frame(mset_month[,1:2], mset_month[,which(mset_month_colinds %notin% c(1, 2, 12))]) # NAs are omitted from which()!
mset_month_nodjf %>% na.omit %>% nrow()
table(na.omit(mset_month_nodjf)$genus)

#####



#####

# 2 week summaries
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as month-year combinations
# Ending at 12-18 to avoid last day of year indexing, last "2 week" period is slightly longer
unq_names <- c(seq(ymd('2018-01-01'), ymd('2018-12-18'), by = "2 weeks"),
              seq(ymd('2019-01-01'), ymd('2019-12-18'), by = "2 weeks"),
              seq(ymd('2020-01-01'), ymd('2020-12-18'), by = "2 weeks"),
              seq(ymd('2021-01-01'), ymd('2021-12-18'), by = "2 weeks"),
              seq(ymd('2022-01-01'), ymd('2022-12-18'), by = "2 weeks"),
              seq(ymd('2023-01-01'), ymd('2023-12-18'), by = "2 weeks"),
              seq(ymd('2024-01-01'), ymd('2024-12-18'), by = "2 weeks"),
              ymd('2025-01-01')) # keep this last one just for indexing, then drop after aggregation loop
#unq_names <- paste0(year(unq_names), "_", yday(unq_names))
# do composite with yday()


mset_2wk <- cbind.data.frame(mset_genus, as.data.frame(matrix(NA, nrow = nrow(mset_genus), ncol = length(unq_names))))
colnames(mset_2wk) <- c(colnames(mset_genus), as.character(unq_names))

for (i in 1:(length(unq_names)-1)){
  #i <- 1
  print(unq_names[i])
  time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i+1])) # an 2 week time span
  file_inds <- which(ymd(names(file_list)) %within% time_span)
  if (length(file_inds) > 0){
    tree_spectra <- purrr::map_df(file_list[file_inds], fread, .id = 'date') %>% filter(Object_ID %in% mset_genus$Poly_ID)
    
    tree_obs_counts <- as.data.frame(table(tree_spectra$Object_ID))
    colnames(tree_obs_counts)[1] <- "Poly_ID"
    tree_obs_counts$Poly_ID <- as.numeric(as.character(tree_obs_counts$Poly_ID))
    mset_joined_span <- full_join(mset_genus, tree_obs_counts) # check if factor indexing doesn't get changed 
    
    mset_2wk[, i+2] <- mset_joined_span$Freq
  } else {
    print("No extracted data in time range")
  }
}

mset_2wk <- mset_2wk %>% select(!`2025-01-01`)

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
write.csv(mset_2wk, "mset_counts_2week.csv", row.names = FALSE)

#####

# 1 week time range
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as month-year combinations
# Ending at 12-27 to avoid last day of year indexing, last "1 week" period is slightly longer
unq_names <- c(seq(ymd('2018-01-01'), ymd('2018-12-27'), by = "1 week"),
               seq(ymd('2019-01-01'), ymd('2019-12-27'), by = "1 week"),
               seq(ymd('2020-01-01'), ymd('2020-12-27'), by = "1 week"),
               seq(ymd('2021-01-01'), ymd('2021-12-27'), by = "1 week"),
               seq(ymd('2022-01-01'), ymd('2022-12-27'), by = "1 week"),
               seq(ymd('2023-01-01'), ymd('2023-12-27'), by = "1 week"),
               seq(ymd('2024-01-01'), ymd('2024-12-27'), by = "1 week"),
               ymd('2025-01-01')) # keep this last one just for indexing, then drop after aggregation loop

mset_1wk <- cbind.data.frame(mset_genus, as.data.frame(matrix(NA, nrow = nrow(mset_genus), ncol = length(unq_names))))
colnames(mset_1wk) <- c(colnames(mset_genus), as.character(unq_names))

for (i in 1:(length(unq_names)-1)){
  #i <- 1
  print(unq_names[i])
  time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i+1])) # an 1 week time span
  file_inds <- which(ymd(names(file_list)) %within% time_span)
  if (length(file_inds) > 0){
    tree_spectra <- purrr::map_df(file_list[file_inds], fread, .id = 'date') %>% filter(Object_ID %in% mset_genus$Poly_ID)
    
    tree_obs_counts <- as.data.frame(table(tree_spectra$Object_ID))
    colnames(tree_obs_counts)[1] <- "Poly_ID"
    tree_obs_counts$Poly_ID <- as.numeric(as.character(tree_obs_counts$Poly_ID))
    mset_joined_span <- full_join(mset_genus, tree_obs_counts) # check if factor indexing doesn't get changed 
    
    mset_1wk[, i+2] <- mset_joined_span$Freq
  } else {
    print("No extracted data in time range")
  }
}

mset_1wk <- mset_1wk %>% select(!`2025-01-01`)

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
write.csv(mset_1wk, "mset_counts_1week.csv", row.names = FALSE)

#####

# Compare counts 
setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
mset_month <- fread("mset_counts_4b_month.csv")
mset_2wk <- fread("mset_counts_4b_2week.csv")
mset_1wk <- fread("mset_counts_4b_1week.csv")

mset_month_1 <- mset_month[,3:ncol(mset_month)]
mset_month_1[mset_month_1 > 1] <- 1

mset_2wk_1 <- mset_2wk[,3:ncol(mset_2wk)]
mset_2wk_1[mset_2wk_1 > 1] <- 1

mset_1wk_1 <- mset_1wk[,3:ncol(mset_1wk)]
mset_1wk_1[mset_1wk_1 > 1] <- 1

mset_month_1 %>% colSums(na.rm = TRUE) %>% sort()
mset_2wk_1 %>% colSums(na.rm = TRUE) %>% sort()
mset_1wk_1 %>% colSums(na.rm = TRUE) %>% sort()

# This will be take some minimum threshold number of trees
# Get all dates with at least some value, or an intersected value

mset_month_1_sums <- mset_month_1 %>% colSums(na.rm = TRUE) %>% as.data.frame()
mset_month_1_sums$date <- rownames(mset_month_1_sums)
mset_month_1_sums$agg_span <- "Month"
colnames(mset_month_1_sums)[1] <- "ntrees"

mset_2wk_1_sums <- mset_2wk_1 %>% colSums(na.rm = TRUE) %>% as.data.frame()
mset_2wk_1_sums$date <- rownames(mset_2wk_1_sums)
mset_2wk_1_sums$agg_span <- "Two_Week"
colnames(mset_2wk_1_sums)[1] <- "ntrees"

mset_1wk_1_sums <- mset_1wk_1 %>% colSums(na.rm = TRUE) %>% as.data.frame()
mset_1wk_1_sums$date <- rownames(mset_1wk_1_sums)
mset_1wk_1_sums$agg_span <- "One_Week"
colnames(mset_1wk_1_sums)[1] <- "ntrees"

mset_compare_sums <- rbind.data.frame(mset_month_1_sums, mset_2wk_1_sums, mset_1wk_1_sums)

mset_compare_sums$date <- ymd(mset_compare_sums$date)

ggplot(mset_compare_sums) +
  geom_point(aes(x = date, y = ntrees, color = agg_span)) +
  scale_x_date(breaks = "1 year", date_labels = "%Y", date_minor_breaks = "3 months", expand = c(0.03,0.03)) +
  facet_wrap(~agg_span)
ggsave("ntree_counts_w_time_spans.png", width = 10, height = 3, units = "in")

# Just need to pick threshold, drop low observations, this remaining intersection is the best that we can do

#####

# 8 band!
# Need to do the same thing for 8 band, likely 2021 through 2024
# Monthly summaries
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as month-year combinations
unq_names <- seq(ymd('2021-01-01'), ymd('2024-12-31'), by = "1 month")

mset_month <- cbind.data.frame(mset_genus, as.data.frame(matrix(NA, nrow = nrow(mset_genus), ncol = length(unq_names))))
colnames(mset_month) <- c(colnames(mset_genus), as.character(unq_names))

for (i in 1:length(unq_names)){
  #i <- 1
  print(unq_names[i])
  time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i]) + months(1) - days(1)) # a one month time span
  file_inds <- which(ymd(names(file_list)) %within% time_span)
  if (length(file_inds) > 0){
    tree_spectra <- purrr::map_df(file_list[file_inds], fread, .id = 'date') %>% filter(Object_ID %in% mset_genus$Poly_ID)
    
    tree_obs_counts <- as.data.frame(table(tree_spectra$Object_ID))
    colnames(tree_obs_counts)[1] <- "Poly_ID"
    tree_obs_counts$Poly_ID <- as.numeric(as.character(tree_obs_counts$Poly_ID))
    mset_joined_span <- full_join(mset_genus, tree_obs_counts) # check if factor indexing doesn't get changed 
    
    mset_month[, i+2] <- mset_joined_span$Freq
  } else {
    print("No extracted data in time range")
  }
}

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
write.csv(mset_month, "mset_counts_8b_month.csv", row.names = FALSE)

# 2 week summaries
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as month-year combinations
# Ending at 12-18 to avoid last day of year indexing, last "2 week" period is slightly longer
unq_names <- c(seq(ymd('2021-01-01'), ymd('2021-12-18'), by = "2 weeks"),
               seq(ymd('2022-01-01'), ymd('2022-12-18'), by = "2 weeks"),
               seq(ymd('2023-01-01'), ymd('2023-12-18'), by = "2 weeks"),
               seq(ymd('2024-01-01'), ymd('2024-12-18'), by = "2 weeks"),
               ymd('2025-01-01')) # keep this last one just for indexing, then drop after aggregation loop
#unq_names <- paste0(year(unq_names), "_", yday(unq_names))
# do composite with yday()


mset_2wk <- cbind.data.frame(mset_genus, as.data.frame(matrix(NA, nrow = nrow(mset_genus), ncol = length(unq_names))))
colnames(mset_2wk) <- c(colnames(mset_genus), as.character(unq_names))

for (i in 1:(length(unq_names)-1)){
  #i <- 1
  print(unq_names[i])
  time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i+1])) # an 2 week time span
  file_inds <- which(ymd(names(file_list)) %within% time_span)
  if (length(file_inds) > 0){
    tree_spectra <- purrr::map_df(file_list[file_inds], fread, .id = 'date') %>% filter(Object_ID %in% mset_genus$Poly_ID)
    
    tree_obs_counts <- as.data.frame(table(tree_spectra$Object_ID))
    colnames(tree_obs_counts)[1] <- "Poly_ID"
    tree_obs_counts$Poly_ID <- as.numeric(as.character(tree_obs_counts$Poly_ID))
    mset_joined_span <- full_join(mset_genus, tree_obs_counts) # check if factor indexing doesn't get changed 
    
    mset_2wk[, i+2] <- mset_joined_span$Freq
  } else {
    print("No extracted data in time range")
  }
}

mset_2wk <- mset_2wk %>% select(!`2025-01-01`)

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
write.csv(mset_2wk, "mset_counts_8b_2week.csv", row.names = FALSE)

# 1 week time range
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
file_list <- list.files()
names(file_list) <- strsplit(file_list, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 26)  # name the list

# setup unique names as month-year combinations
# Ending at 12-27 to avoid last day of year indexing, last "1 week" period is slightly longer
unq_names <- c(seq(ymd('2021-01-01'), ymd('2021-12-27'), by = "1 week"),
               seq(ymd('2022-01-01'), ymd('2022-12-27'), by = "1 week"),
               seq(ymd('2023-01-01'), ymd('2023-12-27'), by = "1 week"),
               seq(ymd('2024-01-01'), ymd('2024-12-27'), by = "1 week"),
               ymd('2025-01-01')) # keep this last one just for indexing, then drop after aggregation loop

mset_1wk <- cbind.data.frame(mset_genus, as.data.frame(matrix(NA, nrow = nrow(mset_genus), ncol = length(unq_names))))
colnames(mset_1wk) <- c(colnames(mset_genus), as.character(unq_names))

for (i in 1:(length(unq_names)-1)){
  #i <- 1
  print(unq_names[i])
  time_span <- interval(ymd(unq_names[i]), ymd(unq_names[i+1])) # an 1 week time span
  file_inds <- which(ymd(names(file_list)) %within% time_span)
  if (length(file_inds) > 0){
    tree_spectra <- purrr::map_df(file_list[file_inds], fread, .id = 'date') %>% filter(Object_ID %in% mset_genus$Poly_ID)
    
    tree_obs_counts <- as.data.frame(table(tree_spectra$Object_ID))
    colnames(tree_obs_counts)[1] <- "Poly_ID"
    tree_obs_counts$Poly_ID <- as.numeric(as.character(tree_obs_counts$Poly_ID))
    mset_joined_span <- full_join(mset_genus, tree_obs_counts) # check if factor indexing doesn't get changed 
    
    mset_1wk[, i+2] <- mset_joined_span$Freq
  } else {
    print("No extracted data in time range")
  }
}

mset_1wk <- mset_1wk %>% select(!`2025-01-01`)

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
write.csv(mset_1wk, "mset_counts_8b_1week.csv", row.names = FALSE)



#####

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
mset_month <- fread("mset_counts_8b_month.csv")
mset_2wk <- fread("mset_counts_8b_2week.csv")
mset_1wk <- fread("mset_counts_8b_1week.csv")

mset_month_1 <- mset_month[,3:ncol(mset_month)]
mset_month_1[mset_month_1 > 1] <- 1

mset_2wk_1 <- mset_2wk[,3:ncol(mset_2wk)]
mset_2wk_1[mset_2wk_1 > 1] <- 1

mset_1wk_1 <- mset_1wk[,3:ncol(mset_1wk)]
mset_1wk_1[mset_1wk_1 > 1] <- 1

mset_month_1 %>% colSums(na.rm = TRUE) %>% sort()
mset_2wk_1 %>% colSums(na.rm = TRUE) %>% sort()
mset_1wk_1 %>% colSums(na.rm = TRUE) %>% sort()

# This will be take some minimum threshold number of trees
# Get all dates with at least some value, or an intersected value

mset_month_1_sums <- mset_month_1 %>% colSums(na.rm = TRUE) %>% as.data.frame()
mset_month_1_sums$date <- rownames(mset_month_1_sums)
mset_month_1_sums$agg_span <- "Month"
colnames(mset_month_1_sums)[1] <- "ntrees"

mset_2wk_1_sums <- mset_2wk_1 %>% colSums(na.rm = TRUE) %>% as.data.frame()
mset_2wk_1_sums$date <- rownames(mset_2wk_1_sums)
mset_2wk_1_sums$agg_span <- "Two_Week"
colnames(mset_2wk_1_sums)[1] <- "ntrees"

mset_1wk_1_sums <- mset_1wk_1 %>% colSums(na.rm = TRUE) %>% as.data.frame()
mset_1wk_1_sums$date <- rownames(mset_1wk_1_sums)
mset_1wk_1_sums$agg_span <- "One_Week"
colnames(mset_1wk_1_sums)[1] <- "ntrees"

mset_compare_sums <- rbind.data.frame(mset_month_1_sums, mset_2wk_1_sums, mset_1wk_1_sums)

mset_compare_sums$date <- ymd(mset_compare_sums$date)

ggplot(mset_compare_sums) +
  geom_point(aes(x = date, y = ntrees, color = agg_span)) +
  scale_x_date(breaks = "1 year", date_labels = "%Y", date_minor_breaks = "3 months", expand = c(0.05,0.05)) +
  facet_wrap(~agg_span)
ggsave("ntree_counts_8b_w_time_spans.png", width = 10, height = 3, units = "in")



#####



#####

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/tree_count_figures")
mset_4b_month <- fread("mset_counts_4b_month.csv")
mset_4b_2wk <- fread("mset_counts_4b_2week.csv")
mset_4b_1wk <- fread("mset_counts_4b_1week.csv")

mset_8b_month <- fread("mset_counts_8b_month.csv")
mset_8b_2wk <- fread("mset_counts_8b_2week.csv")
mset_8b_1wk <- fread("mset_counts_8b_1week.csv")

# Load 4b and 8b time ranges, sort, then sort to drop off smallest and get intersection to see which trees are left

# Set it up so that we retain X% of the maximum possible trees (with some % of time points, threshold?)
countTreesSpan <- function(mset_span, mset_label){
  #mset_span <- mset_8b_month
  mset_span_1 <- mset_span[,3:ncol(mset_span)]
  mset_span_1[mset_span_1 > 1] <- 1
  
  span_order <- mset_span_1 %>% colSums(na.rm = TRUE) %>% order()
  
  mset_span_1 <- as.data.frame(mset_span_1)
  
  tree_counts <- cbind.data.frame(0:(length(span_order)-1), NA)
  colnames(tree_counts) <- c("Time_points_removed", "Num_trees")
  for (i in 0:(length(span_order)-1)){
    selected_span <- span_order[(i+1):length(span_order)] # column order is scrambled, but does not matter here. Just doing row counts
    mset_span_sub <- cbind.data.frame(mset_span[,1:2], mset_span_1[, selected_span])
    mset_span_sub2 <- na.omit(mset_span_sub)
    tree_counts$Num_trees[i+1] <- nrow(mset_span_sub2)
  }
  tree_counts$label <- mset_label
  return(tree_counts)
}

# now need another function that does this for both 4b and 8b, and then this is the retained set? or do this post-hoc?

tree_count_4b_month <- countTreesSpan(mset_4b_month, "4b_month")
tree_count_4b_2wk <- countTreesSpan(mset_4b_2wk, "4b_2wk")
tree_count_4b_1wk <- countTreesSpan(mset_4b_1wk, "4b_1wk")
tree_count_8b_month <- countTreesSpan(mset_8b_month, "8b_month")
tree_count_8b_2wk <- countTreesSpan(mset_8b_2wk, "8b_2wk")
tree_count_8b_1wk <- countTreesSpan(mset_8b_1wk, "8b_1wk")

tree_count_all <- rbind.data.frame(tree_count_4b_month, tree_count_4b_2wk, tree_count_4b_1wk, tree_count_8b_month, tree_count_8b_2wk, tree_count_8b_1wk)


ggplot(tree_count_all) +
  geom_point(aes(x = Time_points_removed, y = Num_trees)) +
  facet_wrap(~label, scales = "free_x")
ggsave("ntree_time_points_removed.png", width = 10, height = 5, units = "in")


#
# Get dates for which intersections of tree set are > 95% of the maximum, for both 4b and 8b
# Take intersection of the set for 4b and 8b, these are the Poly_IDs of interest
# Time windows to do for both 4b and 8b are then known based on this previous thresholding

mset_span <- mset_4b_2wk
mset_span_1 <- mset_span[,3:ncol(mset_span)]
mset_span_1[mset_span_1 > 1] <- 1

span_order <- mset_span_1 %>% colSums(na.rm = TRUE) %>% order()

mset_span_1 <- as.data.frame(mset_span_1)

tree_counts_wdate <- cbind.data.frame(0:(length(span_order)-1), NA, NA)
colnames(tree_counts_wdate) <- c("Time_points_removed", "Num_trees", "Addl_range_removed")
for (i in 0:(length(span_order)-1)){
  selected_span <- span_order[(i+1):length(span_order)] # column order is scrambled, but does not matter here. Just doing row counts
  mset_span_sub <- cbind.data.frame(mset_span[,1:2], mset_span_1[, selected_span])
  mset_span_sub2 <- na.omit(mset_span_sub)
  tree_counts_wdate$Num_trees[i+1] <- nrow(mset_span_sub2)
  
  if (i > 0){
    tree_counts_wdate$Addl_range_removed[i+1] <- colnames(mset_span_1)[span_order[i]]
  }
}
#tree_counts$label <- mset_label
min_tree_count <- quantile(1:243675, 0.95)
#min(which(tree_counts$Num_trees > quantile(1:243675, 0.95)))
tree_counts_wdate$in_ts <- FALSE
tree_counts_wdate$in_ts[which(tree_counts_wdate$Num_trees > min_tree_count)] <- TRUE
tree_counts_wdate$in_ts[min(which(tree_counts_wdate$in_ts == TRUE))] <- FALSE #set this smallest index to be false, since this was removed
tree_counts_wdate$Addl_range_removed <- ymd(tree_counts_wdate$Addl_range_removed)

# add in last remaining date for full plotting
last_remaining_date <- colnames(mset_span_1)[which(colnames(mset_span_1) %notin% tree_counts_wdate$Addl_range_removed)]
tree_counts_wdate_f <- rbind.data.frame(tree_counts_wdate, c(NA, NA, last_remaining_date, TRUE))

ggplot(tree_counts_wdate_f) +
  geom_point(aes(x = yday(Addl_range_removed), y = year(Addl_range_removed), color = in_ts)) +
  scale_color_manual(values = c("gray90", "black")) +
  scale_x_continuous(breaks = c(0, 60, 150, 240, 330), minor_breaks = c(30, 90, 120, 180, 210, 270, 300, 360)) +
  scale_y_continuous(breaks = 2018:2024, minor_breaks = NULL) +
  labs(x = "Day of Year (DOY)", y = "Year", color = "Window kept?", title = "4 band PlanetScope @ 2 week time windows") +
  theme_bw()
ggsave("ps_time_windows_95pct_of_trees_4band.png", width = 8, height = 5, units = "in")

# get selected date ranges
date_start_4b <- tree_counts_wdate_f$Addl_range_removed[which(tree_counts_wdate_f$in_ts == TRUE)]
date_ranges_4b <- cbind.data.frame(sort(date_start_4b), NA)
colnames(date_ranges_4b) <- c("start_date", "end_date")
date_ranges_4b$end_date <- ymd(date_ranges_4b$start_date) + days(13)

year_end_dates <- which(yday(date_ranges_4b$end_date) > 360)
date_ranges_4b$end_date[year_end_dates] <- ymd(paste0(as.character(year(date_ranges_4b$end_date[year_end_dates])), "-12-31"))

# Get tree Poly IDs
selected_span <- span_order[(min(which(tree_counts_wdate$in_ts == TRUE))-1):length(span_order)] # subtract 1 to get correct indexing
mset_span_sub <- cbind.data.frame(mset_span[,1:2], mset_span_1[, selected_span])
mset_span_sub2 <- na.omit(mset_span_sub)

tree_id_list_4b <- mset_span_sub2$Poly_ID
table(mset_span_sub2$genus)


#####
mset_span <- mset_8b_2wk
mset_span_1 <- mset_span[,3:ncol(mset_span)]
mset_span_1[mset_span_1 > 1] <- 1

span_order <- mset_span_1 %>% colSums(na.rm = TRUE) %>% order()

mset_span_1 <- as.data.frame(mset_span_1)

tree_counts_wdate <- cbind.data.frame(0:(length(span_order)-1), NA, NA)
colnames(tree_counts_wdate) <- c("Time_points_removed", "Num_trees", "Addl_range_removed")
for (i in 0:(length(span_order)-1)){
  selected_span <- span_order[(i+1):length(span_order)] # column order is scrambled, but does not matter here. Just doing row counts
  mset_span_sub <- cbind.data.frame(mset_span[,1:2], mset_span_1[, selected_span])
  mset_span_sub2 <- na.omit(mset_span_sub)
  tree_counts_wdate$Num_trees[i+1] <- nrow(mset_span_sub2)
  
  if (i > 0){
    tree_counts_wdate$Addl_range_removed[i+1] <- colnames(mset_span_1)[span_order[i]]
  }
}
#tree_counts$label <- mset_label
min_tree_count <- quantile(1:243675, 0.95)
#min(which(tree_counts$Num_trees > quantile(1:243675, 0.95)))
tree_counts_wdate$in_ts <- FALSE
tree_counts_wdate$in_ts[which(tree_counts_wdate$Num_trees > min_tree_count)] <- TRUE
tree_counts_wdate$in_ts[min(which(tree_counts_wdate$in_ts == TRUE))] <- FALSE #set this smallest index to be false, since this was removed
tree_counts_wdate$Addl_range_removed <- ymd(tree_counts_wdate$Addl_range_removed)

# add in last remaining date for full plotting
last_remaining_date <- colnames(mset_span_1)[which(colnames(mset_span_1) %notin% tree_counts_wdate$Addl_range_removed)]
tree_counts_wdate_f <- rbind.data.frame(tree_counts_wdate, c(NA, NA, last_remaining_date, TRUE))

ggplot(tree_counts_wdate_f) +
  geom_point(aes(x = yday(Addl_range_removed), y = year(Addl_range_removed), color = in_ts)) +
  scale_color_manual(values = c("gray90", "black")) +
  #scale_x_continuous(breaks = c(0, 60, 150, 240, 330), minor_breaks = c(30, 90, 120, 180, 210, 270, 300, 360)) +
  scale_y_continuous(breaks = 2021:2024, minor_breaks = NULL) +
  labs(x = "Day of Year (DOY)", y = "Year", color = "Window kept?", title = "8 band PlanetScope @ 2 week time windows") +
  theme_bw()
ggsave("ps_time_windows_95pct_of_trees_8band.png", width = 8, height = 3, units = "in")

# From these two data sets (4b and 8b), generate list of date ranges for each, and the intersected Poly_ID for both


# get selected date ranges
date_start_8b <- tree_counts_wdate_f$Addl_range_removed[which(tree_counts_wdate_f$in_ts == TRUE)]
date_ranges_8b <- cbind.data.frame(sort(date_start_8b), NA)
colnames(date_ranges_8b) <- c("start_date", "end_date")
date_ranges_8b$end_date <- ymd(date_ranges_8b$start_date) + days(13)

year_end_dates <- which(yday(date_ranges_8b$end_date) > 360)
date_ranges_8b$end_date[year_end_dates] <- ymd(paste0(as.character(year(date_ranges_8b$end_date[year_end_dates])), "-12-31"))

# Get tree Poly IDs
selected_span <- span_order[(min(which(tree_counts_wdate$in_ts == TRUE))-1):length(span_order)] # subtract 1 to get correct indexing
mset_span_sub <- cbind.data.frame(mset_span[,1:2], mset_span_1[, selected_span])
mset_span_sub2 <- na.omit(mset_span_sub)

tree_id_list_8b <- mset_span_sub2$Poly_ID
table(mset_span_sub2$genus)

# Intersect 4b and 8b, this is what we can use everywhere
tree_id_list_intersect <- intersect(tree_id_list_4b, tree_id_list_8b)
write.csv(tree_id_list_intersect, "tree_id_list_95pct_2week_intersect.csv", row.names = FALSE)
write.csv(tree_id_list_4b, "tree_id_list_95pct_2week_4b.csv", row.names = FALSE)
write.csv(tree_id_list_8b, "tree_id_list_95pct_2week_8b.csv", row.names = FALSE)

# and the date_ranges_
write.csv(date_ranges_4b, "tree_id_date_ranges_95pct_2week_4b.csv", row.names = FALSE)
write.csv(date_ranges_8b, "tree_id_date_ranges_95pct_2week_8b.csv", row.names = FALSE)

# Use this list of tree Poly_IDs, and the 4b and 8b image date ranges, to make composites for each tree for each date range.