

library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(rpart)
library(rpart.plot)
library(caret)
library(ranger)
library(paletteer)

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
}

# Set up a classification with a peak, baseline, and change between them (difference or normalized ratio

# Peak
# Nov 23, 2023
# Alternate is Nov 27, 2023 if point is unavailable in Nov 23, 2023 (test this second)

# Baseline
# Sep 14, 2023

tree_pheno_all <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point/")

# could set interval here if doing additional dates
date_strings <-  c("20231123", "20230914")

for (i in 1:length(date_strings)){
  dt <- date_strings[i]
  print(dt)
  dt_filename <- paste0("*_", dt, "_*.csv")
  file_list <- list.files(pattern = glob2rx(dt_filename))
  df_all_spectra_dt <- purrr::map_df(file_list, fread) # load them all in
  df_all_spectra_dt$date <- dt
  if (i == 1) {
    df_all_spectra <- df_all_spectra_dt
  } else {
    df_all_spectra <- rbind.data.frame(df_all_spectra, df_all_spectra_dt)
  }
}

# filter tree pheno input
tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023)

df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_input$Poly_ID)
df_all_spectra_sub <- merge(df_all_spectra_sub, tree_pheno_input[,c("Poly_ID", "genus", "species", "Year", "EOS_50", "R2")], by.x = "Object_ID", by.y = "Poly_ID")

input <- as.data.frame(df_all_spectra_sub)

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
    
    nd_index <- calcNormDif(input[,b2], input[,b1])
    
    if (i == 1 & j == 2){
      omnbr <- nd_index
    } else {
      omnbr <- cbind.data.frame(omnbr, nd_index)
      colnames(omnbr) <- colname_list
    }
    
  }
}

df_spectra_omnbr <- bind_cols(input, omnbr)

# Make df long based on Object_ID, Year, EOS_50, and R2 as row IDs, with variable-value pairs of date and the bands + indices and their values
df_spectra_omnbr_long <- pivot_longer(df_spectra_omnbr, all_of(c(band_list, colname_list)))
df_spectra_omnbr_long <- df_spectra_omnbr_long %>% mutate(date_band = paste0(name, "_", date)) # need this order because of ranger
df_spectra_omnbr_wide <- df_spectra_omnbr_long %>% subset(select = -c(name, date)) %>% pivot_wider(names_from = date_band, values_from = value)

# keep only trees with both dates
df_spectra_omnbr_cc <- df_spectra_omnbr_wide[complete.cases(df_spectra_omnbr_wide),] 

# Need to add in the difference or ratio step here!
var_names <- unique(df_spectra_omnbr_long$name)
date_names <- unique(df_spectra_omnbr_long$date) # could be date_strings, but might be generic and transformed later by this point. Should be length of 2 for two dates

added_bands_list <- c()

for (i in 1:length(var_names)){
  band_ind_names <- paste0(var_names[i], '_', date_names)
  
  # Difference example
  new_name <- paste0(var_names[i], '_dif')
  df_spectra_omnbr_cc[,new_name] <- df_spectra_omnbr_cc[,band_ind_names[1]] - df_spectra_omnbr_cc[,band_ind_names[2]]
  
  added_bands_list <- c(added_bands_list, new_name)
}

# loop indexing with the same names for bands or norm indices

#
full_band_list <- c(unique(df_spectra_omnbr_long$date_band), added_bands_list)
comp_list <- c("genus", full_band_list)

df_cls <- df_spectra_omnbr_cc[, comp_list]

set.seed(14)

df_cls <- df_cls[which(df_cls$genus != "Unknown"), ]
df_cls$genus[which(df_cls$genus != "Ginkgo")] <- "Other"

df_cls$genus <- as.factor(df_cls$genus)

ginkgo <- df_cls[which(df_cls$genus == "Ginkgo"),]
other <- df_cls[which(df_cls$genus == "Other"),]

samp_size <- nrow(ginkgo)
samp_inds <- sample(1:nrow(other), samp_size)
df_cls_sub <- rbind.data.frame(ginkgo, other[samp_inds,])

train_index <- sample(1:nrow(df_cls_sub), round(0.8*nrow(df_cls_sub), 0))
train_data <- df_cls_sub[train_index, ]
test_data <- df_cls_sub[-train_index, ]

tree_model <- rpart(genus ~ ., 
                    data = train_data, 
                    method = "class")

rpart.plot(tree_model, box.palette = "auto", nn = TRUE)

p <- predict(tree_model, test_data[,2:ncol(test_data)], type = 'class')
confusionMatrix(p, reference = test_data$genus)

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)

# roughly the same accuracy with the two dates for 11/23/23 as the only source, so now add in a difference ratio layer too.

#####



