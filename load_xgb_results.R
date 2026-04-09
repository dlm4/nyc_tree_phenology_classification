# compare_xgb_outputs.R

# compare confusion matrices, overall accuracy, etc.

library(tidyverse)
library(xgboost)

setwd("/Users/dlm356/Library/CloudStorage/Box-Box/dlm_box/data/classification/outputs")

cm_file_list <- list.files(pattern = glob2rx("*confmat.rds"))

#file_name <- "everything_xgb_full_confmat.rds"
getConfMatOA <- function(file_name){
  cm <- readRDS(file_name)
  output <- c(file_name, cm$overall["Accuracy"], cm$overall["Kappa"], cm$overall["AccuracyLower"], cm$overall["AccuracyUpper"])
  return(output)
}

#cm_df <- cbind.data.frame(matrix(NA, length(cm_file_list), 5))
cm_df <- sapply(cm_file_list, FUN = getConfMatOA)
cm_df %>% t() %>% as.data.frame()


xgb_model_list <- list.files(pattern = glob2rx("*.model"))

getXGBiter <- function(xgb_file){
  xgb_model <- xgb.load(xgb_file)
  output <- c(xgb_file, xgb_model$best_iteration, xgb_model$best_score)
  return(output)
}

xgb_df <- sapply(xgb_model_list, FUN = getXGBiter)
xgb_df %>% t() %>% as.data.frame()


# make plot showing OA vs date for single dates
cm_df <- cm_df %>% t() %>% as.data.frame()

cm_df_1tp <- cm_df[grep("_1timepoint_confmat", cm_df$V1),]
cm_df_1tp$date <- cm_df_1tp$V1 %>% strsplit("_") %>% lapply('[[', 2) %>% unlist() %>% ymd()

ggplot(cm_df_1tp) +
  geom_hline(yintercept = 0.482) +
  geom_linerange(aes(xmin = yday(date), xmax = yday(date)+13, y = as.numeric(Accuracy), color = as.factor(year(date)))) +
  #geom_point(aes(x = yday(date), y = as.numeric(Accuracy), color = as.factor(year(date)))) +
  #scale_x_continuous(breaks = c(60, ))
  labs(x = "Day of Year", y = "Overall Accuracy", color = "Year") +
  theme_bw()

#####
# best single time points for classification? Will def be with 4+8b data
# Keeping 4 and 8 band separate because of some different acquisitions in overlapping years