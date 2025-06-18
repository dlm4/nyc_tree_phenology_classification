

library(tidyverse)
library(data.table)
library(purrr)

ids_acru <- c(25127, 78876, 386894, 478190, 861512, 901596,
              948007, 1261989, 1350939, 1465337, 1927262, 1256444)

# fread by objset, filter to these ids

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")

for (i in 1:7){
  print(paste0("Objset ", i))
  objset_text <- paste0("*objset", i, "_*.csv")
  file_list_sub <- list.files(pattern = glob2rx(objset_text))
  names(file_list_sub) <- strsplit(file_list_sub, "[.]") %>% lapply('[[', 1) %>% stringr::str_sub(start = 19, end = 27) %>% ymd() # name the list
  tree_df_sub <- purrr::map_df(file_list_sub, fread, .id = 'date') 
  tree_df_sub2 <- tree_df_sub %>% filter(Object_ID %in% ids_acru)
  if (i == 1){
    tree_df_agg <- tree_df_sub2
  } else {
    tree_df_agg <- bind_rows(tree_df_agg, tree_df_sub2)
  }
}

