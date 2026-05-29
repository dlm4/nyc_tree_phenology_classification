

library(data.table)
library(tidyverse)
library(tidytext) # for reorder_within
# https://stackoverflow.com/questions/61340327/how-to-obtain-feature-importance-by-class-using-ranger
#as.data.table(rf.iris$variable.importance.local)[,Species := iris$Species][,lapply(.SD,mean),by=Species]

local_varimp_list_readin <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/genus_importance/everything_local_var_imp_list_mainiter.rds")

# test
#local_varimp_test <- local_varimp_list_readin[[1]]
#local_varimp_test_summary <- as.data.table(local_varimp_test[,2:ncol(local_varimp_test)])[,genus := local_varimp_test$genus][,lapply(.SD,mean),by=genus]


local_varimp_giantdf <- rbind.data.frame(local_varimp_list[[1]], local_varimp_list[[2]], local_varimp_list[[3]], local_varimp_list[[4]], local_varimp_list[[5]],
                                         local_varimp_list[[6]], local_varimp_list[[7]], local_varimp_list[[8]], local_varimp_list[[9]], local_varimp_list[[10]])
local_varimp_giantdf_summary <- as.data.table(local_varimp_giantdf[,2:ncol(local_varimp_giantdf)])[,genus := local_varimp_giantdf$genus][,lapply(.SD,mean),by=genus]

# get means of each test, so we can get SDs of these
for (i in 1:10){
  local_varimp_test <- local_varimp_list[[i]]
  local_varimp_summary <- as.data.table(local_varimp_test[,2:ncol(local_varimp_test)])[,genus := local_varimp_test$genus][,lapply(.SD,mean),by=genus]
  if (i == 1){
    local_varimp_summary_mean_each <- local_varimp_summary
  } else {
    local_varimp_summary_mean_each <- rbind.data.frame(local_varimp_summary_mean_each, local_varimp_summary)
  }
}

#local_varimp_giantdf_summary_SD <- as.data.table(local_varimp_giantdf[,2:ncol(local_varimp_giantdf)])[,genus := local_varimp_giantdf$genus][,lapply(.SD,sd),by=genus]
local_varimp_giantdf_summary_SD <- local_varimp_summary_mean_each %>% group_by(genus) %>% summarize_all(sd)

local_varimp_giantdf_summary_long <- local_varimp_giantdf_summary %>% as_tibble() %>% pivot_longer(cols = !genus, names_to = "time_range_feature", values_to = "value")

local_varimp_giantdf_summary_SD_long <- local_varimp_giantdf_summary_SD %>% as_tibble() %>% pivot_longer(cols = !genus, names_to = "time_range_feature", values_to = "value_sd")

# local_varimp_giantdf_summary_long %>% filter(genus == "Acer") %>%
#   ggplot() +
#   geom_point(aes(x = value, y = time_range_feature))


topx <- local_varimp_giantdf_summary_long %>%
  group_by(genus) %>%
  top_n(10, value)

topx_sd <- inner_join(topx, local_varimp_giantdf_summary_SD_long)

# top1 <- local_varimp_giantdf_summary_long %>%
#   group_by(genus) %>%
#   slice_max(value)

p_imp <- ggplot(topx_sd) + 
  geom_point(aes(x = value, y = reorder_within(time_range_feature, value, genus)), size = 1) +
  geom_linerange(aes(xmin = value - value_sd, xmax = value + value_sd, y = reorder_within(time_range_feature, value, genus)), linewidth = 0.3) +
  facet_wrap(~genus, scales = "free", ncol = 3) +
  scale_y_reordered() +
  labs(x = "Mean Permutation Importance", y = "Most Important Features by Genus") +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        strip.text = element_text(size = 7, color = "black", face = "italic"),
        strip.background = element_rect(fill = "white", color = NA))
setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
ggsave(paste0(feature_set_name, "_local_imp.jpg"), p_imp, width = 190, height = 150, units = "mm")
 