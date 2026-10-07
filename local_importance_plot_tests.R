

library(data.table)
library(tidyverse)
library(tidytext) # for reorder_within
library(ggridges)
# https://stackoverflow.com/questions/61340327/how-to-obtain-feature-importance-by-class-using-ranger
#as.data.table(rf.iris$variable.importance.local)[,Species := iris$Species][,lapply(.SD,mean),by=Species]

local_varimp_list <- readRDS("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/genus_importance/everything_local_var_imp_list_mainiter.rds")

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
#ggsave(paste0(feature_set_name, "_local_imp.jpg"), p_imp, width = 190, height = 150, units = "mm")


# Duration ----------------------------------------------------------------

# parse for duration vector
local_varimp_giantdf_summary_long <- local_varimp_giantdf_summary_long %>% 
  mutate(dur = str_split_i(time_range_feature, "_", -2),
         tr = str_split_i(time_range_feature, "_", -1))

tr_prefix <- local_varimp_giantdf_summary_long$tr %>% str_sub(1, 1)
d_prefix <- which(tr_prefix == "d")
local_varimp_giantdf_summary_long$dur2 <- local_varimp_giantdf_summary_long$dur
local_varimp_giantdf_summary_long$dur2[d_prefix] <- paste0("c_", local_varimp_giantdf_summary_long$dur[d_prefix])
local_varimp_giantdf_summary_long$dur2[-d_prefix] <- paste0("sy_", local_varimp_giantdf_summary_long$dur[-d_prefix])

dur_sums <- local_varimp_giantdf_summary_long %>% 
  select(genus, value, dur2) %>%
  group_by(genus, dur2) %>%
  summarize(imp_sum = sum(value))

# ggplot(dur_sums) +
#   geom_col(aes(x = genus, y = imp_sum, fill = dur2))

# currently normalized by the sum of all
# dur_sums <- dur_sums %>% group_by(genus) %>% mutate(imp_sum_norm = imp_sum/(sum(imp_sum)))
#, but could normalize to max instead
dur_sums <- dur_sums %>% group_by(genus) %>% mutate(imp_sum_norm = imp_sum/(max(imp_sum)))

p_durcomp <- ggplot(dur_sums) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray90") +
  geom_vline(xintercept = seq(1.5, 11.5), linewidth = 0.3, color = "gray80") +
  geom_col(aes(x = genus, y = imp_sum_norm, fill = dur2), position = position_dodge(), color = "gray20", linewidth = 0.2) +
  labs(x = "Genus", y = "Sum of Permutation Importance\n(scaled by genus)", fill = "Time Agg.") +
  scale_fill_discrete(palette = scales::pal_brewer(palette = "Dark2")) +
  theme_bw() + 
  theme(panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 7, color = "black"),
        legend.key.size = unit(2, "mm"),
        strip.text = element_text(size = 7, color = "black", face = "italic"),
        strip.background = element_rect(fill = "white", color = NA))
#setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
#ggsave("fig_duration_comparison.jpg", p_durcomp, width = 190, height = 50, units = "mm")
ggsave("fig_duration_comparison_max.jpg", p_durcomp, width = 190, height = 50, units = "mm")

 #https://www.kaggle.com/code/dansbecker/permutation-importance
# permuation importance can be negative if the randomly shuffled column is better for prediction than the actual data. Means very bad predictors


# Show most important single feature for each duration
dur_max <- local_varimp_giantdf_summary_long %>% group_by(genus, dur2) %>% summarize(imp_max = max(value))
dur_max <- dur_max %>% group_by(genus) %>% mutate(imp_max_norm = imp_max/(max(imp_max)))

p_durmax <- ggplot(dur_max) +
  geom_hline(yintercept = seq(0, 1, 0.25), linewidth = 0.3, color = "gray90", linetype = "dotted") +
  geom_vline(xintercept = seq(1.5, 11.5), linewidth = 0.3, color = "gray80") +
  geom_point(aes(x = genus, y = imp_max_norm, color = dur2, shape = dur2), position = position_dodge(width = 1), size = 1) +
  labs(x = "Genus", y = "Max of Permutation Importance\n(scaled by genus)", color = "Time Agg.", shape = "Time Agg.") +
  scale_color_discrete(palette = scales::pal_brewer(palette = "Dark2")) +
  theme_bw() + 
  theme(#panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_blank(),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        axis.text.x = element_text(size = 5, color = "black", face = "italic"),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 7, color = "black"),
        legend.key.size = unit(2, "mm"),
        strip.text = element_text(size = 7, color = "black", face = "italic"),
        strip.background = element_rect(fill = "white", color = NA))
setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
ggsave("fig_duration_max_only.jpg", p_durmax, width = 190, height = 50, units = "mm")


# Time of year, year ------------------------------------------------------

local_varimp_giantdf_summary_long$date_start <- NA
local_varimp_giantdf_summary_long$date_start[-d_prefix] <- as.character(ymd(local_varimp_giantdf_summary_long$tr[-d_prefix]))

# by year
local_varimp_giantdf_summary_long$yr <- "collapsed"
local_varimp_giantdf_summary_long$yr[-d_prefix] <- year(ymd(local_varimp_giantdf_summary_long$date_start[-d_prefix]))

yr_sums <- local_varimp_giantdf_summary_long %>% 
  select(genus, value, yr) %>%
  group_by(genus, yr) %>%
  summarize(imp_sum = sum(value))
yr_sums <- yr_sums %>% group_by(genus) %>% mutate(imp_sum_norm = imp_sum/(sum(imp_sum)))
#yr_sums <- yr_sums %>% group_by(genus) %>% mutate(imp_sum_norm = imp_sum/(max(imp_sum))) # max
p_yrcomp <- ggplot(yr_sums) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray90") +
  geom_vline(xintercept = seq(1.5, 11.5), linewidth = 0.3, color = "gray80") +
  geom_col(aes(x = genus, y = imp_sum_norm, fill = yr), position = position_dodge(), color = "gray20", linewidth = 0.2) +
  labs(x = "Genus", y = "Sum of Permutation Importance\n(scaled by genus)", fill = "Year") +
  scale_fill_discrete(palette = scales::pal_brewer(palette = "RdYlBu")) +
  theme_bw() + 
  theme(panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 7, color = "black"),
        legend.key.size = unit(2, "mm"),
        strip.text = element_text(size = 7, color = "black", face = "italic"),
        strip.background = element_rect(fill = "white", color = NA))
#setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
ggsave("fig_yr_comparison.jpg", p_yrcomp, width = 190, height = 50, units = "mm")
#ggsave("fig_yr_comparison_max.jpg", p_yrcomp, width = 190, height = 50, units = "mm")


# by time of year
local_varimp_giantdf_summary_long$doy_start <- NA
local_varimp_giantdf_summary_long$doy_start[-d_prefix] <- yday(ymd(local_varimp_giantdf_summary_long$date_start[-d_prefix]))
local_varimp_giantdf_summary_long$doy_start[d_prefix] <- as.numeric(substring(local_varimp_giantdf_summary_long$tr[d_prefix], 4))
local_varimp_giantdf_summary_long$doy_end <- local_varimp_giantdf_summary_long$doy_start + as.numeric(substr(local_varimp_giantdf_summary_long$dur,1,1))*7
# max doy_end is 337, so don't need to correct for end of year wraparound, thankfully

df_summary_lineplot <- local_varimp_giantdf_summary_long %>% 
  group_by(genus) %>% 
  mutate(value_norm = value/(max(value)))

df_summary_lineplot %>% filter(genus %in% c("Ginkgo", "Pyrus", "Platanus")) %>%
  ggplot() +
  geom_hline(yintercept = seq(1.5, 7.5), linewidth = 0.1, color = "gray80") +
  geom_vline(xintercept = seq(120, 330, 30), linewidth = 0.1, linetype = "dotted", color = "gray20") +
  geom_linerange(aes(xmin = doy_start, xmax = doy_end, y = yr, color = dur, alpha = value_norm), position = position_dodge(width = 1), linewidth = 0.5) +
  scale_color_discrete(palette = scales::pal_brewer(palette = "Set1")) +
  scale_alpha(range = c(0,1)) +
  facet_wrap(~genus) +
  scale_y_discrete(limits = rev) +
  labs(x = "Day of Year", y = "Year", alpha = "Scaled Imp.", color = "Time Agg.") +
  theme_bw() +
  theme(#panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_blank(),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 7, color = "black"),
        legend.key.size = unit(3, "mm"),
        strip.text = element_text(size = 7, color = "black", face = "italic"),
        strip.background = element_rect(fill = "white", color = NA))
#setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
#ggsave("fig_time_ranges_genus.jpg", width = 190, height = 120, units = "mm")
ggsave("fig_time_ranges_genus_ex.jpg", width = 190, height = 60, units = "mm")
# repeat for only a few genera, for example. Maybe just set this up as a comparison

df_summary_lineplot <- df_summary_lineplot %>% 
  mutate(doy_mid = ceiling((doy_start + doy_end)/2))

# Try this again with all years together, ridgeline plot
df_summary_lineplot %>% 
ggplot() +
  geom_ridgeline(aes(x = doy_mid, y = genus, height = value_norm))

# this setup didn't work, need to sum the importances another way...

doy_range <- seq(min(df_summary_lineplot$doy_start), max(df_summary_lineplot$doy_end))
gen_list <- unique(df_summary_lineplot$genus)
df_by_doy <- cbind.data.frame(rep(gen_list, each = length(doy_range)), rep(doy_range, length(gen_list)))
colnames(df_by_doy) <- c("Genus", "DOY")

# create daily sum and max of importances...
df_by_doy$value_sum <- 0
df_by_doy$value_max <- 0
for (i in 1:nrow(df_by_doy)){
  g <- df_by_doy$Genus[i]
  d <- df_by_doy$DOY[i]
  df_sub <- df_summary_lineplot %>% filter(genus == g, doy_start <= d, doy_end >= d)
  df_by_doy$value_sum[i] <- sum(df_sub$value)
  df_by_doy$value_max[i] <- max(df_sub$value)
}

df_by_doy <- df_by_doy %>% 
  group_by(Genus) %>% 
  mutate(value_sum_norm = value_sum/(max(value_sum)),
         value_max_norm = value_max/(max(value_max)))

#month_doy_1 <- c(1, 32, 60, 91, 121, 152, 183, 213, 244, 274, 305, 335)
month_doy_1 <- c(121, 152, 183, 213, 244, 274, 305, 335)
#month_doy_mid <- month_doy_1 + 14

month_names <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
month_doy_mid <- c(110, month_doy_1[1:7] + 15)

df_by_doy %>% 
  ggplot() +
  geom_ridgeline(aes(x = DOY, y = Genus, height = value_sum_norm, fill = Genus), min_height = 0, linewidth = 0.3) +
  scale_y_discrete(limits = rev) +
  scale_fill_discrete(palette = scales::pal_brewer(palette = "Set3")) +
  guides(fill = "none") +
  #geom_vline(xintercept = seq(120, 330, 30), linewidth = 0.1, linetype = "dotted", color = "gray20") +
  geom_vline(xintercept = month_doy_1, linewidth = 0.1, linetype = "dotted", color = "gray20") +
  annotate("text", month_doy_mid, y = 0.5, label = month_names, size = 5/.pt) +
  theme_bw() +
  theme(#panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black", face = "italic"),
    #legend.text = element_text(size = 5, color = "black"),
    #legend.title = element_text(size = 7, color = "black"),
    #legend.key.size = unit(3, "mm"),
    strip.text = element_text(size = 7, color = "black", face = "italic"),
    strip.background = element_rect(fill = "white", color = NA))
ggsave("fig_genus_ridgeline_test.jpg", width = 120, height = 90, units = "mm")


df_by_doy %>% 
  ggplot() +
  geom_ridgeline(aes(x = DOY, y = Genus, height = value_max_norm, fill = Genus), min_height = 0, linewidth = 0.3) +
  scale_y_discrete(limits = rev) +
  scale_fill_discrete(palette = scales::pal_brewer(palette = "Set3")) +
  guides(fill = "none") +
  #geom_vline(xintercept = seq(120, 330, 30), linewidth = 0.1, linetype = "dotted", color = "gray20") +
  geom_vline(xintercept = month_doy_1, linewidth = 0.1, linetype = "dotted", color = "gray20") +
  annotate("text", month_doy_mid, y = 0.5, label = month_names, size = 5/.pt) +
  theme_bw() +
  theme(#panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black", face = "italic"),
    #legend.text = element_text(size = 5, color = "black"),
    #legend.title = element_text(size = 7, color = "black"),
    #legend.key.size = unit(3, "mm"),
    strip.text = element_text(size = 7, color = "black", face = "italic"),
    strip.background = element_rect(fill = "white", color = NA))
ggsave("fig_genus_ridgeline_test_v2.jpg", width = 120, height = 90, units = "mm")


# df_by_doy %>% 
#   ggplot() +
#   geom_line(aes(x = DOY, y = value_norm_sum)) +
#   facet_wrap(~Genus) +
#   theme_bw()




# Spectra -----------------------------------------------------------------

local_varimp_giantdf_summary_long <- local_varimp_giantdf_summary_long %>% 
  mutate(band_name = str_split_i(time_range_feature, "_median_", -2))

local_varimp_giantdf_summary_long <- local_varimp_giantdf_summary_long %>% 
  mutate(doy_mid = ceiling((doy_start + doy_end)/2))

df_for_spectra <- local_varimp_giantdf_summary_long %>% 
  group_by(genus) %>% 
  mutate(value_norm = value/(max(value)))

df_for_spectra <- df_for_spectra %>% 
  group_by(genus) %>% 
  mutate(value_norm_rank = rank(value_norm))

df_for_spectra %>% filter(genus %in% c("Ginkgo", "Pyrus", "Platanus")) %>%
ggplot() +
  geom_vline(xintercept = seq(1.5,5.5), color = "gray70", linewidth = 0.3) +
  geom_jitter(aes(x = band_name, y = value_norm, color = doy_mid, alpha = value_norm_rank), shape = 16, size = 1) +#, shape = yr)) +
  scale_color_continuous(palette = rainbow(12)) +
  #scale_shape_manual(values = c(16, 17, 15, 18, 3, 8, 10, 25)) +
  facet_wrap(~genus) + 
  scale_alpha(range = c(0,1)) +
  labs(x = "Band or Band Combination", y = "Permutation Importance\n(scaled to max by genus)", color = "DOY", alpha = "Imp. Rank") +
  theme_bw() +
  theme(panel.grid.minor = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
    panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.text = element_text(size = 5, color = "black"),
    legend.title = element_text(size = 7, color = "black"),
    legend.key.size = unit(3, "mm"),
    strip.text = element_text(size = 7, color = "black", face = "italic"),
    strip.background = element_rect(fill = "white", color = NA))
#setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
#ggsave("fig_spectra_doy.jpg", width = 190, height = 120, units = "mm")
#ggsave("fig_time_ranges_genus_ex.jpg", width = 190, height = 60, units = "mm")
ggsave("fig_spectra_doy_ex.jpg", width = 190, height = 60, units = "mm")
