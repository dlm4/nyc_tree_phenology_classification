# Comparing the accuracies from the different classification outputs

# load each variation, with OOB, OA, kappa, and cilo and cihi for each with the parsimonious iteration

setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/")
#agg_name <- "everything"

getBestIter <- function(agg_name){
  mean_acc_filename <- list.files(pattern = glob2rx(paste0(agg_name, "*_feature_selection_data_mean.csv")))
  mean_acc <- read.csv(mean_acc_filename)
  sel_feat_filename <- list.files(pattern = glob2rx(paste0(agg_name, "*_selected_features_info_cv.rds")))
  sel_feat <- readRDS(sel_feat_filename)
  output_df <- cbind.data.frame(agg_name, mean_acc[sel_feat$iter_round_selected,])
  return(output_df)
}

mean_acc_all <- rbind.data.frame(getBestIter("1week_singleyears"),
                                 getBestIter("2week_singleyears"),
                                 getBestIter("4week_singleyears"),
                                 getBestIter("8week_singleyears"),
                                 getBestIter("1week_collapsed"),
                                 getBestIter("2week_collapsed"),
                                 getBestIter("4week_collapsed"),
                                 getBestIter("8week_collapsed"),
                                 getBestIter("allweek_singleyear"),
                                 getBestIter("allweek_collapsed"),
                                 getBestIter("everything"))

p1 <- ggplot(mean_acc_all) +
  geom_point(aes(x = agg_name, y = oa_out), size = 1) +
  geom_linerange(aes(x = agg_name, ymin = cilo_out, ymax = cihi_out), linewidth = 0.3) +
  labs(x = "Aggregation Time Range", y = "Mean Overall Accuracy (%)") +
  scale_y_continuous(breaks = seq(0.5, 0.8, 0.05), minor_breaks = seq(0.5, 0.8, 0.01)) +
  theme_bw() +
  theme(panel.grid.minor.y = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5, color = "black"))
setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
ggsave("mean_acc_runs.jpg", p1, width = 90, height = 80, units = "mm")
ggsave("mean_acc_runs.eps", p1, width = 90, height = 80, units = "mm")


# Compare accuracies for individual tree genera --------------------------

agg_name <- "everything"
setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/varsel_outputs_cv/")

#getBestIterGen <- function(agg_name)
mean_acc_filename <- list.files(pattern = glob2rx(paste0(agg_name, "*_feature_selection_data_mean.csv")))
mean_acc <- read.csv(mean_acc_filename)
sel_feat_filename <- list.files(pattern = glob2rx(paste0(agg_name, "*_selected_features_info_cv.rds")))
sel_feat <- readRDS(sel_feat_filename)
#output_df <- cbind.data.frame(agg_name, mean_acc[sel_feat$iter_round_selected,])
confmat_filename <- list.files(pattern = glob2rx(paste0(agg_name, "*_conf_mat_biglist.rds")))
confmat_biglist <- readRDS(confmat_filename)
#return(output_df)

# mean of confusion matrix tables
# example for everything version, can re-do for others
x <- iter_round_selected # 3 for everything
confmat_tables <-  list(conf_mat_biglist[[x]][[1]]$table, conf_mat_biglist[[x]][[2]]$table, conf_mat_biglist[[x]][[3]]$table, conf_mat_biglist[[x]][[4]]$table, conf_mat_biglist[[x]][[5]]$table,
                        conf_mat_biglist[[x]][[6]]$table, conf_mat_biglist[[x]][[7]]$table, conf_mat_biglist[[x]][[8]]$table, conf_mat_biglist[[x]][[9]]$table, conf_mat_biglist[[x]][[10]]$table)
confmat_table_mean <- Reduce('+', confmat_tables)/length(confmat_tables)

confmat_overalls <- list(conf_mat_biglist[[x]][[1]]$overall, conf_mat_biglist[[x]][[2]]$overall, conf_mat_biglist[[x]][[3]]$overall, conf_mat_biglist[[x]][[4]]$overall, conf_mat_biglist[[x]][[5]]$overall,
                         conf_mat_biglist[[x]][[6]]$overall, conf_mat_biglist[[x]][[7]]$overall, conf_mat_biglist[[x]][[8]]$overall, conf_mat_biglist[[x]][[9]]$overall, conf_mat_biglist[[x]][[10]]$overall)
confmat_overall_mean <- Reduce('+', confmat_overalls)/length(confmat_overalls)

confmat_byClass <-  list(conf_mat_biglist[[x]][[1]]$byClass, conf_mat_biglist[[x]][[2]]$byClass, conf_mat_biglist[[x]][[3]]$byClass, conf_mat_biglist[[x]][[4]]$byClass, conf_mat_biglist[[x]][[5]]$byClass,
                         conf_mat_biglist[[x]][[6]]$byClass, conf_mat_biglist[[x]][[7]]$byClass, conf_mat_biglist[[x]][[8]]$byClass, conf_mat_biglist[[x]][[9]]$byClass, conf_mat_biglist[[x]][[10]]$byClass)
confmat_byClass_mean <- Reduce('+', confmat_byClass)/length(confmat_byClass)



# Importance scores -------------------------------------------------------

# This is variable importance for the entire runs
var_imp <- readRDS("everything_var_imp_list.rds")
imp_ordered <- order(var_imp[[x]]$imp_mean, decreasing = TRUE)
var_imp[[x]]$imp_sd <- apply(var_imp[[x]][,2:11], 1, sd)
var_imp_mean_ordered <- cbind.data.frame(var_imp[[x]]$var_names[imp_ordered], var_imp[[x]]$imp_mean[imp_ordered], var_imp[[x]]$imp_sd[imp_ordered])
colnames(var_imp_mean_ordered) <- c("var_names", "imp_mean", "imp_sd")
var_imp_mean_ordered$var_names <- factor(var_imp_mean_ordered$var_names, levels = unique(var_imp_mean_ordered$var_names))

p2 <- ggplot(var_imp_mean_ordered[1:30,]) +
  geom_point(aes(x = imp_mean, y = var_names), size = 1) +
  geom_linerange(aes(xmin = imp_mean - imp_sd, xmax = imp_mean + imp_sd, y = var_names), linewidth = 0.3) +
  scale_y_discrete(limits = rev) +
  labs(x = "Mean Permutation Importance", y = "Variable Name") +
  theme_bw() +
  theme(#panel.grid.minor.y = element_line(color = "gray90", linetype = "dotted", linewidth = 0.1),
        panel.grid = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 5, color = "black"))
setwd("/Volumes/NYC_geo/tree_classification/extracted_timeranges/median_wide/output_figs_and_tables")
ggsave("mean_var_imp_everything.jpg", p2, width = 90, height = 100, units = "mm")
ggsave("mean_var_imp_everything.eps", p2, width = 90, height = 100, units = "mm")



# Need to re-run appropriate iterations with seeds to get variable importance for individual species for the specific RF model run rounds
# Just do this for the most parsimonious set for each
# Would not have known which round was it a priori, so need to do it again and retain it, but should not take too long to process