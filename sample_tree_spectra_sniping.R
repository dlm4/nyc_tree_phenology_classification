

library(tidyverse)
library(data.table)
library(purrr)
'%notin%' <- Negate('%in%')
library(rpart)
library(rpart.plot)
library(caret)
library(ranger)
library(paletteer)

#####
# Load in individual tree spectra, organize, and write out in this first part

# Acer rubrum
ids_acru <- c(25127, 78876, 386894, 478190, 861512, 901596,
              948007, 1261989, 1350939, 1465337, 1927262, 1256444)

# Ginkgo biloba
ids_gibi <- c(58875, 80889, 108001, 129315, 160069, 402050,
              570059, 630759, 836013, 1328766, 1559444, 1614591)

# Gleditsia triacanthos
ids_gltr <- c(109918, 167837, 185594, 217624, 413182, 1010339,
              1056336, 1149348, 1318659, 1319566, 1438777, 1799222)

# Acer platanoides
ids_acpl <- c(183103, 183113, 183124, 1056200, 1731090, 2073483,
              1495495, 802344, 63612, 247345)

# Prunus serrulata
ids_prse <- c(164336, 200605, 200790, 489289, 569229, 612685,
              1030178, 1127870, 1127971, 1164209, 1333868, 494589)

# Robinia pseudoacacia
ids_rops <- c(73746, 85797, 344275, 516452, 520121, 523155,
              536633, 1134258, 1167681, 1198795, 1873218, 1000011)

# Combine
ids_all <- c(ids_acru, ids_gibi, ids_gltr, ids_acpl, ids_prse, ids_rops)

# fread by objset, filter to these ids

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

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point")
tree_df_agg_8b <- getTreePlanetSpectra(ids_all)
setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_extract/tree_outputs_point")
tree_df_agg_4b <- getTreePlanetSpectra(ids_all)

# Merge and remove duplicates

tree_df_agg_4b$objdate <- paste(tree_df_agg_4b$Object_ID, tree_df_agg_4b$date, sep = "_")
tree_df_agg_8b$objdate <- paste(tree_df_agg_8b$Object_ID, tree_df_agg_8b$date, sep = "_")

dupes <- which(tree_df_agg_4b$objdate %in% tree_df_agg_8b$objdate)
tree_df_agg_4b_sub <- tree_df_agg_4b[!dupes,]

tree_df_agg <- bind_rows(tree_df_agg_4b_sub, tree_df_agg_8b)

tree_df_agg_new <- tree_df_agg[which(tree_df_agg$objdate %notin% tree_df_agg_4b$objdate),] # some 8 band data wasn't in the 4 band
tree_df_agg_4b_new <- tree_df_agg_4b[which(tree_df_agg_4b$objdate %notin% tree_df_agg$objdate),] # all 4 band is included in the merged df

# Label tree species
tree_df_agg$species <- NA
tree_df_agg$species[which(tree_df_agg$Object_ID %in% ids_acru)] <- "ACRU"
tree_df_agg$species[which(tree_df_agg$Object_ID %in% ids_gibi)] <- "GIBI"
tree_df_agg$species[which(tree_df_agg$Object_ID %in% ids_gltr)] <- "GLTR"
tree_df_agg$species[which(tree_df_agg$Object_ID %in% ids_acpl)] <- "ACPL"
tree_df_agg$species[which(tree_df_agg$Object_ID %in% ids_prse)] <- "PRSE"
tree_df_agg$species[which(tree_df_agg$Object_ID %in% ids_rops)] <- "ROPS"
tree_df_agg$Object_ID %>% unique() %>% length()

tree_df_agg <- subset(tree_df_agg, select = -objdate)
col_order <- c("date", "Object_ID", "species", "coastal_blue", "blue", "green_i", "green", "yellow", "red", "rededge", "nir")
tree_df_agg <- as_tibble(tree_df_agg) # convert to tibble for col reordering
tree_df_agg <- tree_df_agg[, col_order]

# write this out
write.csv(tree_df_agg, "/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/sampled_trees_for_color_spectra.csv", row.names = FALSE)

#####
# get phenology data too, can use aggregated file

tree_pheno <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")
tree_pheno <- tree_pheno %>% filter(Poly_ID %in% ids_all)

#####
# Plotting with new data

tree_df_agg <- read.csv("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/sampled_trees_for_color_spectra.csv")

tree_df_agg$date <- tree_df_agg$date %>% ymd()
#tree_df_agg_sub <- tree_df_agg %>% filter(month(date) %in% c(3,4,5,6) & Object_ID %in% c(494589, 1198795)) # PRSE: 494589, ROSP: 1198795

#
tree_df_agg_sub <- tree_df_agg %>% filter(month(date) %in% c(3,4,5,6) &
                                            year(date) %in% 2017:2024 &
                                            Object_ID %in% ids_acru)

# try min max scaling
col_min <- min(c(min(tree_df_agg_sub$red), min(tree_df_agg_sub$green), min(tree_df_agg_sub$blue)))
col_max <- max(c(max(tree_df_agg_sub$red), max(tree_df_agg_sub$green), max(tree_df_agg_sub$blue)))

# try plotting the actual color
tree_df_agg_sub <- tree_df_agg_sub %>% mutate(rgb_hex = rgb((red - col_min)/(col_max - col_min)*255, 
                                                            (green - col_min)/(col_max - col_min)*255, 
                                                            (blue - col_min)/(col_max - col_min)*255, maxColorValue = 255))

tree_df_agg_sub$ndvi <- (tree_df_agg_sub$nir - tree_df_agg_sub$red)/(tree_df_agg_sub$nir + tree_df_agg_sub$red)

ggplot(tree_df_agg_sub) +
  geom_point(aes(x = yday(date), y = year(date), color = rgb_hex, size = ndvi)) +
  scale_color_identity() +
  facet_wrap(~Object_ID, nrow = 2) +
  scale_x_continuous(breaks = seq(60,180,10))

# need to do these at 10% and 90% rather than 100% (or more), and set everything beyond that to 255
# Need to do this for each tree uniquely? Maybe too dark for some trees rather than others...
r <- tree_df_agg_sub$red
g <- tree_df_agg_sub$green
b <- tree_df_agg_sub$blue

qmax <- 0.90
qmin <- 0.10
col_max <- max(quantile(r, qmax), quantile(g, qmax), quantile(b, qmax))
col_min <- min(quantile(r, qmin), quantile(g, qmin), quantile(b, qmin))

r_s <- (r - col_min)/(col_max - col_min)*255
g_s <- (g - col_min)/(col_max - col_min)*255
b_s <- (b - col_min)/(col_max - col_min)*255

scale255 <- function(x){
  x[x > 255] <- 255
  x[x < 0] <- 0
  return(x)
}

r_s <- scale255(r_s)
g_s <- scale255(g_s)
b_s <- scale255(b_s)


# try plotting the actual color??
tree_df_agg_sub <- tree_df_agg_sub %>% mutate(rgb_hex = rgb(r_s, g_s, b_s, maxColorValue = 255))

tree_df_agg_sub$ndvi <- (tree_df_agg_sub$nir - tree_df_agg_sub$red)/(tree_df_agg_sub$nir + tree_df_agg_sub$red)

tree_pheno_sub <- tree_pheno %>% filter(Poly_ID %in% unique(tree_df_agg_sub$Object_ID))
tree_pheno_sub$species[tree_pheno_sub$species == "Prunus serrulata"] <- "PRSE"
tree_pheno_sub$species[tree_pheno_sub$species == "Robinia pseudoacacia"] <- "ROPS"
tree_pheno_sub$species[tree_pheno_sub$species == "Ginkgo biloba"] <- "GIBI"
tree_pheno_sub$species[tree_pheno_sub$species == "Acer rubrum"] <- "ACRU"
tree_pheno_sub$species[tree_pheno_sub$species == "Acer platanoides"] <- "ACPL"
tree_pheno_sub$species[tree_pheno_sub$species == "Gleditsia triacanthos"] <- "GLTR"

tree_df_agg_sub$unqname <- paste0(tree_df_agg_sub$species, "_", tree_df_agg_sub$Object_ID)
tree_pheno_sub$unqname <- paste0(tree_pheno_sub$species, "_", tree_pheno_sub$Poly_ID)

ggplot(tree_df_agg_sub) +
  geom_hline(aes(yintercept = year(date)), linewidth = 0.5, col = "gray70") +
  geom_point(aes(x = yday(date), y = year(date), fill = rgb_hex, size = ndvi), shape = 21) +
  geom_point(data = tree_pheno_sub, aes(x = SOS_50, y = Year), shape = 3, size = 4) +
  scale_fill_identity() +
  facet_wrap(~unqname) +
  scale_x_continuous(breaks = seq(60,180,10)) +
  scale_y_continuous(breaks = seq(2017, 2024)) +
  labs(x = "Day of Year", y = "Year", size = "NDVI") +
  theme_bw()

ggsave("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/figures/tree_color_compare_spring_acru_all.png",
       width = 20, height = 10, units = "in")

# might need to change the lighting adjustment
#####

# try for fall

tree_df_agg_sub <- tree_df_agg %>% filter(month(date) %in% c(9,10,11,12) & Object_ID %in% c(63612, 1465337, 1319566, 1559444)) 


# try min max scaling
col_min <- min(c(min(tree_df_agg_sub$red), min(tree_df_agg_sub$green), min(tree_df_agg_sub$blue)))
col_max <- max(c(max(tree_df_agg_sub$red), max(tree_df_agg_sub$green), max(tree_df_agg_sub$blue)))

# try plotting the actual color??
tree_df_agg_sub <- tree_df_agg_sub %>% mutate(rgb_hex = rgb((red - col_min)/(col_max - col_min)*255, 
                                                            (green - col_min)/(col_max - col_min)*255, 
                                                            (blue - col_min)/(col_max - col_min)*255, maxColorValue = 255))

tree_df_agg_sub$ndvi <- (tree_df_agg_sub$nir - tree_df_agg_sub$red)/(tree_df_agg_sub$nir + tree_df_agg_sub$red)

ggplot(tree_df_agg_sub) +
  geom_point(aes(x = yday(date), y = year(date), color = rgb_hex, size = ndvi)) +
  scale_color_identity() +
  facet_wrap(~Object_ID, nrow = 2) +
  scale_x_continuous(breaks = seq(60,180,10))

# need to do these at 10% and 90% rather than 100% (or more), and set everything beyond that to 255
# Need to do this for each tree uniquely? Maybe too dark for some trees rather than others...
r <- tree_df_agg_sub$red
g <- tree_df_agg_sub$green
b <- tree_df_agg_sub$blue

qmax <- 0.90
qmin <- 0.10
col_max <- max(quantile(r, qmax), quantile(g, qmax), quantile(b, qmax))
col_min <- min(quantile(r, qmin), quantile(g, qmin), quantile(b, qmin))

r_s <- (r - col_min)/(col_max - col_min)*255
g_s <- (g - col_min)/(col_max - col_min)*255
b_s <- (b - col_min)/(col_max - col_min)*255

scale255 <- function(x){
  x[x > 255] <- 255
  x[x < 0] <- 0
  return(x)
}

r_s <- scale255(r_s)
g_s <- scale255(g_s)
b_s <- scale255(b_s)


# color plotting for fall
tree_df_agg_sub <- tree_df_agg_sub %>% mutate(rgb_hex = rgb(r_s, g_s, b_s, maxColorValue = 255))

tree_df_agg_sub$ndvi <- (tree_df_agg_sub$nir - tree_df_agg_sub$red)/(tree_df_agg_sub$nir + tree_df_agg_sub$red)

tree_pheno_sub <- tree_pheno %>% filter(Poly_ID %in% unique(tree_df_agg_sub$Object_ID))
tree_pheno_sub$species[tree_pheno_sub$species == "Prunus serrulata"] <- "PRSE"
tree_pheno_sub$species[tree_pheno_sub$species == "Robinia pseudoacacia"] <- "ROPS"
tree_pheno_sub$species[tree_pheno_sub$species == "Ginkgo biloba"] <- "GIBI"
tree_pheno_sub$species[tree_pheno_sub$species == "Acer rubrum"] <- "ACRU"
tree_pheno_sub$species[tree_pheno_sub$species == "Acer platanoides"] <- "ACPL"
tree_pheno_sub$species[tree_pheno_sub$species == "Gleditsia triacanthos"] <- "GLTR"

tree_df_agg_sub$unqname <- paste0(tree_df_agg_sub$species, "_", tree_df_agg_sub$Object_ID)
tree_pheno_sub$unqname <- paste0(tree_pheno_sub$species, "_", tree_pheno_sub$Poly_ID)

ggplot(tree_df_agg_sub) +
  geom_hline(aes(yintercept = year(date)), linewidth = 0.5, col = "gray70") +
  geom_point(aes(x = yday(date), y = year(date), fill = rgb_hex, size = ndvi), shape = 21) +
  geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = Year), shape = 3, size = 4) +
  scale_fill_identity() +
  facet_wrap(~unqname, nrow = 3) +
  scale_x_continuous(breaks = seq(240,360,10)) +
  scale_y_continuous(breaks = seq(2017, 2024)) +
  labs(x = "Day of Year", y = "Year", size = "NDVI") +
  theme_bw()
ggsave("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/figures/tree_color_compare_fall_acpl_acru_gibi_gltr.png",
       width = 16, height = 10, units = "in", dpi = 400)


#####
# try setting this up for individual trees to check for consistency across individuals, could try tileplot too


tree_df_agg_sub <- tree_df_agg %>% filter(month(date) %in% c(9,10,11,12) &
                                            year(date) %in% 2017:2024 &
                                            Object_ID %in% ids_gibi) 

# try min max scaling
col_min <- min(c(min(tree_df_agg_sub$red), min(tree_df_agg_sub$green), min(tree_df_agg_sub$blue)))
col_max <- max(c(max(tree_df_agg_sub$red), max(tree_df_agg_sub$green), max(tree_df_agg_sub$blue)))

# try plotting the actual color??
tree_df_agg_sub <- tree_df_agg_sub %>% mutate(rgb_hex = rgb((red - col_min)/(col_max - col_min)*255, 
                                                            (green - col_min)/(col_max - col_min)*255, 
                                                            (blue - col_min)/(col_max - col_min)*255, maxColorValue = 255))

tree_df_agg_sub$ndvi <- (tree_df_agg_sub$nir - tree_df_agg_sub$red)/(tree_df_agg_sub$nir + tree_df_agg_sub$red)

ggplot(tree_df_agg_sub) +
  geom_point(aes(x = yday(date), y = year(date), color = rgb_hex, size = ndvi)) +
  scale_color_identity() +
  facet_wrap(~Object_ID, nrow = 2) +
  scale_x_continuous(breaks = seq(60,180,10))

# need to do these at 10% and 90% rather than 100% (or more), and set everything beyond that to 255
# Need to do this for each tree uniquely? Maybe too dark for some trees rather than others...
r <- tree_df_agg_sub$red
g <- tree_df_agg_sub$green
b <- tree_df_agg_sub$blue

qmax <- 0.90
qmin <- 0.10
col_max <- max(quantile(r, qmax), quantile(g, qmax), quantile(b, qmax))
col_min <- min(quantile(r, qmin), quantile(g, qmin), quantile(b, qmin))

r_s <- (r - col_min)/(col_max - col_min)*255
g_s <- (g - col_min)/(col_max - col_min)*255
b_s <- (b - col_min)/(col_max - col_min)*255

scale255 <- function(x){
  x[x > 255] <- 255
  x[x < 0] <- 0
  return(x)
}

r_s <- scale255(r_s)
g_s <- scale255(g_s)
b_s <- scale255(b_s)


# color plotting for fall
tree_df_agg_sub <- tree_df_agg_sub %>% mutate(rgb_hex = rgb(r_s, g_s, b_s, maxColorValue = 255))

tree_df_agg_sub$ndvi <- (tree_df_agg_sub$nir - tree_df_agg_sub$red)/(tree_df_agg_sub$nir + tree_df_agg_sub$red)

tree_pheno_sub <- tree_pheno %>% filter(Poly_ID %in% unique(tree_df_agg_sub$Object_ID))
tree_pheno_sub$species[tree_pheno_sub$species == "Prunus serrulata"] <- "PRSE"
tree_pheno_sub$species[tree_pheno_sub$species == "Robinia pseudoacacia"] <- "ROPS"
tree_pheno_sub$species[tree_pheno_sub$species == "Ginkgo biloba"] <- "GIBI"
tree_pheno_sub$species[tree_pheno_sub$species == "Acer rubrum"] <- "ACRU"
tree_pheno_sub$species[tree_pheno_sub$species == "Acer platanoides"] <- "ACPL"
tree_pheno_sub$species[tree_pheno_sub$species == "Gleditsia triacanthos"] <- "GLTR"

tree_df_agg_sub$unqname <- paste0(tree_df_agg_sub$species, "_", tree_df_agg_sub$Object_ID)
tree_pheno_sub$unqname <- paste0(tree_pheno_sub$species, "_", tree_pheno_sub$Poly_ID)

ggplot(tree_df_agg_sub) +
  geom_hline(aes(yintercept = year(date)), linewidth = 0.5, col = "gray70") +
  geom_point(aes(x = yday(date), y = year(date), fill = rgb_hex, size = ndvi), shape = 21) +
  geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = Year), shape = 3, size = 4) +
  scale_fill_identity() +
  facet_wrap(~unqname) +
  scale_x_continuous(breaks = seq(240,360,10)) +
  scale_y_continuous(breaks = seq(2017, 2024)) +
  labs(x = "Day of Year", y = "Year", size = "NDVI") +
  theme_bw()

# ggsave("/Users/dlm356/dlm356_files/nyc_trees/tree_color_sniping/figures/tree_color_compare_fall_prse_all.png",
#        width = 16, height = 10, units = "in", dpi = 400)

tree_df_agg_sub$Year <- year(tree_df_agg_sub$date)

ggplot(tree_df_agg_sub) +
  geom_tile(aes(x = yday(date), y = unqname, fill = rgb_hex)) +
  geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = unqname), shape = 3, size = 2) +
  scale_fill_identity() +
  #facet_wrap(~as.factor(year(date))) +
  facet_wrap(~as.factor(Year))
  #scale_x_continuous(breaks = seq(240,360,10)) +
  #scale_y_continuous(breaks = seq(2017, 2024)) +
  labs(x = "Day of Year", y = "Tree ID") +
  theme_bw()

ggplot(tree_df_agg_sub) +
  geom_hline(aes(yintercept = year(date)), linewidth = 0.5, col = "gray70") +
  geom_point(aes(x = yday(date), y = unqname, fill = rgb_hex, size = ndvi), shape = 21) +
  geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = unqname), shape = 3, size = 4) +
  scale_fill_identity() +
  facet_wrap(~as.factor(Year)) +
  #scale_x_continuous(breaks = seq(240,360,10)) +
  #scale_y_continuous(breaks = seq(2017, 2024)) +
  labs(x = "Day of Year", y = "Year", size = "NDVI") +
  theme_bw()

# Do an interpolated color plot for clarity of visualization??


daily_dates <- seq(ymd("2017-01-01"), ymd("2024-12-31"), by = "days")
daily_dates <- daily_dates[which(month(daily_dates) %in% c(9, 10, 11, 12))]

# go between first and last date for a given year within this range
unq_ids <- unique(tree_df_agg_sub$unqname)
for (u in 1:length(unq_ids)){
  tr_id <- unq_ids[u]
  color_interp <- rep(NA, length(daily_dates))
  df_daily_colors <- cbind.data.frame(daily_dates, color_interp)
  for (yr in 2017:2024){
    #print(yr)
    ex_tree <- tree_df_agg_sub %>% filter(unqname == tr_id & year(date) == yr)
    ex_tree <- ex_tree[order(ex_tree$date),]
    for (i in 1:(nrow(ex_tree)-1)){
      #print(i)
      date_interval <- seq(ex_tree$date[i], ex_tree$date[i+1], by = "days")
      cp <- colorRampPalette(c(ex_tree$rgb_hex[i], ex_tree$rgb_hex[i+1]))
      df_daily_colors$color_interp[which(df_daily_colors$daily_dates %in% date_interval)] <- cp(length(date_interval))
    }
  }
  df_daily_colors$unqname <- tr_id
  if (u == 1) {
    df_daily_colors_all <- df_daily_colors
  } else {
    df_daily_colors_all <- bind_rows(df_daily_colors_all, df_daily_colors)
  }
}


df_daily_colors_all$Year <- year(df_daily_colors_all$daily_dates)

# test plot
ggplot(df_daily_colors_all) +
  geom_tile(aes(x = yday(daily_dates), y = unqname, fill = color_interp)) +
  scale_fill_identity() +
  geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = unqname), shape = 3, size = 2, color = "magenta") +
  facet_wrap(~as.factor(Year)) +
  #xlim(320,340) +
  labs(x = "Day of Year", y = "Tree ID")

#####

# yellow peak?
tree_df_agg_sub2 <- tree_df_agg_sub %>% mutate(ypeak = (yellow - (green + red)/2)/(yellow + (green + red)/2))
ggplot(tree_df_agg_sub2) +
  geom_point(aes(x = yday(date), y = ypeak, color = as.factor(Object_ID))) + 
  facet_wrap(~as.factor(year(date)))
# no obvious pattern here, would need other bands

#####

# hue color transform
library(plotwidgets)

hsl_interp <- col2hsl(df_daily_colors_all$color_interp)
df_daily_colors_all$hue <- hsl_interp[1,]
df_daily_colors_all$sat <- hsl_interp[2,]
df_daily_colors_all$lit <- hsl_interp[3,]
# hue of 60 "degrees" is yellow

ggplot(df_daily_colors_all) +
  geom_point(aes(x = yday(daily_dates), y = hue-60, color = color_interp, alpha = lit)) +
  scale_color_identity() +
  #geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = unqname), shape = 3, size = 2, color = "magenta") +
  facet_grid(rows = vars(as.factor(Year)), cols = vars(as.factor(unqname))) +
  geom_hline(yintercept = 0) +
  ylim(-50,50) +
  theme_bw()# +
  #labs(x = "Day of Year", y = "Tree ID")

#####
hsl <- col2hsl(tree_df_agg_sub$rgb_hex)
tree_df_agg_sub$hue <- hsl[1,]
tree_df_agg_sub$sat <- hsl[2,]
tree_df_agg_sub$lit <- hsl[3,]


ggplot(tree_df_agg_sub) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x = yday(date), y = hue-60, color = rgb_hex, alpha = sat)) +
  scale_color_identity() +
  #geom_point(data = tree_pheno_sub, aes(x = EOS_50, y = unqname), shape = 3, size = 2, color = "magenta") +
  facet_grid(rows = vars(as.factor(year(date))), cols = vars(as.factor(unqname))) +
  ylim(-50,50) +
  theme_bw()# +
#labs(x = "Day of Year", y = "Tree ID")

#####

# using indices 

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
}


tree_df_agg_sub2 <- tree_df_agg_sub %>% mutate(cci = calcNormDif(green_i, yellow),
                                               ndvi = calcNormDif(nir, red),
                                               pri = calcNormDif(green_i, green),
                                               ccr = calcNormDif(green_i, red))

tree_df_agg_sub2

ggplot(tree_df_agg_sub2) +
  geom_point(aes(x = yday(date), y = -cci, color = as.factor(Object_ID))) + 
  facet_grid(rows = vars(as.factor(year(date))), cols = vars(as.factor(unqname)))

#####
tree_df_agg2 <- tree_df_agg %>% mutate(cci = calcNormDif(green_i, yellow),
                                               ndvi = calcNormDif(nir, red),
                                               pri = calcNormDif(green_i, green),
                                               ccr = calcNormDif(green_i, red))

ggplot(tree_df_agg2) +
  geom_line(aes(x = yday(date), y = -cci, group = Object_ID)) + 
  xlim(240, 365) + ylim(0, 0.3) +
  facet_grid(rows = vars(as.factor(year(date))), cols = vars(as.factor(species)))


tree_df_agg2 %>% filter(date == ymd("2023-11-23")) %>%
  ggplot() + 
  geom_point(aes(x = -cci, y = ndvi, color = species))
# can identify 7/12 gibi based on cci < -0.15 by itself
#####

# test all available trees with this 

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
}

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point/")
file_list <- list.files(pattern = glob2rx("*_20231123_*.csv"))
df_all_spectra <- purrr::map_df(file_list, fread) # load them all in

tree_pheno_all <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
df_all_spectra_sub <- merge(df_all_spectra_sub, tree_pheno_all[,c("Poly_ID", "genus", "species", "Year", "EOS_50", "R2")], by.x = "Object_ID", by.y = "Poly_ID")
df_all_spectra_sub <- df_all_spectra_sub %>% filter(Year == 2023 & R2 > 0.7)


df_all_spectra_sub <- df_all_spectra_sub[!duplicated(df_all_spectra_sub),] # remove duplicate rows
df_all_spectra_sub <- df_all_spectra_sub %>% mutate(cci = calcNormDif(green_i, yellow),
                                       ndvi = calcNormDif(nir, red),
                                       pri = calcNormDif(green_i, green),
                                       ccr = calcNormDif(green_i, red))
gen_counts <- table(df_all_spectra_sub$genus)
sort(gen_counts, decreasing = TRUE)

df_all_spectra_sub <- df_all_spectra_sub %>% filter(genus %in% names(gen_counts[which(gen_counts > 1000)]))

ggplot(df_all_spectra_sub) +
  geom_boxplot(aes(x = genus, y = cci), outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(df_all_spectra_sub) +
  geom_boxplot(aes(x = genus, y = ndvi), outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(df_all_spectra_sub) +
  geom_boxplot(aes(x = genus, y = ccr), outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(df_all_spectra_sub) +
  geom_boxplot(aes(x = genus, y = -cci*1/ndvi*1/(-ccr)), outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#####
# try all combinations of bands as normalized difference indices

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

ggplot(df_spectra_omnbr) +
  geom_boxplot(aes(x = genus, y = nd_yellow_green_i), outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# comparison vars
comp_list <- c("genus", band_list, colname_list)

df_cls <- df_spectra_omnbr[, comp_list]

# try decision tree
# based on: https://www.geeksforgeeks.org/decision-tree-in-r-programming/

# with all genera first, but will want to split into "Ginkgo" and "Not Ginkgo"



set.seed(14)

df_cls <- transform(df_cls, 
                     coastal_blue = as.numeric(coastal_blue),
                     blue = as.numeric(blue),
                     green_i = as.numeric(green_i),
                     green = as.numeric(green),
                     yellow = as.numeric(yellow),
                     red = as.numeric(red),
                     rededge = as.numeric(rededge),
                     nir = as.numeric(nir))

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


ggplot(df_cls_sub) +
  geom_boxplot(aes(x = genus, y = yellow), outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


p <- predict(tree_model, test_data[,2:ncol(test_data)], type = 'class')
confusionMatrix(p, reference = test_data$genus)
# overall accuracy 72%, kappa = 0.43. Will probably get better with time component too

# Visualize these two indices in the data

ggplot(df_cls_sub) +
  geom_point(aes(x = nd_green_blue, y = nd_nir_yellow, color = genus), shape = 1, alpha = 0.5) +
  facet_wrap(~genus) +
  #scale_y_continuous(breaks = seq(-0.05, 0.65, 0.05)) +
  geom_vline(xintercept = c(0.15, 0.094)) +
  geom_hline(yintercept = c(0.41, 0.33, 0.21)) +
  theme_bw()

# try it with a random forest

rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
confusionMatrix(p_rf$predictions, reference = test_data$genus)


#####
# try different subsets with higher quality, bigger trees

calcNormDif <- function(b1, b2){
  return((b1-b2)/(b1+b2))
}

setwd("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_8b_highsunonly_cal_extract/tree_outputs_point/")
#file_list <- list.files(pattern = glob2rx("*_20231123_*.csv"))
file_list <- list.files(pattern = glob2rx("*_20231127_*.csv")) # try with the 11/27 image that includes Staten Island
df_all_spectra <- purrr::map_df(file_list, fread) # load them all in

tree_pheno_all <- fread("/Volumes/NYC_geo/Planet/tests/nyc_daily_stack_4b_highsunonly_cal_pheno/tree_pheno_pointextract_polyid_all_output_2017_2024.csv")

tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023 & 
                                                R2 > 0.7)# & 
                                                #tpstructur == "Full" &
                                                #tpconditio %in% c("Good", "Excellent") &
                                                #dbh > 18)

df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
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

comp_list <- c("genus", band_list, colname_list)

df_cls <- df_spectra_omnbr[, comp_list]


set.seed(14)

df_cls <- transform(df_cls, 
                    coastal_blue = as.numeric(coastal_blue),
                    blue = as.numeric(blue),
                    green_i = as.numeric(green_i),
                    green = as.numeric(green),
                    yellow = as.numeric(yellow),
                    red = as.numeric(red),
                    rededge = as.numeric(rededge),
                    nir = as.numeric(nir))

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


#####
# try setting this up in a filtering loop

# DBH

vals_range <- 0:36
oa <- rep(0, length(vals_range))
kappa <- rep(0, length(vals_range))
n_ginkgo <- rep(0, length(vals_range))

df_acc <- cbind.data.frame(vals_range, oa, kappa, n_ginkgo)

for (ind in 1:nrow(df_acc)){
  val <- df_acc$vals_range[ind]
  print(val)
  tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023 & 
                                                  R2 > 0.7 & 
                                                  #tpstructur == "Full" &
                                                  #tpconditio %in% c("Good", "Excellent") &
                                                  dbh > val)
                                                  #Height > 20)
  
  df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
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
      #print(band_name)
      
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
  
  comp_list <- c("genus", band_list, colname_list)
  
  df_cls <- df_spectra_omnbr[, comp_list]
  
  
  set.seed(14)
  
  df_cls <- transform(df_cls, 
                      coastal_blue = as.numeric(coastal_blue),
                      blue = as.numeric(blue),
                      green_i = as.numeric(green_i),
                      green = as.numeric(green),
                      yellow = as.numeric(yellow),
                      red = as.numeric(red),
                      rededge = as.numeric(rededge),
                      nir = as.numeric(nir))
  
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
  
  rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
  p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
  cm <- confusionMatrix(p_rf$predictions, reference = test_data$genus)
  
  df_acc$oa[ind] <- cm$overall[1]
  df_acc$kappa[ind] <- cm$overall[2]
  df_acc$n_ginkgo[ind] <- nrow(ginkgo)
}

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = oa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum DBH (in)", y = "Overall Accuracy")

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = kappa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum DBH (in)", y = "Kappa")

#####
# Height

vals_range <- seq(0,100,2)
oa <- rep(0, length(vals_range))
kappa <- rep(0, length(vals_range))
n_ginkgo <- rep(0, length(vals_range))

df_acc <- cbind.data.frame(vals_range, oa, kappa, n_ginkgo)

for (ind in 1:nrow(df_acc)){
  val <- df_acc$vals_range[ind]
  print(val)
  tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023 & 
                                                  R2 > 0.7 & 
                                                  #tpstructur == "Full" &
                                                  #tpconditio %in% c("Good", "Excellent") &
                                                  Height > val)
  
  df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
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
      #print(band_name)
      
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
  
  comp_list <- c("genus", band_list, colname_list)
  
  df_cls <- df_spectra_omnbr[, comp_list]
  
  
  set.seed(14)
  
  df_cls <- transform(df_cls, 
                      coastal_blue = as.numeric(coastal_blue),
                      blue = as.numeric(blue),
                      green_i = as.numeric(green_i),
                      green = as.numeric(green),
                      yellow = as.numeric(yellow),
                      red = as.numeric(red),
                      rededge = as.numeric(rededge),
                      nir = as.numeric(nir))
  
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
  
  rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
  p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
  cm <- confusionMatrix(p_rf$predictions, reference = test_data$genus)
  
  df_acc$oa[ind] <- cm$overall[1]
  df_acc$kappa[ind] <- cm$overall[2]
  df_acc$n_ginkgo[ind] <- nrow(ginkgo)
}

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = oa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum Height (ft)", y = "Overall Accuracy")

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = kappa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum Height (ft)", y = "Kappa")

#####
# Radius

vals_range <- seq(0,40,1)
oa <- rep(0, length(vals_range))
kappa <- rep(0, length(vals_range))
n_ginkgo <- rep(0, length(vals_range))

df_acc <- cbind.data.frame(vals_range, oa, kappa, n_ginkgo)

for (ind in 1:nrow(df_acc)){
  val <- df_acc$vals_range[ind]
  print(val)
  tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023 & 
                                                  R2 > 0.7 & 
                                                  #tpstructur == "Full" &
                                                  #tpconditio %in% c("Good", "Excellent") &
                                                  Radius > val)
  
  df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
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
      #print(band_name)
      
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
  
  comp_list <- c("genus", band_list, colname_list)
  
  df_cls <- df_spectra_omnbr[, comp_list]
  
  
  set.seed(14)
  
  df_cls <- transform(df_cls, 
                      coastal_blue = as.numeric(coastal_blue),
                      blue = as.numeric(blue),
                      green_i = as.numeric(green_i),
                      green = as.numeric(green),
                      yellow = as.numeric(yellow),
                      red = as.numeric(red),
                      rededge = as.numeric(rededge),
                      nir = as.numeric(nir))
  
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
  
  rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
  p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
  cm <- confusionMatrix(p_rf$predictions, reference = test_data$genus)
  
  df_acc$oa[ind] <- cm$overall[1]
  df_acc$kappa[ind] <- cm$overall[2]
  df_acc$n_ginkgo[ind] <- nrow(ginkgo)
}

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = oa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum Crown Radius (ft)", y = "Overall Accuracy")

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = kappa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum Crown Radius (ft)", y = "Kappa")


#####
# Shape area

# Radius

vals_range <- seq(0,5000,50)
oa <- rep(0, length(vals_range))
kappa <- rep(0, length(vals_range))
n_ginkgo <- rep(0, length(vals_range))

df_acc <- cbind.data.frame(vals_range, oa, kappa, n_ginkgo)

for (ind in 1:nrow(df_acc)){
  val <- df_acc$vals_range[ind]
  print(val)
  tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023 & 
                                                  R2 > 0.7 & 
                                                  #tpstructur == "Full" &
                                                  #tpconditio %in% c("Good", "Excellent") &
                                                  SHAPE_Area > val)
  
  df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
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
      #print(band_name)
      
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
  
  comp_list <- c("genus", band_list, colname_list)
  
  df_cls <- df_spectra_omnbr[, comp_list]
  
  
  set.seed(14)
  
  df_cls <- transform(df_cls, 
                      coastal_blue = as.numeric(coastal_blue),
                      blue = as.numeric(blue),
                      green_i = as.numeric(green_i),
                      green = as.numeric(green),
                      yellow = as.numeric(yellow),
                      red = as.numeric(red),
                      rededge = as.numeric(rededge),
                      nir = as.numeric(nir))
  
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
  
  rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
  p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
  cm <- confusionMatrix(p_rf$predictions, reference = test_data$genus)
  
  df_acc$oa[ind] <- cm$overall[1]
  df_acc$kappa[ind] <- cm$overall[2]
  df_acc$n_ginkgo[ind] <- nrow(ginkgo)
}

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = oa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum Shape Area (sq. ft)", y = "Overall Accuracy")

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = kappa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  labs(x = "Minimum Shape Area (sq. ft)", y = "Kappa")

#####
# R2

vals_range <- seq(0, 0.95,0.01)
oa <- rep(0, length(vals_range))
kappa <- rep(0, length(vals_range))
n_ginkgo <- rep(0, length(vals_range))

df_acc <- cbind.data.frame(vals_range, oa, kappa, n_ginkgo)

for (ind in 1:nrow(df_acc)){
  val <- df_acc$vals_range[ind]
  print(val)
  tree_pheno_input <- tree_pheno_all %>% filter(Year == 2023 & 
                                                  R2 > val) #& 
                                                  #tpstructur == "Full" &
                                                  #tpconditio %in% c("Good", "Excellent") &
                                                  #Radius > val)
  
  df_all_spectra_sub <- df_all_spectra %>% filter(Object_ID %in% tree_pheno_all$Poly_ID)
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
      #print(band_name)
      
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
  
  comp_list <- c("genus", band_list, colname_list)
  
  df_cls <- df_spectra_omnbr[, comp_list]
  
  
  set.seed(14)
  
  df_cls <- transform(df_cls, 
                      coastal_blue = as.numeric(coastal_blue),
                      blue = as.numeric(blue),
                      green_i = as.numeric(green_i),
                      green = as.numeric(green),
                      yellow = as.numeric(yellow),
                      red = as.numeric(red),
                      rededge = as.numeric(rededge),
                      nir = as.numeric(nir))
  
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
  
  rf_model <- ranger(genus ~ ., train_data, importance = "impurity")
  p_rf <- predict(rf_model, test_data[,2:ncol(test_data)])
  cm <- confusionMatrix(p_rf$predictions, reference = test_data$genus)
  
  df_acc$oa[ind] <- cm$overall[1]
  df_acc$kappa[ind] <- cm$overall[2]
  df_acc$n_ginkgo[ind] <- nrow(ginkgo)
}

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = oa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  labs(x = "Minimum R2", y = "Overall Accuracy")

ggplot(df_acc) +
  geom_point(aes(x = vals_range, y = kappa, color = log10(n_ginkgo))) +
  paletteer::scale_color_paletteer_c("grDevices::Viridis", direction = -1) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  labs(x = "Minimum R2", y = "Kappa")
