# plot tree_flower_spectra

# Plotting script for overlaying ASD and Planet-resolution spectra

library(tidyverse)
library(reshape2)

# ASD file
asd_res <- read.csv('/Users/dlm356/Library/CloudStorage/Box-Box/dlm_box/ASD_tree_spectra/asd_alldates_meta.csv')

# "Planet" (8-band) file
planet_res <- read.csv('/Users/dlm356/Library/CloudStorage/Box-Box/dlm_box/ASD_tree_spectra/planet_alldates_meta.csv')

# note that read.csv() appends an X to the front of the numbered band center names

# Make long

asd_res_long <- pivot_longer(asd_res, cols = X325:X1075, names_to = "band")
asd_res_long$wvl <- str_split(asd_res_long$band, "X") %>% sapply(function(x) x[[2]]) %>% as.numeric() # I remember there being something more elegant than this


planet_res_long <- pivot_longer(planet_res, cols = X443:X865, names_to = "band")
planet_res_long$wvl <- str_split(planet_res_long$band, "X") %>% sapply(function(x) x[[2]]) %>% as.numeric() 

# ref data is in: tree_asd_spectra_spring2025.xlsx
# ID_2 is the unique identifier, subset
sel_ids_2 <- c("ACPL_20250425_01", # farmhouse (should be flowers)
             "ACPL_20250502_01", # farmhouse flowers with some leaves
             "ACPL_20250502_02") # farmhouse flowers with even more leaves

# specify which one of these samples because we take multiple spectra of samples (most of the time)
sample_ids <- c("LOM00080",
                "LOM00280",
                "LOM00290")

asd_res_long_sub <- asd_res_long %>% filter(ID_2 %in% sel_ids_2 & Sample_ID %in% sample_ids)
planet_res_long_sub <- planet_res_long %>% filter(ID_2 %in% sel_ids_2 & Sample_ID %in% sample_ids)


# plot
band_names <- c("CoastalBlue", "Blue", "GreenI", "Green", "Yellow", "Red", "RedEdge", "NIR") # planet band names for bottom labeling

ggplot() +
  geom_line(data = asd_res_long_sub, aes(x = wvl, y = value, col = ID_2, group = ID_2), linewidth = 0.5, alpha = 0.2) +
  geom_point(data = planet_res_long_sub, aes(x = wvl, y = value, col = ID_2, group = ID_2)) +
  geom_line(data = planet_res_long_sub, aes(x = wvl, y = value, col = ID_2, group = ID_2), linewidth = 0.5) +
  ylim(-0.05, 0.9) +
  geom_text(data = planet_res_long_sub[1:8,], aes(x = wvl, y = -0.04, label = band_names), angle = 45, size = 2) +
  scale_x_continuous(breaks = seq(300, 1100, 100)) +
  labs(x = "Wavelength (nm)", y = "Reflectance", col = "Sample Name") +
  theme_bw()

ggsave('/Users/dlm356/Library/CloudStorage/Box-Box/dlm_box/ASD_tree_spectra/acpl_example_plot.pdf',
       height = 6, width = 8, units = "in")
