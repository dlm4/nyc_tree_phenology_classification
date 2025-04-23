library(tidyverse)
library(reshape2)

setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_flower_spectra/2025-04-22_KatzLab")

filelist <- list.files()

# spec <- read.table(filelist[1], skip = 33, header = TRUE, sep = ",")
# 
# ggplot(spec) +
#   geom_line(aes(Wavelength, `LOM01540.asd`))

# ChatGPT
# Function to read a file and keep only first column + data
read_file <- function(file) {
  df <- read.table(file, header = TRUE, sep = ",", skip = 33)  # adjust as needed
  df
}
# Read all files
data_list <- lapply(filelist, read_file)
# Merge all by the first column (assumes same column name)
merged_data <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), data_list)

merged_data_long <- melt(merged_data, id.vars = "Wavelength")

ggplot(merged_data_long) +
  geom_line(aes(Wavelength, value, col = variable))

#
# Campus trees ID'd with this map:
# https://cugir.library.cornell.edu/catalog/cugir-009100

labs_acpl <- paste0("LOM0", c(1540:1549, 1560:1579),".asd") # Acer platanoides, Norway maple
labs_prse2 <- paste0("LOM0", 1550:1559,".asd") # Prunus serrulata, Kwanzan Cherry, treeid: 5005
labs_maso <- paste0("LOM0", 1580:1609,".asd") # Magnolia x soulangiana, Saucer Magnolia, treeid: 4805

merged_data_long$spname <- NA
merged_data_long$spname[which(merged_data_long$variable %in% labs_acpl)] <- "Acer platanoides"
merged_data_long$spname[which(merged_data_long$variable %in% labs_prse2)] <- "Prunus serrulata"
merged_data_long$spname[which(merged_data_long$variable %in% labs_maso)] <- "Magnolia x soulangeana"

merged_data_long %>% #filter(Wavelength %in% 400:1100) %>%
ggplot() +
  geom_line(aes(x = Wavelength, y = value, col = spname, group = variable), linewidth = 0.5) +
  labs(x = "Wavelength (nm)", y = "Reflectance", col = "Tree Species") +
  theme_bw()

# Average by sample as recommended by AY
merged_data_long$sample_name <- NA

merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1540:1549,".asd"))] <- "ACPL_01"
merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1560:1569,".asd"))] <- "ACPL_02"
merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1570:1579,".asd"))] <- "ACPL_03"

merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1550:1559,".asd"))] <- "PRSE_01"

merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1580:1589,".asd"))] <- "MASO_01"
merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1590:1599,".asd"))] <- "MASO_02"
merged_data_long$sample_name[which(merged_data_long$variable %in% paste0("LOM0", 1600:1609,".asd"))] <- "MASO_03"

merged_data_long_means <- aggregate(merged_data_long$value, by = list(merged_data_long$Wavelength, merged_data_long$sample_name), FUN = "mean")
colnames(merged_data_long_means) <- c("wavelength", "sample_name", "refl")

merged_data_long_means %>% #filter(Wavelength %in% 400:1100) %>%
  ggplot() +
  geom_line(aes(x = wavelength, y = refl, col = sample_name, group = sample_name), linewidth = 0.5) +
  labs(x = "Wavelength (nm)", y = "Reflectance", col = "Sample Name") +
  theme_bw()



# Alex's ASD converter
current_df <- pivot_wider(merged_data_long_means, names_from = sample_name, values_from = refl) # needs wide format

planet_bands <- "/Users/dlm356/dlm356_files/nyc_trees/tree_flower_spectra/PlanetBands_PSBSD.csv"

# read spectral bands of Planet PSB.SD
bands = read.csv(planet_bands)

# Initialize a data frame to store average spectra for each band
avg_spectra_df = data.frame(Band = character(0), stringsAsFactors = FALSE)

# Iterate over each band
for (j in seq_len(nrow(bands))) {
  # Extract band information
  band_name = bands$BandName[j]
  low_end = bands$low_nm[j]
  high_end = bands$hi_nm[j]
  
  # Select spectral data within the band's wavelength range
  band_data = current_df[current_df$wavelength >= low_end & current_df$wavelength <= high_end, ]
  
  # Calculate the average spectra for each band
  avg_spectra = colMeans(band_data[, -1]) # Exclude the first column (wavelength) while calculating the mean
  
  # Create a data frame with band name and average spectra
  avg_spectra_df_band = data.frame(Band = band_name, t(avg_spectra))
  
  # Append to the output data frame
  avg_spectra_df = rbind(avg_spectra_df, avg_spectra_df_band)
}

colnames(avg_spectra_df)[-1] = colnames(current_df)[-1]

# back to my own scripting
avg_spectra_df_wvl <- merge(bands, avg_spectra_df, by.x = "BandName", by.y = "Band")
avg_spectra_df_wvl_long <- melt(avg_spectra_df_wvl, id.vars = c("BandName", "low_nm", "hi_nm", "wvlngth"))


avg_spectra_df_wvl_long %>% #filter(Wavelength %in% 400:1100) %>%
  ggplot(aes(x = wvlngth, y = value, col = variable, group = variable)) +
  geom_point() +
  geom_line(linewidth = 0.5) +
  xlim(300,1100) +
  labs(x = "Wavelength (nm)", y = "Reflectance", col = "Sample Name") +
  theme_bw()

ggplot() +
  geom_line(data = merged_data_long_means, aes(x = wavelength, y = refl, col = sample_name, group = sample_name), linewidth = 0.5, alpha = 0.2) +
  geom_point(data = avg_spectra_df_wvl_long, aes(x = wvlngth, y = value, col = variable, group = variable)) +
  geom_line(data = avg_spectra_df_wvl_long, aes(x = wvlngth, y = value, col = variable, group = variable), linewidth = 0.5) +
  ylim(-0.05, 0.9) +
  geom_text(data = avg_spectra_df_wvl_long[1:8,], aes(x = wvlngth, y = -0.04, label = BandName), angle = 90, size = 2.5) +
  scale_x_continuous(breaks = seq(300, 1100, 100)) +
  labs(x = "Wavelength (nm)", y = "Reflectance", col = "Sample Name") +
  theme_bw()
