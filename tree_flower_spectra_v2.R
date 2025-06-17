# Read tree flower spectra data,
# merge into metadata
# convolve to PlanetScope bands

library(tidyverse)
library(reshape2)
library(readxl)

#setwd("2025-04-25_KatzLab")
#setwd("2025-05-02_KatzLab")
#setwd("2025-05-12_KatzLab")

readASDFile <- function(file) {
  df <- read.table(file, header = TRUE, sep = ",", skip = 33)  # adjust as needed
  return(df)
}

organizeASDData <- function(file_path, date_name){
  setwd(file_path)
  filelist <- list.files()
  
  # file reader
  data_list <- lapply(filelist, readASDFile)
  
  # Merge all by the first column (assumes same column name)
  merged_data <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), data_list)
  
  merged_data_t <- t(merged_data)
  wvl <- as.integer(merged_data_t[1,])
  
  asd_data <- as.data.frame(merged_data_t[2:nrow(merged_data_t),])
  colnames(asd_data) <- as.character(wvl)
  
  block_size <- 10
  blocks_vec <- rep(1:(nrow(asd_data)/block_size), each = block_size)
  asd_data_agg <- aggregate(asd_data, by = list(blocks_vec), FUN = "mean")
  asd_grouping <- rownames(asd_data)[seq(1, nrow(asd_data), block_size)] %>% str_split("[.]") %>% sapply( '[[', 1)
  asd_data_agg$Group.1 <- asd_grouping
  colnames(asd_data_agg)[1] <- "Sample_ID"
  
  asd_data_agg2 <- cbind.data.frame(rep(date_name, nrow(asd_data_agg)), 
                                    asd_data_agg)
  colnames(asd_data_agg2)[1] <- "Date"
  setwd("..") # move back up one directory level to reset, works only if going one level below top!
  return(asd_data_agg2)
}

#file_path <- "2025-04-25_KatzLab"
#date_name <- "20250425"
setwd("/Users/dlm356/dlm356_files/nyc_trees/tree_flower_spectra/")
asd_20250425 <- organizeASDData("2025-04-25_KatzLab", "20250425")
asd_20250502 <- organizeASDData("2025-05-02_KatzLab", "20250502")
asd_20250512 <- organizeASDData("2025-05-12_KatzLab", "20250512")
asd_20250606 <- organizeASDData("2025-06-06_KatzLab", "20250606")

######
# Aggregating to PlanetScope bands
planet_bands <- "/Users/dlm356/dlm356_files/nyc_trees/tree_flower_spectra/PlanetBands_PSBSD.csv"

# read spectral bands of Planet PSB.SD
bands <- read.csv(planet_bands)

asd2planet <- function(raw_df){
  raw_df_wvl <- as.numeric(colnames(raw_df)[3:length(raw_df)])
  
  raw_df_pre <- raw_df[,1:2]
  planet_df <- as.data.frame(matrix(NA, nrow = nrow(raw_df), ncol = length(bands$wvlngth)))
  colnames(planet_df) <- bands$wvlngth
  planet_df <- cbind.data.frame(raw_df_pre, planet_df)
  
  for (i in 1:nrow(bands)){
    band_inds <- which(raw_df_wvl %in% seq(bands$low_nm[i], bands$hi_nm[i])) + 2
    planet_df[, i+2] <- rowMeans(raw_df[,band_inds])
  }
  return(planet_df)
}

planet_20250425 <- asd2planet(asd_20250425)
planet_20250502 <- asd2planet(asd_20250502)
planet_20250512 <- asd2planet(asd_20250512)
planet_20250606 <- asd2planet(asd_20250606)

#####

asd_alldates <- rbind.data.frame(asd_20250425, asd_20250502, asd_20250512, asd_20250606)
planet_alldates <- rbind.data.frame(planet_20250425, planet_20250502, planet_20250512, planet_20250606)

# Merge data together with metadata

meta <- read_excel("/Users/dlm356/dlm356_files/nyc_trees/tree_flower_spectra/tree_asd_spectra_spring2025.xlsx", sheet = 1)

asd_alldates_meta <- cbind.data.frame(meta, asd_alldates[,3:ncol(asd_alldates)])
planet_alldates_meta <- cbind.data.frame(meta, planet_alldates[,3:ncol(planet_alldates)])

write.csv(asd_alldates_meta, "asd_alldates_meta.csv", row.names = FALSE)
write.csv(planet_alldates_meta, "planet_alldates_meta.csv", row.names = FALSE)
