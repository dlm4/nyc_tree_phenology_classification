# Extract lidar intensity values for vegetation points within TNC tree crown polygons

library(terra)
library(sf)
library(lidR)
library(tidyverse)

# Load lidar data
# will need to apply across catalog
# likely will need to chunk this up by tree polygon id

# Load tnc polygon data gdb

# Load las file index poly

# Fix projections
# Intersect to get list of las files needed for a given subset area

# for each tree crown polygon, loop or catalog or alt parallel process
# Subset for just vegetation classified lidar points (low, medium, high)
# intensity
# To start: get mean of all high vegetation points
# Alt: get mean all vegetation points, 95th percentile, SD
# retain these values for each tree id
#

# Export dataframe with these intensity values