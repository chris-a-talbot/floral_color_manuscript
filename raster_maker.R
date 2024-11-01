################################################################################
## This script is used to generate presence/absence rasters from occurrences. ##
##                      Chris Talbot, March 2024                              ##
################################################################################

# Load pacman package for loading the other packages
if (!require("pacman")) {install.packages("pacman")}
library("pacman")

# Load required packages
p_load("data.table", "stringr", "dplyr", "parallel", "tools", "speciesRaster", "terra")

# Contains functions to create rasters
source("helper_raster.R")

# Run multi-core processes with max-1 cores
core_num = detectCores() - 1

# Get all of the occurrence data with only species names, lat, and lons
occ_trimmed = fread("./occurrences/occ_trimmed.csv")[,-1] 

# Get the full flowers data from the file with updated iNaturalist phenology
flowers = fread("./data/better_flowers.csv")

# Some processing is necessary because of the way occurrences were downloaded using inconsistent nomenclature
# Valid species in the occurrences may be named via the species or taxon columns in the main file
valid_species = c(word(flowers$taxon, 1, 2), flowers$species)

# Get occurrences for only the species which exist in the main data file
occ_trimmed = occ_trimmed[species %in% valid_species]

# Update the species column such that all occurrences are named by the taxon, 
# rather than species, column of the main data file
occ_trimmed[flowers, on = .(species = species), species := sapply(strsplit(i.taxon, " "), function(x) paste(x[1:min(2,length(x))], collapse = " ")), by = .EACHI]

# Restrict occurrences to those within the desired extent of analyses
occ_trimmed = occ_trimmed[decimalLatitude >= 38 & decimalLatitude <= 48 & decimalLongitude >= -92.5 & decimalLongitude <= -67.5]

# Get a list of species names
species_list = word(flowers$taxon, 1, 2)

# Get the phenological trait data
traits_pheno = cbind(species_list, flowers[,16:27])
colnames(traits_pheno) = c("species", colnames(traits_pheno)[-1])

# Get, and write, the rasters for each month, feb:nov, using the desired extent
# Resolution 0.09 ~ 10km^2 cells
# Change cores to 1 to run on Windows
rasters = get_monthly_rasters_from_occurrences(occ_trimmed, pheno=traits_pheno, 
                                               ext=c(-92.5, -67.5, 38, 48), 
                                               res=0.09, write=FALSE, cores=1)
