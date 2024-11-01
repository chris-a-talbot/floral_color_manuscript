# Load necessary libraries
library(sf)
library(data.table)
library(dplyr)
library(stringr)

# Load the species ranges from RDS
load("./RData/shapes_flat.RData")  # This should load the object(s) into the workspace

# Assuming the loaded object is named 'shapes_flat' or similar, assign it to 'species_ranges'
species_ranges <- shapes_flat  # Replace 'shapes_flat' with the actual variable name loaded

# Load the species list
flowers <- fread("G:/My Drive/Documents/Research/Weber/honors-thesis/data/flowers_master.csv")
species_list <- word(flowers$taxon, 1, 2)

# Subset the species ranges to only those in the species list
species_ranges <- species_ranges[names(species_ranges) %in% species_list]
species_names <- names(species_ranges)

# Ensure CRS of species ranges are consistent
target_crs <- st_crs(species_ranges[[1]])
species_ranges <- lapply(species_ranges, function(shape) st_transform(shape, target_crs))

# Create a bounding box for the specified coordinates
bbox <- st_as_sfc(st_bbox(c(xmin = -92.5, ymin = 38, xmax = -67.5, ymax = 48), crs = target_crs))

# Filter out species ranges that do not intersect with the bounding box
species_ranges <- lapply(species_ranges, function(shape) st_union(shape))
species_ranges <- species_ranges[sapply(species_ranges, function(shape) length(shape[[1]]) > 0)]
species_ranges <- lapply(species_ranges, function(shape) st_simplify(shape, dTolerance = 0.01)) # Simplify geometries with a small tolerance
species_ranges <- lapply(species_ranges, function(shape) st_make_valid(shape))
species_ranges <- species_ranges[sapply(species_ranges, function(shape) length(shape[[1]]) > 0)]
species_names <- names(species_ranges)

# Initialize data structures
valid_species_ranges <- list()
ssc_dt <- data.table(species1 = character(), species2 = character(), ssc = numeric())

# Function to calculate Szymkiewicz-Simpson coefficient
calculate_ssc <- function(shape1, shape2) {
  intersection <- st_intersection(shape1, shape2)
  if (nrow(intersection) == 0) {
    return(0)
  }
  area_intersection <- st_area(intersection)
  area1 <- st_area(shape1)
  area2 <- st_area(shape2)
  min_area <- pmin(area1, area2)
  ssc <- as.numeric(area_intersection / min_area)
  return(ssc)
}

# Check if a saved state exists
ssc_file <- "sykiewicz_simpson_coefficients.csv"
shapefiles_rds <- "valid_species_ranges.rds"
pairs_file <- "pairs.csv"

if (file.exists(ssc_file)) {
  ssc_dt <- fread(ssc_file)
  existing_pairs <- with(ssc_dt, paste(str_replace(species1, " ", "_"), str_replace(species2, " ", "_")))
} else {
  existing_pairs <- character()
}

if (file.exists(shapefiles_rds)) {
  valid_species_ranges <- readRDS(shapefiles_rds)
}

if (file.exists(pairs_file)) {
  pairs <- readRDS(pairs_file)
}

# Loop through species pairs and calculate SSC, saving progress periodically
save_interval <- 100  # Save progress every 100 comparisons
comparison_count <- 0
num_species <- length(species_ranges)

for (i in 1:(num_species - 1)) {
  for (j in (i + 1):num_species) {
    species_pair <- paste(str_replace(species_names[i], " ", "_"), str_replace(species_names[j], " ", "_"))
    if (!species_pair %in% pairs) {
      if (!species_names[i] %in% names(valid_species_ranges)) {
        valid_species_ranges[[species_names[i]]] <- species_ranges[[i]]
      }
      if (!species_names[j] %in% names(valid_species_ranges)) {
        valid_species_ranges[[species_names[j]]] <- species_ranges[[j]]
      }
      # Add error handling for invalid geometries
      tryCatch({
        ssc <- calculate_ssc(valid_species_ranges[[species_names[i]]], valid_species_ranges[[species_names[j]]])
        ssc_dt <- rbind(ssc_dt, data.table(species1 = str_replace(species_names[i], " ", "_"), species2 = str_replace(species_names[j], " ", "_"), ssc = ssc))
        # Add the symmetric pair
        ssc_dt <- rbind(ssc_dt, data.table(species1 = str_replace(species_names[j], " ", "_"), species2 = str_replace(species_names[i], " ", "_"), ssc = ssc))
        comparison_count <- comparison_count + 1
        if (comparison_count %% save_interval == 0) {
          fwrite(ssc_dt, ssc_file)
          saveRDS(pairs, pairs_file)
          saveRDS(valid_species_ranges, shapefiles_rds)
          cat("Progress saved after", comparison_count, "comparisons\n")
        }
      }, error = function(e) {
        cat("Error with species pair", species_pair, ":", e$message, "\n")
      })
      pairs <- c(pairs, species_pair)
    }
  }
}

# Save the final table and valid species ranges
fwrite(ssc_dt, ssc_file)
saveRDS(valid_species_ranges, shapefiles_rds)
saveRDS(pairs, pairs_file)

library(ape)
library(phytools)
library(vegan)

megatree = read.tree("data/ALLMB.tre")
phylo <- drop.tip(megatree, megatree$tip.label[which(!megatree$tip.label %in% str_replace(species_names, " ", "_"))])
phylo_resolved <- multi2di(phylo)
gen_dist <- cophenetic(phylo_resolved)

ssc_dt_short <- ssc_dt[species1 < species2]
# ssc_dt_short <- ssc_dt_short[ssc < 1]
# ssc_dt_short <- ssc_dt_short[ssc > 0]

# Convert the matrix to long format for easy merging
gen_dist_long <- melt(gen_dist, variable.name = "species2", value.name = "gd", id.vars = "species1")
gen_dist_long$Var1 <- as.character(gen_dist_long$Var1)
gen_dist_long$Var2 <- as.character(gen_dist_long$Var2)
gen_dist_long <- gen_dist_long[gen_dist_long$Var1 < gen_dist_long$Var2, ]

# Merge based on formatted species names
ssc_combined <- merge(ssc_dt_short, gen_dist_long, by.x = c("species1", "species2"), by.y = c("Var1", "Var2"), all.x = TRUE)

plot(data=ssc_combined, ssc~gd)
lin = lm(ssc~gd, data=ssc_combined)
abline(lin, col="red")
summary(lin)

# Create lists of unique species
unique_species <- unique(c(ssc_combined$species1, ssc_combined$species2))

# Prepare distance matrices for the Mantel test
# SSC distance matrix
ssc_dist_matrix <- matrix(NA, nrow = length(unique_species), ncol = length(unique_species))
rownames(ssc_dist_matrix) <- unique_species
colnames(ssc_dist_matrix) <- unique_species

for (i in 1:nrow(ssc_combined)) {
  ssc_dist_matrix[ssc_combined$species1[i], ssc_combined$species2[i]] <- ssc_combined$ssc[i]
  ssc_dist_matrix[ssc_combined$species2[i], ssc_combined$species1[i]] <- ssc_combined$ssc[i]
}

# Genetic distance matrix
gen_dist_matrix <- matrix(NA, nrow = length(unique_species), ncol = length(unique_species))
rownames(gen_dist_matrix) <- unique_species
colnames(gen_dist_matrix) <- unique_species

for (i in 1:nrow(ssc_combined)) {
  gen_dist_matrix[ssc_combined$species1[i], ssc_combined$species2[i]] <- ssc_combined$gd[i]
  gen_dist_matrix[ssc_combined$species2[i], ssc_combined$species1[i]] <- ssc_combined$gd[i]
}

# Convert to distance objects
ssc_dist <- as.dist(ssc_dist_matrix)
gen_dist <- as.dist(gen_dist_matrix)

# Perform the Mantel test
mantel_test <- mantel(ssc_dist, gen_dist, method = "pearson", na.rm=T)

# Print the results
print(mantel_test)










# Load necessary package
library(sf)

# Define the directory containing the shapefiles
shp_dir <- "./shp/pts/"

# List all shapefiles ending with "_0.09.shp"
shapefiles <- list.files(shp_dir, pattern = "_0.09.shp$", full.names = TRUE)
shapefiles <- shapefiles[1:50]
shapefile_names <- sub("_0.09.shp$", "", basename(shapefiles)) 
shapefile_names <- shapefile_names[1:50]

# Initialize lists to store the geometries and areas
shapefile_geoms <- list()
shapefile_areas <- list()

# Read all shapefiles and calculate their union and area
for (i in seq_along(shapefiles)) {
  shp <- st_union(read_sf(shapefiles[i]))
  shapefile_geoms[[i]] <- shp
  shapefile_areas[[i]] <- st_area(shp)
}

# Create an empty matrix to store the overlap coefficients
n <- length(shapefiles)
overlap_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(shapefile_names, shapefile_names))

saveRDS(shapefile_geoms, "./shp/pts/geoms.RDS")

# Function to calculate overlap coefficient between two geometries
calculate_overlap <- function(geom1, geom2, area1, area2) {
  intersection <- st_intersection(geom1, geom2)
  
  if(length(intersection) == 0) {
    return(NA)
  }
  
  # Check if intersection is empty
  if (st_is_empty(intersection)) {
    return(NA)
  }
  
  intersection_area <- st_area(intersection)
  min_area <- pmin(area1, area2)
  overlap_coefficient <- as.numeric(intersection_area / min_area)
  return(overlap_coefficient)
}

# Loop through all pairs of shapefiles
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    overlap_matrix[i, j] <- calculate_overlap(shapefile_geoms[[i]], shapefile_geoms[[j]], shapefile_areas[[i]], shapefile_areas[[j]])
  }
}

# Print the result
print(overlap_matrix)

