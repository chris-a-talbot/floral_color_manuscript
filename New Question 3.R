library(dplyr)
library(stringr)
library(data.table)
library(phytools)

#Function to load HSL values from CSV files in a specified directory
load_hsl_values <- function(directory) {
  # List all files in the directory matching the pattern "_stats.csv"
  files <- list.files(directory, pattern = "_stats.csv$", full.names = TRUE)
  hsl_data <- lapply(files, function(file) {
    # Extract the species name from the file name
    species_name <- tools::file_path_sans_ext(basename(file))
    species_name <- str_remove(species_name, "_stats")
    # Read the CSV file into a data frame
    df <- fread(file)
    if ("HSL" %in% colnames(df)) {
      # Extract and convert the HSL values to numeric
      hsl_values <- str_remove_all(df$HSL[1], "\\[|\\]") %>% str_split(",") %>% unlist() %>% as.numeric()
      return(data.frame(species = str_replace(species_name, " ", "_"), H = hsl_values[1] * 2, S = hsl_values[2] / 255, L = hsl_values[3] / 255))
    }
  }) %>% bind_rows()
  return(hsl_data)
}

# Load HSL data from the specified directory
hsl_data <- load_hsl_values("C:/Users/chris/Documents/pfc/AppData/Output")
rownames(hsl_data) <- hsl_data$species
hsl_data <- hsl_data[, -1]
# Calculate Euclidean distance matrix for HSL data
hsl_dist <- as.matrix(dist(hsl_data, method = "euclidean"))

overlap_matrix <- readRDS("data/overlap_matrix.RDS")

phylogeny_main <- readRDS("data/phylogeny_read.rds")

common_species <- Reduce(intersect, list(phylogeny_main$tip.label, rownames(overlap_matrix), rownames(hsl_data)))

# Trim the phylogeny to include only common species
phylogeny <- drop.tip(phylogeny_main, phylogeny_main$tip.label[which(!phylogeny_main$tip.label %in% common_species)]) %>%
  multi2di() %>% force.ultrametric()

hsl_dist_c <- hsl_dist[common_species, common_species]
overlap_matrix_c <- overlap_matrix[common_species, common_species]
gen_dist_c <- cophenetic(phylogeny)

# Create all unordered pairs of species
species_pairs <- t(combn(phylogeny$tip.label, 2))

# Extract values and create data frame
pairs <- data.frame(
  species1 = species_pairs[,1],
  species2 = species_pairs[,2],
  gen_dist = sapply(1:nrow(species_pairs), function(i) 
    gen_dist_c[species_pairs[i,1], species_pairs[i,2]]),
  hsl_dist = sapply(1:nrow(species_pairs), function(i) 
    hsl_dist_c[species_pairs[i,1], species_pairs[i,2]]),
  overlap = sapply(1:nrow(species_pairs), function(i) 
    overlap_matrix_c[species_pairs[i,1], species_pairs[i,2]])
)

# Order the data frame by genetic distance
pairs <- pairs[order(pairs$gen_dist), ]

# Initialize the close_pairs dataframe
close_pairs <- data.frame()

# Continue until pairs is empty
while(nrow(pairs) > 0) {
  # Add the top row (lowest gen_dist) to close_pairs
  close_pairs <- rbind(close_pairs, pairs[1,])
  
  # Get the species from the selected pair
  species1 <- pairs$species1[1]
  species2 <- pairs$species2[1]
  
  # Remove all rows containing either of these species
  pairs <- pairs[!(pairs$species1 %in% c(species1, species2) | 
                     pairs$species2 %in% c(species1, species2)), ]
}

# Reset row names
rownames(close_pairs) <- NULL

close_pairs <- close_pairs[close_pairs$gen_dist <= 30,]
close_pairs <- close_pairs[close_pairs$overlap > 0,]

plot(data=close_pairs, overlap ~ gen_dist)
plot(data=close_pairs, hsl_dist ~ overlap)

# Create the multiple linear regression model
model <- lm(overlap ~ gen_dist + hsl_dist, data = close_pairs)

# View a summary of the model
summary(model)
