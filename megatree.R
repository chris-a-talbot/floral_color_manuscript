library(ape)
library(phytools)
library(data.table)
library(stringr)
library(tidyverse)
library(phytools)
library(geiger)
library(vegan)
library(phyr)
library(terra)
library(speciesRaster)
library(grDevices)
library(dplyr)

megatree = read.tree("data/ALLMB.tre")

flowers <- fread("data/better_flowers.csv")
sp_names <- str_replace(word(flowers$taxon, 1, 2), " ", "_")
sp_names_backup <- str_replace(word(flowers$species, 1, 2), " ", "_")
sp_names_full <- c(sp_names, sp_names_backup)

phylo <- drop.tip(megatree, megatree$tip.label[which(!megatree$tip.label %in% sp_names_full)])
phylo_resolved <- multi2di(phylo)
gen_dist <- cophenetic(phylo_resolved)

# Function to load RGB values from CSV files
load_rgb_values <- function(directory) {
  files <- list.files(directory, pattern = "_stats.csv$", full.names = TRUE)
  rgb_data <- lapply(files, function(file) {
    species_name <- tools::file_path_sans_ext(basename(file))
    species_name <- str_remove(species_name, "_stats")
    df <- fread(file)
    if ("RGB" %in% colnames(df)) {
      # Strip non-numeric characters and convert to numeric vector
      rgb_values <- str_remove_all(df$RGB[1], "\\[|\\]") %>% str_split(",") %>% unlist() %>% as.numeric()
      return(data.frame(species = str_replace(species_name, " ", "_"), R = rgb_values[1], G = rgb_values[2], B = rgb_values[3]))
    }
  }) %>% bind_rows()
  return(rgb_data)
}

# Function to load RGB values from CSV files
load_hsl_values <- function(directory) {
  files <- list.files(directory, pattern = "_stats.csv$", full.names = TRUE)
  hsl_data <- lapply(files, function(file) {
    species_name <- tools::file_path_sans_ext(basename(file))
    species_name <- str_remove(species_name, "_stats")
    df <- fread(file)
    if ("HSL" %in% colnames(df)) {
      # Strip non-numeric characters and convert to numeric vector
      hsl_values <- str_remove_all(df$HSL[1], "\\[|\\]") %>% str_split(",") %>% unlist() %>% as.numeric()
      return(data.frame(species = str_replace(species_name, " ", "_"), H = hsl_values[1] * 2, S = hsl_values[2] / 255, L = hsl_values[3] / 255))
    }
  }) %>% bind_rows()
  return(hsl_data)
}

directory <- "C:/Users/chris/Documents/pfc/AppData/Output"
rgb_data <- load_rgb_values(directory)
saveRDS(rgb_data, "data/rgb_data.RDS")
hsl_data <- load_hsl_values(directory)
saveRDS(hsl_data, "data/hsl_data.RDS")

phylo_colors <- drop.tip(phylo_resolved, phylo_resolved$tip.label[which((!phylo_resolved$tip.label %in% rgb_data$species))])

rgb_data <- rgb_data[rgb_data$species %in% phylo_colors$tip.label, ]
rgb_data_named <- rgb_data
rownames(rgb_data_named) = rgb_data_named[,1]
rgb_data_named <- rgb_data_named[,-1]

hsl_data <- hsl_data[hsl_data$species %in% phylo_colors$tip.label, ]
hsl_data_named <- hsl_data
rownames(hsl_data_named) = hsl_data_named[,1]
hsl_data_named <- hsl_data_named[,-1]




## PHYLOGENETIC SIGNAL

K_vals = sapply(1:ncol(hsl_data_named), function(i) {
  data = hsl_data_named[,i]
  names(data) = rownames(hsl_data_named)
  phylosig(phylo_colors, data)
})

pagel_lambda_multivariate <- function(traits, phy) {
  lambda_vals <- sapply(1:ncol(traits), function(i) {
    trait <- traits[, i, drop = FALSE]
    fit <- fitContinuous(phy, trait, model = "lambda")
    return(fit$opt$lambda)
  })
  return(lambda_vals)
}

color_lambda <- pagel_lambda_multivariate(hsl_data_named, phylo_colors)
overlap_lambda <- pagel_lambda_multivariate()



## MANTEL TEST

color_dist <- as.matrix(dist(hsl_data_named, method = "euclidean"))
gen_dist_a <- cophenetic(phylo_colors)
color_dist <- color_dist[rownames(gen_dist_a), colnames(gen_dist_a)]
mantel_result <- mantel(color_dist, gen_dist_a)


### GET OCCS


cell_overlap_mat = readRDS("./data/overlap_matrix.RDS")


### MANTEL TEST AGAIN

# Assuming gen_dist_a and color_dist are already defined
# Identify the common species
common_species <- intersect(intersect(rownames(gen_dist_a), rownames(color_dist)), rownames(cell_overlap_mat))

# Filter the matrices to include only the common species
gen_dist_a <- gen_dist_a[common_species, common_species]
color_dist <- color_dist[common_species, common_species]
cell_overlap_mat <- cell_overlap_mat[common_species, common_species]

partial_mantel_result <- mantel.partial(cell_overlap_mat, color_dist, gen_dist_a)

mr = mantel(cell_overlap_mat, gen_dist_a)

## PGLMM

# Prepare a data frame with pairwise distances
species_pairs <- expand.grid(species1 = rownames(color_dist), species2 = rownames(color_dist))
species_pairs <- species_pairs %>% 
  filter(species1 != species2) %>% 
  mutate(color_distance = mapply(function(x, y) color_dist[x, y], species1, species2),
         phylo_distance = mapply(function(x, y) gen_dist_a[x, y], species1, species2),
         overlap = mapply(function(x, y) cell_overlap_mat[x, y], species1, species2))

species_pairs_test <- species_pairs[sample(nrow(species_pairs), 10000),]

pglmm_model_function <- function() {
  pglmm(
    color_distance ~ overlap + (1|species1__) + (1|species2__) + (1|species1__@species2) + (1|species2__@species1),
    data = species_pairs,
    cov_ranef = list(species1 = phylo, species2 = phylo),
    family = "gaussian"
  )
}

# Timing the PGLMM function
time_taken <- system.time(
  pglmm_model <- pglmm_model_function()
)

# Summarize the model
summary(pglmm_model)

# Fit a linear model
lm_model <- lm(overlap ~ color_distance + phylo_distance, data = species_pairs)
summary_lm <- summary(lm_model)
summary_lm

# Extract t-value and p-value
t_value <- as.numeric(summary_lm$coefficients[2, "t value"])
p_value <- as.numeric(summary_lm$coefficients[2, "Pr(>|t|)"])

species_pairs$color_similarity = max(species_pairs$color_distance) - species_pairs$color_distance

ggplot(species_pairs, aes(x = color_similarity, y = overlap, color = phylo_distance)) +
  geom_point(size = 3, alpha=0.2) +  # Increase point size for better visibility
  scale_color_gradientn(colors = c("red", "blue")) +
  labs(title = "Color similarity vs. range overlap across all species",
       x = "Color similarity",
       y = "Overlap coefficient",
       color = "Genetic distance") +
  theme_minimal(base_family = "Atkinson Hyperlegible") +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Atkinson Hyperlegible", size = 20),  # Center the title and set the font and size
    text = element_text(family = "Atkinson Hyperlegible"),  # Set the font for the entire plot
    axis.title = element_text(size = 16),  # Increase axis title size
    axis.text = element_text(size = 14),  # Increase axis text size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 12),  # Increase legend text size
    # plot.margin = margin(10, 50, 10, 10)  # Increase right margin to accommodate the annotation
  )


m = mantel.partial(cell_overlap_mat, color_dist, gen_dist_a)







common_species <- intersect(rownames(gen_dist), rownames(cell_overlap_mat))
gen_dist_overlap <- gen_dist[common_species, common_species]
cell_overlap_gen <- cell_overlap_mat[common_species, common_species]

# Function to convert distance matrices to pairs with distances
matrix_to_pairs <- function(matrix1, matrix2) {
  # Ensure row and column names are identical
  if (!all(rownames(matrix1) == rownames(matrix2)) || !all(colnames(matrix1) == colnames(matrix2))) {
    stop("The row and column names of the two matrices must be identical and in the same order.")
  }
  
  # Get the upper triangle indices
  upper_triangle_indices <- which(upper.tri(matrix1), arr.ind = TRUE)
  
  # Extract the row and column names
  row_names <- rownames(matrix1)[upper_triangle_indices[, 1]]
  col_names <- colnames(matrix1)[upper_triangle_indices[, 2]]
  
  # Create a data frame with the pairs and their distances
  pairs <- data.frame(
    row = row_names,
    col = col_names,
    distance1 = matrix1[upper_triangle_indices],
    distance2 = matrix2[upper_triangle_indices]
  )
  
  return(pairs)
}

# Convert the distance matrices to pairs
pairs <- matrix_to_pairs(gen_dist_overlap, cell_overlap_gen)
colnames(pairs) <- c("sp1", "sp2", "gen_dist", "overlap")

plot(data=pairs, overlap~gen_dist) + abline(lm(data=pairs, overlap~gen_dist), col="red")
summary(lm(data=pairs, overlap~gen_dist))

library(vegan)
m = mantel(cell_overlap_gen, gen_dist_overlap, na.rm=T)
print(m)















# Ensuring colors are correctly matched with the tree tip labels
plot_tree_with_colors <- function(tree, rgb_data) {
  # Create a color vector
  colors <- setNames(rgb(rgb_data$R/255, rgb_data$G/255, rgb_data$B/255, maxColorValue = 1), rgb_data$species)
  
  # Ensure colors are applied to the correct species
  tip_colors <- sapply(tree$tip.label, function(x) colors[x])
  
  # Setting up the plot background color and plotting the tree
  par(bg = "#5E5E5E")
  plot(tree, tip.color = tip_colors, show.node.label = FALSE, cex=1, type="fan")
}

png(filename = "phylogenetic_tree.png", width = 36, height = 36, units = "in", res = 300)
plot_tree_with_colors(phylo_colors, rgb_data)
dev.off()

plot_tree_with_colors(phylo_colors, rgb_data)

print(setdiff(rgb_data$species, phylo_colors$tip.label))
print(setdiff(phylo_colors$tip.label, rgb_data$species))



# SSC is generated from convex hulls
# Overlap is generated from cells
# Point overlap is generated from buffered points


pt_overlap_matrix <- as.matrix(readRDS("./data/pt_overlap_matrix.RDS"))
hull_overlap_matrix <- readRDS("./data/hull_overlap_matrix.RDS")
cell_overlap_mat = readRDS("./data/overlap_matrix.RDS")

common_species <- intersect(intersect(rownames(hull_overlap_matrix), 
                                      rownames(cell_overlap_mat)), 
                            rownames(pt_overlap_matrix))


hull_overlap_min <- as.matrix(hull_overlap_matrix[common_species, common_species])
cell_overlap_min <- as.matrix(cell_overlap_mat[common_species, common_species])
pt_overlap_min <- as.matrix(pt_overlap_matrix[common_species, common_species])

hull_min <- hull_overlap_min[lower.tri(hull_overlap_min)]
cell_min <- cell_overlap_min[lower.tri(cell_overlap_min)]
pt_min <- pt_overlap_min[lower.tri(pt_overlap_min)]

plot(hull_min ~ cell_min, na.rm=T, xlab="Cell overlap coefficient", ylab="Convex hull overlap coefficient") + abline(lm(hull_min ~ cell_min), col="red")
plot(hull_min ~ pt_min, na.rm=T, xlab="Buffered point overlap coefficient", ylab="Convex hull overlap coefficient") + abline(lm(hull_min ~ pt_min), col="red")
plot(cell_min ~ pt_min, na.rm=T, xlab="Buffered point overlap coefficient", ylab="Cell overlap coefficient") + abline(lm(cell_min ~ pt_min), col="red")
summary(lm(cell_min ~ pt_min))

library(vegan)
m = mantel(cell_overlap_min, pt_overlap_min, na.rm=T)
print(m)









