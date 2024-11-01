### STEP 1: Get all necessary data
### Packages, shapefiles, color traits, phenology, and phylogeny

# setwd("G:/My Drive/Documents/Research/Weber/honors-thesis")

# Pre-made R functions for handling rasters & presence/absence matrices
source("helper_raster.R")
source("helper_ecdf.R")

# Load pacman package for loading the other packages
if (!require("pacman")) {install.packages("pacman")}
library("pacman")

# Load required packages
p_load("data.table", "stringr", "dplyr", 
       "parallel", "speciesRaster", "terra")

core_num = detectCores() - 1

# Get the flowers masterfile
flowers = fread("./data/flowers_master.csv")

# Get a list of species names
species_list = word(flowers$taxon, 1, 2)

# Get the phenological trait data
traits_pheno = cbind(species_list, flowers[,16:27])
colnames(traits_pheno) = c("species", colnames(traits_pheno)[-1])

# Get the color trait data
traits_color = as.data.frame(flowers[,7:15])
rownames(traits_color) = species_list

# # Get the phylogenetic data
# tree_path = "./data/flowers_tree.txt"
# tree = read.newick(tree_path)
# tree$tip.label = str_replace(word(tree$tip.label, 1, 2, sep="_"), "_", " ")
# gen_dist_mat = cophenetic.phylo(tree)
# species_names = tree$tip.label

# Get the presence/absence data
months = c("feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov")
pa_matrices_by_month = get_month_pa_matrix_from_tif(months, res=0.09)
pa_dts_by_month = mclapply(pa_matrices_by_month, function(matrix) {
  return(as.data.table(matrix))
}, mc.cores=1)
cell_num = nrow(pa_matrices_by_month[[1]])

better_flowers = fread("./data/better_flowers.csv")

# # Adjusting the approach to safely filter columns
# filtered_pa_dts_by_month <- lapply(pa_dts_by_month, function(dt) {
#   # Find the intersection of column names in dt with species_list_present
#   valid_cols <- intersect(names(dt), species_list_present)
#   # Subset dt to keep only columns that exist in both dt and species_list_present
#   dt[, ..valid_cols]
# })

sums_dt = mclapply(pa_dts_by_month, function(dt) rowSums(dt), 
                   mc.cores=1
                   )
sums <- Reduce("+", lapply(sums_dt, function(x) replace(x, is.na(x), 0)))
names(sums) <- 1:length(sums)
sums[sums < 100] <- NA

# site_qualities = mclapply(1:length(sums), function(index) {
#   if(is.na(sums[index])) { return(NA) }
#   row = pa_matrices_by_month[[4]][index,]
#   species_site = colnames(row)[which(row == 1, arr.ind = TRUE)[, "col"]]
#   nob = sum(better_flowers[word(taxon, 1, 2) %in% species_site]$ob_num)
#   print(index)
#   return(nob*sums[index])
# }, mc.cores=1)
# 
# save.image("./RData/sitepicker.RData")
# 
# site_qualities_list = unlist(site_qualities) 
# site_qualities_list[!(site_qualities_list %in% site_qualities_list[order(site_qualities_list, decreasing=TRUE)][1:2500])] = NA

# # Step 1: Identify non-NA elements
# non_na_indices <- which(!is.na(site_qualities_list))
# 
# # Step 2: Randomly select 100 non-NA elements, if there are more than 100 non-NA elements
# if (length(non_na_indices) > 100) {
#   selected_indices <- sample(non_na_indices, 100)
# } else {
#   selected_indices <- non_na_indices
# }

# # Step 3: Create a new list that includes only the selected non-NA and all NA elements
# best_sites <- rep(NA, length(site_qualities_list)) # Start with all NAs
# best_sites[selected_indices] <- site_qualities_list[selected_indices] # Replace selected indices with their original values
# 
# # Convert NA entries to actual NA values, not as character "NA" if my_list is character type
# best_sites <- lapply(best_sites, function(x) ifelse(is.na(x) || x == "NA", NA, x))
# 
# best_indices <- which(!is.na(best_sites))

set.seed(1025)
non_na_indices <- which(!is.na(sums))
if (length(non_na_indices) > 100) {
  selected_indices <- sample(non_na_indices, 100)
} else {
  selected_indices <- non_na_indices
}
best_sites <- rep(NA, length(sums))
best_sites[selected_indices] <- sums[selected_indices]
best_sites <- lapply(best_sites, function(x) ifelse(is.na(x) || x=="NA", NA, x))
best_indices <- which(!is.na(best_sites))

saveRDS(best_indices, "./data/indices.rds")