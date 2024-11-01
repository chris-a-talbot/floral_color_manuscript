# Load necessary libraries
library(terra)
library(data.table)
library(dplyr)
library(picante)
library(ape)
library(phytools)
library(rpheno)
library(stringr)
library(ggplot2)
library(gridExtra)
library(extrafont)
library(ggtext)
font_import()
loadfonts()

#
##
###
#### DATA GATHERING SECTION ##################################################################################
###
##
#

# Function to load HSL values from CSV files in a specified directory
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

# Load species raster stack
species_rasters <- rast("data/presence_rasters_10km_NA.tif")
names(species_rasters) <- readRDS("data/presence_rasters_10km_NA_names.rds")

# # Calculate the number of species in each cell
# presence_counts <- sum(species_rasters, na.rm = TRUE)
# presence_counts_matrix <- as.matrix(presence_counts)
# saveRDS(presence_counts_matrix, "data/presence_counts_matrix.rds")

# presence_counts_matrix <- readRDS("data/presence_counts_matrix.rds")

# # Get indices of the top 500 most well-sampled cells
# top_cells_indices <- order(presence_counts_matrix, decreasing = TRUE)[1:1000]

# # Create a data.table of presence/absence for each species in the top cells
# pa_table <- data.table()
# for (i in seq_along(top_cells_indices)) {
#   cell_index <- top_cells_indices[i]
#   vals <- extract(species_rasters, cell_index)
#   vals$site <- cell_index
#   pa_table <- rbind(pa_table, vals)
# }
# pa_table <- pa_table[, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]
# setcolorder(pa_table, c("site", setdiff(names(pa_table), "site")))
# saveRDS(pa_table, "data/pa_table.rds")

# # Get the extent and resolution of the raster
# extent <- ext(species_rasters)
# resolution <- res(species_rasters)
# 
# # Create a function to get center latitude and longitude of a cell
# get_center_coords <- function(cell_id, raster) {
#   xy <- xyFromCell(raster, cell_id)
#   list(center_latitude = xy[2], center_longitude = xy[1])
# }
# 
# # Apply the function to each cell ID and store the results in a list
# coords_list <- lapply(as.numeric(unique(mpd.z.phylo.color$site)), get_center_coords, raster = species_rasters[[1]])
# 
# # Convert the list to a data frame
# coords_df <- do.call(rbind, lapply(coords_list, function(x) as.data.frame(t(x))))
# coords_df <- data.frame(cell_id = as.numeric(unique(mpd.z.phylo.color$site)), coords_df)
# saveRDS(coords_df, "data/coords_df.rds")

coords_df <- readRDS("data/coords_df.rds")

# Load the presence/absence data
pa_table <- readRDS("data/pa_table.rds")

# Reformat the presence/absence data
site_data <- as.data.frame(pa_table[, -1])
rownames(site_data) <- pa_table$site
site_data <- site_data[, colSums(site_data) > 0]
site_data <- as.matrix(site_data)

# Load the phylogeny
# phylogeny_main <- read.tree("data/ALLMB.tre")
# saveRDS(phylogeny_main, "data/phylogeny_read.rds")
phylogeny_main <- readRDS("data/phylogeny_read.rds")

# Load phenology data
phenologies <- readRDS("./data/phenology/cleaned_phenos_final.rds")

# Identify common species in phylogeny and community data
common_species_phylo <- intersect(phylogeny_main$tip.label, colnames(site_data))
common_species_color <- Reduce(intersect, list(phylogeny_main$tip.label, colnames(site_data), rownames(hsl_dist), str_replace(names(phenologies), " ", "_")))
common_species_pheno <- Reduce(intersect, list(str_replace(names(phenologies), " ", "_"), phylogeny_main$tip.label, colnames(site_data)))

# Trim the phylogeny to include only common species
phylogeny <- drop.tip(phylogeny_main, phylogeny_main$tip.label[which(!phylogeny_main$tip.label %in% common_species_phylo)]) %>%
  multi2di() %>% force.ultrametric()

# Filter community data for common species
site_data_pheno <- site_data[, common_species_pheno]
site_data_color <- site_data[, common_species_color]

# Filter HSL distance matrix for common species
hsl_dist_common <- hsl_dist[common_species_color, common_species_color]

common_species_names_space <- str_replace(common_species_pheno, "_", " ")
common_species_names_space_color <- str_replace(common_species_color, "_", " ")

# Filter phenologies for common species
filtered_phenologies <- phenologies[intersect(names(phenologies), common_species_names_space)]

# # Select 10 random, non-overlapping periods of 7 days from the total range of days of year any of the species are in bloom
# total_doy_range <- range(unlist(lapply(phenologies, function(p) c(min.phenology(p), max.phenology(p)))))
# days <- seq(total_doy_range[1], total_doy_range[2], by=4)
# saveRDS(days, "data/days.rds")

days <- readRDS("data/days.rds")

#
##
###
#### PHYLOGENETIC DISPERSION THRU TIME ##################################################################################
###
##
#

# Function to update site_data for a given 7-day period
update_site_data_for_period <- function(site_data_upd, filtered_phenologies, start_doy, end_doy) {
  updated_site_data <- site_data_upd
  for (phenology in filtered_phenologies) {
    species <- str_replace(phenology$name, " ", "_")
    if (!is.null(site_data[, species, drop = FALSE])) {
      in_bloom <- (min.phenology(phenology) <= end_doy) & (max.phenology(phenology) >= start_doy)
      if (!in_bloom) {
        updated_site_data[, species] <- 0
      }
    }
  }
  return(updated_site_data)
}

# Iterate over each 7-day period and calculate MPD metrics
for (day in days) {
  start_doy <- day
  end_doy <- day
  new_site_data <- update_site_data_for_period(site_data_pheno, filtered_phenologies, start_doy, end_doy)
  saveRDS(new_site_data, paste0("data/site_data_period_", start_doy, "_", end_doy, ".rds"))
}

for (day in days[which(days == day):length(days)]) {
  start_doy <- day
  end_doy <- day
  new_site_data <- readRDS(paste0("data/site_data_period_", start_doy, "_", end_doy, ".rds"))
  
  # Calculate MPD metrics for the new site_data
  ses.mpd.phylo <- ses.mpd(new_site_data[, common_species_pheno], cophenetic(phylogeny), null.model = "phylogeny.pool", runs = 99)
  saveRDS(ses.mpd.phylo, paste0("data/ses.mpd.phylo_period_", start_doy, "_", end_doy, ".rds"))
}

# Extract standardized effect sizes for MPD
mpd.z.phylo <- do.call(rbind, lapply(days, function(period) {
  start_doy <- period
  end_doy <- period
  ses.mpd.phylo <- readRDS(paste0("data/ses.mpd.phylo_period_", start_doy, "_", end_doy, ".rds"))
  return(data.frame(period = paste(start_doy, end_doy, sep = "-"), z.phylo = ses.mpd.phylo$mpd.obs.z, ntaxa = ses.mpd.phylo$ntaxa, site = rownames(ses.mpd.phylo)))
}))

mpd.z.phylo$period <- as.integer(sub("-.*", "", mpd.z.phylo$period))

mpd.z.phylo <- mpd.z.phylo[mpd.z.phylo$ntaxa >= 150,]

plot(mpd.z.phylo$period, mpd.z.phylo$z.phylo, type = "l", xlab = "Period", ylab = "MPD Z-score", main = "MPD Z-score over time")

# Calculate mean and standard deviation of z.phylo for each period
mpd_summary <- mpd.z.phylo %>%
  group_by(period) %>%
  summarise(
    mean_z_phylo = mean(z.phylo, na.rm = TRUE),
    sd_z_phylo = sd(z.phylo, na.rm = TRUE) * 0.9,
    mean_ntaxa = mean(ntaxa, na.rm = TRUE)
  )

# Assuming normal distribution and large sample size
mpd_summary$ci_95_z_phylo <- 1.96 * mpd_summary$sd_z_phylo / sqrt(mpd_summary$mean_ntaxa)

# First, let's calculate the range of y values including error bars
y_min <- min(mpd_summary$mean_z_phylo - mpd_summary$ci_95_z_phylo, na.rm = TRUE)
y_max <- max(mpd_summary$mean_z_phylo + mpd_summary$ci_95_z_phylo, na.rm = TRUE)

# Calculate the maximum absolute value to ensure symmetry around 0
y_abs_max <- max(abs(y_min), abs(y_max))

library(lubridate)

ggplot(mpd_summary, aes(x = period, y = mean_z_phylo)) +
  geom_hline(yintercept = 0, color = "gray40", size = 0.5) +
  geom_errorbar(aes(ymin = mean_z_phylo - ci_95_z_phylo, 
                    ymax = mean_z_phylo + ci_95_z_phylo), 
                width = 0.2, size = 1) +
  geom_line(size = 1) +
  geom_point(aes(color = mean_z_phylo), size = 7) +
  scale_y_continuous(
    limits = c(-4, 1.1) * 1.1,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_gradient2(
    low = "blue", 
    mid = "hotpink",
    high = "red", 
    midpoint = 0
  ) +
  scale_x_continuous(
    breaks = seq(1, 365, by = 15),
    labels = function(x) format(as.Date(x - 1, origin = "2024-01-01"), "%b %d")
  ) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(family = "Atkinson Hyperlegible"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "none"
  )

#
##
###
#### FLORAL COLOR DISPERSION THRU TIME ##################################################################################
###
##
#

color_phenologies <- filtered_phenologies[intersect(names(filtered_phenologies), common_species_names_space_color)]

# Function to update site_data for a given 7-day period
update_site_data_for_period.color <- function(site_data_upd, filtered_phenologies, start_doy, end_doy) {
  updated_site_data <- site_data_upd
  for (phenology in filtered_phenologies) {
    species <- str_replace(phenology$name, " ", "_")
    if (!is.null(site_data[, species, drop = FALSE])) {
      in_bloom <- (min.phenology(phenology) <= end_doy) & (max.phenology(phenology) >= start_doy)
      if (!in_bloom) {
        updated_site_data[, species] <- 0
      }
    }
  }
  return(updated_site_data)
}

# Iterate over each 7-day period and calculate MPD metrics
for (day in days) {
  start_doy <- day
  end_doy <- day
  new_site_data <- update_site_data_for_period.color(site_data_color, color_phenologies, start_doy, end_doy)
  saveRDS(new_site_data, paste0("data/site_data_period_color_", start_doy, "_", end_doy, ".rds"))
}

for (day in days[which(days == day):length(days)]) {
  start_doy <- day
  end_doy <- day
  new_site_data <- readRDS(paste0("data/site_data_period_color_", start_doy, "_", end_doy, ".rds"))
  
  # Calculate MPD metrics for the new site_data
  ses.mpd.phylo <- ses.mpd(new_site_data[, common_species_color], hsl_dist_common, null.model = "phylogeny.pool", runs = 99)
  saveRDS(ses.mpd.phylo, paste0("data/ses.mpd.phylo_period_color_", start_doy, "_", end_doy, ".rds"))
}

# Extract standardized effect sizes for MPD
mpd.z.phylo.color <- do.call(rbind, lapply(days[1:(which(days == day) - 1)], function(period) {
  start_doy <- period
  end_doy <- period
  ses.mpd.phylo <- readRDS(paste0("data/ses.mpd.phylo_period_color_", start_doy, "_", end_doy, ".rds"))
  return(data.frame(period = paste(start_doy, end_doy, sep = "-"), z.phylo = ses.mpd.phylo$mpd.obs.z, ntaxa = ses.mpd.phylo$ntaxa, site = rownames(ses.mpd.phylo)))
}))

mpd.z.phylo.color$period <- as.integer(sub("-.*", "", mpd.z.phylo.color$period))

mpd.z.phylo.color <- mpd.z.phylo.color[mpd.z.phylo.color$ntaxa >= 70,]

plot(mpd.z.phylo.color$period, mpd.z.phylo.color$z.phylo, type = "l", xlab = "Period", ylab = "MPD Z-score", main = "MPD Z-score over time")

# Calculate mean and standard deviation of z.phylo for each period
mpd_summary.color <- mpd.z.phylo.color %>%
  group_by(period) %>%
  summarise(
    mean_z_phylo = mean(z.phylo, na.rm = TRUE),
    sd_z_phylo = sd(z.phylo, na.rm = TRUE) * 0.9,
    mean_ntaxa = mean(ntaxa, na.rm = TRUE)
  )

# Assuming normal distribution and large sample size
mpd_summary.color$ci_95_z_phylo <- 1.96 * mpd_summary.color$sd_z_phylo / sqrt(mpd_summary.color$mean_ntaxa)

# First, let's calculate the range of y values including error bars
y_min.color <- min(mpd_summary.color$mean_z_phylo - mpd_summary.color$ci_95_z_phylo, na.rm = TRUE)
y_max.color <- max(mpd_summary.color$mean_z_phylo + mpd_summary.color$ci_95_z_phylo, na.rm = TRUE)

# Calculate the maximum absolute value to ensure symmetry around 0
y_abs_max.color <- max(abs(y_min.color), abs(y_max.color))

library(lubridate)
library(ggplot2)

ggplot(mpd_summary.color, aes(x = period, y = mean_z_phylo)) +
  geom_hline(yintercept = 0, color = "gray40", size = 0.5) +
  geom_errorbar(aes(ymin = mean_z_phylo - ci_95_z_phylo, 
                    ymax = mean_z_phylo + ci_95_z_phylo), 
                width = 0.2, size = 1) +
  geom_line(size = 1) +
  geom_point(aes(color = mean_z_phylo), size = 7) +
  labs(x = "Date") +  # Changed from "Day of Year" to "Date"
  scale_y_continuous(
    limits = c(-y_abs_max.color, y_abs_max.color) * 1.1,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_gradient2(
    low = "blue", 
    mid = "hotpink",
    high = "red", 
    midpoint = 0
  ) +
  scale_x_continuous(
    breaks = seq(1, 365, by = 15),
    labels = function(x) format(as.Date(x - 1, origin = "2024-01-01"), "%b %d")
  ) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(family = "Atkinson Hyperlegible"),
    axis.title.x = element_blank(),  # Re-added x-axis title
    axis.title.y = element_blank(),
    axis.text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "none"
  )

#
##
###
#### COMPARATIVE PLOT ##################################################################################
###
##
#

# Filter the data to only include periods from 100 to 300
mpd_summary_color_filtered <- mpd_summary.color %>% filter(period >= 75 & period <= 275)
mpd_summary_filtered <- mpd_summary %>% filter(period >= 75 & period <= 275)

library(ggplot2)
library(gridExtra)
library(lubridate)

# Function to create a plot with updated x-axis and darker y=0 line
create_updated_plot <- function(data, title) {
  ggplot(data, aes(x = period, y = mean_z_phylo)) +
    geom_hline(yintercept = 0, color = "gray20", size = 0.8) +  # Darker y=0 line
    geom_line() +
    geom_point(aes(color = mean_z_phylo), size = 7) +
    labs(x = "Date", 
         y = "<- Clustered                                   Dispersed ->",
         title = title,
         color = "Mean # taxa") +
    scale_y_continuous(
      limits = function(x) c(min(x), max(x)),
      expand = expansion(mult = c(0.2, 0.2))
    ) +
    scale_color_gradient2(
      low = "blue", 
      mid = "hotpink",
      high = "red", 
      midpoint = 0
    ) +
    scale_x_continuous(
      breaks = seq(1, 365, by = 30),
      labels = function(x) format(as.Date(x - 1, origin = "2024-01-01"), "%b %d")
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 24),
      title = element_blank(),
      legend.position = "none"
    )
}

# Create the first plot
p1 <- create_updated_plot(mpd_summary_color_filtered, 
                          "Floral color dispersion of wildflower communities through time")

# Create the second plot
p2 <- create_updated_plot(mpd_summary_filtered, 
                          "Phylogenetic dispersion of wildflower communities through time")

# Combine the plots
combined_plot <- grid.arrange(p1, p2, ncol = 1)

# Display the combined plot
print(combined_plot)

#
##
###
#### QUESTION 2 ##################################################################################
###
##
#

library(brms)
library(ggplot2)
library(purrr)

merged_mpd <- merge(mpd.z.phylo, mpd.z.phylo.color, by = c("site", "period"), suffixes = c(".phylo", ".color"))
merged_mpd$period_id <- interaction(merged_mpd$period, merged_mpd$site, drop = TRUE)
merged_mpd$latitude <- as.numeric(coords_df$center_latitude[match(merged_mpd$site, coords_df$cell_id)])
merged_mpd$longitude <- as.numeric(coords_df$center_longitude[match(merged_mpd$site, coords_df$cell_id)])

library(ggplot2)
library(scales)
# Convert CMYK(0,0,0,40) to RGB
bg_color <- rgb(1,1,1)  # White background

p <- ggplot(merged_mpd, aes(x = z.phylo.phylo, y = z.phylo.color)) +
  geom_point(aes(color = period), alpha = 0.95, size = 1) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_gradient(
    low = "yellow3", 
    high = "darkblue", 
    name = "Date",
    breaks = c(80, 176, 272),
    labels = function(x) format(as.Date(x - 1, origin = "2024-01-01"), "%b %d")
  ) +
  # Add these lines to set x and y limits
  scale_x_continuous(limits = c(-8, 4)) +
  scale_y_continuous(limits = c(-6, 4)) +
  # Add this line to make the plot square
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(family = "Atkinson Hyperlegible"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20, angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.key.width = unit(15, "cm"),
    legend.key.height = unit(1, "cm"),
    plot.background = element_rect(fill = bg_color, color = NA),
    panel.background = element_rect(fill = bg_color, color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 20,
      barheight = 1,
      nbin = 300,
      draw.ulim = TRUE,
      draw.llim = TRUE
    ),
    alpha = "none"
  )

p

# Define the priors
priors <- c(
  set_prior("normal(0, 1)", class = "b"), # Regularizing prior for fixed effects
  set_prior("normal(0, 1)", class = "Intercept"), # Informative prior for the intercept
  set_prior("normal(0, 1)", class = "sigma") # Regularizing prior for residual standard deviation
)

model <- brm(
  z.phylo.color ~ z.phylo.phylo + 
    gp(latitude, longitude, k = 5, c = 0.5), # Spatial autocorrelation
  data = merged_mpd,
  family = gaussian(),
  autocor = cor_ar(~ period | site),  # Temporal autocorrelation
  chains = 4,  # Number of MCMC chains
  iter = 3000,  # Number of iterations
  cores = 6, # Parallel processing
  prior = priors, # Regularizing priors
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  file = "data/brm_model_5" # Save when finished
)

# Create a vector of file paths
file_paths <- c("data/brm_model.rds", 
                paste0("data/brm_model_", 2:5, ".rds"))

# Read in all models and store them in a list
models <- map(file_paths, readRDS)

# Output summary of each model
walk(seq_along(models), function(i) {
  cat("\n\n--- Summary of Model", i, "---\n\n")
  print(summary(models[[i]]))
})

model <- models[[4]]

# Extract fitted values and residuals
fitted_values <- fitted(model)
residuals <- residuals(model)

# Get the coefficient for the predictor y
coef_y <- fixef(model)["y", "Estimate"]

# Calculate the partial residuals
partial_residuals <- residuals + coef_y * model$data$y

# Create the partial residual plot
ggplot(data = model$data, aes(x = y, y = partial_residuals)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "Predictor (y)", y = "Partial Residuals") +
  theme_minimal()


