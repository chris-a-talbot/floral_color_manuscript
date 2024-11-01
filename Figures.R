library(ggplot2)
library(patchwork)
library(extrafont)

# Load the Atkinson Hyperlegible font (make sure it's installed on your system)
# You may need to run font_import() first if you haven't already
loadfonts(device = "win")

# Create dummy data for the plots
set.seed(123)
n <- 100
x <- seq(0, 10, length.out = n)
y1 <- 2 + 0.5 * x + rnorm(n, sd = 0.5)  # Positive relationship
y2 <- rep(0, n)  # Flat relationship at y = 0
y3 <- 8 - 0.5 * x + rnorm(n, sd = 0.5)  # Negative relationship

# Function to create a square plot
create_square_plot <- function(x, y, title, show_x_label = FALSE, show_y_label = FALSE) {
  p <- ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = title) +
    theme_minimal(base_family = "Atkinson Hyperlegible") +
    theme(
      text = element_text(size = 34),
      plot.title = element_text(size = 34, hjust = 0.5),
      aspect.ratio = 1,  # This makes the plot square
      axis.text = element_blank(),  # Remove tick labels
      axis.ticks = element_blank(),  # Remove tick marks
      panel.grid = element_blank(),  # Remove all grid lines
      axis.title = element_blank()  # Remove axis titles by default
    ) +
    coord_fixed(ratio = 1)  # This ensures the aspect ratio is maintained
  
  if (show_x_label) {
    p <- p + labs(x = "Phylogenetic diversity") +
      theme(axis.title.x = element_text(size = 26, margin = margin(t = 10)))
  }
  if (show_y_label) {
    p <- p + labs(y = "Floral color diversity") +
      theme(axis.title.y = element_text(size = 26, margin = margin(r = 10), angle = 90))
  }
  
  return(p)
}

# Create the three plots with correct relationships and labels
plot1 <- create_square_plot(x, y1, "Niche conservatism", show_y_label = TRUE)
plot2 <- create_square_plot(x, y2, "Neither", show_x_label = TRUE)
plot3 <- create_square_plot(x, y3, "Character displacement")

# Combine the plots side by side with vertical lines
combined_plot <- plot1 + plot_spacer() + plot2 + plot_spacer() + plot3 +
  plot_layout(
    ncol = 5, 
    widths = c(1, 0.02, 1, 0.02, 1)
  ) &
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

# Save the combined plot
ggsave("combined_plot.png", combined_plot, width = 15, height = 5, dpi = 300)











library(ggplot2)
library(patchwork)
library(dplyr)

# Determine the common x-axis range
x_min <- max(min(mpd_summary.color$period), min(mpd_summary$period))
x_max <- min(max(mpd_summary.color$period), max(mpd_summary$period))

# Filter the data to include only the common range
mpd_summary.color_filtered <- mpd_summary.color %>%
  filter(period >= x_min & period <= x_max)

mpd_summary_filtered <- mpd_summary %>%
  filter(period >= x_min & period <= x_max)

# Create the first plot (for mpd_summary.color)
plot1 <- ggplot(mpd_summary.color_filtered, aes(x = period, y = mean_z_phylo)) +
  geom_hline(yintercept = 0, color = "gray40", size = 0.5) +
  geom_errorbar(aes(ymin = mean_z_phylo - ci_95_z_phylo, 
                    ymax = mean_z_phylo + ci_95_z_phylo), 
                width = 0.2, size = 1) +
  geom_line(size = 1) +
  geom_point(aes(color = mean_z_phylo), size = 7) +
  scale_y_continuous(
    limits = c(-y_abs_max.color, y_abs_max.color) * 1.1,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_gradient2(
    low = "blue", 
    mid = "hotpink",
    high = "red", 
    midpoint = -0.25
  ) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(family = "Atkinson Hyperlegible"),
    axis.text = element_text(size = 26),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "none",
    axis.title = element_blank()
  )

# Create the second plot (for mpd_summary)
plot2 <- ggplot(mpd_summary_filtered, aes(x = period, y = mean_z_phylo)) +
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
    midpoint = -2
  ) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(family = "Atkinson Hyperlegible"),
    axis.text = element_text(size = 26),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "none",
    axis.title = element_blank()
  )

# Combine the plots using patchwork
combined_plot <- plot1 + plot2 +
  plot_layout(ncol = 1, heights = c(1, 1))

# Print the combined plot
print(combined_plot)



# Load required library
library(ape)

# Set seed for reproducibility
set.seed(123)

# Create a random ultrametric tree with 5 species
tree <- rcoal(5)

# Rename the tip labels
tree$tip.label <- paste("Species", 1:5)

# Function to plot the tree with specified species in bold
plot_tree_with_bold <- function(tree, bold_species, title) {
  plot(tree, show.tip.label = FALSE, main = title, edge.width = 2, 
       x.lim = c(0, 1.2))
  
  # Add labels for non-bold species
  non_bold_indices <- which(!(tree$tip.label %in% bold_species))
  tiplabels(tree$tip.label[non_bold_indices], adj = c(0, 0.5), 
            frame = "none", cex = 2, offset = 0.1, 
            tip = non_bold_indices)
  
  # Add bold labels for specified species
  bold_indices <- which(tree$tip.label %in% bold_species)
  tiplabels(tree$tip.label[bold_indices], adj = c(0, 0.5), 
            frame = "none", font = 2, cex = 2, offset = 0.1, 
            tip = bold_indices)
}

# Set up the plotting area for two side-by-side plots
par(mfrow = c(1, 2), mar = c(1, 1, 2, 5))

# Plot the first phylogeny with two sister species in bold
plot_tree_with_bold(tree, c("Species 5", "Species 4"), 
                    "Phylogeny 1: Sister Species in Bold")

# Plot the second phylogeny with the two most distantly related species in bold
plot_tree_with_bold(tree, c("Species 1", "Species 5"), 
                    "Phylogeny 2: Most Distant Species in Bold")

# Reset the plotting area
par(mfrow = c(1, 1))



# Required libraries
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

# Load map data
north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = "united states of america", returnclass = "sf")
great_lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = "sf")

# Generate random points within the map extent
set.seed(42) # For reproducibility
n_points <- 5000 # Generate more points initially to ensure enough fall on land
random_points <- data.frame(
  lon = runif(n_points, -90, -65),
  lat = runif(n_points, 35, 50)
)

# Convert random points to an sf object
random_points_sf <- st_as_sf(random_points, coords = c("lon", "lat"), crs = st_crs(north_america))

# Filter points to keep only those on land
land_points <- st_intersection(random_points_sf, st_union(north_america))

# Select a sample of 19 points from the filtered land points
set.seed(42) # For reproducibility
selected_land_points <- land_points %>% sample_n(350)

# Add the Pennsylvania point
pennsylvania_point <- data.frame(
  lon = -77.5,
  lat = 41
) %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(north_america))

# Combine the points
all_points <- bind_rows(
  selected_land_points %>% mutate(alpha = 0.25),
  pennsylvania_point %>% mutate(alpha = 1)
)

# Create the map
ggplot() +
  # Add water background
  geom_sf(data = st_bbox(north_america) %>% st_as_sfc() %>% st_set_crs(st_crs(north_america)), 
          fill = "lightblue", color = NA) +
  
  # Add land areas
  geom_sf(data = north_america, fill = "gray50", color = "white") +
  geom_sf(data = states, fill = "gray50", color = "white") +
  
  # Add Great Lakes
  geom_sf(data = great_lakes, fill = "lightblue", color = "lightblue") +
  
  # Add random red dots with varying alpha
  geom_sf(data = all_points, aes(alpha = alpha), color = "red", size = 3) +
  
  # Set the coordinate system and map limits
  coord_sf(xlim = c(-90, -65), ylim = c(35, 50), expand = FALSE) +
  
  # Customize the appearance
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_alpha_identity() # Use identity scale for alpha







library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load your raster
raster <- species_rasters["Aquilegia_canadensis"]

# Get North America map and lakes
na_map <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
lakes <- ne_load(type = "lakes", scale = "medium", category = "physical", returnclass = "sf")

# Convert raster to points and then to data frame
raster_points <- as.points(raster, values = TRUE)
raster_df <- as.data.frame(raster_points, geom = "XY")
raster_df <- raster_df[raster_df$Aquilegia_canadensis == 1, ]

# Create the plot
p <- ggplot() +
  geom_sf(data = na_map, fill = "gray80", color = NA) +
  geom_sf(data = lakes, fill = "lightblue", color = NA) +
  geom_point(data = raster_df, aes(x = x, y = y), color = "red", size = 0.5) +
  coord_sf(xlim = c(min(raster_df$x), max(raster_df$x)),
           ylim = c(min(raster_df$y), max(raster_df$y))) +
  theme_void() +
  theme(panel.background = element_rect(fill = "lightblue"))

# Save the plot as PNG
ggsave("Aquilegia_canadensis_distribution.png", p, width = 10, height = 8, dpi = 300)
