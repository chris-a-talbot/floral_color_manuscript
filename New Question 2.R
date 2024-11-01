install.packages("brms")
install.packages("dplyr")

library(brms)
library(dplyr)

mpd_summary_phylo <- readRDS("data/mpd_summary.rds")
mpd_summary_color <- readRDS("data/mpd_summary_color.rds")

merged_mpd <- merge(mpd.z.phylo, mpd.z.phylo.color, by = c("site", "period"), suffixes = c(".phylo", ".color"))
merged_mpd$period_id <- interaction(merged_mpd$period, merged_mpd$site, drop = TRUE)
merged_mpd$latitude <- as.numeric(coords_df$center_latitude[match(merged_mpd$site, coords_df$cell_id)])
merged_mpd$longitude <- as.numeric(coords_df$center_longitude[match(merged_mpd$site, coords_df$cell_id)])

model <- brm(
  z.phylo.color ~ z.phylo.phylo + 
    (1 | period) +  # Random intercept and slope for site within period
    (1 | site) + # +  # Random intercept and slope for period within site
    gp(latitude, longitude, k = 5, c = 0.5),
  data = merged_mpd,
  family = gaussian(),
  autocor = cor_ar(~ period | site),  # Temporal autocorrelation
  chains = 4,  # Number of MCMC chains
  iter = 200,  # Number of iterations
  cores = 6,
  control = list(adapt_delta = 0.99)
)

summary(model)

saveRDS(model, "data/brm_model.rds")