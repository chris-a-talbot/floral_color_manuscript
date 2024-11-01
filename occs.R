library(pacman)
p_load("parallel", "data.table", "dplyr", "sf", "stringr", "terra")

get_occ_trimmed = function(wd="G:/My Drive/Documents/Research/Weber/honors-thesis") {
  occ_trimmed = fread("./occurrences/occ_trimmed.csv")[,-1]
  return(occ_trimmed)  
}

occ_trimmed = get_occ_trimmed() # Get the full list of species occurences for all species

species_list_occ = unique(occ_trimmed$species) # Get the list of species represented in the occurrences

flowers = fread("./data/flowers_master.csv")

valid_species = c(word(flowers$taxon, 1, 2), flowers$species)

occ_trimmed = occ_trimmed[species %in% valid_species]

occ_trimmed[flowers, on = .(species = species), species := sapply(strsplit(i.taxon, " "), function(x) paste(x[1:min(2,length(x))], collapse = " ")), by = .EACHI]

occ_trimmed = occ_trimmed[decimalLatitude >= 38 & decimalLatitude <= 48 & decimalLongitude >= -92.5 & decimalLongitude <= -67.5]

occ_points <- st_as_sf(occ_trimmed, coords = c("decimalLongitude", "decimalLatitude"), crs = 32618, agr = "constant")

# Parallel loop for buffering and writing shapefiles (adapted for parallel execution)
mclapply(unique(occ_points$species), function(species) {
  if(file.exists(paste0("./shp/pts/", species, ".shp"))) { return(NA) }
  species_points <- occ_points[occ_points$species == species, ]
  occ_buffers_01 <- st_buffer(species_points$geometry, dist = 0.09)
  st_write(occ_buffers_01, paste0("./shp/pts/", species, "_0.09.shp"))
  print(paste0(species, " written"))
}, mc.cores=7
)