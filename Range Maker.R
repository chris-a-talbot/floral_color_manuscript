# # - - - - - - Start Prep work - - - - - - - - # #

get_occ_trimmed = function(wd="G:/My Drive/Documents/Research/Weber/honors-thesis") {
  
  library(data.table) # For data manipulation
  
  # Set working directory with flowers.csv
  setwd(wd)
  
  occ_trimmed = fread("occurrences/occ_trimmed.csv")[,-1]
  return(occ_trimmed)  
}

library(data.table) # For data management
library(rangeBuilder) # For making the ranges
library(dplyr) # For mapping and whatnot
library(sf) # For reading/writing spatial objects

occ_trimmed = get_occ_trimmed() # Get the full list of species occurences for all species

species_list = unique(occ_trimmed$species) # Get the list of species represented in the occurrences

# UPDATE REGULARLY - This is the place i left off on during the previous run
start = 1

# # - - - - - - - End Prep work - - - - - - - - # #
# # # # # # # # # # # # # # # # # # # # # # # # # #
# # - - - - - - Start Range work - - - - - - -  # #

# Run for every species in the list
for(i in start:length(species_list)) {
  
  if(file.exists( paste0("shp/points/", species_name, ".shp"))) next
  
  species_name = species_list[i] # Next species name
  species_occ = occ_trimmed[species == species_name][decimalLatitude >= 17][decimalLongitude <= -50][order(decimalLatitude, decimalLongitude)] %>%
    unique() # Gets occurrences for only the given species, also isolated to North America, then sorts by lat & lon and removes duplicates
  
  species_occ_trimmed = data.table() # Empty data.table for later use
  
  # Indices for the loop
  num_occ = nrow(species_occ) # Number of occurrences for this species
  if(num_occ <= 5) { next } # Skip this species if there's 5 or fewer occurrences within the specifications
  
  # Preliminary work for the while loop
  starting_index = 1
  ending_index = 1000
  
  # If <=2500 occurrences, only run the while() loop contents once
  run_once = FALSE
  if(num_occ <= 2500) { ending_index = num_occ; run_once = TRUE }
  
  # Loops through the occurrences in 1000-data-point chunks
  while(num_occ > 0) {
    
    chunk = species_occ[starting_index:ending_index,] # Get the next chunk
    chunk_trimmed = filterByProximity(chunk[,2:3], 7) # Filter by proximity
    
    num_occ_debug = nrow(chunk_trimmed)
    
    # Cast geometry back to data.table, remove "L" column, and reintroduce the species column (necessary for getDynamicAlphaHull())
    chunk_trimmed = data.table(st_coordinates(st_cast(chunk_trimmed$geometry,"MULTIPOINT")))[,-3][,species:=species_name]
    # Fix the column order (necessary for getDynamicAlphaHull())
    setcolorder(chunk_trimmed, c("species", "X", "Y")) 
    # Append the filtered chunk to the final data.table of filtered occurrences
    species_occ_trimmed = rbind(species_occ_trimmed, chunk_trimmed) 
    
    # If there were less than 2500 occurrences, we only run this loop once and can stop here
    if(run_once) { break }
    
    # Remove 1000 from the remaining number of occurrences to do
    num_occ = num_occ - 1000
    
    # If there's less than 1000 occurrences left, set the ending index to the final occurrence
    # If there's more than 1000 occurrences left, move on to the next 1000
    if(num_occ < 1000) {
      starting_index = ending_index + 1
      ending_index = nrow(species_occ)
    } else {
      starting_index = ending_index + 1
      ending_index = ending_index + 1000
    } # End if/else - ready to do the while loop again if necessary
    
  } # End while() - species_occ_trimmed now contains a reduced list of species occurrences
  
  # Randomly sample down to 4500 points if applicable (helps with running getAlphaDynamicHull)
  if(nrow(species_occ_trimmed)>4500) {species_occ_trimmed = species_occ_trimmed[sample(.N,4500)]}
  
  # Make the range from the trimmed occurrence data
  range = getDynamicAlphaHull(species_occ_trimmed, coordHeaders=c('Y','X'), partCount=5)
  
  # Save the shapefile
  st_write(range[[1]], paste0("shp/", species_name, ".shp"), append=FALSE)
  
} # End for()

# # - - - - - - - End Range work - - - - - - -  # #
# # # # # # # # # # # # # # # # # # # # # # # # # #
# # - - - - - Start Visualization work - - - -  # #

# Get a map of the world
world = rangeBuilder:::loadWorldMap()

# Get a range from shapefile
species_name = "Berberis vulgaris"
range = st_cast(st_read(paste0("shp/", species_name, ".shp"))$geometry, "MULTIPOLYGON")

# Plot the range & world map
plot(range, col=transparentColor('dark green', 0.5), border = NA)
plot(world, add = TRUE, lwd = 0.5)

# Plot the full set of unique occurrence points
species_occ = occ_trimmed[species == species_name][decimalLatitude >= 17][decimalLongitude <= -50][order(decimalLatitude, decimalLongitude)] %>% 
  unique()
points(species_occ[,c('decimalLongitude','decimalLatitude')], cex = 0.5, pch = 3)

# # - - - - - - End Visualization work - - - -  # #