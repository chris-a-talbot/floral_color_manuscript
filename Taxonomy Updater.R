update_names = function(wd="G:/My Drive/Documents/Research/Weber/honors-thesis/data", 
                        input_file="flowers_master.csv", write=TRUE, output_file = "flowers_master_updated.csv",
                        overwrite=FALSE) {
  library(lcvplants)
  
  setwd(wd)
  
  flowers_master = fread(input_file)
  
  for(i in 1:nrow(flowers_master)){
    species = flowers_master$species[i] # fix this line
    taxon = flowers_master$taxon[i] # fix this line
    if((overwrite==TRUE) || (is.na(taxon))) {
      result = lcvp_search(species, max_distance = 2)
      species_data$taxon[i] = result["Output.Taxon"]
      species_data$family[i] = result["Family"]
    }
  }
  
  if(write) write.csv(flowers_master, output_file, row.names=FALSE)
}
