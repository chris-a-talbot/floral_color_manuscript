library(sf) # For working with shapefiles
library(rangeBuilder) # Just for the worldmap
library(dplyr) # For coding things
library(data.table) # For data manipulation
library(stringr) # For working with strings

input_file="data/flowers_master.csv"
write=TRUE
premade_pct_mat_file = "data/distmat_range_pct_lin.csv"
premade_bin_mat_file = "data/distmat_range_bin_lin.csv"
wd="G:/My Drive/Documents/Research/Weber/honors-thesis"

sf_use_s2(FALSE)

setwd(wd)

flowers_master = fread(input_file)[order(species),]
num_species = nrow(flowers_master) 

bin_distmat = data.frame(matrix(NA, nrow = num_species, ncol = num_species))
pct_distmat = data.frame(matrix(NA, nrow = num_species, ncol = num_species))

if(!is.na(premade_pct_mat_file)) {
  pct_distmat = read.csv(premade_pct_mat_file)[,-1]
} else pct_distmat = data.frame()

if(!is.na(premade_bin_mat_file)) {
  bin_distmat = read.csv(premade_bin_mat_file)[,-1]
} else bin_distmat = data.frame()

i_start = 1010

for(i in i_start:num_species) {
  
  if(i %% 5 == 0) {
    
    pct_distmat = pct_distmat %>% as.matrix() %>% 
      Matrix::forceSymmetric(uplo="U") %>% as.matrix()
    rownames(pct_distmat) = flowers_master$species
    colnames(pct_distmat) = flowers_master$species
    pct_distmat = pct_distmat[order(rownames(pct_distmat)), 
                              order(colnames(pct_distmat))]
    
    bin_distmat = bin_distmat %>% as.matrix() %>% 
      Matrix::forceSymmetric(uplo="U") %>% as.matrix()
    rownames(bin_distmat) = flowers_master$species
    colnames(bin_distmat) = flowers_master$species
    bin_distmat = bin_distmat[order(rownames(bin_distmat)), 
                              order(colnames(bin_distmat))]
    
    write.csv(pct_distmat, "data/distmat_range_pct_lin.csv")
    write.csv(bin_distmat, "data/distmat_range_bin_lin.csv")
  }
  
  i_null = FALSE
  
  # Get the name of species i and ensure it has range data
  spc_name_1 = flowers_master$taxon[i] %>% word(1,2, sep=" ")
  spc_path_1 = paste0("shp/", spc_name_1, ".shp")
  if(!file.exists(spc_path_1)) { i_null = TRUE }
  else {
    
    # Get the range and range area for species i
    r1 = st_read(spc_path_1)$geometry %>% 
       st_cast("MULTIPOLYGON") %>% st_make_valid() %>% st_union() %>%
       st_make_valid()
    area_r1 = st_area(r1) %>% sum()
    
    # r1 = shapefiles[i]
    # area_r1 = areas[i]
    
  }
  
  for(j in 1:num_species) {
    
    if(!is.na(bin_distmat[j,i])) next
    if(!is.na(bin_distmat[i,j])) next
    
    # If species i doesn't exist, set values to NA and move on
    if(i_null) { 
      bin_distmat[i,j] = NA
      pct_distmat[i,j] = NA
      next
    }
    
    # Get the name of species j and ensure it has range data
    spc_name_2 = flowers_master$taxon[j] %>% word(1,2, sep=" ")
    spc_path_2 = paste0("shp/", spc_name_2, ".shp")
    if(!file.exists(spc_path_2)) {
      bin_distmat[i,j] = NA
      pct_distmat[i,j] = NA
      next
    }
    
    # Get the range and range area for species j
    r2 = st_read(paste0("shp/", spc_name_2, ".shp"))$geometry %>% 
      st_cast("MULTIPOLYGON") %>% st_make_valid() %>% st_union() %>% 
      st_make_valid()
    area_r2 = st_area(r2) %>% sum()
    
    # r2 = st_polygon(shapefiles[j])
    # area_r2 = areas[j]
    
    # Find the area of intersection between the ranges of species i and j
    area_intsct = st_intersection(r1, r2) %>% st_area() %>% sum()
    
    # Get the % overlap for the intersection and smaller range, save it to the matrix
    pct_overlap = as.double(area_intsct / min(c(area_r1, area_r2)))
    pct_distmat[i,j] = pct_overlap
    
    # If >= 15% overlap, consider these species spatially sympatric 
    if(pct_overlap >= 0.15) bin_distmat[i,j] = 1
    else bin_distmat[i,j] = 0
  }
  
}

pct_distmat = pct_distmat %>% as.matrix() %>% 
  Matrix::forceSymmetric(uplo="U") %>% as.matrix()
rownames(pct_distmat) = flowers_master$species
colnames(pct_distmat) = flowers_master$species
pct_distmat = pct_distmat[order(rownames(pct_distmat)), 
                          order(colnames(pct_distmat))]

bin_distmat = bin_distmat %>% as.matrix() %>% 
  Matrix::forceSymmetric(uplo="U") %>% as.matrix()
rownames(bin_distmat) = flowers_master$species
colnames(bin_distmat) = flowers_master$species
bin_distmat = bin_distmat[order(rownames(bin_distmat)), 
                          order(colnames(bin_distmat))]

if(write) {
  write.csv(pct_distmat, "data/distmat_range_pct_lin.csv")
  write.csv(bin_distmat, "data/distmat_range_bin_lin.csv")
} 

shapefiles = list()
areas = list()
for(i in 1:num_species) {
  
  spc_name = flowers_master$taxon[i] %>% word(1,2, sep=" ")
  spc_path = paste0("shp/", spc_name, ".shp")
  
  i_null = FALSE
  
  if(!file.exists(spc_path)) { 
    i_null = TRUE 
    r = NA
    area_r = NA
  }
  else {
    
    # Get the range and range area for species i
    r = st_read(spc_path)$geometry %>% 
      st_cast("MULTIPOLYGON") %>% st_make_valid() %>% st_union() %>%
      st_make_valid()
    area_r = st_area(r) %>% sum()
    
  }
  
 # for(j in 1:length(r[[1]])) {
 #   shapefiles[i,j] = r[[1]][[j]]
 # }
  
  shapefiles[i] = r
  areas[i] = area_r
}


testframe = data.frame()
















make_range_distmat = function(input_file="data/flowers_master.csv", write=TRUE,
                              premade_pct_mat_file = "data/distmat_range_pct.csv",
                              premade_bin_mat_file = "data/distmat_range_bin.csv", 
                              wd="G:/My Drive/Documents/Research/Weber/honors-thesis") {
  
  setwd(wd)
  
  flowers_master = fread(input_file)[order(species),]
  num_species = nrow(flowers_master) 
  
  bin_distmat = data.frame(matrix(NA, nrow = num_species, ncol = num_species))
  pct_distmat = data.frame(matrix(NA, nrow = num_species, ncol = num_species))
  
  if(!is.na(premade_pct_mat_file)) {
    pct_distmat = read.csv(premade_pct_mat_file)[,-1]
  } else pct_distmat = data.frame()
  
  if(!is.na(premade_bin_mat_file)) {
    bin_distmat = read.csv(premade_bin_mat_file)[,-1]
  } else bin_distmat = data.frame()
  
  for(i in 1:num_species) {
    i_null = FALSE
    
    # Get the name of species i and ensure it has range data
    spc_name_1 = flowers_master$taxon[i] %>% word(1,2, sep=" ")
    spc_path_1 = paste0("shp/", spc_name_1, ".shp")
    if(!file.exists(spc_path_1)) { i_null = TRUE }
    else {
      
      # Get the range and range area for species i
      r1 = st_read(spc_path_1)$geometry %>% 
        st_cast("MULTIPOLYGON") %>% st_make_valid()
      area_r1 = st_area(r1) %>% sum()
      
    }
    
    for(j in 1:num_species) {
      
      if(!is.na(bin_distmat[j,i])) next
      if(!is.na(bin_distmat[i,j])) next

      # If species i doesn't exist, set values to NA and move on
      if(i_null) { 
        bin_distmat[i,j] = NA
        pct_distmat[i,j] = NA
        next
      }
      
      # Get the name of species j and ensure it has range data
      spc_name_2 = flowers_master$taxon[j] %>% word(1,2, sep=" ")
      spc_path_2 = paste0("shp/", spc_name_2, ".shp")
      if(!file.exists(spc_path_2)) {
        bin_distmat[i,j] = NA
        pct_distmat[i,j] = NA
        next
      }
      
      # Get the range and range area for species j
      r2 = st_read(paste0("shp/", spc_name_2, ".shp"))$geometry %>% 
        st_cast("MULTIPOLYGON") %>% st_make_valid()
      area_r2 = st_area(r2) %>% sum()
      
      # Find the area of intersection between the ranges of species i and j
      area_intsct = st_intersection(r1, r2) %>% st_area() %>% sum()
      
      # Get the % overlap for the intersection and smaller range, save it to the matrix
      pct_overlap = as.double(area_intsct / min(c(area_r1, area_r2)))
      pct_distmat[i,j] = pct_overlap
      
      # If >= 15% overlap, consider these species spatially sympatric 
      if(pct_overlap >= 0.15) bin_distmat[i,j] = 1
      else bin_distmat[i,j] = 0
    }
    
  }
  
  pct_distmat = pct_distmat %>% as.matrix() %>% 
    Matrix::forceSymmetric(uplo="U") %>% as.matrix()
  rownames(pct_distmat) = flowers_master$species
  colnames(pct_distmat) = flowers_master$species
  pct_distmat = pct_distmat[order(rownames(pct_distmat)), 
                            order(colnames(pct_distmat))]
  
  bin_distmat = bin_distmat %>% as.matrix() %>% 
    Matrix::forceSymmetric(uplo="U") %>% as.matrix()
  rownames(bin_distmat) = flowers_master$species
  colnames(bin_distmat) = flowers_master$species
  bin_distmat = bin_distmat[order(rownames(bin_distmat)), 
                                  order(colnames(bin_distmat))]
  
  if(write) {
    write.csv(pct_distmat, "data/distmat_range_pct.csv")
    write.csv(bin_distmat, "data/distmat_range_bin.csv")
  } else return(c(pct_distmat, bin_distmat))
}

setwd("G:/My Drive/Documents/Research/Weber/honors-thesis")
write.csv(pct_distmat, "data/distmat_range_pct.csv")
write.csv(bin_distmat, "data/distmat_range_bin.csv")





























library(sf) # For working with shapefiles
library(rangeBuilder) # Just for the worldmap
library(dplyr) # For coding things
library(data.table) # For data manipulation

flowers_master = fread(input_file)
num_species = nrow(flowers_master) 

binary_distmat = data.frame()
pct_distmat = data.frame()

setwd("G:/My Drive/Documents/Research/Weber/honors-thesis")

# TO DO:
# Add tests for cases where there is no range overlap
# Add a way to import existing partial matrices
# Add a way to test for self-intersecting ranges
# Do a matrix with both directional pct overlap comparisons?
for(i in 1:num_species) {
  i_null = FALSE
  
  # Get the name of species i and ensure it has range data
  spc_name_1 = flowers_master$species[i]
  spc_path_1 = paste0("shp/", spc_name_1, ".shp")
  if(!file.exists(spc_path_1)) { i_null = TRUE }
  else {
  
    # Get the range and range area for species i
    r1 = st_read(spc_path_1)$geometry %>% 
      st_cast("MULTIPOLYGON") %>% st_make_valid()
    area_r1 = st_area(r1) %>% sum()
  
  }

  for(j in 1:num_species) {
    
    # If this comparison has already been made, don't do it again
    # if(!is.null(pct_distmat[j,i]) && 
    #   !is.null(binary_distmat[j,i])) next
    # if(!is.null(pct_distmat[i,j]) && 
    #   !is.null(binary_distmat[i,j])) next
    
    # If species i doesn't exist, set values to NA and move on
    if(i_null) { 
      binary_distmat[i,j] = NA
      pct_distmat[i,j] = NA
      next
    }
    
    # Get the name of species j and ensure it has range data
    spc_name_2 = flowers_master$species[j]
    spc_path_2 = paste0("shp/", spc_name_2, ".shp")
    if(!file.exists(spc_path_2)) {
      binary_distmat[i,j] = NA
      pct_distmat[i,j] = NA
      next
    }
    
    # Get the range and range area for species j
    r2 = st_read(paste0("shp/", spc_name_2, ".shp"))$geometry %>% 
      st_cast("MULTIPOLYGON") %>% st_make_valid()
    area_r2 = st_area(r2) %>% sum()

    # Find the area of intersection between the ranges of species i and j
    area_intsct = st_intersection(r1, r2) %>% st_area() %>% sum()
    
    # Get the % overlap for the intersection and smaller range, save it to the matrix
    pct_overlap = as.double(area_intsct / min(c(area_r1, area_r2)))
    pct_distmat[i,j] = pct_overlap
  
    # If >= 15% overlap, consider these species spatially sympatric 
    if(pct_overlap >= 0.15) binary_distmat[i,j] = 1
    else binary_distmat[i,j] = 0
  }
  
}

pct_distmat = pct_distmat %>% as.matrix() %>% 
  Matrix::forceSymmetric(uplo="L") %>% as.matrix()
rownames(pct_distmat) = flowers_master$species
colnames(pct_distmat) = flowers_master$species
pct_distmat = pct_distmat[order(rownames(pct_distmat)), 
                          order(colnames(pct_distmat))]

binary_distmat = binary_distmat %>% as.matrix() %>% 
  Matrix::forceSymmetric(uplo="L") %>% as.matrix()
rownames(binary_distmat) = flowers_master$species
colnames(binary_distmat) = flowers_master$species
binary_distmat = binary_distmat[order(rownames(binary_distmat)), 
                                order(colnames(binary_distmat))]










# Get the names of the current species to be compared
spc_name_1 = flowers_master$species[300]
spc_name_2 = flowers_master$species[5]

# Get their ranges
r1 = st_read(paste0("shp/", spc_name_1, ".shp"))$geometry %>% 
  st_cast("MULTIPOLYGON") %>% st_make_valid()
r2 = st_read(paste0("shp/", spc_name_2, ".shp"))$geometry %>% 
  st_cast("MULTIPOLYGON") %>% st_make_valid()

# Get the area of their ranges and determine which is smaller
area_r1 = st_area(r1) %>% sum()
area_r2 = st_area(r2) %>% sum()

# Find the intersection of the two ranges and its area
intsct = st_intersection(r1, r2)
area_intsct = st_area(intsct) %>% sum()

# Find the % overlap of the smaller range with the intersection of the two
pct_overlap = as.double(area_intsct / min(c(area_r1, area_r2)))
                        
# If the % overlap is greater than 10%, consider them sympatric
if(pct_overlap <= 0.1) sympatry = FALSE else sympatry = TRUE


# mapping
world = rangeBuilder:::loadWorldMap()
plot(r1, col=transparentColor('dark green', 0.5), border = NA)
plot(r2, col=transparentColor('dark blue', 0.5), border = NA, add=TRUE)
plot(intsct, col=transparentColor('dark orange', 0.5), border = NA, add=TRUE)
plot(world, add = TRUE, lwd = 0.5)
