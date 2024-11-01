# remotes::install_github("chris-a-talbot/rpheno")

library(rpheno)
library(data.table)
library(dplyr)
library(furrr)
library(future)
library(stringr)

# Function to find the raw phenology RDS file with the largest suffix number
find_largest_raw_rds <- function(directory) {
  if(file.exists(paste0(directory, "/raw_phenos_final.rds"))){
    return("raw_phenos_final.rds")
  }
  
  rds_files <- list.files(directory, pattern = "^raw_phenos_part_\\d+\\.rds$", full.names = TRUE)
  
  if (length(rds_files) == 0) {
    return(NULL)
  }
  
  # Extract the suffix numbers
  suffix_numbers <- as.numeric(gsub("^raw_phenos_part_(\\d+)\\.rds$", "\\1", basename(rds_files)))
  
  # Find the file with the largest suffix number
  largest_file <- rds_files[which.max(suffix_numbers)]
  
  return(largest_file)
}

# Function to read and combine phenologies
process_phenology_files <- function(directory, save_interval = 20, overwrite = FALSE, overwrite_raw = FALSE) {
  raw_phenos_file <- find_largest_raw_rds(directory)
  
  if (file.exists(paste0(directory, "/cleaned_phenos_final.rds")) && !overwrite) {
    cleaned_phenos <- readRDS(paste0(directory, "/cleaned_phenos_final.rds"))
    message("Loaded cleaned phenology objects from 'cleaned_phenos_final.rds'")
    return(cleaned_phenos)
  }
  
  if (!overwrite_raw && !is.null(raw_phenos_file)) {
    raw_phenos <- readRDS(paste0(directory, raw_phenos_file))
    message(paste("Loaded raw phenology objects from", raw_phenos_file))
  } else {
    # List all CSV files in the directory
    files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
    
    raw_phenos <- list()
    valid_files <- list()
    
    # Read phenology files and store them in a list
    for (i in seq_along(files)) {
      result <- tryCatch({
        pheno <- pheno_load(files[i])
        list(pheno = pheno, file = files[i])
      }, error = function(e) {
        message(paste("Skipping file due to error:", files[i], " - ", e$message))
        NULL
      })
      
      if (!is.null(result)) {
        if(str_equal(result$pheno$name, "Asclepias sullivantii")){
          next
        }
        raw_phenos <- c(raw_phenos, list(result$pheno))
        valid_files <- c(valid_files, result$file)
      }
      
      # Save intermediate raw phenologies every 'save_interval' files
      if (i %% save_interval == 0) {
        saveRDS(raw_phenos, file = paste0(directory, "raw_phenos_part_", i, ".rds"))
      }
    }
    
    # Save the final raw phenologies
    saveRDS(raw_phenos, file = paste0(directory, "raw_phenos_final.rds"))
  }
  
  # Combine phenologies with the same name
  reduced_phenos <- list()
  phenos_by_name <- split(raw_phenos, sapply(raw_phenos, function(pheno) pheno$name))
  
  for (name in names(phenos_by_name)) {
    if (length(phenos_by_name[[name]]) == 1) {
      reduced_phenos[[name]] <- phenos_by_name[[name]]
    } else {
      reduced_phenos[[name]] <- phenos_by_name[[name]][1]
    }
  }
  
  reduced_phenos <- sapply(reduced_phenos, `[[`, 1, simplify = FALSE)
  
  # Clean all phenologies
  cleaned_phenos <- lapply(reduced_phenos, function(pheno) {
    tryCatch({
      pheno_clean(pheno)
    }, error = function(e) {
      message(paste("Error cleaning phenology for:", pheno$name, " - ", e$message))
      NULL
    })
  })
  
  # Remove NULL values from the cleaned phenologies
  cleaned_phenos <- cleaned_phenos[!sapply(cleaned_phenos, is.null)]
  
  # Name each object in the cleaned_phenos list
  names(cleaned_phenos) <- sapply(cleaned_phenos, function(pheno) pheno$name)
  
  # Save the final cleaned phenologies
  saveRDS(cleaned_phenos, file = paste0(directory, "cleaned_phenos_final.rds"))
  
  return(cleaned_phenos)
}

create_overlap_matrix <- function(cleaned_phenos, directory="./data/phenology/", file="pheno_binary_overlap_matrix.rds", write=F, overwrite=F) {
  if(!overwrite) {
    if(file.exists(paste0(directory, "pheno_binary_overlap_matrix.rds"))){
      overlap_matrix <- readRDS(paste0(directory, "pheno_binary_overlap_matrix.rds"))
      message("Loaded phenology overlap matrix from 'pheno_binary_overlap_matrix.rds'")
      return(overlap_matrix)
    }
  }
  
  # Get the number of species
  num_species <- length(cleaned_phenos)
  
  # Initialize an empty matrix to store the overlap values
  overlap_matrix <- matrix(NA, nrow = num_species, ncol = num_species)
  colnames(overlap_matrix) <- names(cleaned_phenos)
  rownames(overlap_matrix) <- names(cleaned_phenos)
  
  # Use parallel processing for overlap calculations
  plan(multisession)
  
  # Define a function to calculate overlap for each pair
  calculate_overlap <- function(i, j) {
    tryCatch({
      pheno_overlap(cleaned_phenos[[i]], cleaned_phenos[[j]])
    }, error = function(e) {
      message(paste("Error calculating overlap for indices", i, j, " - ", e$message))
      NA
    })
  }
  
  # Use future_map2 to calculate overlaps in parallel
  overlap_indices <- expand.grid(1:num_species, 1:num_species) %>%
    filter(Var1 < Var2)  # Only upper triangular
  
  overlaps <- future_map2_dbl(overlap_indices$Var1, overlap_indices$Var2, calculate_overlap, .options = furrr_options(seed = TRUE))
  
  # Fill the upper triangular matrix with the calculated overlaps
  overlap_matrix[upper.tri(overlap_matrix)] <- overlaps
  
  overlap_matrix <- as.data.frame(overlap_matrix)
  
  colnames(overlap_matrix) <- str_replace(colnames(overlap_matrix), " ", "_")
  rownames(overlap_matrix) <- str_replace(rownames(overlap_matrix), " ", "_")
  
  overlap_matrix <- as.matrix(overlap_matrix)
  
  # Replace NA with 0
  overlap_matrix[is.na(overlap_matrix)] <- 0
  
  # Calculate the sum of the matrix and its transpose
  sum_matrix <- overlap_matrix + t(overlap_matrix)
  
  # Set the diagonal elements to NA
  diag(sum_matrix) <- NA
  
  saveRDS(sum_matrix, paste0(directory, file))
  
  return(sum_matrix)
}

# Use the function to process your data
directory <- "./data/phenology/"
cleaned_phenos <- process_phenology_files(directory, overwrite=F, overwrite_raw=F)
overlap_matrix <- create_overlap_matrix(cleaned_phenos, directory, write=F, overwrite=F)
cleaned_phenos_NA <- list()
for(i in 1:985) {
  data <- cleaned_phenos[[i]]$occurrences
  new_data <- data[latitude > 38 & latitude < 48]
  if(nrow(new_data) <= 7) next
  new_data[, scientific_name := cleaned_phenos[[i]]$name]
  phen <- phenology(new_data)
  cleaned_phenos_NA <- c(cleaned_phenos_NA, list(phen))
}
names(cleaned_phenos_NA) <- sapply(cleaned_phenos_NA, function(pheno) pheno$name)
overlap_matrix_min <- create_overlap_matrix(cleaned_phenos_NA, directory, file="pheno_binary_overlap_matrix_minlat.rds", write=T, overwrite=T)

