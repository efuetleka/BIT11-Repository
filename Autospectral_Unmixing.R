if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages("remotes")
devtools::install_github("drcytometer/AutoSpectral")
library(flowCore)
library(ggplot2)
library(pheatmap)
library(AutoSpectral)
# define the main folder containing all the files 

main_dir <- "~/SpectralWorkshopZurich_Unmixing"

# Cytometers
cytometers <- c("A8", "Aurora", "Opteon")

# Create result folder

output_dir <- file.path(
  main_dir,
  "AutoSpectral_results"
)

dir.create(
  output_dir,
  showWarnings = FALSE
)

# Loop through cytometers

for (cyto in cytometers) {
  # print progress message 
  cat("\n====================================\n")
  cat("Processing:", cyto, "\n")
  cat("====================================\n")
  
# create a path to the current cytometer folder 
  cyto_dir <- file.path(
    main_dir,
    cyto
  )
  # create a path to the singles folder for that cytometer 
  singles_dir <- file.path(
    cyto_dir,
    "singles"
  )
  # create an output folder for each cytometer
  cyto_output <- file.path(
    output_dir,
    cyto
  )
  
  dir.create(
    cyto_output,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Find full-stain files
 
  full_files <- list.files(
    cyto_dir,
    pattern = "\\.fcs$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  # print the number of full stain files found 
  cat(
    "Full-stain files found:",
    length(full_files),
    "\n"
  )
  
  # Find single-stain files
 
  single_files <- list.files(
    singles_dir,
    pattern = "\\.fcs$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  # print the number of single stain controlls found 
  cat(
    "Single-stain files found:",
    length(single_files),
    "\n"
  )
  
  # Create AutoSpectral settings
 
  asp <- get.autospectral.param(
    cytometer = tolower(cyto)
  )
  
  # Read all single-stain FCS files into R
  
  single_ff <- lapply(
    single_files,
    function(f) {
      
      read.FCS(
        f,
        transformation = FALSE,
        truncate_max_range = FALSE
      )
    }
  )
  # Assign the FCS filenames as names of the flowFrame list elements
  names(single_ff) <- basename(single_files)
  
  cat(
    "Single-stain files loaded successfully\n"
  )
  
  # Create an empty list that will store the spectral signatures extracted
  # from each single-stain control
  
  spectra_list <- list()
  # loop through every single stain control 
  for (i in seq_along(single_ff)) {
    
    # Extract one single-stain flowFrame from the list
    ff <- single_ff[[i]]
    # Extract the numeric fluorescence intensity matrix from the flowFrame
    expr <- exprs(ff)
    
    # Check which columns contain numeric fluorescence data
    numeric_channels <- apply(
      expr,
      2,
      is.numeric
    )
    # keep only the numeric flourescence channels
    expr <- expr[, numeric_channels]
    
    # Calculate the mean fluorescence intensity for every channel
    #This creates a simplified spectral signature for the fluorophore
    
    mean_spec <- colMeans(
      expr,
      na.rm = TRUE
    )
    # Store the calculated spectral signature inside the spectra list
    spectra_list[[i]] <- mean_spec
  }
  
  
  # Find the smallest spectrum length to ensure all 
  # spectra can be combined consistently.
  min_len <- min(
    sapply(spectra_list, length)
  )
  # Combine all spectral signatures into one matrix
  spectra_matrix <- do.call(
    cbind,
    lapply(
      spectra_list,
      function(x) x[1:min_len]
    )
  )
  # Assign channel names to the spectra matrix rows
  rownames(spectra_matrix) <- names(
    spectra_list[[1]][1:min_len]
  )
  # Assign single-stain filenames as column names
  colnames(spectra_matrix) <- basename(
    single_files
  )
  
  # Save spectra matrix
  
  write.csv(
    spectra_matrix,
    file.path(
      cyto_output,
      paste0(
        cyto,
        "_spectra_matrix.csv"
      )
    )
  )
  
  cat(
    "Spectra matrix saved\n"
  )
  
  # Creates a folder for processed/unmixed FCS files
  
  unmixed_dir <- file.path(
    cyto_output,
    "unmixed_files"
  )
  # Create the unmixed output directory
  dir.create(
    unmixed_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  # Loops through each full-stain sample
  for (f in full_files) {
    
   # Display the current full-stain file being processed 
    cat(
      "Processing full-stain file:\n",
      basename(f),
      "\n"
    )
    # Read the full-stain raw FCS file
    ff <- read.FCS(
      f,
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    
    # save the processed file 
    write.FCS(
      ff,
      filename = file.path(
        unmixed_dir,
        paste0(
          "processed_",
          basename(f)
        )
      )
    )
  }
  # print completion message for the current cytometer
  cat(
    "Finished cytometer:",
    cyto,
    "\n"
  )
}
# print completion message after processing the full workflow 
cat("Workflow completed successfully\n")

