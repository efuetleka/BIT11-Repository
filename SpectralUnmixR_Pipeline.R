# Load packages
library(SpectralUnmixR)
library(flowCore)
library(Biobase)

# Root directory
root.dir <- "~/SpectralWorkshopZurich_Unmixing"

# Cytometers

cytometers <- c(
  "Aurora",
  "Opteon",
  "A8"
)

# Loop through cytometers

for (cyto in cytometers) {
  
  cat("\n=========================================\n")
  cat("Processing cytometer:", cyto, "\n")
  cat("=========================================\n\n")
  
  # Cytometer folders
  
  cyto.dir <- file.path(root.dir, cyto)
  
  singles.dir <- file.path(cyto.dir, "Singles")
  
  # Results folder
  
  results.dir <- file.path(
    cyto.dir,
    "SpectralUnmixR_Results"
  )
  
  dir.create(
    results.dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # List single stains
  
  single.files <- list.files(
    singles.dir,
    pattern = "\\.fcs$",
    full.names = TRUE
  )
  
  # Identify negative controls 
  
  if (cyto == "Aurora") {
    
    neg.bead <- single.files[
      grepl(
        "Unstained \\(Beads\\)",
        basename(single.files)
      )
    ]
    
    neg.cell <- single.files[
      grepl(
        "Neg cells",
        basename(single.files)
      )
    ]
    
  } else if (cyto == "Opteon") {
    
    neg.bead <- single.files[
      grepl(
        "Unstained beads",
        basename(single.files),
        ignore.case = TRUE
      )
    ]
    
    neg.cell <- single.files[
      grepl(
        "Unstained cells",
        basename(single.files),
        ignore.case = TRUE
      )
    ]
    
  } else if (cyto == "A8") {
    
    neg.bead <- single.files[
      grepl(
        "C01-Beads",
        basename(single.files)
      )
    ]
    
    neg.cell <- single.files[
      grepl(
        "C02-Cells",
        basename(single.files)
      )
    ]
    
  }
  
  # Positive controls
  positive.files <- setdiff(
    single.files,
    c(neg.bead, neg.cell)
  )
  
  cat("Positive controls:", length(positive.files), "\n")
  
  # Build input list
  
  input.list <- list()
  
  for (f in positive.files) {
    
    file.name <- basename(f)
    
    # Cell controls (Zombie / LD)
    
    if (grepl("Zombie|LD", file.name, ignore.case = TRUE)) {
      
      input.list[[length(input.list) + 1]] <- list(
        
        pos = f,
        
        neg = neg.cell
        
      )
      
    } else {
    
      # Bead controls
      
      input.list[[length(input.list) + 1]] <- list(
        
        pos = f,
        
        neg = neg.bead
        
      )
      
    }
    
  }
  
  # Check
  
  cat("Input list size:", length(input.list), "\n")
  
  for (i in seq_along(input.list)) {
    
    cat(
      sprintf(
        "%2d  %s\n",
        i,
        basename(input.list[[i]]$pos)
      )
    )
    
  }
  
  cat("\n")

# Build reference spectra

cat("\nBuilding reference spectra...\n")

spec <- MedianSpectra(
  input.list,
  transformation = FALSE,
  truncate_max_range = FALSE,
  emptyValue = FALSE
)
cat("\nColumns in spec:\n")
print(colnames(spec))
cat("Reference spectra built successfully.\n")

# Remove non-fluorescence channels

if (cyto == "Aurora") {
  
  remove.channels <- c(
    "Time",
    "FSC-H",
    "FSC-A",
    "SSC-H",
    "SSC-A",
    "SSC-B-H",
    "SSC-B-A",
    "file"
  )
  
} else if (cyto == "Opteon") {
  
  remove.channels <- c(
    "Time",
    "FSC-H",
    "FSC-A",
    "USSC-H","USSC-A",
    "VSSC-H","VSSC-A",
    "BSSC-H","BSSC-A",
    "YSSC-H","YSSC-A",
    "RSSC-H","RSSC-A",
    "IRSSC-H","IRSSC-A",
    "FSC-Width",
    "BSSC-Width",
    "file"
  )
  
} else if (cyto == "A8") {
  
  remove.channels <- c(
    "Time",
    "FSC-H",
    "FSC-A",
    "SSC (Violet)-H",
    "SSC (Violet)-A",
    "LightLoss (Imaging)-H",
    "LightLoss (Imaging)-A",
    "file"
  )
  
}
detector.channels <- setdiff(
  colnames(spec),
  remove.channels
)

# Build spectral matrix

M <- t(
  as.matrix(
    spec[, detector.channels]
  )
)

# Extract fluorophore names

fluor.names <- basename(spec$file)

if (cyto == "Aurora") {
  
  fluor.names <- sub("\\.fcs$", "", fluor.names)
  fluor.names <- sub(".*? ", "", fluor.names)
  fluor.names <- sub(" \\(Beads\\).*", "", fluor.names)
  fluor.names <- sub(" \\(Cells\\).*", "", fluor.names)
  
} else if (cyto == "Opteon") {
  
  fluor.names <- sub("_raw\\.fcs$", "", fluor.names)
} else if (cyto == "A8") {
  
  fluor.names <- basename(spec$file)
  
  fluor.names <- sub("\\.fcs$", "", fluor.names)
  
  fluor.names <- sub("^Single stains-[A-Z][0-9]+-", "", fluor.names)
  
  fluor.names <- vapply(
    strsplit(fluor.names, "_"),
    function(x) {
      
      if (length(x) >= 2) {
        
        paste(x[-c(1, length(x))], collapse = "_")
        
      } else {
        
        x
        
      }
      
    },
    character(1)
  )
  
  fluor.names <- gsub("_Fire", "-Fire", fluor.names)
  cat("\nFluorophore names:\n")
  print(fluor.names)
}



colnames(M) <- fluor.names

# Save reference matrix
write.csv(
  M,
  file.path(
    results.dir,
    "ReferenceSpectra.csv"
  ),
  row.names = TRUE
)
# Check

cat("\nReference matrix dimensions:\n")
print(dim(M))

cat("\nDetector channels:\n")
print(head(rownames(M)))

cat("\nFluorophores:\n")
print(colnames(M))

# Find full stain FCS files

full.files <- list.files(
  cyto.dir,
  pattern = "\\.fcs$",
  recursive = FALSE,
  full.names = TRUE
)

# Remove all files inside Singles
full.files <- full.files[
  !grepl(
    paste0(.Platform$file.sep, "Singles", .Platform$file.sep),
    full.files,
    fixed = TRUE
  )
]

cat("\nFull stain files found:", length(full.files), "\n")
print(basename(full.files))

# Loop through full stain samples

for(full.file in full.files){
  
  cat("\n-----------------------------------------\n")
  cat("Processing:", basename(full.file), "\n")
  cat("-----------------------------------------\n")

  # Sample results folder
  
  sample.name <- tools::file_path_sans_ext(
    basename(full.file)
  )
  
  sample.dir <- file.path(
    results.dir,
    sample.name
  )
  
  dir.create(
    sample.dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Read fluorescence data
  
  dat <- GetData(
    full.file,
    transformation = FALSE,
    truncate_max_range = FALSE,
    emptyValue = FALSE
  )
  
  # Run SpectralUnmixR
  
  res <- UnmixFile(
    data = dat,
    M = M,
    unmixing = "OLS",
    error = "SoS",
    write_FCS = FALSE,
    fullOutput = TRUE
  )
  
  cat("Unmixing completed.\n")
    # Keep scatter channels
    
  if (cyto == "Aurora") {
    
    scatter.channels <- intersect(
      c(
        "Time",
        "FSC-H",
        "FSC-A",
        "SSC-H",
        "SSC-A",
        "SSC-B-H",
        "SSC-B-A"
      ),
      colnames(res$unmixed)
    )
    
  } else if (cyto == "Opteon") {
    
    scatter.channels <- intersect(
      c(
        "Time",
        "FSC-H",
        "FSC-A",
        "USSC-H",
        "USSC-A",
        "VSSC-H",
        "VSSC-A",
        "BSSC-H",
        "BSSC-A",
        "YSSC-H",
        "YSSC-A",
        "RSSC-H",
        "RSSC-A",
        "IRSSC-H",
        "IRSSC-A"
      ),
      colnames(res$unmixed)
    )
    
  } else if (cyto == "A8") {
    
    scatter.channels <- intersect(
      c(
        "Time",
        "FSC-H",
        "FSC-A",
        "SSC (Violet)-H",
        "SSC (Violet)-A"
      ),
      colnames(res$unmixed)
    )
    
  }
  cat("\nColumns in res$unmixed:\n")
  print(colnames(res$unmixed))
    scatter <- res$unmixed[, scatter.channels, drop = FALSE]
    
    # Keep unmixed fluorophores
    
    fluorescence <- res$unmixed[
      ,
      colnames(M),
      drop = FALSE
    ]
    
    # Combine scatter + fluorescence
    
    final.data <- cbind(
      
      scatter,
      
      fluorescence
      
    )
    
    # Create flowFrame
    
    ff.final <- flowFrame(
      
      as.matrix(final.data)
      
    )
    
    # Update parameter information
    
    pd <- pData(parameters(ff.final))
    
    pd$name <- colnames(final.data)
    
    pd$desc <- colnames(final.data)
    
    parameters(ff.final) <- AnnotatedDataFrame(pd)
    
    # Output filename
    
    outfile <- file.path(
      
      sample.dir,
      
      paste0(
        
        sample.name,
        
        "_SpectralUnmixR.fcs"
        
      )
      
    )
    
    # Write FCS
    
    cat("\nFinal channels:\n")
    print(colnames(final.data))
    write.FCS(
      
      ff.final,
      
      outfile
      
    )
    
    cat("Saved:", basename(outfile), "\n")
  }   # end full.file loop
  
  cat("\nFinished", cyto, "\n")
  
}   # end cytometer loop

