library(AutoSpectral)
library(flowCore)

# ============================================================
# Root directory
# ============================================================

root.dir <- "~/SpectralWorkshopZurich_Unmixing"

# ============================================================
# Cytometers
# ============================================================

cytometers <- c(
  "A8",
  "Aurora",
  "Opteon"
)

# ============================================================
# Process each cytometer
# ============================================================

for (cyto in cytometers) {
  
  tryCatch({
    
    cat("\n====================================\n")
    cat("Processing:", cyto, "\n")
    cat("====================================\n")
    
    cyto.dir <- file.path(root.dir, cyto)
    
    control.dir <- file.path(
      cyto.dir,
      "Singles"
    )
    
    result.dir <- file.path(
      cyto.dir,
      "Results"
    )
    
    dir.create(
      result.dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # --------------------------------------------------------
    # AutoSpectral parameters
    # --------------------------------------------------------
    
    asp <- get.autospectral.param(
      cytometer = tolower(cyto),
      figures = TRUE
    )
    
    # --------------------------------------------------------
    # Opteon scatter fix
    # --------------------------------------------------------
    
    if (cyto == "Opteon") {
      
      cat("Applying Opteon scatter fix\n")
      
      asp$default.scatter.parameter <- c(
        "FSC-A",
        "BSSC-A"
      )
    }
    
    # --------------------------------------------------------
    # Control file
    # --------------------------------------------------------
    
    control.file <- file.path(
      cyto.dir,
      "fcs_control_file.csv"
    )
    
    setwd(cyto.dir)
    
    # Create control file if it does not exist
    
    if (!file.exists(control.file)) {
      
      cat("Creating control file...\n")
      
      create.control.file(
        control.dir,
        asp
      )
      
      cat(
        "\nControl file created for ",
        cyto,
        "\nEdit the file and rerun the script.\n"
      )
      
      next
    }
    
    cat("Using existing control file\n")
    
    # --------------------------------------------------------
    # Create flow.control
    # --------------------------------------------------------
    
    flow.control <- define.flow.control(
      control.dir,
      control.file,
      asp
    )
    
    # --------------------------------------------------------
    # Clean controls
    # --------------------------------------------------------
    
    flow.control <- clean.controls(
      flow.control,
      asp
    )
    
    # --------------------------------------------------------
    # Generate fluorophore spectra
    # --------------------------------------------------------
    
    cat("Generating fluorophore spectra...\n")
    
    spectra <- get.fluorophore.spectra(
      flow.control,
      asp
    )
    
    # --------------------------------------------------------
    # Generate AF spectra
    # --------------------------------------------------------
    
    ctrl <- read.csv(
      control.file,
      stringsAsFactors = FALSE
    )
    
    af.file <- file.path(
      control.dir,
      ctrl$filename[
        ctrl$fluorophore == "AF"
      ]
    )
    file.exists(af.file)
    cat(
      "Generating AF spectra from:",
      basename(af.file),
      "\n"
    )
    
    af.spectra <- get.af.spectra(
      unstained.sample = af.file,
      asp = asp,
      spectra = spectra
    )
    
    # --------------------------------------------------------
    # QC plots
    # --------------------------------------------------------
    
    cat("Generating QC plots...\n")
    
    spectral.heatmap(
      spectra = spectra,
      title = paste(cyto, "Spectral Heatmap"),
      plot.dir = result.dir,
      save = TRUE
    )
    
    spectral.trace(
      spectral.matrix = spectra,
      asp = asp,
      title = paste(cyto, "Spectral Signatures"),
      plot.dir = result.dir,
      save = TRUE
    )
    
    # --------------------------------------------------------
    # Unmix full stain files
    # --------------------------------------------------------
    
    cat("Unmixing files...\n")
    
    unmix.folder(
      fcs.dir = cyto.dir,
      spectra = spectra,
      af.spectra = af.spectra,
      asp = asp,
      flow.control = flow.control,
      method = "OLS",
      output.dir = result.dir,
      verbose = TRUE
    )
    
    cat("\nFinished:", cyto, "\n")
    
  }, error = function(e) {
    
    cat(
      "\nFAILED:",
      cyto,
      "\n",
      conditionMessage(e),
      "\n"
    )
    
  })
  
}
