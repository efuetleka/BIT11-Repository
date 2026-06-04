library(AutoSpectral)
library(flowCore)

# PATHS

root.dir <- "~/SpectralWorkshopZurich_Unmixing/Aurora"
setwd(root.dir)
# create directory for the control files 
control.dir <- file.path(root.dir, "Singles")
# create a folder to store the unmixing results

result.dir <- file.path(root.dir, "Results")
dir.create(result.dir, showWarnings = FALSE)

# define the cytometer parameters 
asp <- get.autospectral.param(
  cytometer = "aurora",
  figures = TRUE
)

# create a control file for autospectral

control.file <- file.path(getwd(), "fcs_control_file.csv")
# this prevents overwriting the file if code is runned second time 
if (!file.exists(control.file)) {
  create.control.file(control.dir, asp)
} else {
  message("Control file already exists. Using existing file.")
}
# generate a flow.control obeject

flow.control <- define.flow.control(
  control.dir,
  control.file,
  asp
)
# clean the controls to avoid contaminating spectral signatures 
flow.control <- clean.controls(
  flow.control,
  asp
)
# generate a spectral matrix that stores each flourophores intensity per detector
spectra <- get.fluorophore.spectra(
  flow.control,
  asp
)
# create an auto flouresence spectral matrix
af.file <- file.path(
  control.dir,
  "C1 Neg cells BSB (Cells)_plate FM.fcs"
)
# check if the file exist
file.exists(af.file)

af.spectra <- get.af.spectra(
  unstained.sample = af.file,
  asp = asp,
  spectra = spectra
)
# generate a spectral heatmap 
spectral.heatmap(
  spectra = spectra,
  title = "Aurora Spectral Heatmap",
  save = FALSE
)
# plot the spectra traces
spectral.trace(
  spectral.matrix = spectra,
  asp = asp,
  title = "Aurora Spectral Signatures",
  save = FALSE
)
# unmix the full stain files 
unmix.folder(
  fcs.dir = root.dir,
  spectra = spectra,
  af.spectra = af.spectra,
  asp = asp,
  flow.control = flow.control,
  method = "OLS",
  output.dir = result.dir,
  verbose = TRUE
)


