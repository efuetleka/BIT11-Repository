# Load packages
library(flowCore)
library(CompensAID)
library(ggplot2)

# Read one cleaned FCS file
clean_file <- clean_file <- "~/SpectralWorkshopZurich_Unmixing/AutoSpectral_Preprocessing/Aurora/H1 Well_001_plate FM OLS/H1 Well_001_plate FM OLS_cleaned.fcs"
ff_clean <- read.FCS(
  clean_file,
  transformation = FALSE,
  truncate_max_range = FALSE,
  emptyValue = FALSE
)

# Keep only marker channels, not FSC/SSC/Time
# extract the metadata
meta <- ff_clean@parameters@data
# replace channels with missing description 
meta$desc[is.na(meta$desc)] <- ""

marker_channels <- meta$name[
  !grepl("FSC|SSC|Time|Original_ID", meta$name, ignore.case = TRUE) &
    !grepl("FSC|SSC|Time|Original_ID", meta$desc, ignore.case = TRUE)
]
# create a new flow frame containing  only the markers of interest
marker_channels <- marker_channels[
  marker_channels != "AF-A"
]

ff_markers <- ff_clean[, marker_channels]
## Give every channel a unique marker names

pd <- pData(parameters(ff_markers))

pd$desc <- pd$name

parameters(ff_markers) <- AnnotatedDataFrame(pd)

flowCore::markernames(ff_markers)
cat("Marker names:\n")
print(flowCore::markernames(ff_markers))

cat("\nChannel names:\n")
print(colnames(exprs(ff_markers)))
# Run CompensAID
res_compensaid <- CompensAID::CompensAID(
  ff = ff_markers
)

# Plot SSI matrix
ssi_matrix_plot <- CompensAID::PlotMatrix(
  output = res_compensaid
)

print(ssi_matrix_plot)

# Find suspicious marker combinations and return combinations less than -1
index <- which(res_compensaid[["matrix"]] < -1, arr.ind = TRUE)
# create a table for flagged paires from the matrix using channel names 
flagged_combinations <- data.frame(
  primary_channel = rownames(res_compensaid[["matrix"]])[index[, "col"]],
  secondary_channel = rownames(res_compensaid[["matrix"]])[index[, "row"]]
)
# display the flagged combinations 
flagged_combinations

# visualize the first flagged pair:
df_pair <- data.frame(
  BUV661_A = exprs(ff_markers)[, "BUV661-A"],
  APC_A = exprs(ff_markers)[, "APC-A"]
)
# same randon sample every run 
set.seed(123)
df_pair <- df_pair[sample(nrow(df_pair), min(10000, nrow(df_pair))), ]
# create scatter plot for the flagged pair 
ggplot(df_pair, aes(x = BUV661_A, y = APC_A)) +
  geom_point(color = "blue", size = 0.4, alpha = 0.5) +
  xlab("BUV661-A: CD38") +
  ylab("APC-A: CD127") +
  ggtitle("Flagged CompensAID pair: BUV661-A vs APC-A") +
  theme_minimal()
