#  Load the flowCore package
library(flowCore)

# define the path to the Main folder containing all cytometer folders

root.dir <- "~/SpectralWorkshopZurich_Unmixing"
# create a Vector containing the names of the cytometers
cytometers <- c(
  "A8",
  "Aurora",
  "Opteon"
)

# Event count report
# Loop through each cytometer
for (cyto in cytometers) {
  
  cat("\n====================================\n")
  cat("Checking:", cyto, "\n")
  cat("====================================\n")
# Build path to the Singles folder
  
  singles.dir <- file.path(
    root.dir,
    cyto,
    "Singles"
  )
# Find all FCS files inside the Singles folder
  fcs.files <- list.files(
    singles.dir,
    pattern = "\\.fcs$",
    full.names = TRUE
  )
  
  event.report <- data.frame()
  
  for (f in fcs.files) {
    
    ff <- read.FCS(
      f,
      transformation = FALSE
    )
    
    n.events <- nrow(
      exprs(ff)
    )
    
    event.report <- rbind(
      event.report,
      data.frame(
        File = basename(f),
        Events = n.events,
        Low_Event_Count = n.events < 1000
      )
    )
  }
  
  event.report <- event.report[
    order(event.report$Events),
  ]
  
  print(event.report)
  
  cat(
    "\nFiles below 1000 events:",
    sum(event.report$Low_Event_Count),
    "\n"
  )
  
  write.csv(
    event.report,
    file.path(
      root.dir,
      cyto,
      paste0(
        cyto,
        "_event_count_report.csv"
      )
    ),
    row.names = FALSE
  )
}
#########################################################
library(flowCore)
library(ggplot2)

pos <- read.FCS(
  file.path(
    control.dir,
    "FITC_raw.fcs"
  ),
  transformation = FALSE
)

neg <- read.FCS(
  file.path(
    control.dir,
    "Unstained beads_raw.fcs"
  ),
  transformation = FALSE
)

df <- rbind(
  data.frame(
    Signal = exprs(pos)[, "B525-A"],
    Group = "FITC"
  ),
  data.frame(
    Signal = exprs(neg)[, "B525-A"],
    Group = "Negative"
  )
)

ggplot(
  df,
  aes(
    x = Group,
    y = Signal,
    colour = Group
  )
) +
  geom_jitter(
    width = 0.2,
    alpha = 0.3
  ) +
  theme_minimal() +
  ggtitle("FITC bead positivity check")
