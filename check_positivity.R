library(flowCore)
library(ggplot2)

root.dir <- "~/SpectralWorkshopZurich_Unmixing"

cytometers <- c(
  "Aurora",
  "Opteon"
)

for(cyto in cytometers){
  
  cat(
    "\n====================================\n",
    "Processing:",
    cyto,
    "\n====================================\n"
  )
  
  control.dir <- file.path(
    root.dir,
    cyto,
    "Singles"
  )
  
  control.file <- file.path(
    root.dir,
    cyto,
    "fcs_control_file.csv"
  )
  
  plot.dir <- file.path(
    root.dir,
    cyto,
    "Bead_Positivity_QC"
  )
  
  dir.create(
    plot.dir,
    showWarnings = FALSE
  )
  
  ctrl <- read.csv(
    control.file,
    stringsAsFactors = FALSE
  )
  
  beads <- subset(
    ctrl,
    control.type == "beads" &
      fluorophore != "Negative"
  )
  
  for(i in seq_len(nrow(beads))){
    
    ff <- read.FCS(
      file.path(
        control.dir,
        beads$filename[i]
      ),
      transformation = FALSE
    )
    
    detector <- beads$channel[i]
    
    if(is.na(detector) || detector == ""){
      next
    }
    
    df <- data.frame(
      Event = seq_len(
        nrow(exprs(ff))
      ),
      Signal = exprs(ff)[, detector]
    )
    
    p <- ggplot(
      df,
      aes(
        x = Event,
        y = Signal
      )
    ) +
      geom_point(
        size = 0.3,
        alpha = 0.4
      ) +
      theme_minimal() +
      labs(
        title = paste(
          cyto,
          "-",
          beads$fluorophore[i]
        ),
        subtitle = paste(
          "Main detector:",
          detector
        ),
        x = "Event",
        y = detector
      )
    
    ggsave(
      file.path(
        plot.dir,
        paste0(
          beads$fluorophore[i],
          "_positivity.png"
        )
      ),
      p,
      width = 7,
      height = 5
    )
  }
}
