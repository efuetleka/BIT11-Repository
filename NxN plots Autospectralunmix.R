
# Load packages

library(flowCore)
library(ggcyto)
library(ggplot2)
library(ggpubr)

# Root directory

root.dir <- "~/SpectralWorkshopZurich_Unmixing"

# Cytometers

cytometers <- c(
  "Aurora",
  "Opteon",
  "A8"
)
# function for NxN plot

PlotNxN <- function(ff,
                    channels,
                    lims = NULL,
                    plotFile = "PlotNxN.png"){
  plots <- list()
  
  message("Creating plots")
  pb_i <- 0
  pb <- utils::txtProgressBar(min = 0, max = length(channels)-1, 
                              initial = 0, char = "+", style = 3, width = 50)
  utils::setTxtProgressBar(pb, pb_i)
  
  for(i in 1:(length(channels) - 1)){
    for(j in (i+1):length(channels)){
      suppressMessages({
        p <- ggcyto::autoplot(ff, 
                              channels[j], channels[i], 
                              bins = 128) +
          theme_minimal() +
          
          theme(legend.position = "none")
        
        if(!is.null(lims)){
          p <- p +
            ggcyto::ggcyto_par_set(limits = lims)
        }
        
        p <- ggcyto::as.ggplot(p)
      })
      
      plots[[paste0(channels[i], "_", channels[j])]] <- p
    }
    pb_i <- pb_i+1
    utils::setTxtProgressBar(pb, pb_i)
  }
  
  message("\nArranging plots")
  pb_i <- 0
  pb <- utils::txtProgressBar(min = 0, max = length(channels)-1, 
                              initial = 0, char = "+", style = 3, width = 50)
  utils::setTxtProgressBar(pb, pb_i)
  
  blank_plot <- ggplot() + theme_void()
  
  all_plots <- list()
  k <- 1
  for(row in 1:(length(channels) - 1)){
    n_plots <- length(channels) - row
    plot_row <- plots[k:(k + n_plots - 1)]
    k <- k + n_plots
    
    n_blank <- row - 1
    row_with_blanks <- c(rep(list(blank_plot), n_blank), plot_row)
    
    plots_row <- ggpubr::ggarrange(plotlist = row_with_blanks,
                                   ncol = length(channels)-1,
                                   nrow = 1)
    
    all_plots[[as.character(row)]] <- plots_row
    
    pb_i <- pb_i+1
    utils::setTxtProgressBar(pb, pb_i)
  }
  
  message("\nPrinting")
  
  grDevices::png(
                 filename = plotFile,
                 width = min(500 * length(channels)-1, 20000), 
                 height = min(500 * length(channels)-1, 20000),
                 res = 150
                 )
  plot <- ggpubr::ggarrange(plotlist = all_plots,
                            ncol = 1, 
                            nrow = length(channels)-1)
  print(plot)
  grDevices::dev.off()
}
# Loop through all cytometers
for(cyto in cytometers)
{
  
  cat("\n====================================\n")
  cat("Processing:", cyto, "\n")
  cat("====================================\n")
  
  tryCatch({
    # AutoSpectral results folder
    preprocess.dir <- file.path(
      root.dir,
      "AutoSpectral_Preprocessing",
      cyto
    )
    # Find all unmixed OLS FCS files
    files <- list.files(
      preprocess.dir,
      pattern = "OLS_cleaned\\.fcs$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    
    cat(
      "Found",
      length(files),
      "AutoSpectral cleaned unmixed FCS files\n"
    )
    # Loop through every file
    for(f in files)
    {
      
      cat(
        "\nCreating NxN plot for:",
        basename(f),
        "\n"
      )
      # Read FCS
      ff <- flowCore::read.FCS(
        f,
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      # Keep only fluorescence channels
      channels <- colnames(ff)
      
      channels <- channels[
        grepl(
          "BUV|BV|FITC|RB|PE|APC|R718|Zombie|AF",
          channels
        )
      ]
      # Remove Autofluorescence channel
      channels <- channels[
        channels != "AF-A"
      ]
      # Show channels being plotted
      print(channels)
      # Generate NxN plot
      PlotNxN(
        ff = ff,
        channels = channels,
        plotFile = file.path(
          dirname(f),
          paste0(
            tools::file_path_sans_ext(
              basename(f)
            ),
            "_NxN.png"
          )
        )
      )
      
      cat(
        "Saved:",
        paste0(
          tools::file_path_sans_ext(
            basename(f)
          ),
          "_NxN.png"
        ),
        "\n"
      )
    }
    # Finished processing all files for this cytometer
    cat(
      "\nCompleted:",
      cyto,
      "\n"
    )
    
  },
  error = function(e){
    
    cat(
      "\nFAILED:",
      cyto,
      "\n"
    )
    
    cat(
      conditionMessage(e),
      "\n"
    )
    
  })
  
}
