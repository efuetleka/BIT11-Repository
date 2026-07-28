
# Load required packages

library(flowCore)
library(ggcyto)
library(ggplot2)
library(ggpubr)

# Root directory containing all six cytometer folders
root.dir <- root.dir <- "~/FullStains_FM/Preprocessing_results"

# Cytometers to analyse

cytometers <- c(
  
  "AuroraEvo_FM",
  
  "Opteon_FM",
  
  "A8_FM",
  
  "ID7000_FM",
  
  "Mosaic_FM",
  
  "Xenith_FM"
  
)

# Function to generate NxN scatter plots

PlotNxN <- function(ff,
                    channels,
                    lims = NULL,
                    plotFile = "PlotNxN.png")
{
  
  plots <- list()
  
  # Create all pairwise scatter plots
  
  message("Creating plots")
  
  pb_i <- 0
  
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(channels)-1,
    initial = 0,
    char = "+",
    style = 3,
    width = 50
  )
  
  utils::setTxtProgressBar(pb, pb_i)
  
  for(i in 1:(length(channels)-1))
  {
    
    for(j in (i+1):length(channels))
    {
      
      suppressMessages({
        
        p <- ggcyto::autoplot(
          
          ff,
          
          channels[j],
          
          channels[i],
          
          bins = 128
          
        
        ) +
          
          theme_minimal() +
          
          theme(
            
            legend.position = "none"
            
          )
        
        if(!is.null(lims))
        {
          
          p <- p +
            
            ggcyto::ggcyto_par_set(
              
              limits = lims
              
            )
          
        }
        
        p <- ggcyto::as.ggplot(p)
        
      })
      
      plots[[paste0(
        
        channels[i],
        
        "_",
        
        channels[j]
        
      )]] <- p
      
    }
    
    pb_i <- pb_i + 1
    
    utils::setTxtProgressBar(pb, pb_i)
    
  }
  
  # Arrange plots into NxN layout
  
  message("\nArranging plots")
  
  pb_i <- 0
  
  pb <- utils::txtProgressBar(
    
    min = 0,
    
    max = length(channels)-1,
    
    initial = 0,
    
    char = "+",
    
    style = 3,
    
    width = 50
    
  )
  
  utils::setTxtProgressBar(pb, pb_i)
  
  blank_plot <- ggplot() +
    
    theme_void()
  
  all_plots <- list()
  
  k <- 1
  
  for(row in 1:(length(channels)-1))
  {
    
    n_plots <- length(channels) - row
    
    plot_row <- plots[
      
      k:(k + n_plots - 1)
      
    ]
    
    k <- k + n_plots
    
    n_blank <- row - 1
    
    row_with_blanks <- c(
      
      rep(list(blank_plot), n_blank),
      
      plot_row
      
    )
    
    plots_row <- ggpubr::ggarrange(
      
      plotlist = row_with_blanks,
      
      ncol = length(channels)-1,
      
      nrow = 1
      
    )
    
    all_plots[[
      
      as.character(row)
      
    ]] <- plots_row
    
    pb_i <- pb_i + 1
    
    utils::setTxtProgressBar(pb, pb_i)
    
  }
  
  # Save NxN plot
  
  message("\nSaving plot")
  
  grDevices::png(
    
    filename = plotFile,
    
    width = min(
      
      500 * (length(channels)-1),
      
      20000
      
    ),
    
    height = min(
      
      500 * (length(channels)-1),
      
      20000
      
    ),
    
    res = 150
    
  )
  
  plot <- ggpubr::ggarrange(
    
    plotlist = all_plots,
    
    ncol = 1,
    
    nrow = length(channels)-1
    
  )
  
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
    
    # Cytometer folder
    
    cyto.dir <- file.path(
      root.dir,
      cyto
    )
    
    # Find all FCS files
    
    files <- list.files(
      cyto.dir,
      pattern = "\\.fcs$",
      recursive = TRUE,
      full.names = TRUE
    )
    cat(
      "Found",
      length(files),
      "FCS files\n"
    )
    
    # Loop through every sample
    
    for(f in files)
    {
      
      cat(
        "\nProcessing:",
        basename(f),
        "\n"
      )
      
      # Read FCS
      
      ff <- flowCore::read.FCS(
        f,
        transformation = FALSE,
        truncate_max_range = FALSE
      )
      
      # Select fluorescence channels
      
      if(cyto %in% c(
        "AuroraEvo_FM",
        "A8_FM",
        "Opteon_FM"
      ))
      {
        
        channels <- colnames(ff)
        
        channels <- channels[
          grepl(
            "BUV|BV|FITC|RB|PE|APC|R718|Zombie|AF",
            channels
          )
        ]
        # remove the autoflouresence channel
        channels <- channels[
          channels != "AF-A"
        ]
        
      }
      
      else if(cyto == "ID7000_FM")
      {
        
        channels <- colnames(ff)
        
        channels <- channels[
          grepl(
            "BUV|BV|FITC|RB|PE|APC|R718|Zombie",
            channels
          )
        ]
        
        channels <- channels[
          !grepl("\\[AF", channels)
        ]
        
      }
      
      else if(cyto == "Mosaic_FM")
      {
        
        channels <- colnames(ff)
        
        channels <- channels[
          grepl(
            "BUV|BV|FITC|RB|PE|APC|R718|Zombie|APCF810",
            channels
          )
        ]
        
        channels <- channels[
          grepl("-A$", channels)
        ]
        
      }
      
      else if(cyto == "Xenith_FM")
      {
        
        channels <- colnames(ff)
        
        channels <- channels[
          grepl("^FL[0-9]+-Comp$", channels)
        ]
        
      }
      
      # Show selected channels
      
      cat(
        "Using",
        length(channels),
        "channels\n"
      )
      
      print(channels)
      
      # Skip file if no channels found
      
      if(length(channels) == 0)
      {
        
        warning(
          paste(
            "No fluorescence channels found for",
            basename(f)
          )
        )
        
        next
        
      }
      
      # Create NxN plot
      
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
    cat(
      "\nCompleted:",
      cyto,
      "\n"
    )
    
  },
  
  error = function(e)
  {
    
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
