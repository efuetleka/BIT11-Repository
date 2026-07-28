# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
# Install CompensAID from GitHub
devtools::install_github("Olsman/CompensAID")
library(flowCore)
library(CytoExploreR)
library(PeacoQC)
library(ggplot2)
library(ggpubr)
library(Biobase)
library(CompensAID)
library(pheatmap)
library(grid)

# Root directory
root.dir <- "~/SpectralWorkshopZurich_Unmixing"
# Cytometers
cytometers <- c(
  "Aurora",
  "Opteon",
  "A8"
)
   # Find AutoSpectral unmixed FCS files
fcs_files <- c()

for(cyto in cytometers)
{
  
  result.dir <- file.path(
    root.dir,
    cyto,
    "Results"
  )
  
  files <- list.files(
    result.dir,
    pattern = "\\.fcs$",
    full.names = TRUE
  )
  
  fcs_files <- c(
    fcs_files,
    files
  )
  
}

cat(
  "Found",
  length(fcs_files),
  "AutoSpectral unmixed FCS files\n"
)
# create a folder where the results will be saved 
# Main output folder
output.dir <- file.path(
  root.dir,
  "AutoSpectral_Preprocessing"
)
dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# function to add  Original_ID to each cell/event to know which were kept or removed 

add_original_id <- function(ff) {
  # extract numeric event from original data 
  expr <- flowCore::exprs(ff)
  
  # add a new column original_ID to identify each event/cell
  expr <- cbind(expr, Original_ID = seq_len(nrow(expr)))
  
  # extract metadata describing each chanell 
  pd <- Biobase::pData(flowCore::parameters(ff))
  
  # copy first metadata row as template
  new_param <- pd[1, , drop = FALSE]
  
  # name the new parameter original_ID
  new_param$name <- "Original_ID"
  new_param$desc <- "Original_ID"
  
  # define the range of the original_ID collumn 
  new_param$range <- max(expr[, "Original_ID"])
  new_param$minRange <- 1
  new_param$maxRange <- max(expr[, "Original_ID"])
  # add the new metadata row and make sure the metadata names match the data columns.
  pd <- rbind(pd, new_param)
  rownames(pd) <- colnames(expr)
  # rebuild the flowFrame with the new Original_ID column.
  # description = list () used to avoid problem with fcs 3.2 metadata
  flowCore::flowFrame(
    exprs = expr,
    parameters = Biobase::AnnotatedDataFrame(pd),
    description = list()
  )
}
# function to select channels of interest 
pick_channel <- function(meta, preferred, fallback_pattern) {
  # check if the channel exist 
  hit <- preferred[preferred %in% meta$name]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  # prepare the description column for searching
  desc <- meta$desc
  desc[is.na(desc)] <- ""
  hits <- meta$name[
    grepl(fallback_pattern, meta$name, ignore.case = TRUE) |
      grepl(fallback_pattern, desc, ignore.case = TRUE)
  ]
  if (length(hits) == 0) {
    return(NA)
  }
  
  hits[1]
}

# General channel finder
# find other channels such as Time or Zombie NIR
find_channel <- function(meta, pattern) {
  # search both the channel name and description 
  desc <- meta$desc
  desc[is.na(desc)] <- ""
  
  hits <- meta$name[
    grepl(pattern, meta$name, ignore.case = TRUE) |
      grepl(pattern, desc, ignore.case = TRUE)
  ]
  
  if (length(hits) == 0) {
    return(NA)
  }
  
  hits[1]
}

# Main processing function
# process one fcs file 

cat("\nProcessing file:\n")
print(file)

process_fcs <- function(fcs_file) {
  # extract the sample name from the file name
  sample_name <- tools::file_path_sans_ext(basename(fcs_file))
  # extract the cytometer folder name
  cytometer <- basename(dirname( dirname(fcs_file)))
  # create a separate results folder for each sample
  sample_dir <- file.path(output.dir, cytometer, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Processing: ", sample_name)
  
  # Read FCS file
  
  fcs_data <- flowCore::read.FCS(
    fcs_file,
    transformation = FALSE,
    truncate_max_range = FALSE,
    emptyValue = FALSE
  )
  
  # Add Original_ID before any preprocessing
  fcs_data <- add_original_id(fcs_data)
  
  # extract Metadata
  # missing descriptions are replaced with empty text
  
  meta <- fcs_data@parameters@data
  meta$desc[is.na(meta$desc)] <- ""
  
  # Detect FSC, SSC, Time, and Zombie channels
  # find the forward scatter area channel 
  fsc_a <- pick_channel(
    meta,
    preferred = c("FSC-A", "FSC51-A", "FSC53-A"),
    fallback_pattern = "FSC.*-A"
  )
  # find the foward scater hight channel 
  fsc_h <- pick_channel(
    meta,
    preferred = c("FSC-H", "FSC51-H", "FSC53-H"),
    fallback_pattern = "FSC.*-H"
  )
  # find the side scatter area channel with the highest variance 
  ssc_candidates <- meta$name[
    grepl("SSC.*-A", meta$name, ignore.case=TRUE)
  ]
  
  if(length(ssc_candidates)==0)
    stop("No SSC channel found")
  
  vars <- sapply(
    ssc_candidates,
    function(ch)
      var(exprs(fcs_data)[, ch], na.rm=TRUE)
  )
  
  ssc_a <- ssc_candidates[which.max(vars)]
  # find the time and the viability channel 
  
  time_ch <- find_channel(meta, "^Time$|TIME|Time")
  zombie <- find_channel(meta, "Zombie|dead|LD|viability")
  
  # stop if any key channel is missing 
  
  if (any(is.na(c(fsc_a, fsc_h, ssc_a, time_ch, zombie)))) {
    stop(paste("Important channel missing in", sample_name))
  }
  # First try using the marker descriptions (raw files)
  # Detect marker channels automatically
  marker_words <- c(
    "CD3","CD4","CD8","CD11c","CD14","CD16","CD19",
    "CD25","CD27","CD38","CD39","CD45RA","CD56","CD57",
    "CD95","CD127","CCR7","CXCR5","HLA","TCR","PD",
    "TIGIT","BTLA"
  )
  
  # First try using marker descriptions
  markers <- meta$name[
    grepl(
      paste(marker_words, collapse="|"),
      meta$desc,
      ignore.case = TRUE
    )
  ]
  # If that finds very few markers,
  # assume this is an unmixed file
  if(length(markers) < 10){
    
    markers <- meta$name[
      grepl(
        "BUV|BV|BB|FITC|PE|PerCP|APC|AF647|AF700|R718|Spark|CF|Zombie",
        meta$name,
        ignore.case = TRUE
      )
    ]
    
  }
  # Remove channels that are not biological markers
    markers <- markers[
      !grepl(
        "AF|Autofluorescence|Residual|Error|error",
        markers,
        ignore.case = TRUE
      )
    ]
  
  # Remove duplicate Zombie channel
  markers <- setdiff(markers, zombie)
  # Add Zombie channel back once
  markers <- c(markers, zombie)
  # Keep unique channels
  markers <- unique(markers)
  # Show detected markers
  cat("\nDetected marker channels:\n")
  print(markers)
  
  cat("\nNumber of markers:", length(markers), "\n") 
  
  
  # Include live/dead marker for transformation 
  markers <- unique(c(markers, zombie))
  
  if (length(markers) == 0) {
    stop(paste("No marker channels found in", sample_name))
  }
  
  # Remove margins
  
  ff_m <- PeacoQC::RemoveMargins(
    fcs_data,
    channels = markers,
    remove_min = NULL,
  )
  
  # Logicle transformation
  
  trans <- tryCatch(
    estimateLogicle(ff_m, channels = markers),
    # if logicletransformation fails use arcsinh transformation 
    error = function(e) {
      message("Logicle failed for ", sample_name, ". Using arcsinh instead.")
      transformList(markers, arcsinhTransform(a = 0, b = 1/500, c = 0))
    }
  )
  #apply the transformation
  ff_t <- flowCore::transform(ff_m, trans)
  
  cat("\nAfter logicle transformation.\n")
  
  # 7. Linear transformation of SSC-A
  # using reference marker
  # select a reference marker
  reference_marker <- meta$name[
    grepl(
      "CD3",
      meta$desc,
      ignore.case = TRUE
    )
  ][1]
  # If CD3 is missing the script uses the first available marker
  if (is.na(reference_marker)) {
    reference_marker <- markers[1]
  }
  # Calculate the 5th and 95th percentile of the reference marker
  q5_goal <- quantile(exprs(ff_t)[, reference_marker], 0.05, na.rm = TRUE)
  q95_goal <- quantile(exprs(ff_t)[, reference_marker], 0.95, na.rm = TRUE)
  # calculate the 5th and 95th percentile of SSC-A.
  q5_ssc <- quantile(exprs(ff_t)[, ssc_a], 0.05, na.rm = TRUE)
  q95_ssc <- quantile(exprs(ff_t)[, ssc_a], 0.95, na.rm = TRUE)
  # calculate the slope
  ssc_a_slope <- (q95_goal - q5_goal) / (q95_ssc - q5_ssc)
  # calculate the intercept 
  ssc_a_intercept <- q5_goal - q5_ssc * ssc_a_slope
  # apply linear transformation 
  trans <- c(
    trans,
    transformList(
      ssc_a,
      flowCore::linearTransform(
        a = ssc_a_slope,
        b = ssc_a_intercept
      )
    )
  )
  ff_t <- flowCore::transform(ff_m, trans)
  # create density plots after transformation for all markers of interest
  
  df <- as.data.frame(exprs(ff_t)[, markers, drop = FALSE])
  # convert to long format for ploting 
  df_long <- stack(df)
  # use marker descriptions as plot labels 
  labels <- meta$desc[
    match(df_long$ind, meta$name)
  ]
  
  labels[is.na(labels)] <- ""
  
  # If descriptions are missing or all called "other",
  # use the channel names instead
  if(length(unique(labels)) <= 2 ||
     all(labels == "") ||
     all(labels == "other"))
  {
    labels <- df_long$ind
  }
  
  df_long$desc <- labels
  # create the density plot for all transformed markers 
  p_density <- ggplot(df_long, aes(x = values)) +
    geom_density(fill = "blue", alpha = 0.4) +
    facet_wrap(~ desc, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_text(size = 7))
  # save the density plot in the result folder 
  ggsave(
    file.path(sample_dir, "density_plots.png"),
    p_density,
    width = 16,
    height = 10
  )
  
  # Remove doublets
  
  ff_s <- PeacoQC::RemoveDoublets(
    ff_t,
    channel1 = fsc_a,
    channel2 = fsc_h
  )
  
  # Live/dead gating using cytoexploreR 
  
  message("Draw live/dead gate for: ", sample_name)
  
  gate_list <- cyto_gate_draw(
    ff_s,
    alias = "Live",
    channels = c(fsc_a, zombie),
    type = "polygon",
    display = 50000,
    axes_limits = "data"
  )
  # extract the drawn gate
  live_gate <- gate_list$Live
  # apply the gate and keep only live cells.
  live_filter <- flowCore::filter(ff_s, live_gate)
  ff_l <- ff_s[live_filter@subSet, ]
  
  # PeacoQC quality control on the live cells 
  
  PQC <- PeacoQC::PeacoQC(
    ff = ff_l,
    channels = markers,
    determine_good_cells = "all",
    plot = TRUE,
    save_fcs = TRUE,
    output_directory = sample_dir
  )
  # extract the final cleaned flowFrame
  ff_qc <- PQC$FinalFF
  cat("\nAfter PeacoQC:\n")
  cat("\nAfter PeacoQC:\n")
  print(summary(exprs(ff_qc)[, markers[1]]))
  # CompensAID analysis
  
  message("Running CompensAID for: ", sample_name)
  
  compensaid_dir <- file.path(
    sample_dir,
    "CompensAID_results"
  )
  
  dir.create(
    compensaid_dir,
    showWarnings = FALSE
  )
  
  compensaid_result <- tryCatch({
    
    ############################################################
    ## Prepare CompensAID input
    ############################################################
    
    meta_comp <- ff_qc@parameters@data
    
    meta_comp$desc[is.na(meta_comp$desc)] <- ""
    
    marker_channels <- meta_comp$name[
      !grepl(
        "FSC|SSC|Time|Original_ID",
        meta_comp$name,
        ignore.case = TRUE
      ) &
        !grepl(
          "FSC|SSC|Time|Original_ID",
          meta_comp$desc,
          ignore.case = TRUE
        )
    ]
    
    # Remove autofluorescence channel
    marker_channels <- marker_channels[
      marker_channels != "AF-A"
    ]
    
    # Create CompensAID flowFrame
    ff_markers <- ff_qc[, marker_channels]
    
    # Give every channel a unique marker name
    pd <- pData(parameters(ff_markers))
    
    pd$desc <- pd$name
    
    parameters(ff_markers) <- AnnotatedDataFrame(pd)
    
    cat("\nMarker names:\n")
    print(flowCore::markernames(ff_markers))
    
    cat("\nChannel names:\n")
    print(colnames(exprs(ff_markers)))
    
    cat("\nNumber of channels:\n")
    print(ncol(exprs(ff_markers)))
    
    # Run CompensAID
    res_compensaid <- CompensAID::CompensAID(
      ff = ff_markers,
    )
    # Save SSI matrix
    write.csv(
      res_compensaid$matrix,
      file.path(
        compensaid_dir,
        "CompensAID_matrix.csv"
      )
    )
    ssi_matrix_plot <- CompensAID::PlotMatrix(
      output = res_compensaid
    )
    
    ggsave(
      filename = file.path(
        compensaid_dir,
        "CompensAID_heatmap.png"
      ),
      plot = ssi_matrix_plot,
      width = 10,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    # Identify suspicious interactions
    # Using SSI < -1 as described in the paper
    
    flagged_index <- which(
      res_compensaid$matrix < -1,
      arr.ind = TRUE
    )
    
    if (nrow(flagged_index) > 0) {
      
      flagged_pairs <- data.frame(
        secondary_channel =
          rownames(res_compensaid$matrix)[
            flagged_index[, "row"]
          ],
        primary_channel =
          colnames(res_compensaid$matrix)[
            flagged_index[, "col"]
          ],
        SSI =
          res_compensaid$matrix[
            flagged_index
          ]
      )
      
      flagged_pairs <- unique(
        flagged_pairs
      )
      
      write.csv(
        flagged_pairs,
        file.path(
          compensaid_dir,
          "flagged_compensation_pairs.csv"
        ),
        row.names = FALSE
      )
      
      # Plot strongest problematic pair
      
      strongest <- flagged_pairs[
        which.min(flagged_pairs$SSI),
      ]
      
      df_pair <- data.frame(
        primary =
          exprs(ff_markers)[,
                            strongest$primary_channel],
        secondary =
          exprs(ff_markers)[,
                            strongest$secondary_channel]
      )
      
      set.seed(123)
      
      df_pair <- df_pair[
        sample(
          nrow(df_pair),
          min(
            10000,
            nrow(df_pair)
          )
        ),
      ]
      
      p_pair <- ggplot(
        df_pair,
        aes(
          x = primary,
          y = secondary
        )
      ) +
        geom_point(
          color = "blue",
          size = 0.4,
          alpha = 0.5
        ) +
        xlab(
          strongest$primary_channel
        ) +
        ylab(
          strongest$secondary_channel
        ) +
        ggtitle(
          paste(
            "Strongest flagged pair:",
            strongest$primary_channel,
            "vs",
            strongest$secondary_channel,
            "\nSSI =",
            round(
              strongest$SSI,
              2
            )
          )
        ) +
        theme_minimal()
      
      ggsave(
        file.path(
          compensaid_dir,
          "strongest_flagged_pair.png"
        ),
        p_pair,
        width = 8,
        height = 6
      )
      
    } else {
      
      writeLines(
        "No suspicious marker combinations detected.",
        file.path(
          compensaid_dir,
          "CompensAID_summary.txt"
        )
      )
      
    }
    
    "CompensAID completed"
    
  },
  error = function(e) {
    
    writeLines(
      paste(
        "CompensAID failed:",
        e$message
      ),
      file.path(
        compensaid_dir,
        "CompensAID_error.txt"
      )
    )
    
    "CompensAID failed"
    
  })
  # Create population plots
  # make a random selection of FSC-A channel for easy ploting  
  set.seed(123)
  fcs_data_filtered <- fcs_data[
    exprs(fcs_data)[, fsc_a] <
      quantile(exprs(fcs_data)[, fsc_a], 0.995, na.rm = TRUE),
  ]
  # create a filter ploting function 
  filter_plot <- function(ff_pre, ff_post, title, channel_x, channel_y) {
    
    df <- data.frame(
      x = exprs(ff_pre)[, channel_x],
      y = exprs(ff_pre)[, channel_y]
    )
    
    i <- sample(nrow(df), min(10000, nrow(df)))
    
    p <- ggplot(df[i, ], aes(x = x, y = y)) +
      geom_point(
        size = 0.5,
        color = ifelse(
          exprs(ff_pre)[i, "Original_ID"] %in% exprs(ff_post)[, "Original_ID"],
          "blue",
          "red"
        )
      ) +
      xlab(channel_x) +
      ylab(channel_y) +
      theme_minimal() +
      theme(legend.position = "none") +
      ggtitle(title)
    
    return(p)
  }
  # this shows margin removal 
  p1 <- filter_plot(
    fcs_data_filtered,
    ff_m,
    "FSC/SSC",
    fsc_a,
    ssc_a
  )
  # this shows doublet removal 
  p2 <- filter_plot(
    ff_t,
    ff_s,
    "Singlets",
    fsc_a,
    fsc_h
  )
  #this shows live/dead gating 
  p3 <- filter_plot(
    ff_s,
    ff_l,
    "LD",
    fsc_a,
    zombie
  )
  # this shows peacoQC filtering 
  p4 <- filter_plot(
    ff_l,
    ff_qc,
    "PeacoQC",
    time_ch,
    fsc_a
  )
  # combine all four plots as one row 
  combined <- ggpubr::ggarrange(
    p1, p2, p3, p4,
    nrow = 1
  )
  # save to result folder 
  ggsave(
    file.path(sample_dir, "population_plots.png"),
    combined,
    width = 20,
    height = 6
  )
  
  # save the final cleaned fcs file 
  cat("\nBefore saving cleaned file:\n")
  cat("\nAfter PeacoQC:\n")
  print(summary(exprs(ff_qc)[, markers[1]]))
  write.FCS(
    ff_qc,
    filename = file.path(sample_dir, paste0(sample_name, "_cleaned.fcs"))
  )
  
  # create and save summary table for the sample 
  summary <- data.frame(
    Sample = sample_name,
    Cytometer = cytometer,
    Raw_events = nrow(exprs(fcs_data)),
    After_margins = nrow(exprs(ff_m)),
    After_singlets = nrow(exprs(ff_s)),
    After_live_gate = nrow(exprs(ff_l)),
    CompensAID_status = compensaid_result,
    Final_after_PeacoQC = nrow(exprs(ff_qc)),
    FSC_A = fsc_a,
    FSC_H = fsc_h,
    SSC_A = ssc_a,
    Zombie = zombie,
    Reference_marker = reference_marker
  )
  
  write.csv(
    summary,
    file.path(sample_dir, "summary.csv"),
    row.names = FALSE
  )
  
  return(summary)
}

# run the pre processing function on all fcs files 
all_results <- list()

for (file in fcs_files) {
  all_results[[file]] <- process_fcs(file)
}
# combine all sample summaries into one table 
final_summary <- do.call(rbind, all_results)
# save the final summary table for all samples 
write.csv(
  final_summary,
  file.path(
    output.dir,
    "final_summary.csv"
  ),
  row.names = FALSE
)
# display the final summary
final_summary


