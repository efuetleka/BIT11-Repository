
library(flowCore)
library(CytoExploreR)
library(PeacoQC)
library(ggplot2)
library(ggpubr)
library(Biobase)

# Set working directory

setwd("~/FullStains_FM")

# Find FCS files
#This searches for every .fcs file inside the working directory and all subfolders
fcs_files <- list.files(
  pattern = "\\.fcs$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# Exclude the singles subfolder and previous output result
fcs_files <- fcs_files[
  !grepl("Singles|Preprocessing_results|_cleaned", fcs_files, ignore.case = TRUE)
]
# create a main folder where all results will be saved 
dir.create("Preprocessing_results", showWarnings = FALSE)

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

find_channel <- function(meta, pattern) {
  
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


# =========================
# Your filter plot function
# =========================



# =========================
# Main processing function
# =========================

process_fcs <- function(fcs_file) {
  
  sample_name <- tools::file_path_sans_ext(basename(fcs_file))
  cytometer <- basename(dirname(fcs_file))
  
  sample_dir <- file.path("Preprocessing_results", cytometer, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Processing: ", sample_name)
  
  
  # =========================
  # 1. Read FCS file
  # =========================
  
  fcs_data <- read.FCS(
    fcs_file,
    transformation = FALSE,
    truncate_max_range = FALSE,
    emptyValue = FALSE
  )
  
  # Add Original_ID before any preprocessing
  fcs_data <- add_original_id(fcs_data)
  
  
  # =========================
  # 2. Metadata
  # =========================
  
  meta <- fcs_data@parameters@data
  meta$desc[is.na(meta$desc)] <- ""
  
  
  # =========================
  # 3. Detect channels
  # =========================
  
  fsc_a <- pick_channel(
    meta,
    preferred = c("FSC-A", "FSC51-A", "FSC53-A"),
    fallback_pattern = "FSC.*-A"
  )
  
  fsc_h <- pick_channel(
    meta,
    preferred = c("FSC-H", "FSC51-H", "FSC53-H"),
    fallback_pattern = "FSC.*-H"
  )
  
  ssc_a <- pick_channel(
    meta,
    preferred = c("SSC-A", "SSC52-A", "VSSC-A", "BSSC-A"),
    fallback_pattern = "SSC.*-A"
  )
  
  time_ch <- find_channel(meta, "^Time$|TIME|Time")
  zombie <- find_channel(meta, "Zombie|dead|LD|viability")
  
  if (any(is.na(c(fsc_a, fsc_h, ssc_a, time_ch, zombie)))) {
    stop(paste("Important channel missing in", sample_name))
  }
  
  
  # =========================
  # 4. Detect markers of interest
  # =========================
  
  marker_words <- c(
    "CD3", "CD4", "CD8", "CD11c", "CD14", "CD16", "CD19",
    "CD25", "CD27", "CD38", "CD39", "CD45RA", "CD56", "CD57",
    "CD95", "CD127", "CCR7", "CXCR5", "HLA", "TCR", "PD",
    "TIGIT", "BTLA"
  )
  
  markers <- meta$name[
    grepl(paste(marker_words, collapse = "|"), meta$desc, ignore.case = TRUE) &
      grepl("-A$|Comp$", meta$name, ignore.case = TRUE)
  ]
  
  # Include live/dead marker
  markers <- unique(c(markers, zombie))
  
  if (length(markers) == 0) {
    stop(paste("No marker channels found in", sample_name))
  }
  
  
  # =========================
  # 5. Remove margins
  # =========================
  
  ff_m <- PeacoQC::RemoveMargins(
    fcs_data,
    channels = markers,
    remove_min = FALSE,
    remove_max = TRUE
  )
  
  
  # =========================
  # 6. Logicle transformation
  # =========================
  
  trans <- tryCatch(
    estimateLogicle(ff_m, channels = markers),
    error = function(e) {
      message("Logicle failed for ", sample_name, ". Using arcsinh instead.")
      transformList(markers, arcsinhTransform(a = 0, b = 1/500, c = 0))
    }
  )
  
  ff_t <- flowCore::transform(ff_m, trans)
  
  
  # =========================
  # 7. Linear transformation of SSC-A
  # using reference marker
  # =========================
  
  reference_marker <- meta$name[
    grepl("CD3", meta$desc, ignore.case = TRUE) &
      grepl("-A$|Comp$", meta$name, ignore.case = TRUE)
  ][1]
  
  if (is.na(reference_marker)) {
    reference_marker <- markers[1]
  }
  
  q5_goal <- quantile(exprs(ff_t)[, reference_marker], 0.05, na.rm = TRUE)
  q95_goal <- quantile(exprs(ff_t)[, reference_marker], 0.95, na.rm = TRUE)
  
  q5_ssc <- quantile(exprs(ff_t)[, ssc_a], 0.05, na.rm = TRUE)
  q95_ssc <- quantile(exprs(ff_t)[, ssc_a], 0.95, na.rm = TRUE)
  
  ssc_a_slope <- (q95_goal - q5_goal) / (q95_ssc - q5_ssc)
  ssc_a_intercept <- q5_goal - q5_ssc * ssc_a_slope
  
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
  
 
  
  
  # =========================
  # 8. Density plots after transformation
  # =========================
  
  df <- as.data.frame(exprs(ff_t)[, markers, drop = FALSE])
  df_long <- stack(df)
  df_long$desc <- meta$desc[match(df_long$ind, meta$name)]
  df_long$desc <- gsub("Spectral ", "", df_long$desc)
  
  p_density <- ggplot(df_long, aes(x = values)) +
    geom_density(fill = "blue", alpha = 0.4) +
    facet_wrap(~ desc, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_text(size = 7))
  
  ggsave(
    file.path(sample_dir, "density_plots.png"),
    p_density,
    width = 16,
    height = 10
  )
  
  
  # =========================
  # 9. Remove doublets
  # =========================
  
  ff_s <- PeacoQC::RemoveDoublets(
    ff_t,
    channel1 = fsc_a,
    channel2 = fsc_h
  )
  
  
  # =========================
  # 10. Live/dead gating
  # =========================
  
  message("Draw live/dead gate for: ", sample_name)
  
  gate_list <- cyto_gate_draw(
    ff_s,
    alias = "Live",
    channels = c(fsc_a, zombie),
    type = "polygon",
    display = 50000,
    axes_limits = "data"
  )
  
  live_gate <- gate_list$Live
  live_filter <- flowCore::filter(ff_s, live_gate)
  ff_l <- ff_s[live_filter@subSet, ]
  
  
  # =========================
  # 11. PeacoQC
  
  PQC <- PeacoQC::PeacoQC(
    ff = ff_l,
    channels = markers,
    determine_good_cells = "all",
    plot = TRUE,
    save_fcs = FALSE,
    output_directory = sample_dir
  )
  
  ff_qc <- PQC$FinalFF
  
  # 12. Create population plots

  
  set.seed(123)
  
  fcs_data_filtered <- fcs_data[
    exprs(fcs_data)[, fsc_a] <
      quantile(exprs(fcs_data)[, fsc_a], 0.995, na.rm = TRUE),
  ]
  # create function for filter plot 
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
  
  p1 <- filter_plot(
    fcs_data_filtered,
    ff_m,
    "FSC/SSC",
    fsc_a,
    ssc_a
  )
  
  p2 <- filter_plot(
    ff_t,
    ff_s,
    "Singlets",
    fsc_a,
    fsc_h
  )
  
  p3 <- filter_plot(
    ff_s,
    ff_l,
    "LD",
    fsc_a,
    zombie
  )
  
  p4 <- filter_plot(
    ff_l,
    ff_qc,
    "PeacoQC",
    time_ch,
    fsc_a
  )
  
  combined <- ggpubr::ggarrange(
    p1, p2, p3, p4,
    nrow = 1
  )
  
  ggsave(
    file.path(sample_dir, "population_plots.png"),
    combined,
    width = 20,
    height = 6
  )
  
  # 13. Save cleaned FCS
  
  write.FCS(
    ff_qc,
    filename = file.path(sample_dir, paste0(sample_name, "_cleaned.fcs"))
  )
  
  # 14. Save summary
  
  summary <- data.frame(
    Sample = sample_name,
    Cytometer = cytometer,
    Raw_events = nrow(exprs(fcs_data)),
    After_margins = nrow(exprs(ff_m)),
    After_singlets = nrow(exprs(ff_s)),
    After_live_gate = nrow(exprs(ff_l)),
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

# Run all FCS files

all_results <- list()

for (file in fcs_files) {
  all_results[[file]] <- process_fcs(file)
}

final_summary <- do.call(rbind, all_results)

write.csv(
  final_summary,
  "Preprocessing_results/final_summary.csv",
  row.names = FALSE
)

final_summary
