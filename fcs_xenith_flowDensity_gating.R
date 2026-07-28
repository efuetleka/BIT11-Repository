# install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("ggplot2")
BiocManager::install("ggpubr")
BiocManager::install("pheatmap")
BiocManager::install("tidyr")
BiocManager::install("FlowSOM")
BiocManager::install("ggcyto")
#load packages
library(flowCore)
library(flowDensity)
library(CytoExploreR)
library(FlowSOM)
library(ggplot2)
library(PeacoQC)
library(ggcyto)
#set the working directory
setwd("C:/Users/IT SOLLUTION Sarl/Documents/shared/FullStains_FM/Xenith_FM")
#read in the FCS file
fcs_file <- "C:/Users/IT SOLLUTION Sarl/Documents/shared/FullStains_FM/Xenith_FM/Sample H2 Run 1 20260129180729.fcs"
fcs_data <- read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
#check structure of the data
summary(fcs_data)
#check dimentions
dim(fcs_data)
#check column names 
colnames(exprs(fcs_data))
# check metadata of the columns 
fcs_data@parameters@data
# select a reference marker as a standard to align or rescale other channels
reference_marker <- "FL01-Comp"
# select markers of interest 
# Automatically select spectral marker channels
# Extract metadata
meta <- fcs_data@parameters@data

# Select channels based on their descriptions
markers_of_interest <- meta$name[
  grepl("^Spectral", meta$desc) &
    grepl("-Comp$", meta$name) &
    !grepl("Unstained|Scatter", meta$desc)
]
# Make sure Zombie NIR is included
zombie_channel <- meta$name[
  grepl("Zombie NIR", meta$desc, ignore.case = TRUE)
]

markers_of_interest <- unique(c(markers_of_interest, zombie_channel))
# Check selected channel names and descriptions
data.frame(
  name = markers_of_interest,
  desc = meta$desc[match(markers_of_interest, meta$name)]
)


# select channels of interest 
channels_of_interest <- GetChannels(object = fcs_data,
                                    markers = markers_of_interest,
                                    exact = FALSE)
# remove margins from the channels of interest in the dataset 
ff_m <- PeacoQC::RemoveMargins(fcs_data, channels_of_interest, remove_min = NULL)
# do a logicle transformation on all the markers of interest including the viability marker
#for easy visualization and interpretation of plots 
translist <- estimateLogicle(ff_m, markers_of_interest,)
ff_t <- flowCore::transform(ff_m, translist)
# do a density plot for all markers after tranformation 
# Extract transformed data
df <- as.data.frame(exprs(ff_t))

# Keep ONLY automatically selected marker channels
df <- df[, markers_of_interest]

# Convert to long format for ggplot faceting
df_long <- stack(df)
# Add descriptive labels
df_long$desc <- meta$desc[match(df_long$ind, meta$name)]
# Plot density plots for all markers
ggplot(df_long, aes(x = values)) +
  geom_density(fill = "blue", alpha = 0.4) +
  facet_wrap(~ desc, scales = "free") +
  theme_minimal()

#calculate the 5th and 95th percentiles of your reference marker
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05)
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95)
#calculate the 5th and 95th percentiles of your scatter channel
q5_FSCA <- quantile(exprs(ff_t)[,"FSC51-A"], 0.05)
q95_FSCA <- quantile(exprs(ff_t)[,"FSC51-A"], 0.95)
# calculate the slope 
FSCA_a <- (q95_goal- q5_goal) / (q95_FSCA- q5_FSCA)
#calculate the intercept
FSCA_b <- q5_goal- q5_FSCA * (q95_goal- q5_goal) / (q95_FSCA- q5_FSCA)
# create a linear transformation list
translist <- c(translist,
               transformList("FSC51-A", flowCore::linearTransform(a = FSCA_a,
                                                                  b =FSCA_b)))
# create a data frame 
df <- as.data.frame(exprs(fcs_data))
# remove doublet cells from the transformed data 
ff_s <- PeacoQC::RemoveDoublets(
  ff_t,
  channel1 = "FSC51-A",
  channel2 = "FSC51-H"
)

#set gates to filter live cells from dead cells 

# Find the Zombie NIR channel automatically
zombie_channel <- meta$name[
  grepl("Zombie NIR", meta$desc, ignore.case = TRUE)
]

# Automatically estimate live/dead cutoff
zombie_cutoff <- flowDensity::deGate(
  obj = ff_s,
  channel = zombie_channel
)

# Live cells have LOW Zombie NIR signal
live_cells <- exprs(ff_s)[, zombie_channel] < zombie_cutoff

# Keep only live cells
ff_l <- ff_s[live_cells, ]

# run quality control on filtered data and selected channels 
PQC <- PeacoQC::PeacoQC(ff = ff_l,
                        channels = channels_of_interest,
                        plot = TRUE, save_fcs = FALSE)
# randomise the selection and select a subset of the data before preprocessing
# to plot with the data after preprocessing 
set.seed(123)
fcs_data_filtered <- fcs_data[exprs(fcs_data)[, "FSC51-A"] < 2e6, ]
# create a plot function 
filter_plot <- function(ff_pre, ff_post, title, channel_x, channel_y){
  df <- data.frame(x = exprs(ff_pre)[,channel_x],
                   y =exprs(ff_pre)[,channel_y])
  i <- sample(nrow(df), min(10000, nrow(df)))
  if (!"Original_ID" %in% colnames(exprs(ff_pre))) {
    ff_pre@exprs <- cbind(ff_pre@exprs,
                          Original_ID = seq_len(nrow(ff_pre@exprs)))
  }
  p <- ggplot(df[i,], aes(x = x, y = y)) +
    geom_point(size = 0.5,
               color = ifelse(exprs(ff_pre)[i,"Original_ID"] %in%
                                exprs(ff_post)[,"Original_ID"], "blue", "red")) +
    xlab(GetMarkers(ff_pre, channel_x)) +
    ylab(GetMarkers(ff_pre, channel_y)) +
    theme_minimal() + theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}
# create a list of all the plot using the plot function 
to_plot <- list(list(ff_pre = fcs_data_filtered,
                     ff_post = ff_m,
                     title = "FSC51-A FSC51-H",
                     channel_x = "FSC51-A",
                     channel_y = "SSC52-A"),
                list(ff_pre = ff_t,
                     ff_post = ff_s,
                     title = "Singlets",
                     channel_x = "FSC51-A",
                     channel_y = "FSC51-H"),
                list(ff_pre = ff_s,
                     ff_post = ff_l,
                     title = "LD",
                     channel_x = "FSC51-A",
                     channel_y = "FL22-Comp"),
                list(ff_pre = ff_l,
                     ff_post = PQC$FinalFF,
                     title = "PeacoQC",
                     channel_x = "Time",
                     channel_y = "FSC51-A"))
plot_list <- list()
for (plot in to_plot) {
  plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre,
                                                    ff_post = plot$ff_post,
                                                    title = plot$title,
                                                    channel_x = plot$channel_x,
                                                    channel_y = plot$channel_y)
}
# prints all plots on same row 
print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
# use CytoExploreR 
# have a look at flowDensity
# do a density plot for the markers after transformation
