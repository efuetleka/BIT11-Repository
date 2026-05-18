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
setwd("C:/Users/IT SOLLUTION Sarl/Documents/shared/FullStains_FM/AuroraEvo_FM")
#read in the FCS file
fcs_file <- "C:/Users/IT SOLLUTION Sarl/Documents/shared/FullStains_FM/AuroraEvo_FM/H1 Well_001_plate FM.fcs"
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
reference_marker <- "BUV496-A"
# select markers of interest 
markers_of_interest <- c("BUV395-A", "BUV496-A", "BUV563-A", "BUV661-A", "BUV737-A", "BUV805-A", 
                         "BV421-A", "BV570-A", "BV605-A", "BV650-A", "BV711-A",
                         "BV750-A", "BV785-A", "FITC-A", "RB545-A", "RB613-A", "RB705-A", "PE-A", 
                         "PE-Cy5-A", "PE-Cy7-A", "APC-A", "R718-A", "APC-Fire 810-A", "AF-A")



# create a polygon gate to separate live cells from dead cells 
# test several gates 
live_gate1 <- flowCore::polygonGate(
  filterId = "Live_gate_1",
  .gate = matrix(
    data = c(
      400000, 800000, 1600000, 2200000, 2200000, 400000,
      0.2,    0.1,    0.2,     0.8,     1.8,     1.6
    ),
    ncol = 2,
    dimnames = list(c(), c("FSC-A", "Zombie NIR-A"))
  )
)
live_gate2 <- flowCore::polygonGate(
  filterId = "Live",
  .gate = matrix(
    data = c(
      400000, 800000, 1600000, 2200000, 2200000, 400000,
      0.2,    0.1,    0.2,     0.7,     1.5,     1.2
    ),
    ncol = 2,
    dimnames = list(NULL, c("FSC-A", "Zombie NIR-A"))
  )
)

# select channels of interest 
channels_of_interest <- GetChannels(object = fcs_data,
                                    markers = markers_of_interest,
                                    exact = FALSE)
# remove margins from the channels of interest in the dataset 
ff_m <- PeacoQC::RemoveMargins(fcs_data, channels_of_interest, remove_min = NULL)
# do a logicle transformation on all the markers of interest including the viability marker
#for easy visualization and interpretation of plots 
translist <- estimateLogicle(ff_m, c("BUV395-A", "BUV496-A", "BUV563-A", "BUV661-A", "BUV737-A", "BUV805-A", 
                                     "BV421-A", "BV570-A", "BV605-A", "BV650-A", "BV711-A",
                                     "BV750-A", "BV785-A", "FITC-A", "RB545-A", "RB613-A", "RB705-A", "PE-A", 
                                     "PE-Cy5-A", "PE-Cy7-A", "APC-A", "R718-A", "Zombie NIR-A", "APC-Fire 810-A", "AF-A"))
ff_t <- flowCore::transform(ff_m, translist)
# do a density plot for all markers after tranformation 
# Extract data
df <- as.data.frame(exprs(ff_t))
# Remove non-fluorescence channels
exclude <- c("FSC-A","FSC-H","SSC-A","SSC-H","Time","Original_ID")
df <- df[, !(colnames(df) %in% exclude)]
# Convert to long format (base R)
df_long <- stack(df)
# Plot all densities in one figure
ggplot(df_long, aes(x = values)) +
  geom_density(fill = "blue", alpha = 0.4) +
  facet_wrap(~ ind, scales = "free") +
  theme_minimal()

#calculate the 5th and 95th percentiles of your reference marker
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05)
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95)
#calculate the 5th and 95th percentiles of your scatter channel
q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
# calculate the slope 
SSCA_a <- (q95_goal- q5_goal) / (q95_SSCA- q5_SSCA)
#calculate the intercept
SSCA_b <- q5_goal- q5_SSCA * (q95_goal- q5_goal) / (q95_SSCA- q5_SSCA)
# create a linear transformation list
translist <- c(translist,
               transformList("SSC-A", flowCore::linearTransform(a = SSCA_a,
                                                                b =SSCA_b)))
# create a data frame 
df <- as.data.frame(exprs(fcs_data))
# remove double cells
ff_s <- PeacoQC::RemoveDoublets(ff_t)

#using CytoexploreR to manually draw the gates
live_gate_list <- cyto_gate_draw(
  ff_s,
  alias = "Live",
  channels = c("FSC-A", "Zombie NIR-A"),
  type = "polygon",
  display = 50000,
  axes_limits = "data"
)
live_gate3 <- live_gate_list$Live
# class needs to be polygonegate
class(live_gate3)
#apply the polygone gate to your data after doublet removal
selected_live <- filter(ff_s, live_gate3)
# keep only the live cells 
ff_l <- ff_s[selected_live@subSet,]
# run quality control on filtered data and selected channels 
PQC <- PeacoQC::PeacoQC(ff = ff_l,
                        channels = channels_of_interest,
                        plot = TRUE, save_fcs = FALSE)
# randomise the selection and select a subset of the data before preprocessing
# to plot with the data after preprocessing 
set.seed(123)
fcs_data_filtered <- fcs_data[exprs(fcs_data)[, "FSC-A"] < 2e6, ]
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
                     title = "FSC-A SSC-A",
                     channel_x = "FSC-A",
                     channel_y = "SSC-A"),
                list(ff_pre = ff_t,
                     ff_post = ff_s,
                     title = "Singlets",
                     channel_x = "FSC-A",
                     channel_y = "FSC-H"),
                list(ff_pre = ff_s,
                     ff_post = ff_l,
                     title = "LD",
                     channel_x = "FSC-A",
                     channel_y = "Zombie NIR-A"),
                list(ff_pre = ff_l,
                     ff_post = PQC$FinalFF,
                     title = "PeacoQC",
                     channel_x = "Time",
                     channel_y = "FSC-A"))
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
# explore autospectral unmixing tool
# explore SpectralunmixR tool 
# have a look at the compensAID library 
