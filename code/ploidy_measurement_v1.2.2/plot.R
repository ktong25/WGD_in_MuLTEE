library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20230630/20231225-1_Donut_ploidy/")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Themes
add_ridgeplot_theme <- function(p) {
  p + 
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none", 
          panel.grid.major = element_line(size = 0.7, colour = "grey80"),
          panel.grid.minor = element_line(size = 0.35, colour = "grey90"))
}

# Parameters
in_dir <- "results_csvs"
out_fig_path <- "."
control_strains <- c("2N", "4N")  # first 2N then 4N

# Load data
# Input file name: e.g. "<sample>_001_Results.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1)
  measure_vars <- colnames(df)
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  df$Sample <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Strain <- paste(metadata[-length(metadata)], collapse = "_")
  df$ImageID <- as.integer(metadata[length(metadata)])
  df$Strain_ImageID <- paste(df$Strain, df$ImageID, sep = "_")
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = in_dir, pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv)
summary(data)

# Process data
strains <- unique(data$Strain)
data$Strain <- factor(data$Strain, levels = strains)
data <- data %>%
  dplyr::arrange(factor(Strain, levels = strains), ImageID)
data$Strain_ImageID <- factor(data$Strain_ImageID, levels = unique(data$Strain_ImageID))
summary(data)
id_vars <- c("Strain", "ImageID", "Strain_ImageID")

# Plot nuclei intensity
ggplot(data, aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) %>% add_ridgeplot_theme()
ggplot(data, aes(x = RawIntDen_bgs, y = Strain_ImageID, fill = Strain)) %>% add_ridgeplot_theme() +
  theme(axis.text.y = element_text(size = 5))

# Plot nuclei intensity calculated using non-cell background subtraction
# ggplot(data, aes(x = RawIntDen_bgs2, y = Strain, fill = Strain)) %>%
#   add_ridgeplot_theme()
# ggplot(data, aes(x = RawIntDen_bgs2, y = Strain_ImageID, fill = Strain)) %>%
#   add_ridgeplot_theme()

##### Remove small objects (segmentation artifacts) #####

# Plot nuclei area
ggplot(data, aes(x = Area, y = Strain, fill = Strain)) %>% add_ridgeplot_theme()
ggplot(data, aes(x = Area, y = Strain_ImageID, fill = Strain)) %>% add_ridgeplot_theme()

# Determine nuclei area thresholds
# Mean and SD
# ggplot(data, aes(x = Area, y = Strain, fill = Strain)) %>%
#   add_ridgeplot_theme() +
#   geom_segment(data = data %>% 
#                  dplyr::group_by(Strain) %>%
#                  dplyr::summarise(Cutoff = mean(Area) - sd(Area) * 2) %>%
#                  dplyr::mutate(y = 1:length(strains)), 
#              mapping = aes(x = Cutoff, y = y, xend = Cutoff, yend = y+0.9), 
#              color = "black", size = 1)
# Median and MAD (usually work better than mean and SD)
ggplot(data, aes(x = Area, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme() +
  geom_segment(data = data %>% 
                 dplyr::group_by(Strain) %>%
                 dplyr::summarise(Cutoff = median(Area) - mad(Area) * 2) %>%  ###
                 dplyr::mutate(y = 1:length(strains)), 
               mapping = aes(x = Cutoff, y = y, xend = Cutoff, yend = y+0.9), 
               color = "black", size = 1)

# Filter data by nuclei area thresholds to remove little left hump
data_filter <- data %>%
  #dplyr::filter(Area > 2) #%>%  # this allows keeping 2N cells (minority) in all populations
  dplyr::group_by(Strain) %>%
  dplyr::filter(Area > median(Area) - mad(Area) * 2) %>%
  dplyr::filter(!( (!Strain %in% "2N") & Area < 2.6 ))
# dplyr::filter(!( (Strain %in% c("PA_1_t1000_S_rep2", "PA_1_t1000_S_rep3") & Area < 3) | 
#                    (Strain %in% c("PA_1_t1000_S_rep1", "PA_3_t1000_S_rep3") & Area < 2.5) #|
#                  # # #                  (Strain == "PA_4_t400_lacol" & Area < 2.5)
# ))

# Plot nuclei area
ggplot(data_filter, aes(x = Area, y = Strain, fill = Strain)) %>% add_ridgeplot_theme()
ggplot(data_filter, aes(x = Area, y = Strain_ImageID, fill = Strain)) %>% add_ridgeplot_theme()

# Plot nuclei intensity
ggplot(data_filter, aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) %>% add_ridgeplot_theme()
ggplot(data_filter, aes(x = RawIntDen_bgs, y = Strain_ImageID, fill = Strain)) %>% add_ridgeplot_theme() +
  theme(axis.text.y = element_text(size = 5))
save_ggplot("all_images_nuclei_intensity_filtered_1_nuclei_area", width = 10, height = length(unique(data_filter$Strain_ImageID)) * 0.08, format = "pdf")

##### Remove outlier images #####

# Remove obviously outlier images
data_filter2 <- data_filter %>% 
  dplyr::filter(!Strain_ImageID %in% c(paste0("PA_4_t1000_D_rep3_", c(5,7,8)),
                                       paste0("PA_4_t600_D_rep4_", c(6,8))
  ))

# Plot nuclei intensity
ggplot(data_filter2, aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) %>% add_ridgeplot_theme()
ggplot(data_filter2, aes(x = RawIntDen_bgs, y = Strain_ImageID, fill = Strain)) %>% add_ridgeplot_theme() +
  theme(axis.text.y = element_text(size = 5))
save_ggplot("all_images_nuclei_intensity_filtered_2_images", width = 10, height = length(unique(data_filter2$Strain_ImageID)) * 0.08, format = "pdf")

# Assign to data_final
data_final <- data_filter2

##### Data summary #####

data_summary <- data_final %>%
  dplyr::group_by(Strain) %>%
  dplyr::summarize(Count = n())

##### Normalize nuclear intensity by ploidy controls #####

# Get G1 peak of ploidy control strains
control <- data_final %>%
  dplyr::filter(Strain %in% control_strains) %>%
  dplyr::arrange(factor(Strain, levels = control_strains))
control_g1peak <- control %>%
  dplyr::group_by(Strain) %>%
  dplyr::summarise(G1Peak = density(RawIntDen_bgs)$x[which.max(density(RawIntDen_bgs)$y)]) # G1 peak is usually the highest peak
# Plot to check
ggplot(control, aes(x = RawIntDen_bgs, y = Strain, fill = Strain, color = Strain)) %>%
  add_ridgeplot_theme() +
  geom_segment(data = control_g1peak %>%
                 dplyr::mutate(y = 1:length(control_strains)), 
               mapping = aes(x = G1Peak, y = y, xend = G1Peak, yend = y+0.9), color = "black", size = 1)

# Calculate fluorescent intensity of 1N genome
intensity_1N <- (control_g1peak$G1Peak[1] + control_g1peak$G1Peak[2] / 2) / 2 / 2
# Normalize all fluorescent intensity values by it
data_final$DNA_content <- data_final$RawIntDen_bgs / intensity_1N

##### Plot with better formatting #####

# Specify strain order and label
strains <- strains
data_final$Strain <- factor(data_final$Strain, levels = strains)
strains_label <- strains
# strains_label <- strains %>%
#   # gsub("PM_", "PM", ., fixed = TRUE) %>%
#   # gsub("PA_", "PA", ., fixed = TRUE) %>%
#   # gsub("PO_", "PO", ., fixed = TRUE) %>%
#   gsub("_t", " t", ., fixed = TRUE)

# Plot nuclei intensity
ggplot(data_final %>% dplyr::mutate(Strain = factor(Strain, levels = strains, labels = strains_label)), 
       aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) +
  geom_density_ridges(scale = 2) +
  scale_x_continuous(limits = c(0, 1.5e6), 
                     breaks = seq(0, 1.5e6, 0.5e6),
                     minor_breaks = seq(0, 1.5e6, 0.1e6), 
                     labels = c("0", paste0(seq(0.5, 1.5, 0.5), "e6"))) +
  xlab("Total PI intensity of nucleus (a.u.)") +
  geom_vline(xintercept = control_g1peak$G1Peak[1], size = 0.8, linetype = "dashed", color = "#1B9E77") +  # 2Nm G1 peak
  geom_vline(xintercept = control_g1peak$G1Peak[2], size = 0.8, linetype = "dashed", color = "#D95F02") +  # 4Nm G1 peak  
  ggplot_custom_theme +
  theme(legend.position = "none", 
        panel.grid.major = element_line(size = 0.6, colour = "grey80"),
        panel.grid.minor = element_line(size = 0.3, colour = "grey90"), 
        axis.title.x = element_text(vjust = 0), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(vjust = 0), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank())
save_ggplot("nuclear_intensity", width = 6, height = 0.35 * length(strains) + 1.5)

# Plot normalized nuclei intensity (DNA content)
ggplot(data_final %>% dplyr::mutate(Strain = factor(Strain, levels = strains, labels = strains_label)), 
       aes(x = DNA_content, y = Strain, fill = Strain)) +
  geom_density_ridges(scale = 2) +
  scale_x_continuous(limits = c(0, 11), 
                     breaks = seq(0, 10, 2),
                     minor_breaks = seq(0, 10, 1), 
                     labels = c("0", paste0(seq(2, 10, 2), "C"))) +
  xlab("DNA content") +
  ggplot_custom_theme +
  theme(legend.position = "none", 
        panel.grid.major = element_line(size = 0.6, colour = "grey80"),
        panel.grid.minor = element_line(size = 0.3, colour = "grey90"), 
        axis.title.x = element_text(vjust = 0), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(vjust = 0), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank())
save_ggplot("DNA_content", width = 6, height = 0.35 * length(strains) + 1.5)

# Save data
write.csv(data, file = "data_raw.csv", row.names = FALSE)
write.csv(data_final, file = "data_final.csv", row.names = FALSE)
# Save data summary
write.csv(data_summary, file = "data_summary_final.csv", row.names = FALSE)
