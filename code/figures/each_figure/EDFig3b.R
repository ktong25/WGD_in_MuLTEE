library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)
library(scales)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20220923_trial")
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
strains <- c("2Nm", "4Nm", "4Nm2", "t0", paste0("PA", 1:5, "_t1000"))
out_fig_path <- "."

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
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Strain, levels = strains), ImageID)
summary(data)
data$Strain <- factor(data$Strain, levels = strains)
data$ImageID <- factor(data$ImageID, levels = sort(unique(data$ImageID)))
data$Strain_ImageID <- factor(data$Strain_ImageID, levels = unique(data$Strain_ImageID))
summary(data)
id_vars <- c("Strain", "ImageID", "Strain_ImageID")

# Plot nuclei intensity
ggplot(data, aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data, aes(x = RawIntDen_bgs, y = Strain_ImageID, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data, aes(x = RawIntDen_bgs2, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data, aes(x = RawIntDen_bgs2, y = Strain_ImageID, fill = Strain)) %>%
  add_ridgeplot_theme()

# Plot nuclei area
ggplot(data, aes(x = Area, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data, aes(x = Area, y = Strain_ImageID, fill = Strain)) %>%
  add_ridgeplot_theme()

# Determine nuclei area thresholds
for (strain in strains) {
  areas <- data %>% dplyr::filter(Strain == strain) %>% dplyr::pull(Area)
  print(glue("{strain}: {mean(areas) - sd(areas) * 2}"))
}
# 2Nm: 1.44461095046166
# 4Nm: 1.99815800303165
# 4Nm2: 1.69869808249695
# t0: 1.32697519423555
# PA1_t1000: 1.64752013870017
# PA2_t1000: 1.82140082295053
# PA3_t1000: 1.22579155358727
# PA4_t1000: 0.836874084102952
# PA5_t1000: 0.708700641898632
# PA4/5-t1000 have quite many small objects that distorts normal distribution
for (strain in strains) {
  areas <- data %>% dplyr::filter(Strain == strain) %>% dplyr::pull(Area)
  print(glue("{strain}: {median(areas) - mad(areas) * 2}"))
}
# 2Nm: 1.4120055256
# 4Nm: 1.9884425072
# 4Nm2: 1.8370363068
# t0: 1.33199418
# PA1_t1000: 2.0161584888
# PA2_t1000: 2.3050135072
# PA3_t1000: 1.777392542
# PA4_t1000: 1.6377403968
# PA5_t1000: 1.5591959592
# median and MAD makes more sense

# Filter data by nuclei area thresholds to remove little left hump
data_filter <- data %>%
  dplyr::group_by(Strain) %>%
  dplyr::filter(Area > median(Area) - mad(Area) * 2)

# Plot nuclei area
ggplot(data_filter, aes(x = Area, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data_filter, aes(x = Area, y = Strain_ImageID, fill = Strain)) %>%
  add_ridgeplot_theme()

# Plot nuclei intensity
ggplot(data_filter, aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data_filter, aes(x = RawIntDen_bgs, y = Strain_ImageID, fill = Strain)) %>%
  add_ridgeplot_theme()

# Remove obviously outlier images
data_filter2 <- data_filter

# Plot nuclei intensity
ggplot(data_filter2, aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) %>%
  add_ridgeplot_theme()
ggplot(data_filter2, aes(x = RawIntDen_bgs, y = Strain_ImageID, fill = Strain)) %>%
  add_ridgeplot_theme()

# Plot nuclei intensity
# Remove 4Nm and convert 4Nm2 to 4Nm
data_final <- data_filter2
data_final$Strain <- as.character(data_final$Strain)
# data_final <- data_final %>%
#   dplyr::filter(Strain != "4Nm") %>%
#   dplyr::mutate(Strain = ifelse(Strain == "4Nm2", "4Nm", Strain))
# Strain label
strains <- c("2Nm", "4Nm", "t0", paste0("PA", 1:5, "_t1000"))
data_final$Strain <- factor(data_final$Strain, levels = strains)
strains_label <- c("2N multi", "4N multi", "PA t0", paste0("PA", 1:5, " t1000"))
# Plot
ggplot(data_final, # %>% dplyr::filter(Strain %in% c("2Nm", "4Nm")), 
       aes(x = RawIntDen_bgs, 
           y = factor(Strain, levels = strains, labels = strains_label), 
           fill = factor(Strain, levels = strains, labels = strains_label))) +
  geom_density_ridges(scale = 2) +
  scale_x_continuous(limits = c(0, 1.5e6), 
                     breaks = seq(0, 1.5e6, 0.5e6),
                     minor_breaks = seq(0, 1.5e6, 0.1e6), 
                     labels = c("0", paste0(seq(0.5, 1.5, 0.5), "e6"))) + #function(x) format(x, scientific = TRUE)) +
  xlab("Total PI intensity of nucleus (a.u.)") +
  geom_vline(xintercept = 0.238e6, size = 0.8, linetype = "dashed", color = "#1B9E77") +  # 2Nm G1 peak
  geom_vline(xintercept = 0.47e6, size = 0.8, linetype = "dashed", color = "#D95F02") +  # 4Nm G1 peak  
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
save_ggplot("ploidy", width = 6, height = 5)
# Plot ploidy controls
ggplot(data_final %>%
         dplyr::filter(Strain %in% c("2Nm", "4Nm")) %>%
         dplyr::mutate(Strain = factor(Strain, levels = c("4Nm", "2Nm"), labels = c("4N", "2N"))), 
       aes(x = RawIntDen_bgs, y = Strain, fill = Strain)) +
  geom_density_ridges(scale = 2, size = 0.5) +
  scale_fill_manual(values = c("4N" = "#D95F02", "2N" = "#1B9E77") %>% colorspace::lighten(0.4)) + 
  scale_x_continuous(limits = c(0, 1.5e6), 
                     breaks = seq(0, 1.5e6, 0.5e6),
                     minor_breaks = seq(0, 1.5e6, 0.1e6), 
                     labels = c("0", paste0(seq(0.5, 1.5, 0.5), "e6"))) + #function(x) format(x, scientific = TRUE)) +
  xlab("Total PI intensity of nucleus (a.u.)") +
  geom_vline(xintercept = 0.235e6, size = 0.75, linetype = "dashed", color = "black") +  # 2Nm G1 peak
  geom_vline(xintercept = 0.468e6, size = 0.75, linetype = "dashed", color = "black") +  # 4Nm G1 peak  
  ggplot_custom_theme +
  theme(legend.position = "none", 
        panel.grid.major = element_line(size = 0.5, colour = "grey80"),
        panel.grid.minor = element_line(size = 0.25, colour = "grey90"), 
        axis.title.x = element_text(vjust = 0), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(vjust = 0, size = rel(1.25), color = "black"), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank())
save_ggplot("ploidy_controls", width = 6, height = 4)
save_ggplot("ploidy_controls", width = 6, height = 4, mode = "paper")

# Data summary
data_summary <- data_final %>%
  dplyr::group_by(Strain) %>%
  dplyr::summarize(Count = n())
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
data_summary_controls <- data_final %>%
  dplyr::filter(Strain %in% c("2Nm", "4Nm")) %>%
  dplyr::group_by(Strain) %>%
  dplyr::summarize(Count = n())
write.csv(data_summary_controls, file = "data_summary_controls.csv", row.names = FALSE)

# Save data
write.csv(data_final, file = "data.csv", row.names = FALSE)
write.csv(data_final %>% dplyr::filter(Strain %in% c("2Nm", "4Nm")), 
          file = "data_controls.csv", row.names = FALSE)

