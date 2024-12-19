library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggforce)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)
library(introdataviz)  # devtools::install_github("psyteachr/introdataviz")
library(ungeviz)  # devtools::install_github("wilkelab/ungeviz")

setwd("~/Documents/projects/R_data_analysis/cluster_size/20230926_all")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
times <- c("4h", "24h")
times_color <- c("#FE9929", "#41B6C4")  # light orange, light blue
names(times_color) <- times
#show_col(times_color, labels = TRUE)
conditions <- c("PM", "PA")
#conditions_color <- get_pal_colors("Set1")[c(2,1)]  # consistent with Bozdag2023, blue, red # "#377EB8" "#E41A1C"
#names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
#conditions_label <- c("Mixotrophic", "Anaerobic")
lines <- as.character(1:5)
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- paste0("t", c(0, 200, 400, 600, 1000))
evotimes2 <- gsub("t", "", evotimes, fixed = TRUE) %>% as.numeric()
#replicates <- paste0("Rep", 1:3)
wells <- paste0("Well", 1:2)

# Load data
# Input file name: e.g., "4h_PA_1_t200_Results.csv", "4h_PM_NA_t0_Results.csv", "4h_Rep1_PA_1_t600_Well1_Results.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1) %>% dplyr::select(BF.Area)
  measure_vars <- colnames(df)
  # Extract metadata
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  df$Sample <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Time <- metadata[1]
  if (length(metadata) == 4) {
    df$Condition <- metadata[2]
    df$Line <- metadata[3]
    df$EvoTime <- metadata[4]
    df$Replicate <- NA
    df$Well <- NA
  } else if (length(metadata) == 6) {
    df$Condition <- metadata[3]
    df$Line <- metadata[4]
    df$EvoTime <- metadata[5]
    df$Replicate <- metadata[2]
    df$Well <- metadata[6]
  }
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Time, levels = times), 
                 factor(Condition, levels = conditions), 
                 Line, 
                 factor(EvoTime, levels = evotimes), 
                 Replicate, 
                 Well)
colnames(data)[colnames(data) == "BF.Area"] <- "Area"
summary(data)
# Check if data contains objects < 50um^2
# dplyr::filter(data, Area < 50) %>% View()
# None
# nrow(data)
# data <- dplyr::filter(data, Area > 50)
# nrow(data)
# Check if data contains outlier objects like very large area
# dplyr::slice_max(data, order_by = Area, n = 10) %>% View()
# Looks fine
# nrow(data)
# data <- dplyr::filter(data, Area < 2.8e8)
# nrow(data)

# Data formatting
#samples <- unique(data$Sample)
#data$Sample <- factor(data$Sample, levels = samples)
data$Time <- factor(data$Time, levels = times)
data$Condition <- factor(data$Condition, levels = conditions)
data$Line <- factor(data$Line, levels = lines)
data$EvoTime <- factor(data$EvoTime, levels = evotimes)
#data$Replicate <- factor(data$Replicate, levels = replicates)
data$Well <- factor(data$Well, levels = wells)
summary(data)

# Check consistency across wells
id_vars <- c("Sample", "Time", "Condition", "Line", "EvoTime", "Replicate", "Well")
# Data summary: cluster count
data_summary_well <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = n())
# Consistent

# Convert area to radius and volume
data$Radius <- sqrt(data$Area / pi)
data$Volume <- 4/3 * pi * (data$Radius)^3

# Aggregate wells (here only one replicate per strain-time)
# Calculate cluster radius: weighted mean radius, median radius, mean of top 5% radius (proxy for radius at reproduction)
id_vars <- c("Time", "Condition", "Line", "EvoTime")
data_summary <- data %>%
    dplyr::group_by_at(id_vars) %>%
    dplyr::summarize(Count = n(), 
                     Weighted_mean_radius = (weighted.mean(Volume, Volume) / (4/3 * pi)) ^ (1/3), 
                     Median_radius = median(Radius),
                     Mean_top_radius = Radius[Radius > quantile(Radius, prob = 0.95)] %>% mean())

# Prepare data for plotting
data_summary_t0 <- dplyr::filter(data_summary, EvoTime == "t0")
data_summary_plot <- dplyr::filter(data_summary, EvoTime != "t0")
for (line in c(lines, NA)) {  # final NA is for plotting gray color at the top of overlapping points at t0
  data_summary_plot <- rbind(data_summary_plot, 
                             dplyr::mutate(data_summary_t0, Line = line))
}
data_t0 <- dplyr::filter(data, EvoTime == "t0")
data_plot <- dplyr::filter(data, EvoTime != "t0")
for (line in lines) {
  data_plot <- rbind(data_plot, 
                     dplyr::mutate(data_t0, Line = line))
}

# Plot evolution of weighted mean cluster radius (4h and 24h)
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()), 
       aes(x = EvoTime2, y = Weighted_mean_radius, color = Line)) +
  facet_grid(Time ~ Condition) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(trans = "log2", n.breaks = 6, limits = c(16,NA)) +
  labs(x = "Days of evolution", y = expression(Cluster~radius~(mu*m)), color = "Line") +
  ggplot_custom_theme4
save_ggplot("cluster_size_evo_times", width = 7, height = 7)
save_ggplot("cluster_size_evo_times", width = 7, height = 7, mode = "paper")

# Plot evolution of weighted mean cluster radius (only 24h)
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()) %>%
         dplyr::filter(Time == "24h"),  ###
       aes(x = EvoTime2, y = Weighted_mean_radius, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +  ###
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(trans = "log2", n.breaks = 6, limits = c(16,NA)) +
  labs(x = "Days of evolution", y = expression(Cluster~radius~(mu*m)), color = "Line") +
  ggplot_custom_theme4
save_ggplot("cluster_size_evo_24h", width = 7, height = 4)
save_ggplot("cluster_size_evo_24h", width = 7, height = 4, mode = "paper")

# Plot evolution of weighted mean cluster radius (only 4h)
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()) %>%
         dplyr::filter(Time == "4h"),  ###
       aes(x = EvoTime2, y = Weighted_mean_radius, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(trans = "log2", n.breaks = 6, limits = c(16,NA)) +
  labs(x = "Days of evolution", y = expression(Cluster~radius~(mu*m)), color = "Line") +
  ggplot_custom_theme4
save_ggplot("cluster_size_evo_4h", width = 7, height = 4)
save_ggplot("cluster_size_evo_4h", width = 7, height = 4, mode = "paper")

# Plot evolution of weighted mean cluster radius (only 24h), add artificial ploidy data
# Load data
art_data_summary <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/cluster_size/20230926_artificial/data_summary.rds") %>% 
  dplyr::filter(Time == "24h")
summary(art_data_summary)
# Parameters
ploidys <- c("2N", "4N")
ploidys_linetype <- c("12", "42")  # dotted, dashed
names(ploidys_linetype) <- ploidys
# Plot
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()) %>%
         dplyr::filter(Time == "24h"),
       aes(x = EvoTime2, y = Weighted_mean_radius, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(trans = "log2", n.breaks = 6, limits = c(16,NA)) +
  labs(x = "Days of evolution", y = expression(Cluster~radius~(mu*m)), color = "Line") +
  geom_hline(data = art_data_summary,
             mapping = aes(yintercept = Weighted_mean_radius, linetype = Ploidy),
             color = "gray60", size = 0.75) +
  scale_linetype_manual(values = ploidys_linetype) +
  ggplot_custom_theme4 + 
  guides(color = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))
save_ggplot("cluster_size_evo_24h_artificial", width = 7, height = 4)
save_ggplot("cluster_size_evo_24h_artificial", width = 7, height = 4, mode = "paper")

# Plot distribution of cluster radius, 4h/24h as lefe/right violin, not weighted by cluster volume
ggplot(data_plot, aes(x = EvoTime, y = Radius, fill = Time, color = Time)) +
  facet_grid(Line ~ Condition) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, size = 0.25) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, size = 0.25) +
  scale_fill_manual(values = times_color) +
  scale_y_continuous(trans = "log2", limits = c(3.5, 1050), n.breaks = 7) +
  labs(y = expression(Cluster~radius~(mu*m))) +
  ggplot_custom_theme4 +
  theme(axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cluster_size_evo_distribution", width = 10, height = 12)
save_ggplot("cluster_size_evo_distribution", width = 10, height = 12, mode = "paper")

# Plot distribution of cluster radius, 4h/24h as lefe/right violin, weighted by cluster volume
# Plot weighted mean radius
ggplot(data_plot, aes(x = EvoTime, y = Radius, fill = Time, color = Time, weight = Volume)) +
  facet_grid(Line ~ Condition) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, size = 0.25) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, size = 0.25) +
  geom_point(data = data_summary_plot %>% dplyr::filter(!is.na(Line)),
             mapping = aes(x = EvoTime, y = Weighted_mean_radius, fill = Time), inherit.aes = FALSE,
             shape = "circle filled", size = 2, color = "black", alpha = 0.5, 
             position = position_dodge(width = 0.4)) +
  scale_fill_manual(values = times_color) +
  scale_y_continuous(trans = "log2", limits = c(3.5, 1050), n.breaks = 7) +
  labs(y = expression(Cluster~radius~(mu*m))) +
  ggplot_custom_theme4 +
  theme(axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cluster_size_evo_distribution_weighted", width = 10, height = 12)
save_ggplot("cluster_size_evo_distribution_weighted", width = 10, height = 12, mode = "paper")

# Save data
write.csv(data %>% dplyr::select(Time, Condition, Line, EvoTime, Area, Radius, Volume),
          file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")

