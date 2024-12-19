library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggforce)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/cell_size_aspect_ratio/20230926_all")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
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

# Load data
# Input file name: e.g., "PA_1_t200_001_Results.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1)
  colnames(df) <- gsub("BF.", "", colnames(df), fixed = TRUE)
  measure_vars <- colnames(df)
  # Extract metadata
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  df$Sample <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Condition <- metadata[1]
  df$Line <- metadata[2]
  df$EvoTime <- metadata[3]
  df$ImageID <- as.integer(metadata[4])
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Condition, levels = conditions), 
                 Line, 
                 factor(EvoTime, levels = evotimes), 
                 ImageID)
summary(data)
# Check if data contains objects < 1um^2
# dplyr::filter(data, Area < 1) %>% View()
# None
# nrow(data)
# data <- dplyr::filter(data, Area > 1)
# nrow(data)
# Check if data contains outlier objects like very large area
# dplyr::slice_max(data, order_by = Area, n = 10) %>% View()
# Looks fine
# nrow(data)
# data <- dplyr::filter(data, Area < 2.8e8)
# nrow(data)

# Data formatting
data$Condition <- factor(data$Condition, levels = conditions)
data$Line <- factor(data$Line, levels = lines)
data$EvoTime <- factor(data$EvoTime, levels = evotimes)
imageids <- unique(data$ImageID)
data$ImageID <- factor(data$ImageID, levels = imageids)
summary(data)

# Calculate cell volume
data$Volume <- 4/3 * pi * (data$Major/2) * (data$Minor/2)^2
summary(data)

# Check all images
id_vars <- c("Sample", "Condition", "Line", "EvoTime", "ImageID")
data_summary_image <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = n())

# Aggregate images of each strain
# Calculate mean/median/sd of cell volume/AR
id_vars <- c("Condition", "Line", "EvoTime")
data_summary <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarise(Count = n(), 
                   Mean_volume = mean(Volume), 
                   Mean_AR = mean(AR), 
                   Median_volume = median(Volume), 
                   Median_AR = median(AR), 
                   SD_volume = sd(Volume), 
                   SD_AR = sd(AR))

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

# Plot evolution of mean cell volume
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()),
       aes(x = EvoTime2, y = Mean_volume, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0,NA), n.breaks = 6, 
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days of evolution", y = expression(Cell~volume~(mu*m^3)), color = "Line") +
  ggplot_custom_theme4
save_ggplot("cell_volume_evo", width = 7, height = 4)
save_ggplot("cell_volume_evo", width = 7, height = 4, mode = "paper")

# Plot evolution of mean cell AR
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()),
       aes(x = EvoTime2, y = Mean_AR, color = Line)) +  ###
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(1,NA), n.breaks = 6,  ###
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days of evolution", y = "Cell aspect ratio", color = "Line") +
  ggplot_custom_theme4
save_ggplot("cell_AR_evo", width = 7, height = 4)
save_ggplot("cell_AR_evo", width = 7, height = 4, mode = "paper")

# Plot evolution of mean cell volume/AR, add artificial ploidy data
# Load data
art_data_summary <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/cell_size_aspect_ratio/20230926_artificial/data_summary.rds")
summary(art_data_summary)
# Parameters
ploidys <- c("2N", "4N")
ploidys_linetype <- c("12", "42")  # dotted, dashed
names(ploidys_linetype) <- ploidys
# Plot for mean cell volume
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()),
       aes(x = EvoTime2, y = Mean_volume, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0,NA), n.breaks = 6, 
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days of evolution", y = expression(Cell~volume~(mu*m^3)), color = "Line") +
  geom_hline(data = art_data_summary,
             mapping = aes(yintercept = Mean_volume, linetype = Ploidy),
             color = "gray60", size = 0.75) +
  scale_linetype_manual(values = ploidys_linetype) +
  ggplot_custom_theme4 + 
  guides(color = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))
save_ggplot("cell_volume_evo_artificial", width = 7, height = 4)
save_ggplot("cell_volume_evo_artificial", width = 7, height = 4, mode = "paper")
# Plot for mean cell AR
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()),
       aes(x = EvoTime2, y = Mean_AR, color = Line)) +  ###
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(1,NA), n.breaks = 6,  ###
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days of evolution", y = "Cell aspect ratio", color = "Line") +
  geom_hline(data = art_data_summary,
             mapping = aes(yintercept = Mean_AR, linetype = Ploidy),
             color = "gray60", size = 0.75) +
  scale_linetype_manual(values = ploidys_linetype) +
  ggplot_custom_theme4 + 
  guides(color = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))
save_ggplot("cell_AR_evo_artificial", width = 7, height = 4)
save_ggplot("cell_AR_evo_artificial", width = 7, height = 4, mode = "paper")

# Plot distribution of cell volume
ggplot(data_plot, aes(x = EvoTime, y = Volume)) +
  facet_grid(Line ~ Condition) +
  geom_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) + 
  geom_violin(scale = "width", trim = TRUE, fill = "#B284BE", color = "black", alpha = 0.3, linewidth = 0.25) + 
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0,0.05)), n.breaks = 5) +
  labs(x = NULL, y = expression(Cell~volume~(mu*m^3))) +
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cell_volume_distribution", width = 7, height = 10)
save_ggplot("cell_volume_distribution", width = 7, height = 10, mode = "paper")

# Plot distribution of cell AR
ggplot(data_plot, aes(x = EvoTime, y = AR)) +  ###
  facet_grid(Line ~ Condition) +
  geom_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) + 
  geom_violin(scale = "width", trim = TRUE, fill = "#C2B280", color = "black", alpha = 0.3, linewidth = 0.25) + ###
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_y_continuous(limits = c(1,NA), expand = expansion(mult = c(0,0.05)), n.breaks = 5) +  ###
  labs(x = NULL, y = "Cell aspect ratio") +  ###
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cell_AR_distribution", width = 7, height = 10)
save_ggplot("cell_AR_distribution", width = 7, height = 10, mode = "paper")

# Save data
write.csv(data %>% dplyr::select(!Sample), file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")
