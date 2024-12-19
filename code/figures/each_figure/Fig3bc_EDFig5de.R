library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggforce)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)
library(rstatix)
library(ggpubr)

setwd("~/Documents/projects/R_data_analysis/cell_size_aspect_ratio/20230926_artificial")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
#conditions_color <- get_pal_colors("Set1")[c(2,1)]  # consistent with Bozdag2023, blue, red # "#377EB8" "#E41A1C"
#names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
conditions_label <- c("Mixotrophic", "Anaerobic")
ploidys <- c("2N", "4N")
ploidys_color <- get_pal_colors("Dark2")[c(1,2)]  # green, orange # "#1B9E77" "#D95F02"
names(ploidys_color) <- ploidys
#show_col(ploidys_color, labels = TRUE)
replicates <- paste0("Rep", 1:4)

# Load data
# Input file name: e.g., "PA_2N_Rep1_001_Results.csv"
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
  df$Ploidy <- metadata[2]
  df$Replicate <- metadata[3]
  df$ImageID <- as.integer(metadata[4])
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Condition, levels = conditions), 
                 Ploidy, 
                 Replicate,  
                 ImageID)
summary(data)
# Check if data contains objects < 1um^2
# dplyr::filter(data, Area < 1) %>% View()
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
# Check if data contains outlier objects like very large AR
# dplyr::slice_max(data, order_by = AR, n = 10) %>% View()
# nrow(data)
# data <- dplyr::filter(data, AR < 2.3)
# nrow(data)

# Data formatting
data$Condition <- factor(data$Condition, levels = conditions)
data$Ploidy <- factor(data$Ploidy, levels = ploidys)
data$Replicate <- factor(data$Replicate, levels = replicates)
imageids <- unique(data$ImageID)
data$ImageID <- factor(data$ImageID, levels = imageids)
summary(data)

# Calculate cell volume
data$Volume <- 4/3 * pi * (data$Major/2) * (data$Minor/2)^2
summary(data)

# This dataset is not manually corrected (due to re-running using newest cellpose version)
# Thus roughly filter based on the upper cell AR limit observed in the previous manually-corrected version of this dataset (albeit running using a older cellpose version)
# Cell volume looks fine
data <- dplyr::filter(data, 
                      (Condition == "PM" & Ploidy == "2N" & AR < 2.1) |
                      (Condition == "PM" & Ploidy == "4N" & AR < 2.4) |
                      (Condition == "PA" & Ploidy == "2N" & AR < 2) |
                      (Condition == "PA" & Ploidy == "4N" & AR < 2.2))

# Check all replicates and images
id_vars <- c("Sample", "Condition", "Ploidy", "Replicate", "ImageID")
data_summary_rep_image <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = n())

# Check consistency across replicates
# Aggregate images of each replicate
# Calculate mean/median/sd of cell volume/AR
id_vars <- c("Condition", "Ploidy", "Replicate")
data_summary_rep <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarise(Count = n(), 
                   Mean_volume = mean(Volume), 
                   Mean_AR = mean(AR), 
                   Median_volume = median(Volume), 
                   Median_AR = median(AR), 
                   SD_volume = sd(Volume), 
                   SD_AR = sd(AR))

# Aggregate replicates of each strain
id_vars <- c("Condition", "Ploidy")
data_summary <- data_summary_rep %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = mean(Count), 
                   Mean_volume = mean(Mean_volume), 
                   Mean_AR = mean(Mean_AR), 
                   Median_volume = mean(Median_volume), 
                   Median_AR = mean(Median_AR), 
                   SD_volume = mean(SD_volume), 
                   SD_AR = mean(SD_AR))

# Perform statistical tests
data_summary_rep_stat <- list()
# Cell volume
data_summary_rep_stat[["Volume"]] <- data_summary_rep %>%
  dplyr::group_by(Condition) %>%
  rstatix::t_test(Mean_volume ~ Ploidy, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  #adjust_pvalue(method = "bonferroni") %>%
  #add_significance("p.adj") %>%
  add_xy_position(x = "Condition", fun = "max", dodge = 0.8)
# Cell AR
data_summary_rep_stat[["AR"]] <- data_summary_rep %>%
  dplyr::group_by(Condition) %>%
  rstatix::t_test(Mean_AR ~ Ploidy, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  #adjust_pvalue(method = "bonferroni") %>%
  #add_significance("p.adj") %>%
  add_xy_position(x = "Condition", fun = "max", dodge = 0.8)

# Plot mean cell volume
ggplot(data_summary_rep %>%
         dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = conditions_label)),
       aes(x = Condition, y = Mean_volume)) +
  stat_summary(mapping = aes(fill = Ploidy), 
               geom = "bar", fun = mean, color = NA, alpha = 0.3, width = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_point(mapping = aes(color = Ploidy), 
             alpha = 0.6, size = 2.5, 
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8, seed = 0)) +
  stat_summary(mapping = aes(color = Ploidy), 
               geom = "errorbar", fun.data = mean_se, width = 0.4, linewidth = 0.75,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = ploidys_color) +
  scale_color_manual(values = ploidys_color) +
  stat_pvalue_manual(data_summary_rep_stat[["Volume"]], label = "p = {p}",  # Ploidy should not be in ggplot(aes()) but in separate geom_XXX(aes())
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0,0.2)), n.breaks = 5) +
  labs(x = NULL, y = expression(Cell~volume~(mu*m^3))) +
  ggplot_custom_theme4 +
  theme(legend.position = c(0.9, 0.95),
        legend.title = element_blank())
save_ggplot("cell_volume_ploidy", width = 4, height = 4)
save_ggplot("cell_volume_ploidy", width = 4, height = 4, mode = "paper")

# Plot mean cell AR
# Let bar plot baseline be 1, by shifting y values -1 and y axis ticks +1
y_baseline <- 1
ggplot(data_summary_rep %>%
         dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = conditions_label)),
       aes(x = Condition, y = Mean_AR - y_baseline)) +  ###
  stat_summary(mapping = aes(fill = Ploidy), 
               geom = "bar", fun = mean, color = NA, alpha = 0.3, width = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_point(mapping = aes(color = Ploidy), 
             alpha = 0.6, size = 2.5, 
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8, seed = 0)) +
  stat_summary(mapping = aes(color = Ploidy), 
               geom = "errorbar", fun.data = mean_se, width = 0.4, linewidth = 0.75,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = ploidys_color) +
  scale_color_manual(values = ploidys_color) +
  stat_pvalue_manual(data_summary_rep_stat[["AR"]] %>% 
                       dplyr::mutate(y.position = y.position-1),   ###
                     label = "p = {p}",  # Ploidy should not be in ggplot(aes()) but in separate geom_XXX(aes())
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_y_continuous(limits = c(0,NA), labels = function(y) y + y_baseline, ###
                     expand = expansion(mult = c(0,0.15)), n.breaks = 5) +  ###
  labs(x = NULL, y = "Cell aspect ratio") +  ###
  ggplot_custom_theme4 +
  theme(legend.position = c(0.9, 0.95),
        legend.title = element_blank())
save_ggplot("cell_AR_ploidy", width = 4, height = 4)
save_ggplot("cell_AR_ploidy", width = 4, height = 4, mode = "paper")

# Plot distribution of cell volume
ggplot(data %>% dplyr::mutate(Replicate = gsub("Rep", "", Replicate, fixed = TRUE), 
                              Condition = factor(Condition, levels = conditions, labels = conditions_label)), 
       aes(x = Replicate, y = Volume, fill = Ploidy, color = Ploidy)) +
  facet_grid(Condition ~ Ploidy) +
  geom_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) + 
  geom_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, linewidth = 0.25) + 
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = ploidys_color) +
  scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0,0.05)), n.breaks = 5) +
  labs(y = expression(Cell~volume~(mu*m^3))) +
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cell_volume_ploidy_distribution", width = 6, height = 5)
save_ggplot("cell_volume_ploidy_distribution", width = 6, height = 5, mode = "paper")

# Plot distribution of cell AR
ggplot(data %>% dplyr::mutate(Replicate = gsub("Rep", "", Replicate, fixed = TRUE), 
                              Condition = factor(Condition, levels = conditions, labels = conditions_label)), 
       aes(x = Replicate, y = AR, fill = Ploidy, color = Ploidy)) +  ###
  facet_grid(Condition ~ Ploidy) +
  geom_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) + 
  geom_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, linewidth = 0.25) + 
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = ploidys_color) +
  scale_y_continuous(limits = c(1,NA), expand = expansion(mult = c(0,0.05)), n.breaks = 6) +  ###
  labs(y = "Cell aspect ratio") +   ###
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cell_AR_ploidy_distribution", width = 6, height = 5)
save_ggplot("cell_AR_ploidy_distribution", width = 6, height = 5, mode = "paper")

# Save data
write.csv(data %>% dplyr::select(!Sample), file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
write.csv(data_summary_rep, file = "data_summary_rep.csv", row.names = FALSE)
write.csv(data_summary_rep_stat[["Volume"]] %>% apply(2, as.character), 
          file = "data_summary_rep_volume_stat.csv", row.names = FALSE)
write.csv(data_summary_rep_stat[["AR"]] %>% apply(2, as.character), 
          file = "data_summary_rep_AR_stat.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")
saveRDS(data_summary_rep, file = "data_summary_rep.rds")
saveRDS(data_summary_rep_stat, file = "data_summary_rep_stat.rds")
