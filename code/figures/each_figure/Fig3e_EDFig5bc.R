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
library(rstatix)
library(ggpubr)

setwd("~/Documents/projects/R_data_analysis/cluster_size/20230926_artificial")
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
conditions_label <- c("Mixotrophic", "Anaerobic")
ploidys <- c("2N", "4N")
ploidys_color <- get_pal_colors("Dark2")[c(1,2)]  # green, orange # "#1B9E77" "#D95F02"
names(ploidys_color) <- ploidys
#show_col(ploidys_color, labels = TRUE)
#replicates <- paste0("Rep", 1:5)
replicates <- paste0("Rep", 1:4)

# Load data
# Input file name: e.g., "4h_PA_2N_Rep1_Results.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1) %>% dplyr::select(BF.Area)
  measure_vars <- colnames(df)
  # Extract metadata
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  df$Sample <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Time <- metadata[1]
  df$Condition <- metadata[2]
  df$Ploidy <- metadata[3]
  df$Replicate <- metadata[4]
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Time, levels = times), 
                 factor(Condition, levels = conditions), 
                 Ploidy, 
                 Replicate)
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
data$Ploidy <- factor(data$Ploidy, levels = ploidys)
data$Replicate <- factor(data$Replicate, levels = replicates)
summary(data)

# Convert area to radius and volume
data$Radius <- sqrt(data$Area / pi)
data$Volume <- 4/3 * pi * (data$Radius)^3

# Check consistency across replicates
# Calculate cluster radius: weighted mean radius, median radius, mean of top 5% radius (proxy for radius at reproduction)
id_vars <- c("Time", "Condition", "Ploidy", "Replicate")
data_summary_rep <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = n(), 
                   Weighted_mean_radius = (weighted.mean(Volume, Volume) / (4/3 * pi)) ^ (1/3), 
                   Median_radius = median(Radius),
                   Mean_top_radius = Radius[Radius > quantile(Radius, prob = 0.95)] %>% mean())
#write.csv(data_summary_rep, file = "data_summary_all_reps.csv", row.names = FALSE)

# Aggregate replicates
id_vars <- c("Time", "Condition", "Ploidy")
data_summary <- data_summary_rep %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = mean(Count), 
                   Weighted_mean_radius = mean(Weighted_mean_radius), 
                   Median_radius = mean(Median_radius), 
                   Mean_top_radius = mean(Mean_top_radius))

# Perform statistical tests (24h only)
data_summary_rep_24h_stat <- data_summary_rep %>% dplyr::filter(Time == "24h") %>%
  dplyr::group_by(Condition) %>%
  rstatix::t_test(Weighted_mean_radius ~ Ploidy, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  #adjust_pvalue(method = "bonferroni") %>%
  #add_significance("p.adj") %>%
  add_xy_position(x = "Condition", fun = "max", dodge = 0.8)

# Plot 24h weighted mean radius
ggplot(data_summary_rep %>% 
         dplyr::filter(Time == "24h") %>%
         dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = conditions_label)),
       aes(x = Condition, y = Weighted_mean_radius)) +
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
  stat_pvalue_manual(data_summary_rep_24h_stat, label = "p = {p}",  # Ploidy should not be in ggplot(aes()) but in separate geom_XXX(aes())
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), n.breaks = 5) +
  labs(x = NULL, y = expression(Cluster~radius~(mu*m))) +
  ggplot_custom_theme4 +
  theme(legend.position = c(0.9, 0.95),
        legend.title = element_blank())
save_ggplot("cluster_size_ploidy", width = 4, height = 4)
save_ggplot("cluster_size_ploidy", width = 4, height = 4, mode = "paper")

# Plot distribution of cluster radius, 4h/24h as lefe/right violin, not weighted by cluster volume
ggplot(data %>% dplyr::mutate(Replicate = gsub("Rep", "", Replicate, fixed = TRUE), 
                              Condition = factor(Condition, levels = conditions, labels = conditions_label)), 
       aes(x = Replicate, y = Radius, fill = Time, color = Time)) +
  facet_grid(Condition ~ Ploidy) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, linewidth = 0.25) +
  scale_fill_manual(values = times_color) +
  labs(y = expression(Cluster~radius~(mu*m))) +
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cluster_size_ploidy_distribution", width = 6, height = 5)
save_ggplot("cluster_size_ploidy_distribution", width = 6, height = 5, mode = "paper")

# Plot distribution of cluster radius, 4h/24h as lefe/right violin, weighted by cluster volume
# Plot weighted mean radius
ggplot(data %>% dplyr::mutate(Replicate = gsub("Rep", "", Replicate, fixed = TRUE), 
                              Condition = factor(Condition, levels = conditions, labels = conditions_label)), 
       aes(x = Replicate, y = Radius, fill = Time, color = Time, weight = Volume)) +
  facet_grid(Condition ~ Ploidy) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, size = 0.25) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, size = 0.25) +
  geom_point(data = data_summary_rep %>%
               dplyr::mutate(Replicate = gsub("Rep", "", Replicate, fixed = TRUE),
                             Condition = factor(Condition, levels = conditions, labels = conditions_label)),
             mapping = aes(x = Replicate, y = Weighted_mean_radius, fill = Time), inherit.aes = FALSE,
             shape = "circle filled", size = 2, color = "black", alpha = 0.5,
             position = position_dodge(width = 0.4)) +
  scale_fill_manual(values = times_color) +
  labs(y = expression(Cluster~radius~(mu*m))) +
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25))
save_ggplot("cluster_size_ploidy_distribution_weighted", width = 6, height = 5)
save_ggplot("cluster_size_ploidy_distribution_weighted", width = 6, height = 5, mode = "paper")

# Save data
write.csv(data %>% dplyr::select(!Sample), file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
write.csv(data_summary_rep, file = "data_summary_rep.csv", row.names = FALSE)
write.csv(data_summary_rep_24h_stat %>% apply(2, as.character), 
          file = "data_summary_rep_24h_stat.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")
saveRDS(data_summary_rep, file = "data_summary_rep.rds")
saveRDS(data_summary_rep_24h_stat, file = "data_summary_rep_24h_stat.rds")
