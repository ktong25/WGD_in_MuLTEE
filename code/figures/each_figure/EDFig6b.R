library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/competition_assay/20220909_2N4N_corrected_remove_PO")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
ploidys <- c("2N", "4N")
ploidys_color <- get_pal_colors("Dark2")[c(1,2)]  # green, orange # "#1B9E77" "#D95F02"
names(ploidys_color) <- ploidys
#show_col(ploidys_color, labels = TRUE)
conditions_color <- get_pal_colors("Set1")[c(2,1)]  # like in Bozdag2021, blue, red, green # "#377EB8" "#E41A1C" "#4DAF4A"
names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
conditions_label <- c("Mixotrophic", "Anaerobic")
conditions_label_n <- c("Mixotrophic", "Anaerobic")
conditions_label_n_name <- conditions_label_n  # for modifying strip text
names(conditions_label_n_name) <- conditions

# Load data
# Input file name: e.g., "cs_2Ng_YPD_Rep1_Results_corrected.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1)
  measure_vars <- colnames(df)
  # Extract metadata
  # e.g., "cs_2Ng_YPD_Rep1_Results.csv" -> "2N_PM_Rep1"
  filename <- basename(file)
  metadata <- substr(filename, stri_length("cs_")+1, stri_length(filename)-stri_length("_Results_corrected.csv"))
  metadata <- gsub("g_YPD", "_PM", metadata, fixed = TRUE)
  metadata <- gsub("g_YPGly", "_PO", metadata, fixed = TRUE)
  metadata <- gsub("p_YPD", "_PA", metadata, fixed = TRUE)
  df$Sample <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Condition <- metadata[2]
  df$Ploidy <- metadata[1]
  df$Replicate <- metadata[3]
  df$Ploidy_Replicate <- glue("{df$Ploidy} ({df$Replicate})")
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "04_measurement_cs", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::filter(Condition != "PO") %>%
  dplyr::arrange(factor(Condition, levels = conditions), Ploidy, Replicate)
summary(data)
# Check if data contains objects < 50um^2
# dplyr::filter(data, Area < 50) %>% View()
# This is due to 04_correction, filter out these objects
data <- data %>%
  dplyr::filter(Area > 50)
# Data formatting
samples <- unique(data$Sample)
data$Sample <- factor(data$Sample, levels = samples)
data$Condition <- factor(data$Condition, levels = conditions)
data$Ploidy <- factor(data$Ploidy, levels = ploidys)
replicates <- sort(unique(data$Replicate))
reps <- gsub("Rep", "", replicates, fixed = TRUE)
data$Replicate <- factor(data$Replicate, levels = replicates)
ploidy_replicates <- unique(data$Ploidy_Replicate)
data$Ploidy_Replicate <- factor(data$Ploidy_Replicate, levels = ploidy_replicates)
summary(data)
id_vars <- c("Sample", "Condition", "Ploidy", "Replicate", "Ploidy_Replicate")

# Data summary: cluster count
data_summary <- data %>%
    dplyr::group_by_at(id_vars) %>%
    dplyr::summarize(Count = n())

# Convert area to radius and volume
data$Radius <- sqrt(data$Area / pi)
data$Volume <- 4/3 * pi * (data$Radius)^3

##### Distinguish 2N vs 4N clusters #####

# For each sample (condition x [ploidy-replicate]), plot cluster area vs cell_top5_mean_area
ggplot(data, aes(x = Area, y = cell_top5_mean_area, color = Ploidy)) +
  facet_grid(Condition ~ Ploidy_Replicate) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = ploidys_color) +
  labs(x = expression(Cluster~area~(mu*m^2)), 
       y = expression(Mean~area~of~the~five~largest~detected~cells~(mu*m^2))) + 
  ggplot_custom_theme +
  theme(legend.position = "none", 
        panel.grid = element_line(size = 0.1, color = "gray50"), 
        axis.title = element_text(size = rel(1.75)), 
        strip.text = element_text(size = rel(1.75)))
save_ggplot("cs_area_vs_top5", 
            width = length(ploidy_replicates) * 3, height = length(conditions) * 3)

# Choose cell area cutoff for each sample (no longer use cluster area cutoff)
# Manually determine cutoff
write.csv(data_summary, file = "cs_data_summary_for_cutoff.csv", row.names = FALSE)
# Copy file as "data_summary_for_cutoff_add_cutoff.csv", manually input cutoff values into a new Cutoff column
data_summary$Cutoff <- read.csv("cs_data_summary_for_cutoff_add_cutoff.csv", row.names = NULL)$Cutoff

# Plot cutoff for each sample to double check
ggplot(data, aes(x = Area, y = cell_top5_mean_area, color = Ploidy)) +
  facet_grid(Condition ~ Ploidy_Replicate) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = ploidys_color) +
  labs(x = expression(Cluster~area~(mu*m^2)), 
       y = expression(Mean~area~of~the~five~largest~detected~cells~(mu*m^2))) + 
  geom_hline(data = data_summary,  ###
             mapping = aes(yintercept = Cutoff),
             size = 0.2, linetype = "dashed", color = "gray50") +
  ggplot_custom_theme +
  theme(legend.position = "none", 
        #panel.grid = element_line(size = 0.1, color = "gray50"),  ###
        axis.title = element_text(size = rel(1.75)), 
        strip.text = element_text(size = rel(1.75)))
save_ggplot("cs_area_vs_top5_cutoff", 
            width = length(ploidy_replicates) * 3, height = length(conditions) * 3)

# Plot with two ploidies merged into one plot
# Use Rep1 of each condition
data_subset <- data %>% dplyr::filter(Replicate == "Rep1")
set.seed(42)
data_subset <- data_subset[sample(1:nrow(data_subset)), ]  # randomize rows for plotting data points
ggplot(data_subset, 
       aes(x = Area, y = cell_top5_mean_area, color = Ploidy)) +
  facet_wrap(~Condition, nrow = 1, 
             labeller = labeller(Condition = conditions_label_n_name)) +
  geom_point(alpha = 0.5, size = 0.75) +
  scale_color_manual(values = ploidys_color) +
  labs(x = expression(Cluster~area~(mu*m^2)), 
       y = expression(atop(Average~area~of~the~five, 
                           largest~cells~detected~(mu*m^2)))) +
  geom_hline(data = data_summary,  ###
             mapping = aes(yintercept = Cutoff), 
             size = 0.75, linetype = "dashed", color = "gray30") +
  ggplot_custom_theme +
  theme(legend.position = "right", 
        legend.title = element_blank(),
        strip.background = element_blank(), 
        axis.title = element_text(size = rel(1.75)), 
        strip.text = element_text(size = rel(2)), 
        legend.text = element_text(size = rel(1.75)))
save_ggplot("cs_area_vs_top5_cutoff_merge", width = 10, height = 6)
save_ggplot("cs_area_vs_top5_cutoff_merge", width = 10, height = 6, mode = "paper")

# Plot cluster size (cluster radius)
# Plot four replicates of each condition
#condition <- "PA"
for (condition in conditions) {
    ggplot(data %>% 
             dplyr::filter(Condition == condition) %>%
             dplyr::mutate(Replicate = factor(Replicate, levels = replicates, labels = reps)), 
         aes(x = Replicate, 
             y = Radius, fill = Ploidy)) +  ###
    facet_wrap(~Ploidy, nrow = 1) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.4, size = 0.75) + 
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white") +
    scale_fill_manual(values = ploidys_color) +
    {if (condition == "PA") xlab("Petite mutant") else xlab("Replicate")} +
    ylab(expression(Cluster~radius~(mu*m))) +  ###
    theme_classic(base_size = 12) + 
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = rel(1.75)), 
          axis.title = element_text(size = rel(1.5)), 
          axis.text = element_text(size = rel(1.25)), 
          legend.position = "none")
  save_ggplot(glue("{condition}_all_reps_cluster_radius"), width = 5, height = 4)  ###
  if (condition == "PA") save_ggplot(glue("{condition}_all_reps_cluster_radius"), width = 7, height = 5, mode = "paper")  ###
}
# Plot one replicate (rep1) of all conditions
ggplot(data %>% 
         dplyr::filter(Replicate == "Rep1") %>%
         dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = conditions_label_n)), 
       aes(x = Ploidy, y = Radius, fill = Ploidy)) +  ###
  facet_wrap(~Condition, 
             nrow = 1, strip.position = "bottom") +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.4, size = 0.75) + 
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white") +
  scale_fill_manual(values = ploidys_color) +
  ylab(expression(Cluster~radius~(mu*m))) +   ###
  theme_classic(base_size = 12) + 
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = rel(1.5)), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)), 
        axis.title.x = element_blank(), 
        legend.position = "none") +
  ggpubr::stat_compare_means(
    comparisons = list( c("2N", "4N") ), 
    method = "t.test",  # default var.equal = FALSE
    label = "p.signif",   # "p.format"
    bracket.size = 0.5, size = 6, vjust = 0.4)
save_ggplot(glue("all_conditions_rep1_cluster_radius"), width = 5, height = 4)  ###
save_ggplot(glue("all_conditions_rep1_cluster_radius"), width = 7, height = 5, mode = "paper")  ###

# Save data
write.csv(data, file = "cs_data.csv", row.names = FALSE)
write.csv(data %>% dplyr::filter( Condition == "PA" | Replicate == "Rep1" ), 
          file = "cs_data_used.csv", row.names = FALSE)
write.csv(data_summary, file = "cs_data_summary.csv", row.names = FALSE)
write.csv(data_summary %>% dplyr::filter( Condition == "PA" | Replicate == "Rep1" ), 
          file = "cs_data_summary_used.csv", row.names = FALSE)



