library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(readxl)
library(reshape2)

setwd("~/Documents/projects/R_data_analysis/cluster_size/20231211_Bozdag2023_vs_2N4N")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
lines <- paste0("PA", as.character(1:5))
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
evotimes <- seq(0,200,50)
ploidys <- c("2N", "4N")
ploidys_linetype <- c("12", "42")  # dotted, dashed
names(ploidys_linetype) <- ploidys

# Cluster size

# Load data
data_cluster <- read_excel("01_Source_data.xlsx", sheet = "Fig1e", skip = 1)
colnames(data_cluster)[1] <- "EvoTime"
data_cluster <- data_cluster %>% 
  dplyr::filter(EvoTime %in% evotimes) %>%
  reshape2::melt(id.vars = "EvoTime",
                 variable.name = "Line",
                 value.name = "Weighted_mean_radius") %>%
  dplyr::mutate(Line = factor(Line, levels = lines))
summary(data_cluster)
data_ploidy_cluster <- read.csv("cluster_size_20230926_artificial_data_summary.csv", 
                                row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::filter(Time == "24h" & Condition == "PA" & Ploidy == "4N") %>%
  dplyr::select(Ploidy, Weighted_mean_radius)

# Plot evolution of weighted mean cluster radius, add artificial ploidy data
ggplot(data_cluster,
       aes(x = EvoTime, y = Weighted_mean_radius, color = Line)) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes, 
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(trans = "log2", n.breaks = 6, limits = c(16,NA)) +
  labs(x = "Days of evolution", y = expression(Cluster~radius~(mu*m)), color = "Line") +
  geom_hline(data = data_ploidy_cluster,
             mapping = aes(yintercept = Weighted_mean_radius, linetype = Ploidy),
             color = "gray60", size = 0.75) +
  scale_linetype_manual(values = ploidys_linetype) +
  ggplot_custom_theme4 + 
  guides(color = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))
save_ggplot("cluster_size_evo_pop_artificial", width = 5, height = 4)
save_ggplot("cluster_size_evo_pop_artificial", width = 5, height = 4, mode = "paper")

### Cell aspect ratio

# Load data
data_cell <- read_excel("01_Source_data.xlsx", sheet = "Fig2d", skip = 1)
colnames(data_cell)[1] <- "EvoTime"
data_cell <- data_cell %>% 
  dplyr::filter(EvoTime %in% evotimes) %>%
  reshape2::melt(id.vars = "EvoTime",
                 variable.name = "Line",
                 value.name = "Mean_AR") %>%
  dplyr::mutate(Line = factor(Line, levels = lines))
summary(data_cell)
data_ploidy_cell <- read.csv("cell_size_aspect_ratio_20230926_artificial_data_summary.csv", 
                             row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::filter(Condition == "PA" & Ploidy == "4N") %>%
  dplyr::select(Ploidy, Mean_AR)

# Plot evolution of weighted mean cluster radius, add artificial ploidy data
ggplot(data_cell,
       aes(x = EvoTime, y = Mean_AR, color = Line)) +  ###
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes, 
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(limits = c(1,NA), n.breaks = 4,  ###
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days of evolution", y = "Cell aspect ratio", color = "Line") +
  geom_hline(data = data_ploidy_cell,
             mapping = aes(yintercept = Mean_AR, linetype = Ploidy),
             color = "gray60", size = 0.75) +
  scale_linetype_manual(values = ploidys_linetype) +
  ggplot_custom_theme4 + 
  guides(color = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))
save_ggplot("cell_AR_evo_pop_artificial", width = 5, height = 4)
save_ggplot("cell_AR_evo_pop_artificial", width = 5, height = 4, mode = "paper")
