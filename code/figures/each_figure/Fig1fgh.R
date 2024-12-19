library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("~/Documents/projects/R_data_analysis/snowflake_physics_modeling/20231211")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."

# Load data
data <- read.csv("ForKai_all_sim_measurements_PM_overlap.csv", row.names = NULL, 
                 col.names = c("Cell_volume", "Cell_AR", "Cluster_volume", "Cluster_ncells", "Cluster_packfrac"))
unique(data$Cell_volume)
length(unique(data$Cell_volume))
unique(data$Cell_AR)
length(unique(data$Cell_AR))
data$Cluster_radius <- (data$Cluster_volume / (4/3*pi)) ^ (1/3)

# Get mean and standard error
se <- function(x) sd(x)/sqrt(length(x))  # Function for calculating standard error
data_summary <- data %>%
  dplyr::group_by(Cell_volume, Cell_AR) %>%
  dplyr::summarise(Count = n(), 
                   Cluster_radius_mean = mean(Cluster_radius), 
                   Cluster_radius_se = se(Cluster_radius), 
                   Cluster_ncells_mean = mean(Cluster_ncells), 
                   Cluster_ncells_se = se(Cluster_ncells), 
                   Cluster_packfrac_mean = mean(Cluster_packfrac), 
                   Cluster_packfrac_se = se(Cluster_packfrac))

# Plot cell volume/AR vs cluster radius/ncells/packfrac
clabel_map <- c("Cluster_radius" = expression(Cluster~radius~at~division~(mu*m)),
                "Cluster_ncells" = "Cluster cell number at division", 
                "Cluster_packfrac" = "Cluster packing fraction at division")
for (cvar in c("Cluster_radius", "Cluster_ncells", "Cluster_packfrac")) {
  # cvar <- "Cluster_radius"
  ggplot(data_summary,
         aes(x = factor(Cell_volume), y = factor(Cell_AR), fill = .data[[paste0(cvar, "_mean")]])) +
    geom_tile(color = NA) + 
    scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "Blues"))(100)[90:1]) +
    #coord_fixed(ratio = 1) +
    labs(title = clabel_map[cvar], 
         x = expression(Cell~volume~(mu*m^3)), 
         y = "Cell aspect ratio", 
         fill = NULL) +
    ggplot_custom_theme3 +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(), 
          legend.title = element_blank(), 
          legend.text = element_text(size = rel(0.75)), 
          axis.title = element_text(size = rel(1.25)),
          axis.text.x = element_text(size = rel(0.75), angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = rel(0.75)),
          plot.title = element_text(size = rel(1.5), hjust = 0.5))
  save_ggplot(cvar, width = 4.5, height = 4.5)
  save_ggplot(cvar, width = 4.5, height = 4.5, mode = "paper") 
}

# Save data
write.csv(data, file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
