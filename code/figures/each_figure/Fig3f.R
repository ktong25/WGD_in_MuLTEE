library(dplyr)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/cluster_cell_correlation/20230926")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
conditions_color <- get_pal_colors("Set1")[c(2,1)]  # consistent with Bozdag2021, blue, red # "#377EB8" "#E41A1C"
names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
conditions_label <- c("Mixotrophic", "Anaerobic")
conditions_label_color <- conditions_color
names(conditions_label_color) <- conditions_label
lines <- as.character(1:5)
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- paste0("t", c(0, 200, 400, 600, 1000))
#evotimes2 <- gsub("t", "", evotimes, fixed = TRUE) %>% as.numeric()
ploidys <- c("2N", "4N")
ploidys_color <- get_pal_colors("Dark2")[c(1,2)]  # green, orange # "#1B9E77" "#D95F02"
names(ploidys_color) <- ploidys
#show_col(ploidys_color, labels = TRUE)

# Prepare data
# Load data summary files for evolved isolates and artificially-constructed ploidy strains
cluster_evo <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/cluster_size/20230926_all/data_summary.rds")
cluster_art <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/cluster_size/20230926_artificial/data_summary.rds")
cell_evo <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/cell_size_aspect_ratio/20230926_all/data_summary.rds")
cell_art <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/cell_size_aspect_ratio/20230926_artificial/data_summary.rds")
# Prefix column names (starting from "Count" column) with cluster or cell
colnames(cluster_evo)[which(colnames(cluster_evo) == "Count"):ncol(cluster_evo)] <- paste0("Cluster.", colnames(cluster_evo)[which(colnames(cluster_evo) == "Count"):ncol(cluster_evo)])
colnames(cluster_art)[which(colnames(cluster_art) == "Count"):ncol(cluster_art)] <- paste0("Cluster.", colnames(cluster_art)[which(colnames(cluster_art) == "Count"):ncol(cluster_art)])
colnames(cell_evo)[which(colnames(cell_evo) == "Count"):ncol(cell_evo)] <- paste0("Cell.", colnames(cell_evo)[which(colnames(cell_evo) == "Count"):ncol(cell_evo)])
colnames(cell_art)[which(colnames(cell_art) == "Count"):ncol(cell_art)] <- paste0("Cell.", colnames(cell_art)[which(colnames(cell_art) == "Count"):ncol(cell_art)])
# Combine data (use only cluster 24h data)
data_evo <- dplyr::full_join(x = cluster_evo %>% dplyr::filter(Time == "24h") %>% dplyr::ungroup() %>% dplyr::select(!c(Time, Cluster.Count)), 
                             y = cell_evo %>% dplyr::select(!Cell.Count),
                             by = c("Condition", "Line", "EvoTime"))
data_art <- dplyr::full_join(x = cluster_art %>% dplyr::filter(Time == "24h") %>% dplyr::ungroup() %>% dplyr::select(!c(Time, Cluster.Count)), 
                             y = cell_art %>% dplyr::select(!Cell.Count),
                             by = c("Condition", "Ploidy"))

# Plot cluster radius and cell volume/AR for artificial ploidy strains
# Show slope of 2N->4N cell volume/AR change that is similar between two conditions
ggplot(data_art %>% dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = conditions_label)), 
       aes(x = Cell.Mean_volume, y = Cell.Mean_AR, size = Cluster.Weighted_mean_radius, fill = Ploidy, color = Ploidy)) +
  geom_line(mapping = aes(group = Condition, color = Condition), alpha = 0.7, size = 1.25) +
  geom_point(shape = "circle filled", fill = "white", alpha = 1, stroke = 0.75) +
  geom_point(shape = "circle filled", color = "black", alpha = 0.5, stroke = 0.75) +
  scale_fill_manual(values = ploidys_color) +
  scale_color_manual(values = conditions_label_color, 
                     breaks = conditions_label) +  # somehow this is necessary to set the right order for legend texts
  scale_radius(range = c(4, 12), breaks = c(15,30,45), limits = c(15,45)) + 
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(x = expression(Cell~volume~(mu*m^3)), 
       y = "Cell aspect ratio", 
       fill = "Ploidy", 
       color = "Condition", 
       size = expression(Cluster~radius~(mu*m))) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 4)), 
         color = guide_legend(order = 2), 
         size = guide_legend(order = 3)) +
  ggplot_custom_theme2 +
  theme(legend.title = element_text(size = rel(1.25)))
save_ggplot("cell_volume_vs_AR_artificial", width = 6, height = 5)
save_ggplot("cell_volume_vs_AR_artificial", width = 6, height = 5, mode = "paper")

# Plot x vs y for evolved and artificial

# Parameters
evotimes_color <- viridis_pal(option = "viridis", begin = 0.2, end = 1, direction = -1)(4)  # "#FDE725FF" "#54C568FF" "#23888EFF" "#414487FF"
evotimes_color <- c("gray", evotimes_color)
names(evotimes_color) <- evotimes
#show_col(evotimes_color, labels = TRUE)
lines_shape <- paste(c("circle", "square", "diamond", "triangle", "triangle down"), "open")  # use open points with small size and thick stroke to get the effect of filled points
lines_shape <- c("circle open", lines_shape)
names(lines_shape) <- c("0", lines)  # change t0's Line value from NA to 0
ploidys_size <- c(1.5,3)
names(ploidys_size) <- ploidys

# Plot cell volume vs cell AR (example of the function below, may not be the most up-to-date)
# ggplot(data_evo %>% dplyr::mutate(Line = ifelse(EvoTime == "t0", "0", Line)), 
#        aes(x = Cell.Mean_volume, y = Cell.Mean_AR, color = EvoTime, shape = Line)) +
#   facet_wrap(~Condition, nrow = 1) +
#   geom_point(size = 0.01, stroke = 3, alpha = 0.7) +
#   scale_color_manual(values = evotimes_color) +
#   scale_shape_manual(values = lines_shape, breaks = lines) + 
#   geom_point(data = data_art, mapping = aes(size = Ploidy),
#              color = "red", shape = "asterisk", stroke = 0.75, alpha = 0.7) +
#   scale_size_manual(values = ploidys_size) +
#   labs(x = expression(Cell~volume~(mu*m^3)), y = "Cell aspect ratio", color = "Evo time") + 
#   guides(color = guide_legend(order = 1), shape = guide_legend(order = 2), size = guide_legend(order = 3)) +
#   ggplot_custom_theme4 +
#   theme(legend.title = element_text(size = rel(1.25)), 
#         panel.grid.major = element_line(color = "gray80", size = 0.2))

# Function for plotting
plot_xy_evo <- function(xvar, yvar, xlab, ylab, out_filename, plot_art = FALSE) {
  p <- ggplot(data_evo %>% dplyr::mutate(Line = ifelse(EvoTime == "t0", "0", Line)), 
         aes(x = .data[[xvar]], y = .data[[yvar]], color = EvoTime, shape = Line)) +
    facet_wrap(~Condition, nrow = 1) +
    geom_point(size = 0.01, stroke = 3, alpha = 0.7) +
    scale_color_manual(values = evotimes_color) +
    scale_shape_manual(values = lines_shape, breaks = lines) + 
    labs(x = xlab, y = ylab, color = "Evo time") + 
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2), size = guide_legend(order = 3)) +
    ggplot_custom_theme4 +
    theme(#panel.grid.major = element_line(color = "gray80", size = 0.2), 
          legend.title = element_text(size = rel(1.25)))
  if (plot_art) {
    p <- p + geom_point(data = data_art, mapping = aes(size = Ploidy),
                        color = "red", shape = "asterisk", stroke = 0.75, alpha = 0.7) +
      scale_size_manual(values = ploidys_size) 
  }
  save_ggplot(glue("{out_filename}_evo{ifelse(plot_art, '_artificial', '')}"), width = 8, height = 5)
  save_ggplot(glue("{out_filename}_evo{ifelse(plot_art, '_artificial', '')}"), width = 8, height = 5, mode = "paper")
  return(p)
}

# Plot cell volume vs cell AR
for (plot_art in c(F,T)) {
  plot_xy_evo(xvar = "Cell.Mean_volume", yvar = "Cell.Mean_AR", 
              xlab = expression(Cell~volume~(mu*m^3)), ylab = "Cell aspect ratio", 
              out_filename = "cell_volume_vs_AR", 
              plot_art = plot_art)
}

# Plot cell volume vs cluster radius
for (plot_art in c(F,T)) {
  plot_xy_evo(xvar = "Cell.Mean_volume", yvar = "Cluster.Weighted_mean_radius", 
              xlab = expression(Cell~volume~(mu*m^3)), ylab = expression(Cluster~radius~(mu*m)), 
              out_filename = "cell_volume_vs_cluster_radius", 
              plot_art = plot_art)
}

# Plot cell AR vs cluster radius
for (plot_art in c(F,T)) {
  plot_xy_evo(xvar = "Cell.Mean_AR", yvar = "Cluster.Weighted_mean_radius", 
              xlab = "Cell aspect ratio", ylab = expression(Cluster~radius~(mu*m)), 
              out_filename = "cell_AR_vs_cluster_radius", 
              plot_art = plot_art)
}

# Save data
write.csv(data_evo, file = "data_evo.csv", row.names = FALSE)
write.csv(data_art, file = "data_artificial.csv", row.names = FALSE)
saveRDS(data_evo, file = "data_evo.rds")
saveRDS(data_art, file = "data_artificial.rds")

# # Plot cell volume vs mean top cluster radius (for microscopic strains, thus include PA3 t600)
# ggplot(data_evo %>% dplyr::mutate(Line = ifelse(EvoTime == "t0", "0", Line)) %>%
#          dplyr::filter(!Cluster.Weighted_mean_radius > 200),  ###
#        aes(x = Cell.Mean_volume, y = 4/3 * pi * Cluster.Mean_top_radius^3, color = EvoTime, shape = Line)) +  ###
#   facet_wrap(~Condition, nrow = 1) +
#   geom_point(size = 0.01, stroke = 3, alpha = 0.7) +
#   scale_color_manual(values = evotimes_color) +
#   scale_shape_manual(values = lines_shape, breaks = lines) + 
#   geom_point(data = data_art, mapping = aes(size = Ploidy),
#              color = "red", shape = "asterisk", stroke = 0.75, alpha = 0.7) +
#   scale_size_manual(values = ploidys_size) +
#   labs(x = expression(Cell~volume~(mu*m^3)), y = expression(Cluster~volume~at~division~(mu*m^3)), color = "Evo time") +   ###
#   guides(color = guide_legend(order = 1), shape = guide_legend(order = 2), size = guide_legend(order = 3)) +
#   ggplot_custom_theme4 +
#   theme(legend.title = element_text(size = rel(1.25)), 
#         panel.grid.major = element_line(color = "gray80", size = 0.2))
# save_ggplot("cell_volume_vs_cluster_volume_evo_artificial_micro", width = 8, height = 5)
# save_ggplot("cell_volume_vs_cluster_volume_evo_artificial_micro", width = 8, height = 5, mode = "paper")

# Trend line?
# geom_smooth(method = "loess", formula = y ~ x, 
#             #method = "gam", formula = y ~ s(x, bs = "cs"), 
#             mapping = aes(x = Cell.Mean_volume, y = Cell.Mean_AR), inherit.aes = FALSE, 
#             color = "gray", linetype = "solid", size = 0.5, alpha = 0.2) +
