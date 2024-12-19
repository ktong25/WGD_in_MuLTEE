library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(glue)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggridges)
library(ggforce)
library(ggh4x)
library(ggbreak)
library(ggpubr)
library(cowplot)

setwd("~/Documents/projects/R_data_analysis/papers/WGD")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Paramaters
# Input
in_folder_global <- "~/Documents/projects/R_data_analysis"
read_data_rds <- function(rds_path) {
  readRDS(file.path(in_folder_global, rds_path)) 
}
read_data_csv <- function(csv_path) {
  read.csv(file.path(in_folder_global, csv_path), row.names = NULL, stringsAsFactors = FALSE)
}
# Output
out_fig_path <- "figures"
if (!file.exists(out_fig_path)) {dir.create(out_fig_path, recursive = TRUE)}
save_plot <- function(filename, sup = FALSE, ...) {
  formats <- c("pdf", "png")
  if (sup) {formats <- c(formats, "eps", "tiff")}
  for (format in formats) {
    if (format == "eps") {
      ggsave(paste(filename, format, sep = "."), device = cairo_ps, path = out_fig_path, ...)
    } else {
      ggsave(paste(filename, format, sep = "."), path = out_fig_path, ...)
    }
  }
}
mainfig_name <- function(n) {
  paste0("Fig", n)
}
supfig_name <- function(n) {
  paste0("Ratcliff_EDfig", n)
}
# Figure formatting
lwpt <- 0.352/0.75  # line width: size 1 = 0.75mm, print 1pt = 0.352mm
pspt <- 0.352  # point size/stroke: size 1 = 1mm
adjust_theme <- function(p) {
  p +
    theme(text = element_text(family = "Helvetica"), # Arial requires additional code to work in Mac
          plot.title = element_text(size = 7, hjust = 0.5), 
          axis.title = element_text(size = 7), 
          axis.text = element_text(size = 6), 
          axis.text.x = element_text(color = "black"), 
          axis.text.y = element_text(color = "black"), 
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6), 
          strip.text = element_text(size = 7),
          axis.line = element_line(linewidth = lwpt * 0.5), 
          axis.ticks = element_line(linewidth = lwpt * 0.5), 
          axis.ticks.length = unit(0.75, units = "mm"), 
          legend.spacing = unit(0, 'points'),
          legend.box.spacing = unit(0, units = "points"),
    )
}
panel_labels <- "auto"
panel_labels_trans <- function(labels) {
  if (panel_labels == "auto") {
    return(tolower(labels))
  } else if (panel_labels == "AUTO") {
    return(toupper(labels))
  }
}
panel_label_size <- 8  # points
cluster_radius_text <- bquote("Cluster radius (" * mu * m * ")")
cell_volume_text <- bquote("Cell volume (" * mu * m^3 * ")")
cell_AR_text <- "Cell aspect ratio"
lines <- as.character(1:5)
lines_color <- c(get_pal_colors("Dark2")[4], get_pal_colors("Paired")[c(2,4,10,8)])  # consistent with Bozdag2023, change red to magenta for colorblind-proof
names(lines_color) <- lines
strainbackgrounds_color <- c("#F0AFC8", "#E04692", "#B8D8E9", "#4D93C3", "#5CB356", 
                             "#D5C1DE", "#8764AE", "#FDCC8C", "#FF9933")
neutral_copy_num_color <- "gray95"
chr_copy_nums_color <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(9)[2:9]  #PuOr
chr_copy_nums_color[4] <- neutral_copy_num_color
names(chr_copy_nums_color) <- 1:8
two_signs_color2 <- rep("#71ADF9", 4)

##### Main figures

### Figure 1
# 1ab
p1ab <- plot_grid(NULL, NULL,
                  nrow = 1, rel_widths = c(77, 106), 
                  labels = panel_labels_trans(c("A", "B")), label_size = panel_label_size)
# 1cde
temp_theme_1cde <- function(p) {
  p$layers[[1]]$aes_params$size <- lwpt * 0.75
  p$layers[[2]]$aes_params$size <- pspt * 1.5
  p$layers[[3]]$aes_params$size <- lwpt * 0.75
  p +
    theme(legend.key.width = unit(10, "points"),
          legend.key.height = unit(5, "points"), 
          plot.margin = margin(t = 3, b = 3, r = 6, l = 6, unit = "pt")
          )
}
p1c <- read_data_rds("cluster_size/20230926_all/paper/cluster_size_evo_24h_artificial.rds") %>% adjust_theme() %>% temp_theme_1cde() +
  ylab(cluster_radius_text) +
  theme(legend.position = "none") +
  scale_color_manual(values = lines_color)
#save_plot("p1c", width = 61, height = 40, units = "mm")
p1d <- read_data_rds("cell_size_aspect_ratio/20230926_all/paper/cell_volume_evo_artificial.rds") %>% adjust_theme() %>% temp_theme_1cde() +
  ylab(cell_volume_text) +
  theme(legend.position = "none") +
  scale_color_manual(values = lines_color)
#save_plot("p1d", width = 61, height = 40, units = "mm")
p1e <- read_data_rds("cell_size_aspect_ratio/20230926_all/paper/cell_AR_evo_artificial.rds") %>% adjust_theme() %>% temp_theme_1cde()  +
  theme(legend.position = "none") +
  scale_color_manual(values = lines_color)
#save_plot("p1e", width = 61, height = 40, units = "mm")
p1cde <- plot_grid(p1c, p1d, p1e, 
                   nrow = 1, rel_widths = 1, 
                   align = "hv", axis = "lr", 
                   labels = panel_labels_trans(c("C", "D", "E")), label_size = panel_label_size)
p1cde_legend <- get_legend(p1c + 
                             labs(linetype = "Engineered ploidy") +
                             theme(legend.position = "bottom", legend.background = element_blank()))
p1cde <- plot_grid(p1cde, p1cde_legend, 
                   ncol = 1, rel_heights = c(40,5))
#save_plot("p1cde", width = 183, height = 45, units = "mm")
# 1fgh
temp_theme_1fgh <- function(p, colorbar_breaks) {
  p + 
    xlab(cell_volume_text) +
    coord_fixed(ratio = 1) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4, direction = "vertical")) +
    theme(legend.position = "right", 
          axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.x = element_text(size = 6, margin = margin(t = -2, b = 0)),
          axis.text.y = element_text(size = 6, margin = margin(l = 0, r = -2)), 
          legend.box.margin = margin(t = 0, b = 0, l = 5, r = 0, unit = "pt")
          ) +
    scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "Blues"))(100)[90:1], 
                         breaks = colorbar_breaks)
}
p1f <- read_data_rds("snowflake_physics_modeling/20231211/paper/Cluster_radius.rds") %>% adjust_theme() %>% 
  temp_theme_1fgh(seq(20,60,10)) +
  ggtitle(cluster_radius_text)
#save_plot("p1f", width = 61, height = 50, units = "mm")
p1g <- read_data_rds("snowflake_physics_modeling/20231211/paper/Cluster_ncells.rds") %>% adjust_theme() %>% 
  temp_theme_1fgh(seq(300,800,100)) +
  ggtitle("Cluster cell number")
#save_plot("p1g", width = 61, height = 50, units = "mm")
p1h <- read_data_rds("snowflake_physics_modeling/20231211/paper/Cluster_packfrac.rds") %>% adjust_theme() %>% 
  temp_theme_1fgh(seq(0.2,0.4,0.05)) +
  ggtitle("Cluster packing fraction")
#save_plot("p1h", width = 61, height = 50, units = "mm")
# 1fgh
p1fgh <- plot_grid(p1f, p1g, p1h, 
                   nrow = 1, rel_widths = 1, 
                   labels = panel_labels_trans(c("F", "G", "H")), label_size = panel_label_size)
#save_plot("p1fgh", width = 183, height = 50, units = "mm")
# combine
p1 <- plot_grid(p1ab, p1cde, p1fgh, 
                ncol = 1, rel_heights = c(75, 45, 50))
save_plot(paste0(mainfig_name(1), "_incomplete"), width = 183, height = 170, units = "mm")

### Figure 2
p2a <- read_data_rds("allele_frequency/20231223_MuLTEE_summary/paper/allele_frequency.rds") %>% adjust_theme() +
  theme(panel.border = element_rect(linewidth = lwpt * 0.5), 
        axis.line = element_blank(), 
        axis.ticks = element_line(linewidth = lwpt * 0.25), 
        strip.text.y = element_text(angle = 0))
p2a$layers[[1]]$aes_params$size <- lwpt * 0.25
p2a$layers[[2]]$aes_params$size <- lwpt * 0.5
p2a$layers[[3]]$aes_params$size <- lwpt * 0.5
p2a$layers[[4]]$aes_params$size <- pspt * 0.001
p2a$layers[[4]]$aes_params$stroke <- pspt * 0.8
#save_plot("p2a", width = 183/2, height = 50, units = "mm")
p2b <- read_data_rds("ploidy_measurement/20230606_PMPA_evo_combine/paper/ploidy_evo.rds") %>% adjust_theme()
p2b$layers[[1]]$aes_params$size <- lwpt * 0.75
p2b$layers[[2]]$aes_params$size <- pspt * 1.5
p2b <- p2b +
  theme(legend.key.width = unit(10, "points"),
        legend.key.height = unit(5, "points"), 
        panel.grid.major.y = element_line(linewidth = lwpt * 0.25), 
        axis.line = element_line(linewidth = lwpt * 0.25), 
        axis.ticks = element_line(linewidth = lwpt * 0.25)) +
  scale_color_manual(values = lines_color)
#save_plot("p2b", width = 183/2, height = 50, units = "mm")
temp_theme_2c <- function(p) {  # also used for several other ridge plots
  p$layers[[1]]$aes_params$size <- lwpt * 0.4
  p +
    theme(legend.key.size = unit(8, "points"),
          panel.grid.major.x = element_line(linewidth = lwpt * 0.25), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.x = element_text(margin = margin(r = 1)),
          axis.text.y = element_text(vjust = 0, margin = margin(r = 1))
          )
}
p2c <- read_data_rds("ploidy_measurement/20230606_PMPA_early_combine/paper/PMPA_early.rds") %>% adjust_theme() %>% temp_theme_2c() +
  theme(strip.text.y = element_text(angle = 0))
#save_plot("p2c", width = 183/2, height = 70, units = "mm")
p2d <- read_data_rds("chromosome_cnv/20231214_summary_MuLTEE/paper/copy_number_chr_heatmap_by_line.rds") %>% adjust_theme()
p2d$layers[[1]]$aes_params$size <- lwpt * 0.25
p2d <- p2d +
  theme(legend.key.size = unit(8, "points"),
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)
        ) +
  scale_fill_manual(name = "Copy\nnumber",
                    values = chr_copy_nums_color)
#save_plot("p2d", width = 183/2, height = 70, units = "mm")
# combine
p2 <- plot_grid(p2a, p2b, p2c, p2d, 
                ncol = 2, rel_widths = 1, rel_heights = c(50, 70), 
                labels = panel_labels, label_size = panel_label_size)
save_plot(mainfig_name(2), width = 183, height = 120, units = "mm")

### Figure 3
# 3ad
p3ad <- plot_grid(NULL, NULL,
                  ncol = 1, rel_heights = 1,
                  labels = panel_labels_trans(c("A", "D")), label_size = panel_label_size)
# 3bcef
temp_theme_3bceg <- function(p) {
  p$layers[[2]]$aes_params$size <- pspt * 1
  p$layers[[3]]$aes_params$linewidth <- lwpt * 0.5
  p$layers[[4]]$aes_params$size <- lwpt * 0.5
  p$layers[[4]]$aes_params$label.size <- pspt * 5
  p +
    theme(legend.position = "bottom",
          legend.key.size = unit(8, "points"), 
          plot.margin = margin(t = 6, b = 0, r = 3, l = 3, unit = "pt")
          ) +
    labs(fill = NULL, color = NULL)
}
p3b <- read_data_rds("cell_size_aspect_ratio/20230926_artificial/paper/cell_volume_ploidy.rds") %>% adjust_theme() %>% temp_theme_3bceg() +
  ylab(cell_volume_text)
#save_plot("p3b", width = 37.5, height = 40, units = "mm")
y_baseline <- 1
p3c <- read_data_rds("cell_size_aspect_ratio/20230926_artificial/paper/cell_AR_ploidy.rds") %>% adjust_theme() %>% temp_theme_3bceg()
#save_plot("p3c", width = 37.5, height = 40, units = "mm")
p3e <- read_data_rds("cluster_size/20230926_artificial/paper/cluster_size_ploidy.rds") %>% adjust_theme() %>% temp_theme_3bceg() +
  ylab(cluster_radius_text)
#save_plot("p3e", width = 37.5, height = 40, units = "mm")
p3f <- read_data_rds("cluster_cell_correlation/20230926/paper/cell_volume_vs_AR_artificial.rds") %>% adjust_theme() +
  xlab(cell_volume_text)
p3f$layers[[1]]$aes_params$size <- lwpt * 1
p3f$layers[[2]]$aes_params$stroke <- pspt * 0.75
p3f$layers[[2]]$aes_params$size <- pspt * 5
p3f$layers[[3]]$aes_params$stroke <- pspt * 0.75
p3f$layers[[3]]$aes_params$size <- pspt * 5
p3f <- p3f +
  annotate("text", x = 55, y = 1.35, label = "Mixotrophic", size = pspt * 6, hjust = 0, color = "#377EB8") +
  annotate("text", x = 105, y = 1.28, label = "Anaerobic", size = pspt * 6, hjust = 0, color = "#E41A1C") +
  guides(fill = guide_legend(order = 1), color = "none") + #guide_legend(order = 2, keywidth = unit(10, "points"), override.aes = list(shape = NA))
  labs(fill = NULL) +
  theme(legend.position = c(0.875,0.15), 
        legend.background = element_blank(), 
        legend.key.size = unit(5, "points")) +
  scale_x_continuous(breaks = c(50,100,150), limits = c(45,160)) +
  scale_y_continuous(breaks = c(1.2,1.3,1.4), limits = c(1.2,1.42))
#save_plot("p3f", width = 37.5, height = 40, units = "mm")
p3be <- plot_grid(p3b, p3e, 
                  ncol = 1, rel_heights = 1, 
                  align = "hv", axis = "tblr", 
                  labels = panel_labels_trans(c("B", "E")), label_size = panel_label_size)
#save_plot("p3be", width = 37.5, height = 80, units = "mm")
p3cf <- plot_grid(p3c, p3f, 
                  ncol = 1, rel_heights = 1, 
                  align = "hv", axis = "tblr", 
                  labels = panel_labels_trans(c("C", "F")), label_size = panel_label_size)
#save_plot("p3cf", width = 37.5, height = 80, units = "mm")
p3abcdef <- plot_grid(p3ad, p3be, p3cf,
                      nrow = 1, rel_widths = c(40,37.5,37.5))
#save_plot("p3abcdef", width = 110, height = 80, units = "mm")
# 3gj
set.seed(11)
p3g <- read_data_rds("competition_assay/20220909_2N4N_corrected_remove_PO_polish_20230928/paper/selection_rate.rds") %>% adjust_theme() %>% temp_theme_3bceg() +
  labs(fill = "Settling selection", color = "Settling selection") +
  theme(legend.justification = c(0.8,0), 
        plot.margin = margin(t = 6, b = 0, r = 6, l = 6, unit = "pt"))
#save_plot("p3g", width = 40, height = 50, units = "mm")
p3j <- read_data_rds("ploidy_measurement/20241021_MA_Sayantan/paper/MA_diff.rds") %>% adjust_theme() +
  ylab("Ploidy change\n(re-evolve without selection)") +
  theme(legend.key.size = unit(8, "points"), 
        plot.margin = margin(t = 6, b = 0, r = 6, l = 6, unit = "pt"))
#save_plot("p3j", width = 75, height = 50, units = "mm")
p3gj <- plot_grid(p3g, p3j,
                  nrow = 1, rel_widths = c(40,75), 
                  align = "hv", axis = "tb", 
                  labels = panel_labels_trans(c("G", "J")), label_size = panel_label_size)
#save_plot("p3gj", width = 110, height = 50, units = "mm")
# 3hi
temp_theme_3hi <- function(p) {
  p +
    theme(legend.margin = margin(t = 4, b = 2, l = 0, r = 0, unit = "pt"), 
          legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt"),
          strip.text.x.top = element_text(margin = margin(t = 0, b = 3, l = 0, r = 0, unit = "pt")))
}
p3h <- read_data_rds("ploidy_measurement/20230920_agar_rev_combine/paper/reversion_t0.rds") %>% adjust_theme() %>% temp_theme_2c() %>% temp_theme_3hi()
#save_plot("p3h", width = 73, height = 32, units = "mm")
p3i <- read_data_rds("ploidy_measurement/20230920_agar_rev_combine/paper/reversion_t1000.rds") %>% adjust_theme() %>% temp_theme_2c() %>% temp_theme_3hi() +
  theme(strip.text.y = element_text(angle = 0))
#save_plot("p3i", width = 73, height = 90, units = "mm")
p3hi <- plot_grid(NULL, p3h, p3i,
                  ncol = 1, rel_heights = c(8,32,90), 
                  align = "v", axis = "lr", 
                  labels = panel_labels_trans(c(NA, "H", "I")), label_size = panel_label_size)
#save_plot("p3hi", width = 73, height = 130, units = "mm")
# combine
p3_left <- plot_grid(p3abcdef, p3gj, 
                     ncol = 1, rel_heights = c(80,50))
#save_plot("p3_left", width = 110, height = 130, units = "mm")
p3 <- plot_grid(p3_left, p3hi, 
                nrow = 1, rel_widths = c(110,73))
save_plot(paste0(mainfig_name(3), "_incomplete"), width = 183, height = 130, units = "mm")

### Figure 4
# 4bcd
temp_theme_4bcd <- function(p) {
  p$layers[[1]]$aes_params$size <- lwpt * 0.5
  p$layers[[2]]$aes_params$size <- pspt * 1
  p$layers[[3]]$aes_params$size <- lwpt * 0.5
  p$layers[[4]]$aes_params$size <- lwpt * 0.5
  p$layers[[4]]$aes_params$label.size <- pspt * 5
  p
}
p4b <- read_data_rds("cluster_size/20231227_donutspread/paper/cluster_size.rds") %>% adjust_theme() %>% temp_theme_4bcd() + 
  theme(plot.margin = margin(t = 6, b = 6, l = 3, r = 6, unit = "pt"))
#save_plot("p4b", width = 80/3, height = 35, units = "mm")
p4c <- read_data_rds("cell_size_aspect_ratio/20231227_donutspread/paper/cell_volume.rds") %>% adjust_theme() %>% temp_theme_4bcd()
#save_plot("p4c", width = 80/3, height = 35, units = "mm")
p4d <- read_data_rds("cell_size_aspect_ratio/20231227_donutspread/paper/cell_AR.rds") %>% adjust_theme() %>% temp_theme_4bcd()
#save_plot("p4d", width = 80/3, height = 35, units = "mm")
p4bcd <- plot_grid(p4b, p4c, p4d,
                   nrow = 1, rel_widths = 1,
                   align = "hv", axis = "tblr", 
                   labels = panel_labels_trans(c("B", "C", "D")), label_size = panel_label_size)
#save_plot("p4bcd", width = 80, height = 35, units = "mm")
# 4ef
temp_theme_4ef <- function(p) {
  p + 
    xlab(cell_volume_text) +
    theme(legend.text = element_text(size = 5),
          legend.key.size = unit(5, "points"), 
          legend.background = element_blank())
}
colonymorphs <- c("D", "S") 
colonymorphs_size <- c(2,3)
names(colonymorphs_size) <- colonymorphs
p4e <- read_data_rds("cell_size_aspect_ratio/20231227_donutspread/paper/cell_volume_vs_AR_onlypoints_svm.rds") %>% adjust_theme() %>% temp_theme_4ef() + 
  theme(legend.position = c(0.15,0.95)) +
  scale_size_manual(values = colonymorphs_size * pspt) +
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(size = 6), 
        axis.title.y.left = element_text(margin = margin(t = 6, b = 6, l = 2.4, r = 6, unit = "pt"))
        ) + 
  labs(color = NULL, shape = NULL, size = NULL)
p4e$layers[[2]]$aes_params$size <- lwpt * 0.75
#save_plot("p4e", width = 40, height = 40, units = "mm")
p4f <- read_data_rds("cell_size_aspect_ratio/20231227_donutspread/paper/cell_volume_vs_AR_onlylines.rds") %>% adjust_theme() %>% temp_theme_4ef() + 
  theme(legend.position = "none") +
  scale_color_manual(values = two_signs_color2) 
p4f$layers[[1]]$aes_params$size <- lwpt * 0.4
#save_plot("p4f", width = 40, height = 40, units = "mm")
p4ef <- plot_grid(p4e, p4f,
                  nrow = 1, rel_widths = 1, 
                  align = "hv", axis = "tblr", 
                  labels = panel_labels_trans(c("E", "F")), label_size = panel_label_size)
#save_plot("p4ef", width = 80, height = 40, units = "mm")
# combine
p4bcdef <- plot_grid(p4bcd, p4ef,
                     ncol = 1, rel_heights = c(35,40), 
                     align = "hv", axis = "tblr")
#save_plot("p4bcdef", width = 80, height = 75, units = "mm")
p4 <- plot_grid(NULL, p4bcdef,
                nrow = 1, rel_widths = c(50,80), 
                labels = panel_labels_trans(c("A")), label_size = panel_label_size)
save_plot(paste0(mainfig_name(4), "_incomplete"), width = 130, height = 75, units = "mm")

### Figure 5
# 5abc
temp_theme_5abc <- function(p) {
  p + 
    theme(axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          legend.key.size = unit(5, "points"),
          legend.box.spacing = unit(5, units = "points"),
          axis.title.x = element_text(size = 7),  ###
          axis.text.x = element_text(size = 5),  ###
          axis.text.y = element_text(size = 5, margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")),  ###
          strip.text.y = element_text(size = 6, margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")), ###
          legend.title = element_text(size = 5),  ###
          legend.text = element_text(size = 5),  ###
    )
}
p5a <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/copy_number_chr_heatmap.rds") %>% adjust_theme() %>% temp_theme_5abc() +
  theme(legend.justification = c(0.8,0))
p5a$layers[[1]]$aes_params$size <- lwpt * 0.25
#save_plot("p5a", width = 36, height = 120, units = "mm")
p5b <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/copy_number_chr_heatmap_diff.rds") %>% adjust_theme() %>% temp_theme_5abc() +
  theme(strip.text.y = element_blank(),
        legend.justification = c(0.8,0)) +
  ggtitle("S-D difference")
p5b$layers[[1]]$aes_params$size <- lwpt * 0.25
p5b$layers[[2]]$aes_params$size <- pspt * 4  ###
#save_plot("p5b", width = 29.5, height = 120, units = "mm")
p5c <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/copy_number_chr_heatmap_summary.rds") %>% adjust_theme() %>% temp_theme_5abc() +
  ylab("# of D-S pairs") +
  theme(strip.text.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = lwpt * 0.25), 
        axis.text.y = element_text(size = 5, margin = margin(t = 0, b = 0, l = 0, r = 1, unit = "pt")),
        legend.justification = c(0.8,0)
  )
p5c$layers[[1]]$aes_params$linewidth <- lwpt * 0.5
p5c$layers[[2]]$aes_params$size <- pspt * 3.2
p5c$layers[[4]]$aes_params$linewidth <- lwpt * 0.25
p5c$layers[[5]]$aes_params$linewidth <- lwpt * 0.25
#save_plot("p5c", width = 34.5, height = 120, units = "mm")
p5abc <- plot_grid(p5a, p5b, p5c,
                   nrow = 1, rel_widths = c(36,29.5,34.5), 
                   align = "h", axis = "tb", 
                   labels = panel_labels_trans(c("A", "B", "C")), label_size = panel_label_size)
#save_plot("p5abc", width = 100, height = 120, units = "mm")
# 5d
temp_theme_5dfghi <- function(p) {
  p + 
    theme(legend.text = element_text(size = 5),
          legend.key.size = unit(5, "points"),
          legend.margin = margin(t = 2, b = 2, l = 2, r = 2, unit = "pt"), 
          panel.grid.major.y = element_blank()
    )
}
p5d <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/barplot_dChr_num.rds") %>% adjust_theme() %>% temp_theme_5dfghi() +
  labs(fill = NULL) +
  theme(legend.title = element_text(size = 5), 
        legend.position = c(0.75,0.69),
        axis.title.x = element_text(hjust = 0.9)) +
  scale_fill_manual(values = strainbackgrounds_color)
#save_plot("p5d", width = 33, height = 40, units = "mm")
# 5e
dKar_modes <- c("Same background", "Same line", "Different lines", "Other")
dKar_modes_label <- c("Same\nbackground",  "Same line", "Different lines", "Other")
dKar_modes_color <- c("#FF0090", "#33B050", "#00B0F0", "gray90")  # magenta, green, blue
names(dKar_modes_color) <- dKar_modes
p5e <- read_data_rds("cell_size_aspect_ratio/20231227_donutspread/paper/cell_volume_vs_AR_dKar.rds") %>% adjust_theme() %>% temp_theme_4ef() + 
  theme(legend.title = element_text(size = 5),
        legend.key.width = unit(10, "points"), 
        legend.key.height = unit(7, "points"), 
        legend.margin = margin(t = 6, b = 0, l = 2, r = 0, unit = "pt"), 
        legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")) +
  scale_color_manual(values = dKar_modes_color, breaks = dKar_modes[1:3], labels = dKar_modes_label[1:3]) +
  labs(color = "Convergent\nkaryotype changes") +
  guides(color = guide_legend(order = 1, override.aes = list(shape = NA)))
p5e$layers[[1]]$aes_params$size <- lwpt * 0.5
p5e$layers[[2]]$aes_params$size <- pspt * 1.25
p5e$layers[[3]]$aes_params$size <- lwpt * 0.5
p5e$layers[[4]]$aes_params$size <- pspt * 1.25
#save_plot("p5e", width = 50, height = 40, units = "mm")
# 5de
p5de <- plot_grid(p5d, p5e,
                  nrow = 1, rel_widths = c(33,50), 
                  align = "h", 
                  labels = panel_labels_trans(c("D", "E")), label_size = panel_label_size)
#save_plot("p5de", width = 83, height = 40, units = "mm")
# 5fghi
temp_theme_5fghi <- function(p) {
  p + theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                       margin = margin(t = 0, r = 0, b = 0, l = 0)), 
            axis.ticks.x = element_blank())
}
p5f <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/barplot_chr_strainbkg.rds") %>% adjust_theme() %>% temp_theme_5dfghi() %>% temp_theme_5fghi() +
  theme(legend.position = "none") +
  scale_fill_manual(values = strainbackgrounds_color)
#save_plot("p5f", width = 41.5, height = 45, units = "mm")
p5g <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/barplot_chr_sign.rds") %>% adjust_theme() %>% temp_theme_5dfghi() %>% temp_theme_5fghi() +
  theme(legend.position = "top")  # c(0.135,0.9)
#save_plot("p5g", width = 41.5, height = 45, units = "mm")
p5h <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/barplot_macro_chr_strainbkg.rds") %>% adjust_theme() %>% temp_theme_5dfghi() %>% temp_theme_5fghi() +
  theme(legend.position = "none") +
  scale_fill_manual(values = strainbackgrounds_color) +
  scale_y_continuous(limits = c(0,7), breaks = seq(0,7,2), expand = expansion(mult = c(0,0.05)))
#save_plot("p5h", width = 41.5, height = 35, units = "mm")
p5i <- read_data_rds("chromosome_cnv/20231230_donutspread/paper/barplot_macro_chr_sign.rds") %>% adjust_theme() %>% temp_theme_5dfghi() %>% temp_theme_5fghi() +
  theme(legend.position = "top") + # c(0.275,0.95)
  scale_y_continuous(limits = c(0,7), breaks = seq(0,7,2), expand = expansion(mult = c(0,0.05)))
#save_plot("p5i", width = 41.5, height = 35, units = "mm")
p5fh <- plot_grid(p5f + theme(plot.margin = margin(t = 15, r = 6, b = 6, l = 6, unit = "pt")), 
                  p5h + theme(plot.margin = margin(t = 15, r = 6, b = 6, l = 6, unit = "pt")), 
                  ncol = 1, rel_heights = c(45,35), 
                  align = "v", axis = "lr",
                  labels = panel_labels_trans(c("F", "H")), label_size = panel_label_size)
#save_plot("p5fh", width = 41.5, height = 80, units = "mm")
p5gi <- plot_grid(p5g, p5i, 
                  ncol = 1, rel_heights = c(45,35), 
                  align = "v", axis = "lr",
                  labels = panel_labels_trans(c("G", "I")), label_size = panel_label_size)
#save_plot("p5gi", width = 41.5, height = 80, units = "mm")
p5fghi <- plot_grid(p5fh, p5gi, 
                    nrow = 1, rel_widths = c(41.5,41.5), 
                    align = "h", axis = "tb")
#save_plot("p5fghi", width = 83, height = 80, units = "mm")
# combine
p5_right <- plot_grid(p5de, p5fghi,
                      ncol = 1, rel_heights = c(40,80))
#save_plot("p5_right", width = 83, height = 120, units = "mm")
p5 <- plot_grid(p5abc, p5_right,
                nrow = 1, rel_widths = c(100,83))
save_plot(paste0(mainfig_name(5), "_incomplete"), width = 183, height = 120, units = "mm")


##### Sup figures

### Sup figure 1
sp1 <- plot_grid(NULL, NULL, 
                 ncol = 1, rel_heights = 1, 
                 labels = panel_labels_trans(c("A", "B")), label_size = panel_label_size)
save_plot(paste0(supfig_name(1), "_incomplete"), width = 183, height = 210, units = "mm")

### Sup figure 2
# ab
temp_theme_sp2ab <- function(p) {
  p$layers[[1]]$aes_params$size <- lwpt * 0.25
  p$layers[[2]]$aes_params$size <- lwpt * 0.25
  p +
    ylab(cluster_radius_text) +
    theme(legend.position = "bottom", 
          legend.key.size = unit(8, units = "points"), 
          panel.grid.major.y = element_line(linewidth = lwpt * 0.25), 
          strip.text.y = element_text(angle = 0))
}
sp2a <- read_data_rds("cluster_size/20230926_all/paper/cluster_size_evo_distribution.rds") %>% adjust_theme() %>% temp_theme_sp2ab()
#save_plot("sp2a", width = 183/2, height = 140, units = "mm")
sp2b <- read_data_rds("cluster_size/20230926_all/paper/cluster_size_evo_distribution_weighted.rds") %>% adjust_theme() %>% temp_theme_sp2ab()
sp2b$layers[[3]]$aes_params$size <- pspt * 2
sp2b$layers[[3]]$aes_params$stroke <- pspt * 0.75
#save_plot("sp2b", width = 183/2, height = 140, units = "mm")
# cd
temp_theme_sp2cd <- function(p) {
  p$layers[[1]]$aes_params$linewidth <- lwpt * 0.25
  p$layers[[2]]$aes_params$linewidth <- lwpt * 0.25
  p$layers[[3]]$aes_params$linewidth <- lwpt * 0.25
  p +
    theme(panel.grid.major.y = element_line(linewidth = lwpt * 0.25), 
          strip.text.y = element_text(angle = 0))
}
sp2c <- read_data_rds("cell_size_aspect_ratio/20230926_all/paper/cell_volume_distribution.rds") %>% adjust_theme() %>% temp_theme_sp2cd() +
  ylab(cell_volume_text)
#save_plot("sp2c", width = 183/2, height = 100, units = "mm")
sp2d <- read_data_rds("cell_size_aspect_ratio/20230926_all/paper/cell_AR_distribution.rds") %>% adjust_theme() %>% temp_theme_sp2cd()
#save_plot("sp2d", width = 183/2, height = 100, units = "mm")
# combine
sp2ac <- plot_grid(sp2a, sp2c, 
                   ncol = 1, rel_heights = c(140,100), 
                   align = "v", axis = "lr",
                   labels = panel_labels_trans(c("A", "C")), label_size = panel_label_size)
#save_plot("sp2ac", width = 183/2, height = 240, units = "mm")
sp2bd <- plot_grid(sp2b, sp2d, 
                   ncol = 1, rel_heights = c(140,100), 
                   align = "v", axis = "lr",
                   labels = panel_labels_trans(c("B", "D")), label_size = panel_label_size)
#save_plot("sp2bd", width = 183/2, height = 240, units = "mm")
sp2 <- plot_grid(sp2ac, sp2bd, 
                 nrow = 1, rel_widths = 1, 
                 align = "h", axis = "tb")
save_plot(supfig_name(2), width = 183, height = 240, units = "mm")

### Sup figure 3
sp3b <- read_data_rds("ploidy_measurement/20220923_trial/paper/ploidy_controls.rds") %>% adjust_theme() +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid.major.x = element_line(linewidth = lwpt * 0.5), 
        panel.grid.minor.x = element_line(linewidth = lwpt * 0.25))
sp3b$layers[[1]]$aes_params$size <- lwpt * 0.5
sp3b$layers[[2]]$aes_params$size <- lwpt * 0.75
sp3b$layers[[3]]$aes_params$size <- lwpt * 0.75
#save_plot("sp3b", width = 50, height = 75, units = "mm")
sp3 <- plot_grid(NULL, sp3b, 
                 nrow = 1, rel_widths = c(100, 50), 
                 labels = panel_labels_trans(c("A", "B")), label_size = panel_label_size)
save_plot(paste0(supfig_name(3), "_incomplete"), width = 150, height = 75, units = "mm")

### Sup figure 4
temp_theme_sp4 <- function(p) {
  p$layers[[1]]$aes_params$size <- pspt * 0.001
  p$layers[[1]]$aes_params$stroke <- pspt * 0.4
  p$layers[[2]]$aes_params$size <- lwpt * 0.25
  p$layers[[3]]$aes_params$size <- lwpt * 0.25
  p$layers[[4]]$aes_params$size <- lwpt * 0.75
  p + 
    theme(strip.text.x = element_text(size = 5),
          strip.text.y = element_text(hjust = 0), ###
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(linewidth = lwpt * 0.25), 
          axis.line = element_blank(),
          axis.ticks = element_line(size = lwpt * 0.25))
}
sp4a <- read_data_rds("chromosome_cnv/20231214_summary_MuLTEE/paper/copy_num_bin_PM_by_line.rds") %>% adjust_theme() %>% temp_theme_sp4() +
  ggtitle("PM")
#save_plot("sp4a", width = 183/2, height = 240, units = "mm")
sp4b <- read_data_rds("chromosome_cnv/20231214_summary_MuLTEE/paper/copy_num_bin_PA_by_line.rds") %>% adjust_theme() %>% temp_theme_sp4() +
  ggtitle("PA")
#save_plot("sp4b", width = 183/2, height = 240, units = "mm")
sp4 <- plot_grid(sp4a, sp4b, 
                 nrow = 1, rel_widths = 1, 
                 align = "h", axis = "tb", 
                 labels = panel_labels_trans(c("A", "B")), label_size = panel_label_size)
save_plot(paste0(supfig_name(4), "_incomplete"), width = 183, height = 240, units = "mm")

### Sup figure 5
# bc
temp_theme_sp5bc <- function(p) {
  p$layers[[1]]$aes_params$size <- lwpt * 0.25
  p$layers[[2]]$aes_params$size <- lwpt * 0.25
  p +
    ylab(cluster_radius_text) +
    theme(legend.key.size = unit(8, units = "points"), 
          panel.grid.major.y = element_line(linewidth = lwpt * 0.25))
}
sp5b <- read_data_rds("cluster_size/20230926_artificial/paper/cluster_size_ploidy_distribution.rds") %>% adjust_theme() %>% temp_theme_sp5bc()
sp5b$layers[[1]]$aes_params$linewidth <-  lwpt * 0.25
sp5b$layers[[2]]$aes_params$linewidth <-  lwpt * 0.25
#save_plot("sp5b", width = 150/2, height = 50, units = "mm")
sp5c <- read_data_rds("cluster_size/20230926_artificial/paper/cluster_size_ploidy_distribution_weighted.rds") %>% adjust_theme() %>% temp_theme_sp5bc()
sp5c$layers[[3]]$aes_params$size <- pspt * 2
sp5c$layers[[3]]$aes_params$stroke <- pspt * 0.75
#save_plot("sp5c", width = 150/2, height = 50, units = "mm")
# de
temp_theme_sp5de <- function(p) {
  p$layers[[1]]$aes_params$linewidth <- lwpt * 0.25
  p$layers[[2]]$aes_params$linewidth <- lwpt * 0.25
  p$layers[[3]]$aes_params$linewidth <- lwpt * 0.25
  p +
    theme(legend.key.size = unit(8, units = "points"), 
          panel.grid.major.y = element_line(linewidth = lwpt * 0.25))
}
sp5d <- read_data_rds("cell_size_aspect_ratio/20230926_artificial/paper/cell_volume_ploidy_distribution.rds") %>% adjust_theme() %>% temp_theme_sp5de() +
  ylab(cell_volume_text)
#save_plot("sp5d", width = 150/2, height = 50, units = "mm")
sp5e <- read_data_rds("cell_size_aspect_ratio/20230926_artificial/paper/cell_AR_ploidy_distribution.rds") %>% adjust_theme() %>% temp_theme_sp5de()
#save_plot("sp5e", width = 150/2, height = 50, units = "mm")
# fg
temp_theme_sp5fg <- function(p) {
  p$layers[[1]]$aes_params$size <- lwpt * 0.75
  p$layers[[2]]$aes_params$size <- pspt * 1.5
  p$layers[[3]]$aes_params$size <- lwpt * 0.75
  p +
    theme(legend.key.width = unit(10, "points"),
          legend.key.height = unit(5, "points")
    ) +
  labs(linetype = "Engineered\nploidy")
}
sp5f <- read_data_rds("cluster_size/20231211_Bozdag2023_vs_2N4N/paper/cluster_size_evo_pop_artificial.rds") %>% adjust_theme() %>% temp_theme_sp5fg() +
  ylab(cluster_radius_text) +
  scale_color_manual(values = as.character(lines_color))
#save_plot("sp5f", width = 150/2, height = 40, units = "mm")
sp5g <- read_data_rds("cluster_size/20231211_Bozdag2023_vs_2N4N/paper/cell_AR_evo_pop_artificial.rds") %>% adjust_theme() %>% temp_theme_sp5fg() +
  scale_color_manual(values = as.character(lines_color))
#save_plot("sp5g", width = 150/2, height = 40, units = "mm")
# combine
sp5bdf <- plot_grid(sp5b, sp5d, sp5f,
                    ncol = 1, rel_heights = c(50,50,40), 
                    align = "v", axis = "lr",
                    labels = panel_labels_trans(c("B", "D", "F")), label_size = panel_label_size)
#save_plot("sp5bdf", width = 150/2, height = 140, units = "mm")
sp5ceg <- plot_grid(sp5c, sp5e, sp5g,
                    ncol = 1, rel_heights = c(50,50,40), 
                    align = "v", axis = "lr",
                    labels = panel_labels_trans(c("C", "E", "G")), label_size = panel_label_size)
#save_plot("sp5ceg", width = 150/2, height = 140, units = "mm")
sp5_bottom <- plot_grid(sp5bdf, sp5ceg, 
                        nrow = 1, rel_widths = 1)
#save_plot("sp5_bottom", width = 150, height = 140, units = "mm")
sp5 <- plot_grid(NULL, sp5_bottom, 
                 ncol = 1, rel_heights = c(75,140), 
                 labels = panel_labels_trans(c("A")), label_size = panel_label_size)
save_plot(paste0(supfig_name(5), "_incomplete"), width = 150, height = 215, units = "mm")

### Sup figure 6
# b
sp6b <- read_data_rds("competition_assay/20220909_2N4N_corrected_remove_PO/paper/cs_area_vs_top5_cutoff_merge.rds") %>% adjust_theme() +
  ylab(bquote("Mean area of the five largest cells detected (" * mu * m^2 * ")")) +
  theme(legend.key.size = unit(8, units = "points"), 
        panel.border = element_rect(linewidth = lwpt * 0.5), 
        axis.line = element_blank())
sp6b$layers[[1]]$aes_params$size <- pspt * 0.25
sp6b$layers[[2]]$aes_params$size <- lwpt * 0.75
#save_plot("sp6b", width = 90, height = 70, units = "mm")
# combine
sp6 <- plot_grid(NULL, sp6b, 
                 nrow = 1, rel_widths = c(30,90), 
                 labels = panel_labels_trans(c("A", "B")), label_size = panel_label_size)
save_plot(paste0(supfig_name(6), "_incomplete"), width = 120, height = 70, units = "mm")

### Sup figure 7
# b
convert_sp7b <- function(evotime) {
  sp7b <- read_data_rds(glue("cluster_size/20231201_agar/paper/reversion_size_{evotime}.rds")) %>% adjust_theme() %>% temp_theme_2c() +
    xlab(cluster_radius_text) +
    theme(strip.text.y = element_text(angle = 0), 
          axis.text.y = element_text(margin = margin(r = -4)))
  sp7b$layers[[2]]$aes_params$linewidth <- lwpt * 1
  sp7b$layers[[3]]$aes_params$size <- lwpt * 0.5
  sp7b
}
sp7b <- plot_grid(convert_sp7b("t0") + theme(plot.margin = margin(t = 6, b = 3, l = 6, r = 6, unit = "pt")), 
                  convert_sp7b("t1000") + theme(plot.margin = margin(t = 3, b = 6, l = 6, r = 6, unit = "pt")), 
                  ncol = 1, rel_heights = c(34,86), 
                  align = "v", axis = "lr",
                  labels = panel_labels_trans(c("B")), label_size = panel_label_size)
#save_plot("sp7b", width = 100, height = 120, units = "mm")
# combine
sp7 <- plot_grid(NULL, sp7b, 
                 ncol = 1, rel_heights = c(55,120), 
                 labels = panel_labels_trans(c("A")), label_size = panel_label_size)
save_plot(paste0(supfig_name(7), "_incomplete"), width = 100, height = 175, units = "mm")

### Sup figure 8
sp8a <- read_data_rds("ploidy_measurement/20241021_MA_Sayantan/paper/MA.rds") %>% adjust_theme() %>% temp_theme_2c() +
  theme(legend.key.size = unit(8, "points"))
sp8 <- plot_grid(sp8a, 
                 ncol = 1, rel_heights = 1, 
                 labels = panel_labels_trans(c("A")), label_size = panel_label_size)
save_plot(paste0(supfig_name(8), "_incomplete"), width = 100, height = 120, units = "mm")

### Sup figure 9
# ab
sp9a <- readRDS("paper/barplot_mut-mode_allele-copy-num-change_nolegend.rds") +
  xlab("Change of corrected\nallele frequency") +
  theme(axis.title.y = element_text(margin = margin(r = 0, unit = "pt"), vjust = -1), 
        axis.title.x = element_text(margin = margin(r = 0, unit = "pt"), vjust = 1))
#save_plot("sp9a", width = 40, height = 70, units = "mm")
dACNs_range <- 2:-2
sp9b <- read_data_rds("mutations/20231230_donutspread/paper/barplot_mut-mode_chr-sign_allele-copy-num-change.rds") %>% adjust_theme() +
  xlab("Change of chromosome\ncopy number") +
  theme(legend.key.size = unit(8, "points"),
        legend.margin = margin(t = 0, b = 0, l = 3, r = 1, unit = "pt"))
#save_plot("sp9b", width = 73, height = 70, units = "mm")
sp9ab <- plot_grid(sp9a, sp9b,
                   nrow = 1, rel_widths = c(40,73), 
                   #align = "h", axis = "tb", 
                   labels = panel_labels_trans(c("A", "B")), label_size = panel_label_size)
#save_plot("sp9ab", width = 113, height = 70, units = "mm")
# cd
temp_theme_sp9cd <- function(p) {
  p +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(size = 5), 
          legend.title = element_text(size = 6),
          legend.key.size = unit(6, "points"))
}
sp9c <- read_data_rds("mutations/20231230_donutspread/paper/barplot_strainrep_mut-mode_impact.rds") %>% adjust_theme() %>% temp_theme_sp9cd()
#save_plot("sp9c", width = 70, height = 100, units = "mm")
sp9d <- read_data_rds("mutations/20231230_donutspread/paper/barplot_donut_mut_impact_perc.rds") %>% adjust_theme() %>% temp_theme_sp9cd()
#save_plot("sp9d", width = 70, height = 45, units = "mm")
sp9cd <- plot_grid(sp9c + theme(legend.position = "none"), sp9d,
                   ncol = 1, rel_heights = c(95,45), 
                   align = "v", axis = "lr", 
                   labels = panel_labels_trans(c("C", "D")), label_size = panel_label_size)
#save_plot("sp9cd", width = 70, height = 140, units = "mm")
# ef
sp9e <- read_data_rds("mutations/20231230_donutspread/paper/barplot_mut-mode_impact-perc.rds") %>% adjust_theme() +
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(5, "points")) +
  labs(fill = "Mutation\nimpact")
sp9e$layers[[2]]$aes_params$size <- lwpt * 0.5
sp9e$layers[[2]]$aes_params$label.size <- pspt * 5
sp9e$layers[[3]]$aes_params$size <- lwpt * 0.5
sp9e$layers[[3]]$aes_params$label.size <- pspt * 5
sp9e$layers[[4]]$aes_params$linewidth <- lwpt * 1
#save_plot("sp9e", width = 51, height = 70, units = "mm")
sp9f <- read_data_rds("mutations/20231230_donutspread/paper/mutation-highmoderate-perc_incdec-vs-donut_sim.rds") %>% adjust_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        legend.text = element_text(size = 5), 
        legend.margin = margin(t = 0, b = 0, l = 3, r = 1, unit = "pt"),
        legend.key.size = unit(5, "points")) +
  scale_color_manual(values = strainbackgrounds_color)
sp9f$layers[[1]]$aes_params$size <- pspt * 1.5
sp9f$layers[[2]]$aes_params$size <- lwpt * 1
sp9f$layers[[3]]$aes_params$size <- lwpt * 0.75
#save_plot("sp9f", width = 62, height = 70, units = "mm")
sp9ef <- plot_grid(sp9e, sp9f,
                   nrow = 1, rel_widths = c(51,62), 
                   labels = panel_labels_trans(c("E", "F")), label_size = panel_label_size)
#save_plot("sp9ef", width = 113, height = 70, units = "mm")
# combine
sp9_left <- plot_grid(sp9ab, sp9ef,
                      ncol = 1, rel_heights = 1)
#save_plot("sp9_left", width = 113, height = 140, units = "mm")
sp9 <- plot_grid(sp9_left, sp9cd,
                 nrow = 1, rel_widths = c(113,70))
save_plot(paste0(supfig_name(9), "_incomplete"), width = 183, height = 140, units = "mm")

