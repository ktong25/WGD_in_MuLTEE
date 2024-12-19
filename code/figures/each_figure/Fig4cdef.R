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

setwd("~/Documents/projects/R_data_analysis/cell_size_aspect_ratio/20231227_donutspread")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PA")
lines <- as.character(1:5)
#lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
#names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- c("t600", "t1000")
#evotimes2 <- gsub("t", "", evotimes, fixed = TRUE) %>% as.numeric()
reps <- paste0("rep", 1:3)
reps_label <- gsub("rep", "", reps, fixed = TRUE)
reps_shape <- c("circle", "square", "triangle")
names(reps_shape) <- reps_label
colonymorphs <- c("D", "S") 
colonymorphs_color <- c("#08CAD1", "#9D94FF") #cyan-purple #green-purple c("#2F9698", "#A7226B")  #blue-red c("#47ABD8", "#FF4242"), light c("#9BD2F2", "#ECA6A6") 
names(colonymorphs_color) <- colonymorphs
#show_col(colonymorphs_color, labels = TRUE)
colonymorphs_shape <- c("circle", "diamond")
names(colonymorphs_shape) <- colonymorphs
colonymorphs_size <- c(2,3)
names(colonymorphs_size) <- colonymorphs

# Load data
# Input file name: e.g., "PA_1_t600_D_rep1_001_crop1_Results.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1)
  colnames(df) <- gsub("BF1.", "", colnames(df), fixed = TRUE)
  measure_vars <- colnames(df)
  # Extract metadata
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  df$File <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Strain <- paste(metadata[1:5], collapse = "_")
  df$Condition <- metadata[1]
  df$Line <- metadata[2]
  df$EvoTime <- metadata[3]
  df$ColonyMorph <- metadata[4]
  df$Replicate <- metadata[5]
  df$ImageID <- as.integer(metadata[6])
  df$Crop <- as.integer(substr(metadata[7], 5, stri_length(metadata[7])))
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv)
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
# Format data
data <- data %>%
  tidyr::unite(col = "CondLine", Condition, Line, sep = "", remove = FALSE) %>%
  tidyr::unite(col = "StrainBackground", Condition, Line, EvoTime, sep = "_", remove = FALSE) %>%
  tidyr::unite(col = "StrainBackgroundRep", StrainBackground, Replicate, sep = "_", remove = FALSE) %>%
  dplyr::mutate(StrainBackground = StrainBackground %>%
                  gsub("_t", " t", ., fixed = TRUE) %>%
                  gsub("_", "", ., fixed = TRUE)) %>%
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph, ImageID, Crop)
condlines <- unique(data$CondLine)
strainbackgrounds <- unique(data$StrainBackground)
strainbackgroundreps <- unique(data$StrainBackgroundRep)
data <- data %>%
  dplyr::mutate(CondLine = factor(CondLine, levels = condlines), 
                StrainBackground = factor(StrainBackground, levels = strainbackgrounds),
                StrainBackgroundRep = factor(StrainBackgroundRep, levels = strainbackgroundreps), 
                Line = factor(Line, levels = lines), 
                EvoTime = factor(EvoTime, levels = evotimes), 
                Replicate = factor(Replicate, levels = reps), 
                ColonyMorph = factor(ColonyMorph, levels = colonymorphs)
  )
# Check
length(unique(data$Strain))
unique(data$Strain)
summary(data)

# Calculate cell volume
data$Volume <- 4/3 * pi * (data$Major/2) * (data$Minor/2)^2
summary(data)

# Check all files (image-crops)
data_summary_file <- data %>%
  dplyr::group_by(File) %>%
  dplyr::summarize(Count = n()) #%>%
  #dplyr::arrange(Count)

# Aggregate files of each strain
# Calculate mean/median/sd of cell volume/AR
id_vars <- c("Strain", "CondLine", "StrainBackground", "StrainBackgroundRep", 
             "Condition", "Line", "EvoTime", "Replicate", "ColonyMorph")
data_summary <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarise(Count = n(), 
                   Mean_volume = mean(Volume), 
                   Mean_AR = mean(AR), 
                   Median_volume = median(Volume), 
                   Median_AR = median(AR), 
                   SD_volume = sd(Volume), 
                   SD_AR = sd(AR)) %>%
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph)
summary(data_summary)

# Calculate mean of D and S
data_summary_mean <- data_summary %>%
  dplyr::group_by(ColonyMorph) %>%
  dplyr::summarize(Mean_volume.Mean = mean(Mean_volume), 
                   Mean_AR.Mean = mean(Mean_AR))

##### Plot (show the sample information for each data point)

# Set color palette for significance labels
sign_padjs <- c(paste0("+ ", c("****", "***", "**", "*")), 
                "ns", 
                paste0("- ", c("*", "**", "***", "****")))
sign_padjs_color <- c(colorRampPalette(brewer.pal(9, "Reds"))(100)[seq(90,45,-15)],
                      "gray",
                      colorRampPalette(brewer.pal(9, "Blues"))(100)[seq(45,90,15)])
names(sign_padjs_color) <- sign_padjs
#show_col(sign_padjs_color, labels = TRUE)

# Perform statistical tests (D vs S for each replicate)
# Welch's t-test with Benjamini-Hochberg ("BH" or its alias "fdr") correction
data_stat <- list()
# Cell volume
data_stat[["Volume"]] <- data %>%
  dplyr::group_by(StrainBackgroundRep) %>%
  rstatix::t_test(Volume ~ ColonyMorph, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  dplyr::mutate(SignPadj = ifelse(p.adj.signif == "ns", "ns", 
                                  paste(ifelse(estimate2 > estimate1, "+", "-"), p.adj.signif, sep = " ")), 
                SignPadj = factor(SignPadj, levels = sign_padjs))
data_summary <- dplyr::left_join(data_summary, 
                                 data_stat[["Volume"]] %>% 
                                   dplyr::select(StrainBackgroundRep, SignPadj) %>%
                                   dplyr::rename(Mean_volume.SignPadj = SignPadj), 
                                 by = "StrainBackgroundRep")
# Cell AR
data_stat[["AR"]] <- data %>%
  dplyr::group_by(StrainBackgroundRep) %>%
  rstatix::t_test(AR ~ ColonyMorph, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  dplyr::mutate(SignPadj = ifelse(p.adj.signif == "ns", "ns", 
                                  paste(ifelse(estimate2 > estimate1, "+", "-"), p.adj.signif, sep = " ")), 
                SignPadj = factor(SignPadj, levels = sign_padjs))
data_summary <- dplyr::left_join(data_summary, 
                                 data_stat[["AR"]] %>% 
                                   dplyr::select(StrainBackgroundRep, SignPadj) %>%
                                   dplyr::rename(Mean_AR.SignPadj = SignPadj), 
                                 by = "StrainBackgroundRep")

# Summarize the changes of both cell volume and cell AR
# Only distinguish ns, + significant, - significant, thus 9 possible cases
data_summary <- data_summary %>%
  dplyr::mutate(Mean_volume_AR.Signs = paste(
    ifelse(Mean_volume.SignPadj == "ns", "ns", substr(Mean_volume.SignPadj, 1, 1)), 
    ifelse(Mean_AR.SignPadj == "ns", "ns", substr(Mean_AR.SignPadj, 1, 1)), 
    sep = " "
  ))
data_summary_signs <- data_summary %>%
  dplyr::group_by(Mean_volume_AR.Signs) %>%
  dplyr::summarize(Count = n() / 2)
# - - 17 ~ red
# ns - 3 ~ orange
# + - 4 ~ gold
# - + 1 ~ green
two_signs <- c("- -", "ns -", "+ -", "- +")
data_summary_signs <- data_summary_signs %>% 
  dplyr::arrange(factor(Mean_volume_AR.Signs, levels = two_signs)) %>%
  dplyr::mutate(Mean_volume_AR.Signs = factor(Mean_volume_AR.Signs, levels = two_signs))
# Set color palette for V-AR signs
two_signs_color <- c("#FF0100", "#EF6A18", "#F1CF40", "#78C43E") #c("#fc2847", "#ffa343", "#fdfc74", "#71bc78")  #c("#FF4242", "#EF6A18", "#F1CF40", "#78C43E")
names(two_signs_color) <- two_signs
#show_col(two_signs_color, labels = TRUE)
two_signs_count <- glue("{data_summary_signs$Mean_volume_AR.Signs} ({data_summary_signs$Count})") %>% as.character()

# Plot
# Cell volume
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = ColonyMorph, y = Mean_volume, shape = Replicate2, group = Replicate2)) + ###
  facet_wrap(~StrainBackground, nrow = 3) +
  geom_line(mapping = aes(color = Mean_volume.SignPadj),  ### 
            size = 0.75, alpha = 0.8) +
  geom_point(size = 2, color = "gray", alpha = 1) +
  scale_color_manual(values = sign_padjs_color, labels = sign_padjs, drop = FALSE) +
  scale_shape_manual(values = reps_shape) + 
  scale_y_continuous(n.breaks = 5) +  ###
  labs(x = NULL, y = expression(Cell~volume~(mu*m^3)), ###
       color = expression("Sign P"[adj]), shape = "Replicate") +  
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cell_volume_showinfo", width = 8, height = 7)
save_ggplot("cell_volume_showinfo", width = 8, height = 7, mode = "paper")
# Cell AR
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = ColonyMorph, y = Mean_AR, shape = Replicate2, group = Replicate2)) + ###
  facet_wrap(~StrainBackground, nrow = 3) +
  geom_line(mapping = aes(color = Mean_AR.SignPadj),  ### 
            size = 0.75, alpha = 0.8) +
  geom_point(size = 2, color = "gray", alpha = 1) +
  scale_color_manual(values = sign_padjs_color, labels = sign_padjs, drop = FALSE) +
  scale_shape_manual(values = reps_shape) + 
  scale_y_continuous(n.breaks = 5) +  ###
  labs(x = NULL, y = "Cell aspect ratio",  ###
       color = expression("Sign P"[adj]), shape = "Replicate") + 
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cell_AR_showinfo", width = 8, height = 7)
save_ggplot("cell_AR_showinfo", width = 8, height = 7, mode = "paper")
# Cell volume vs cell AR
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = Mean_volume, y = Mean_AR, fill = ColonyMorph, shape = Replicate2)) +
  facet_wrap(~StrainBackground, nrow = 3) +
  geom_line(mapping = aes(group = Replicate2, color = Mean_volume_AR.Signs),  ###
            size = 0.75, alpha = 0.8) +
  geom_point(size = 2, alpha = 1, color = "black", stroke = 0.25) +
  scale_fill_manual(values = colonymorphs_color) +
  scale_shape_manual(values = paste0(reps_shape, " filled")) + 
  scale_color_manual(values = two_signs_color, breaks = two_signs, labels = two_signs) +
  scale_x_continuous(n.breaks = 5) +  ###
  scale_y_continuous(n.breaks = 5) +  ###
  labs(x = expression(Cell~volume~(mu*m^3)), y = "Cell aspect ratio", ###
       fill = "Colony", color = expression(Delta * V ~ Delta * AR), shape = "Replicate") +  ###
  guides(fill = guide_legend(order = 1, override.aes = list(shape=21)), 
         shape = guide_legend(order = 2), 
         color = guide_legend(order = 3)) +
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cell_volume_vs_AR_showinfo", width = 8, height = 7)
save_ggplot("cell_volume_vs_AR_showinfo", width = 8, height = 7, mode = "paper")

##### Plot (not show the sample information for each data point)

# Perform statistical tests (D vs S for all strain backgrounds combined)
# Paired t-test
data_summary_stat <- list()
# Cell volume
data_summary_stat[["Volume"]] <- data_summary %>%
  dplyr::ungroup() %>%
  rstatix::t_test(Mean_volume ~ ColonyMorph, paired = TRUE, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  add_xy_position(x = "ColonyMorph", fun = "max") %>%
  dplyr::mutate(y.position = y.position * 1.05)
# Cell AR
data_summary_stat[["AR"]] <- data_summary %>%
  dplyr::ungroup() %>%
  rstatix::t_test(Mean_AR ~ ColonyMorph, paired = TRUE, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  add_xy_position(x = "ColonyMorph", fun = "max") %>%
  dplyr::mutate(y.position = y.position * 1.05)

# Plot
# Cell volume
ggplot(data_summary, aes(x = ColonyMorph, y = Mean_volume, color = ColonyMorph)) +  ###
  geom_line(mapping = aes(group = StrainBackgroundRep), 
            size = 0.5, alpha = 0.8, color = "gray") +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_boxplot(width = 0.5, size = 0.5, outlier.shape = NA, fill = NA) +
  stat_pvalue_manual(data_summary_stat[["Volume"]], label = "p = {p}",  ###
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_color_manual(values = colonymorphs_color) +
  scale_y_continuous(n.breaks = 6, expand = expansion(mult = c(0.05, 0.1))) +  ###
  labs(x = NULL, y = expression(Cell~volume~(mu*m^3))) +  ###
  ggplot_custom_theme4 +
  theme(legend.position = "none")
save_ggplot("cell_volume", width = 4, height = 4)
save_ggplot("cell_volume", width = 4, height = 4, mode = "paper")
# Cell AR
ggplot(data_summary, aes(x = ColonyMorph, y = Mean_AR, color = ColonyMorph)) +  ###
  geom_line(mapping = aes(group = StrainBackgroundRep), 
            size = 0.5, alpha = 0.8, color = "gray") +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_boxplot(width = 0.5, size = 0.5, outlier.shape = NA, fill = NA) +
  stat_pvalue_manual(data_summary_stat[["AR"]], label = "p = {p}",  ###
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_color_manual(values = colonymorphs_color) +
  scale_y_continuous(n.breaks = 6, expand = expansion(mult = c(0.05, 0.1))) +  ###
  labs(x = NULL, y = "Cell aspect ratio") +  ###
  ggplot_custom_theme4 +
  theme(legend.position = "none")
save_ggplot("cell_AR", width = 4, height = 4)
save_ggplot("cell_AR", width = 4, height = 4, mode = "paper")
# Cell volume vs AR
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = Mean_volume, y = Mean_AR, fill = ColonyMorph, shape = ColonyMorph)) +
  geom_line(mapping = aes(group = StrainBackgroundRep, color = Mean_volume_AR.Signs),  ###
            size = 0.6, alpha = 0.8) +
  geom_point(size = 2, alpha = 1, color = "black", stroke = 0.25) +
  scale_fill_manual(values = colonymorphs_color) +
  scale_shape_manual(values = paste0(colonymorphs_shape, " filled")) + 
  scale_color_manual(values = two_signs_color, breaks = two_signs, labels = two_signs_count) +
  scale_x_continuous(n.breaks = 5) +  ###
  scale_y_continuous(n.breaks = 5) +  ###
  labs(x = expression(Cell~volume~(mu*m^3)), y = "Cell aspect ratio", ###
       fill = "Colony", shape = "Colony", color = expression(Delta * V ~ Delta * AR)) +  ### ∆V ∆AR
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cell_volume_vs_AR", width = 5, height = 4)
save_ggplot("cell_volume_vs_AR", width = 5, height = 4, mode = "paper")

# Cell volume vs AR (lines only)
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = Mean_volume, y = Mean_AR)) +
  geom_line(mapping = aes(group = StrainBackgroundRep, color = Mean_volume_AR.Signs),  ###
            size = 0.6, alpha = 0.8) +
  scale_color_manual(values = two_signs_color, breaks = two_signs, labels = two_signs_count) +
  scale_x_continuous(n.breaks = 5) +  ###
  scale_y_continuous(n.breaks = 5) +  ###
  labs(x = expression(Cell~volume~(mu*m^3)), y = "Cell aspect ratio", ###
       color = expression(Delta * V ~ Delta * AR)) +  ### ∆V ∆AR
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cell_volume_vs_AR_onlylines", width = 5, height = 4)
save_ggplot("cell_volume_vs_AR_onlylines", width = 5, height = 4, mode = "paper")

# Cell volume vs AR (points only)
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = Mean_volume, y = Mean_AR, color = ColonyMorph, shape = ColonyMorph, size = ColonyMorph)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = colonymorphs_color) +
  scale_shape_manual(values = colonymorphs_shape) + 
  scale_size_manual(values = colonymorphs_size) + 
  scale_x_continuous(n.breaks = 5) +  ###
  scale_y_continuous(n.breaks = 5) +  ###
  labs(x = expression(Cell~volume~(mu*m^3)), y = "Cell aspect ratio", ###
       color = "Colony", shape = "Colony", size = "Colony") +  ###
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cell_volume_vs_AR_onlypoints", width = 5, height = 4)
save_ggplot("cell_volume_vs_AR_onlypoints", width = 5, height = 4, mode = "paper")

# Cell volume vs AR (points only, add decision boundary predicted by SVM with linear kernel by sklearn.svm in Python)
# Run SVM: GPT4 runs the SVM in Python, draws the plot, generates the slope and intercept values, then I re-draw the plot in R
# GPT4 prompts and codes are copied and stored in the file "svm_GPT_20231229.doc"
# Note: this is based on the current data_summary.csv file, GPT4 were run twice with slightly different prompts and generated the same slope/intercept values
# Create end points of linear SVM decision boundary
slope <- -0.01025
intercept <- 3.6963
svm_xs <- c(min(data_summary$Mean_volume), max(data_summary$Mean_volume))
svm_ys <- slope * svm_xs + intercept
svm_xys <- data.frame(Mean_volume = svm_xs, Mean_AR = svm_ys)
# Plot
last_plot() + 
  geom_line(data = svm_xys, mapping = aes(x = Mean_volume, y = Mean_AR), inherit.aes = FALSE, 
            color = "gray", linetype = "dashed")
save_ggplot("cell_volume_vs_AR_onlypoints_svm", width = 5, height = 4)
save_ggplot("cell_volume_vs_AR_onlypoints_svm", width = 5, height = 4, mode = "paper")

##### Plot distribution

# Cell volume
ggplot(data %>%
         dplyr::mutate(ColonyRep = paste0(ColonyMorph, gsub("rep", "", Replicate, fixed = TRUE), sep = ""),
                       ColonyRep = factor(ColonyRep, levels = c("D1","S1","D2","S2","D3","S3"))),
       aes(x = ColonyRep, y = Volume, fill = ColonyMorph, color = ColonyMorph)) + ###
  facet_wrap(~StrainBackground, nrow = 3) +
  geom_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) +
  geom_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, linewidth = 0.25) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = colonymorphs_color) +
  scale_color_manual(values = colonymorphs_color) +
  scale_y_continuous(limits = c(0,NA), n.breaks = 5, expand = expansion(mult = c(0,0.05))) +
  labs(x = "Replicate", y = expression(Cell~volume~(mu*m^3)), fill = "Colony", color = "Colony") +  ###
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25),
        legend.position = "bottom")
save_ggplot("cell_volume_distribution", width = 7, height = 7)
save_ggplot("cell_volume_distribution", width = 7, height = 7, mode = "paper")

# Cell AR
ggplot(data %>%
         dplyr::mutate(ColonyRep = paste0(ColonyMorph, gsub("rep", "", Replicate, fixed = TRUE), sep = ""),
                       ColonyRep = factor(ColonyRep, levels = c("D1","S1","D2","S2","D3","S3"))),
       aes(x = ColonyRep, y = AR, fill = ColonyMorph, color = ColonyMorph)) + ###
  facet_wrap(~StrainBackground, nrow = 3) +
  geom_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, linewidth = 0.25) +
  geom_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, linewidth = 0.25) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = colonymorphs_color) +
  scale_color_manual(values = colonymorphs_color) +
  scale_y_continuous(limits = c(1,NA), n.breaks = 5, expand = expansion(mult = c(0,0.05))) +  ###
  labs(x = "Replicate", y = "Cell aspect ratio", fill = "Colony", color = "Colony") +  ###
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25),
        legend.position = "bottom")
save_ggplot("cell_AR_distribution", width = 7, height = 7)
save_ggplot("cell_AR_distribution", width = 7, height = 7, mode = "paper")

##### Save data
write.csv(data, file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
write.csv(data_summary_mean, file = "data_summary_mean.csv", row.names = FALSE)
write.csv(data_summary_file, file = "data_summary_file.csv", row.names = FALSE)
write.csv(data_stat[["Volume"]] %>% apply(2, as.character), 
          file = "data_volume_stat.csv", row.names = FALSE)
write.csv(data_stat[["AR"]] %>% apply(2, as.character), 
          file = "data_AR_stat.csv", row.names = FALSE)
write.csv(data_summary_stat[["Volume"]] %>% apply(2, as.character) %>% t() %>% data.frame(), 
          file = "data_summary_volume_stat.csv", row.names = FALSE)
write.csv(data_summary_stat[["AR"]] %>% apply(2, as.character) %>% t() %>% data.frame(), 
          file = "data_summary_AR_stat.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")
saveRDS(data_stat, file = "data_stat.rds")
saveRDS(data_summary_stat, file = "data_summary_stat.rds")
