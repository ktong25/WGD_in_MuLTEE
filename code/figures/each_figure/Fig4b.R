library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggforce)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)
library(purrr)
library(rstatix)
library(ggpubr)
library(weights)
library(introdataviz)

setwd("~/Documents/projects/R_data_analysis/cluster_size/20231227_donutspread")
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
# Input file name: e.g., "PA_1_t600_D_rep1_Results.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1) %>% dplyr::select(BF.Area)
  colnames(df) <- gsub("BF.", "", colnames(df), fixed = TRUE)
  measure_vars <- colnames(df)
  # Extract metadata
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  df$Strain <- metadata
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  df$Condition <- metadata[1]
  df$Line <- metadata[2]
  df$EvoTime <- metadata[3]
  df$ColonyMorph <- metadata[4]
  df$Replicate <- metadata[5]
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv)
summary(data)
# Check if data contains objects < 40um^2
# dplyr::filter(data, Area < 40) %>% View()
# None
# nrow(data)
# data <- dplyr::filter(data, Area > 40)
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
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph)
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

# Convert area to radius and volume
data$Radius <- sqrt(data$Area / pi)
data$Volume <- 4/3 * pi * (data$Radius)^3

# Data summary
# Calculate cluster radius: weighted mean radius, median radius, mean of top 5% radius (proxy for radius at reproduction)
id_vars <- c("Strain", "CondLine", "StrainBackground", "StrainBackgroundRep", 
             "Condition", "Line", "EvoTime", "Replicate", "ColonyMorph")
data_summary <- data %>%
    dplyr::group_by_at(id_vars) %>%
    dplyr::summarize(Count = n(), 
                     Weighted_mean_radius = (weighted.mean(Volume, Volume) / (4/3 * pi)) ^ (1/3), 
                     Median_radius = median(Radius),
                     Mean_top_radius = Radius[Radius > quantile(Radius, prob = 0.95)] %>% mean()) %>%
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph)
summary(data_summary)

# Calculate mean of D and S
data_summary_mean <- data_summary %>%
  dplyr::group_by(ColonyMorph) %>%
  dplyr::summarize(Weighted_mean_radius.Mean = mean(Weighted_mean_radius))

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
# Weighted Welch's t-test with Benjamini-Hochberg ("BH" or its alias "fdr") correction
# Weighted t-test is performed on Volume instead of Radius, which will produce slightly different results but should be fine
# Doing such a test on Radius can be hard because there is a cube root from Volume, and I do not know how to incorporate that into a statistical test
data_stat <- map_dfr(strainbackgroundreps, function(s) {
  subdata <- data[data$StrainBackgroundRep == s, ]
  xs <- subdata %>% dplyr::filter(ColonyMorph == colonymorphs[1]) %>% dplyr::pull(Volume)
  ys <- subdata %>% dplyr::filter(ColonyMorph == colonymorphs[2]) %>% dplyr::pull(Volume)
  result <- wtd.t.test(x = xs, y = ys, weight = xs, weighty = ys,
                       samedata = FALSE, mean1 = FALSE)
  df <- data.frame(
    StrainBackgroundRep = s,
    estimate = as.numeric(result$additional["Difference"]),
    estimate1 = as.numeric(result$additional["Mean.x"]),
    estimate2 = as.numeric(result$additional["Mean.y"]),
    .y. = "Volume", 
    group1 = colonymorphs[1],
    group2 = colonymorphs[2],
    n1 = length(xs), 
    n2 = length(ys), 
    statistic = as.numeric(result$coefficients["t.value"]), 
    p = as.numeric(result$coefficients["p.value"]), 
    df = as.numeric(result$coefficients["df"]), 
    method = "Weighted t-test", 
    alternative = "two.sided"
  )
  return(df)
}) %>% data.frame() %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  dplyr::mutate(SignPadj = ifelse(p.adj.signif == "ns", "ns", 
                                  paste(ifelse(estimate2 > estimate1, "+", "-"), p.adj.signif, sep = " ")), 
                SignPadj = factor(SignPadj, levels = sign_padjs))
data_summary <- dplyr::left_join(data_summary, 
                                 data_stat %>% 
                                   dplyr::select(StrainBackgroundRep, SignPadj) %>%
                                   dplyr::rename(Weighted_mean_volume.SignPadj = SignPadj), 
                                 by = "StrainBackgroundRep")

# Plot
ggplot(data_summary %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = ColonyMorph, y = Weighted_mean_radius, shape = Replicate2, group = Replicate2)) + ###
  facet_wrap(~StrainBackground, nrow = 3) +
  geom_line(mapping = aes(color = Weighted_mean_volume.SignPadj),  ### 
            size = 0.75, alpha = 0.8) +
  geom_point(size = 2, color = "gray", alpha = 1) +
  scale_color_manual(values = sign_padjs_color, labels = sign_padjs, drop = FALSE) +
  scale_shape_manual(values = reps_shape) + 
  scale_y_continuous(trans = "log2", n.breaks = 6) +  ###
  labs(x = NULL, y = expression(Cluster~radius~(mu*m)), ###
       color = expression("Sign P"[adj]), shape = "Replicate") +  
  ggplot_custom_theme4 +
  theme(legend.position = "right")
save_ggplot("cluster_size_showinfo", width = 8, height = 7)
save_ggplot("cluster_size_showinfo", width = 8, height = 7, mode = "paper")

##### Plot (not show the sample information for each data point)

# Perform statistical tests (D vs S for all strain backgrounds combined)
# Paired t-test
data_summary_stat <- data_summary %>%
  dplyr::ungroup() %>%
  rstatix::t_test(Weighted_mean_radius ~ ColonyMorph, paired = TRUE, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  add_xy_position(x = "ColonyMorph", fun = "max") %>%
  dplyr::mutate(y.position = log2(y.position) * 1.05)

# Plot
ggplot(data_summary, aes(x = ColonyMorph, y = Weighted_mean_radius, color = ColonyMorph)) +  ###
  geom_line(mapping = aes(group = StrainBackgroundRep), 
            size = 0.5, alpha = 0.8, color = "gray") +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_boxplot(width = 0.5, size = 0.5, outlier.shape = NA, fill = NA) +
  stat_pvalue_manual(data_summary_stat, label = "p = {p}",  ###
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_color_manual(values = colonymorphs_color) +
  scale_y_continuous(trans = "log2", n.breaks = 6, expand = expansion(mult = c(0.05, 0.1))) +  ###
  labs(x = NULL, y = expression(Cluster~radius~(mu*m))) +  ###
  ggplot_custom_theme4 +
  theme(legend.position = "none")
save_ggplot("cluster_size", width = 4, height = 4)
save_ggplot("cluster_size", width = 4, height = 4, mode = "paper")

##### Plot distribution

# Plot distribution of cluster radius, donut/spread as left/right violin, weighted by cluster volume
# Plot weighted mean radius
ggplot(data %>% 
         dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)), 
       aes(x = Replicate2, y = Radius, fill = ColonyMorph, color = ColonyMorph, weight = Volume)) +
  facet_wrap(~StrainBackground, nrow = 3) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, fill = "white", alpha = 1, size = 0.25) +
  introdataviz::geom_split_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.3, size = 0.25) +
  geom_point(data = data_summary %>% 
               dplyr::mutate(Replicate2 = factor(Replicate, levels = reps, labels = reps_label)),
             mapping = aes(x = Replicate2, y = Weighted_mean_radius, fill = ColonyMorph), inherit.aes = FALSE,
             shape = "circle filled", size = 2, color = "black", alpha = 0.5, 
             position = position_dodge(width = 0.4)) +
  scale_fill_manual(values = colonymorphs_color) +
  scale_y_continuous(trans = "log2", n.breaks = 6) +
  labs(x = "Replicate", y = expression(Cluster~radius~(mu*m)), fill = "Colony", color = "Colony") +
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(color = "gray70", size = 0.25), 
        legend.position = "bottom")

save_ggplot("cluster_size_distribution_weighted", width = 7, height = 7)
save_ggplot("cluster_size_evo_distribution_weighted", width = 7, height = 7, mode = "paper")

##### Save data
write.csv(data, file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
write.csv(data_summary_mean, file = "data_summary_mean.csv", row.names = FALSE)
write.csv(data_stat, file = "data_stat.csv", row.names = FALSE)
write.csv(data_summary_stat %>% apply(2, as.character) %>% t() %>% data.frame(), 
          file = "data_summary_stat.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")
saveRDS(data_stat, file = "data_stat.rds")
saveRDS(data_summary_stat, file = "data_summary_stat.rds")

