library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/cluster_size/20231201_agar")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
lines <- as.character(1:5)
evotimes <- c("t0", "t1000")
replicates <- c("anc", "A", "B", "C", "D")

# Load data
# Input file names: 
# "t0_GOB8_28JAN22_a_BF_Results.csv"
# "t70_GOB8A_28JAN22_a_BF_Results.csv"
# "GOB2149_4x_BF_2_16NOV22_a_BF_Results.csv" GOB2149-2153 t0
# "GOB2151D_t70_agarexp_13OCT23_a_BF_Results.csv" x3
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1) %>% dplyr::select(BF.Area)
  colnames(df) <- "Area"
  # Extract metadata
  filename <- basename(file)
  metadata <- substr(filename, 1, stri_length(filename)-stri_length("_Results.csv"))
  metadata <- strsplit(metadata, split = "_", fixed = TRUE)[[1]]
  # Specify Condition / Line / EvoTime / Replicate (anc, A-D)
  if (metadata[1] %in% c("t0", "t70")) {
    df$Time <- metadata[1]  # t0, t70
    df$StrainRep <- metadata[2]  # GOB8, GOB8A etc
  }
  else if (metadata[2] == "4x") {
    df$Time <- "t0"
    df$StrainRep <- metadata[1]
  }
  else if (metadata[2] == "t70") {
    df$Time <- "t70"
    df$StrainRep <- metadata[1]
  }
  df <- df %>%
    dplyr::mutate(Strain = ifelse(sum(endsWith(StrainRep, c("A", "B", "C", "D"))) > 0, 
                                  as.integer(substr(StrainRep, 4, stri_length(StrainRep)-1)), 
                                  as.integer(substr(StrainRep, 4, stri_length(StrainRep)))), 
                  Replicate = ifelse(sum(endsWith(StrainRep, c("A", "B", "C", "D"))) > 0, 
                                     substr(StrainRep, stri_length(StrainRep), stri_length(StrainRep)), 
                                     "anc")
                  ) %>%
    dplyr::mutate(Condition = ifelse(Strain %in% c(8, 2144:2148), "PM", "PA"), 
                  Line = case_when(
                    Strain %in% 2144:2148 ~ as.character(Strain - 2143), 
                    Strain %in% 2149:2153 ~ as.character(Strain - 2148), 
                    Strain %in% c(8, 21) ~ "NA", 
                    TRUE ~ as.character(Strain), 
                  ),
                  EvoTime = ifelse(Strain %in% c(8, 21), "t0", "t1000")
                  ) %>%
    dplyr::mutate(Line = ifelse(Line == "NA", NA, Line))
  # Reorder output column names
  df <- df %>% dplyr::select(Condition, Line, EvoTime, Replicate, Area)
  return(df)
}
data <- ldply(.data = list.files(path = "03_measurement", pattern = "*_Results.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Condition, levels = conditions), 
                 factor(EvoTime, levels = evotimes),
                 Line, 
                 factor(Replicate, levels = replicates))
# Check
unique(data$Condition)
unique(data$Line)
unique(data$EvoTime)
unique(data$Replicate)
unique(paste(data$Condition, data$Line, data$EvoTime, data$Replicate))
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
data$Condition <- factor(data$Condition, levels = conditions)
data$Line <- factor(data$Line, levels = lines)
data$EvoTime <- factor(data$EvoTime, levels = evotimes)
data$Replicate <- factor(data$Replicate, levels = replicates)
summary(data)

# Convert area to radius and volume
data$Radius <- sqrt(data$Area / pi)
data$Volume <- 4/3 * pi * (data$Radius)^3

# Calculate weighted mean radius
id_vars <- colnames(data)[1:4]
data_summary <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarize(Count = n(), 
                   Weighted_mean_radius = (weighted.mean(Volume, Volume) / (4/3 * pi)) ^ (1/3))  %>%
  dplyr::arrange(factor(Condition, levels = conditions), 
                 factor(EvoTime, levels = evotimes),
                 Line, 
                 factor(Replicate, levels = replicates))

# Assign ploidy reduction level
data_summary$PloidyReduction <- 
  c(NA,3,3,3,3, # PM t0 
    NA,2,1,1,1, 
    NA,1,3,3,3,
    NA,2,2,2,2,
    NA,3,3,3,3,
    NA,1,1,1,1,
    NA,3,3,3,3, # PA t0
    NA,3,3,3,3, 
    NA,2,2,2,3,
    NA,2,1,3,2,
    NA,2,3,3,2,
    NA,3,2,3,3
  ) %>%
  factor(., levels = 1:3, labels = c("Complete", "Partial", "No"))
data$PloidyReduction <- rep(data_summary$PloidyReduction, data_summary$Count)
prs_color <- c("#1B9E77", "lightgreen", "#D95F02")  # ploidy reduction colors
names(prs_color) <- c("Complete", "Partial", "No")

# Plot with ploidy reduction level assigned
PM_t0_radius <- data_summary %>% dplyr::filter(Condition == "PM" & EvoTime == "t0" & Replicate == "anc") %>% dplyr::pull(Weighted_mean_radius)
PA_t0_radius <- data_summary %>% dplyr::filter(Condition == "PA" & EvoTime == "t0" & Replicate == "anc") %>% dplyr::pull(Weighted_mean_radius)
# t0
ggplot(data %>% dplyr::filter(is.na(Line)) %>%
         dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = c("PM t0", "PA t0"))), 
       aes(x = Radius, y = Replicate, fill = PloidyReduction)) +
  facet_wrap(~Condition) + 
  geom_density_ridges(mapping = aes(height=after_stat(density), weight = Volume),  # weighted, see https://github.com/wilkelab/ggridges/issues/5
                      stat="density", 
                      scale = 1, size = 0.4, alpha = 0.3) +
  # Add weighted mean radius for each population
  geom_segment(data = data_summary %>%
                 dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = c("PM t0", "PA t0"))) %>% 
                 dplyr::filter(is.na(Line)) %>%
                 dplyr::group_by(Condition, Line, EvoTime) %>%
                 dplyr::mutate(y = length(replicates):1), 
               mapping = aes(x = Weighted_mean_radius, y = y, xend = Weighted_mean_radius, yend = y+0.7), 
               color = "gray20", linewidth = 1, alpha = 0.6, linetype = "solid") +
  # Add a reference line showing t0 weighted mean radius
  geom_vline(data = data_summary %>%
               dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = c("PM t0", "PA t0"))) %>% 
               dplyr::mutate(Weighted_mean_radius = ifelse(Condition == "PM t0", PM_t0_radius, PA_t0_radius)) %>%
               dplyr::filter(is.na(Line)), 
             mapping = aes(xintercept = Weighted_mean_radius), 
             size = 0.5, color = "gray", alpha = 1, linetype = "dashed") +
  scale_x_continuous(trans = "log2", n.breaks = 7, limits = c(6,600)) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = rep("#1B9E77", 3), breaks = c("No"), labels = c("2N or mostly 2N"), drop = TRUE) +
  labs(x = expression(Cluster~radius~(mu*m)), 
       y = "Replicate\npopulation", fill = "Evolved ploidy") +
  ggplot_custom_theme4 +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "gray80"),
        axis.text.y = element_text(vjust = 0, margin = margin(r = 5)), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        axis.title.y = element_text(margin = margin(r = 5)))
save_ggplot(glue("reversion_size_t0"), width = 8, height = 3)
save_ggplot(glue("reversion_size_t0"), width = 8, height = 3, mode = "paper")
# t1000
ggplot(data %>% dplyr::filter(!is.na(Line)) %>%
         dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = c("PM t1000", "PA t1000"))), 
       aes(x = Radius, y = Replicate, fill = PloidyReduction)) +
  facet_grid(Line ~ Condition) + 
  geom_density_ridges(mapping = aes(height=after_stat(density), weight = Volume),  # weighted, see https://github.com/wilkelab/ggridges/issues/5
                      stat="density", 
                      scale = 1, size = 0.4, alpha = 0.3) +
  # Add weighted mean radius for each population
  geom_segment(data = data_summary %>%
                 dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = c("PM t1000", "PA t1000"))) %>% 
                 dplyr::filter(!is.na(Line)) %>%
                 dplyr::group_by(Condition, Line, EvoTime) %>%
                 dplyr::mutate(y = length(replicates):1), 
               mapping = aes(x = Weighted_mean_radius, y = y, xend = Weighted_mean_radius, yend = y+0.7), 
               color = "gray20", linewidth = 1, alpha = 0.6, linetype = "solid") +
  # Add a reference line showing t0 weighted mean radius
  geom_vline(data = data_summary %>%
               dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = c("PM t1000", "PA t1000"))) %>% 
               dplyr::mutate(Weighted_mean_radius = ifelse(Condition == "PM t1000", PM_t0_radius, PA_t0_radius)) %>%
               dplyr::filter(!is.na(Line)), 
             mapping = aes(xintercept = Weighted_mean_radius), 
             size = 0.5, color = "gray", alpha = 1, linetype = "dashed") +
  scale_x_continuous(trans = "log2", n.breaks = 7, limits = c(6,600)) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = prs_color) +
  labs(x = expression(Cluster~radius~(mu*m)), 
       y = "Replicate population", fill = "Ploidy reduction") +
  ggplot_custom_theme4 +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "gray80"),
        axis.text.y = element_text(vjust = 0, margin = margin(r = 5)), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        axis.title.y = element_text(margin = margin(r = 5)))
save_ggplot(glue("reversion_size_t1000"), width = 8, height = 8)
save_ggplot(glue("reversion_size_t1000"), width = 8, height = 8, mode = "paper")

# Save data
write.csv(data, file = "data.csv", row.names = FALSE)
write.csv(data_summary, file = "data_summary.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")
saveRDS(data_summary, file = "data_summary.rds")

