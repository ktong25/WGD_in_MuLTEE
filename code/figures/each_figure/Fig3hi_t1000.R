library(dplyr)
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(scales)
library(glue)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20230920_agar_rev_combine")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
lines <- as.character(1:5)
lines_label <- paste0("L", lines)
names(lines_label) <- lines
evotimes <- "t1000"  # used in several places below
replicates <- c("anc", "A", "B", "C", "D") ###

# Load data
data <- read.csv(glue("data_{evotimes}.csv"), row.names = NULL, stringsAsFactors = FALSE)
data_anc <- read.csv("~/Documents/projects/R_data_analysis/ploidy_measurement/20230606_PMPA_evo_combine/data_ploidy.csv", row.names = NULL, stringsAsFactors = FALSE) %>% ###
  dplyr::filter(EvoTime == "t1000") %>%
  dplyr::mutate(Replicate = "anc") %>%
  dplyr::relocate(Replicate, .after = EvoTime)
data <- rbind(data, data_anc) ###
data <- data %>% dplyr::mutate(
  Condition = factor(Condition, levels = conditions), 
  Line = factor(Line, levels = lines), 
  Replicate = factor(Replicate, levels = replicates)
)
summary(data)

# Assign ploidy reduction level
data_summary <- data %>%
  dplyr::group_by(Batch, Strain, Condition, Line, EvoTime, Replicate) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::arrange(factor(Condition, levels = conditions), Line, EvoTime, Replicate) %>%
  dplyr::mutate(
    Condition = factor(Condition, levels = conditions), 
    Line = factor(Line, levels = lines), 
    Replicate = factor(Replicate, levels = replicates)
  )
data_summary$PloidyReduction <- 
  c(NA,2,1,1,1, # PM
    NA,1,3,3,3,
    NA,2,2,2,2,
    NA,3,3,3,3,
    NA,1,1,1,1,
    NA,3,3,3,3, # PA
    NA,2,2,2,3,
    NA,2,1,3,2,
    NA,2,3,3,2,
    NA,3,2,3,3
  ) %>% 
  factor(., levels = c(1,2,3), labels = c("Complete", "Partial", "No"))
data <- dplyr::left_join(x = data, 
                         y = data_summary %>% dplyr::select(!Count), 
                         by = c("Batch", "Strain", "Condition", "Line", "EvoTime", "Replicate"))
prs_color <- c("#1B9E77", "lightgreen", "#D95F02")  # ploidy reduction colors
names(prs_color) <- c("Complete", "Partial", "No")

# Plot with ploidy reduction level assigned
ggplot(data %>% dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = paste(conditions, evotimes))), ###
       aes(x = DNA_content, y = Replicate, fill = PloidyReduction)) +
  facet_grid(Line ~ Condition) +  #labeller = labeller(Line = lines_label)
  geom_density_ridges(scale = 2, size = 0.4, alpha = 0.3) +
  scale_x_continuous(limits = c(0, 11.5), 
                     breaks = seq(0, 11, 2),
                     #minor_breaks = seq(0, 11, 1), 
                     labels = c("0", paste0(seq(2, 11, 2), "N")), 
                     expand = expansion(mult = c(0,0.05))) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = prs_color) +
  labs(#title = evotimes, ###
       x = "DNA content", y = "Replicate population", fill = "Ploidy reduction") +
  ggplot_custom_theme4 +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_line(size = 0.3, colour = "gray80"),
        axis.text.y = element_text(vjust = 0, margin = margin(r = 5)), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        axis.title.y = element_text(margin = margin(r = 5)),
        #plot.title = element_text(size = rel(1.75), hjust = 0.5, vjust = -1) ###
        )
save_ggplot(glue("reversion_{evotimes}"), width = 8, height = 8)
save_ggplot(glue("reversion_{evotimes}"), width = 8, height = 8, mode = "paper")

# Save data
write.csv(data, file = glue("data_plus_anc_{evotimes}.csv"), row.names = FALSE) ###
write.csv(data_summary, file = glue("data_plus_anc_summary_{evotimes}.csv"), row.names = FALSE) ###
