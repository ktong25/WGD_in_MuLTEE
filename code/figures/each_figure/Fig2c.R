library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)
library(ggridges)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20230606_PMPA_early_combine")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
in_dir <- "data_final_csvs"
out_fig_path <- "."
control_strains <- c("2N", "4N")
conditions <- c("PM", "PA")
lines <- 1:5 %>% as.character()
evotimes <- c("t0", "t50", "t100")
ploidys_color <- get_pal_colors("Dark2")[c(1,2)]  # green, orange # "#1B9E77" "#D95F02"

##### Prepare data #####

# Load data
# Input file name: e.g. "<date>.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = NULL, stringsAsFactors = FALSE) %>%
    dplyr::filter(!Strain %in% c(control_strains)) %>%
    dplyr::select(Strain, RawIntDen_bgs, DNA_content) %>%
    dplyr::mutate(Batch = substr(basename(file), 1, stri_length(basename(file))-stri_length(".csv"))) %>%
    dplyr::select(Batch, Strain, RawIntDen_bgs, DNA_content)
  return(df)
}
data <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::filter(Strain != "PA_4_t1000") %>%
  dplyr::mutate(Strain = gsub("_t0", "_NA_t0", Strain, fixed = TRUE))
# Check each strain comes from only one batch
length(unique(data$Strain))
length(unique(paste(data$Batch, data$Strain)))
# Format data
data <- data %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line", "EvoTime"), sep = "_", 
                  remove = FALSE, convert = FALSE)

# Format data
data <- data %>%
  dplyr::arrange(factor(Condition, levels = conditions), Line, 
                 factor(EvoTime, levels = evotimes))
data$Condition <- factor(data$Condition, levels = conditions)
data$Line <- factor(data$Line, levels = lines)
data$EvoTime <- factor(data$EvoTime, levels = evotimes)
summary(data)

# Assign ploidy
data <- dplyr::mutate(data, Ploidy = case_when(EvoTime == "t0" ~ "2N",
                                               TRUE ~ "4N or mostly 4N"))

# Plot DNA content distribution
ggplot(data, aes(x = DNA_content, y = Line, fill = Ploidy)) +
  facet_grid(EvoTime ~ Condition, scales = "free_y", space = "free_y") + 
  geom_density_ridges(scale = 2, size = 0.4, alpha = 0.3) +
  scale_x_continuous(limits = c(0, 11), 
                     breaks = seq(0, 11, 2),
                     #minor_breaks = seq(0, 11, 1), 
                     labels = c("0", paste0(seq(2, 11, 2), "N")), 
                     expand = expansion(mult = c(0,0.05))) +
  scale_y_discrete(limits = rev, breaks = lines) +  # not show NA in t0
  scale_fill_manual(values = ploidys_color) + 
  labs(x = "DNA content", y = "Line") +
  ggplot_custom_theme4 +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_line(size = 0.3, colour = "gray80"),
        #panel.grid.minor.x = element_line(size = 0.2, colour = "grey90"), 
        axis.text.y = element_text(vjust = 0, margin = margin(r = 5)), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        axis.title.y = element_text(margin = margin(r = 5), hjust = 0.37))
save_ggplot("PMPA_early", width = 6, height = 5)
save_ggplot("PMPA_early", width = 6, height = 5, mode = "paper")

# Save data
write.csv(data, file = "data.csv", row.names = FALSE)
