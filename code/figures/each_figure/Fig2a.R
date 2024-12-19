library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(stringi)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/allele_frequency/20231223_MuLTEE_summary")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
in_file_suffix <- ".tab"
conditions <- c("PM", "PA")
#conditions_color <- get_pal_colors("Set1")[c(2,1,3)]  # like in Bozdag2021, blue, red, green # "#377EB8" "#E41A1C" "#4DAF4A"
#names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
#conditions_label <- c("Mixotrophic", "Anaerobic", "Aerobic")
#ploidys <- c("2N", "3N", "4N")
#ploidys_color <- c("#1B9E77", "#6A3D9A", "#D95F02") # green, purple, orange (Dark2-1, Paired-10, Dark2-2)  
#names(ploidys_color) <- ploidys
#show_col(ploidys_color, labels = TRUE)
lines <- as.character(1:5)
#lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
#names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- paste0("t", c(200, 400, 600, 1000))

# Load data
# Input file name: e.g., "PA1_t200_XXX.tab"
load_one_csv <- function(file) {
  # Read and process file
  df <- read.table(file, header = FALSE, row.names = NULL) %>%
    dplyr::select(V1, V2, V4, V5, V11) %>%
    dplyr::mutate(V1 = gsub("chr", "", V1, fixed = TRUE), 
                  V1 = ifelse(startsWith(V1, "ref|NC_00"),  # some files may have format like "ref|NC_001134|"
                              V1 %>%
                                substr(., stri_length("ref|NC_00")+1, stri_length(.)-1) %>% 
                                factor(., levels = as.character(1:16 + 1132), labels = as.character(as.roman(1:16))) %>%
                                as.character(), 
                              V1), 
                  V11 = as.numeric(gsub("%", "", V11, fixed = TRUE)) / 100)
  colnames(df) <- c("Chr", "Pos", "Nuc_ref", "Nuc_alt", "Alt_freq")
  measure_vars <- colnames(df)
  # Extract metadata
  metadata <- strsplit(basename(file), split = "_", fixed = TRUE)[[1]]
  df$Condition <- substr(metadata[1], 1, 2)
  df$Line <- substr(metadata[1], 3, 3)
  df$EvoTime <- metadata[2]
  # Reorder output column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = "ozan_data", 
                                 pattern = in_file_suffix, full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::arrange(factor(Condition, levels = conditions), Line, 
                 factor(EvoTime, levels = evotimes), 
                 factor(Chr, levels = as.character(as.roman(1:16))),
                 Pos) %>%
  dplyr::mutate(Condition = factor(Condition, levels = conditions), 
                Line = factor(Line, levels = lines), 
                EvoTime = factor(EvoTime, levels = evotimes))
summary(data)
id_vars <- c("Condition", "Line", "EvoTime")

# Plot allele frequency
ggplot(data, aes(x = Line, y = Alt_freq)) +
  facet_grid(Condition ~ EvoTime) +
  geom_hline(yintercept = seq(0.25,0.75,0.25), size = 0.25, color = "gray") +
  geom_violin(scale = "width", trim = TRUE, size = 0.5, color = "black", fill = "white", 
              position = position_dodge(width = 1)) + 
  geom_violin(scale = "width", trim = TRUE, size = 0.5, color = "black", fill = "skyblue3", alpha = 0.2, 
              position = position_dodge(width = 1)) + 
  geom_sina(scale = "width", seed = 42, size = 0.1,  color = "skyblue3", 
            position = position_dodge(width = 1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c(0L, 0.25, 0.5, 0.75, 1L), 
                     expand = expansion(mult = 0.025)) +
  labs(x = "Line", y = "Allele frequency") +
  ggplot_custom_theme3
save_ggplot("allele_frequency", width = 7, height = 4)
save_ggplot("allele_frequency", width = 7, height = 4, mode = "paper")

# Save data
write.csv(data, "data.csv", row.names = FALSE)
saveRDS(data, file = "data.rds")

##### Plot mutation number #####

# Parameters
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- c("t0", "t200", "t400", "t600", "t1000")
evotimes2 <- gsub("t", "", evotimes, fixed = TRUE) %>% as.numeric()

# Get mutation number
data_summary <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarise(Mut_num = n())

# Plot mutation number over time
data_summary_t0 <- data.frame(
  Condition = factor(c("PM", "PA"), levels = conditions), 
  Line = NA, 
  EvoTime = "t0", 
  Mut_num = 0
)
data_summary_plot <- data_summary
for (line in c(lines, NA)) {  # final NA is for plotting gray color at the top of overlapping points at t0
  data_summary_plot <- rbind(data_summary_plot, 
                             dplyr::mutate(data_summary_t0, Line = line))
}
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()), 
       aes(x = EvoTime2, y = Mut_num, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0, NA), n.breaks = 6, expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Days of evolution", y = "Mutation number", color = "Line") +
  ggplot_custom_theme4 #+
  #theme(panel.grid.major.y = element_line(size = 0.25, color = "gray"))
save_ggplot("nmut_evo", width = 7, height = 4)
save_ggplot("nmut_evo", width = 7, height = 4, mode = "paper")

# Save data
write.csv(data_summary, "data_summary.csv", row.names = FALSE)
saveRDS(data, file = "data_summary.rds")
