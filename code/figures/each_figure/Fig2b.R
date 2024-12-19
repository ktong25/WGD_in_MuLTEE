library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20230606_PMPA_evo_combine")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
in_dir <- "data_final_csvs"
out_fig_path <- "."
control_strains <- c("2N", "4N")  # first 2N then 4N
conditions <- c("PM", "PA")
lines <- 1:5 %>% as.character()
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- c("t0", "t200", "t400", "t600", "t1000")
evotimes2 <- gsub("t", "", evotimes, fixed = TRUE) %>% as.numeric()

##### Prepare data #####

# Load PM/PA t0 and PA4 t1000 data
# Input file name: e.g. "<date>.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = NULL, stringsAsFactors = FALSE) %>%
    dplyr::filter(!Strain %in% control_strains) %>%
    dplyr::select(Strain, RawIntDen_bgs, DNA_content) %>%
    dplyr::mutate(Batch = substr(basename(file), 1, stri_length(basename(file))-stri_length(".csv"))) %>%
    dplyr::select(Batch, Strain, RawIntDen_bgs, DNA_content)
  return(df)
}
data <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
              .fun = load_one_csv) %>%
  dplyr::filter(Strain %in% c("PA_4_t1000", "PA_t0", "PM_t0")) %>%
  dplyr::mutate(Strain = gsub("_t0", "_NA_t0", Strain, fixed = TRUE))
# Check each strain comes from only one batch
length(unique(data$Strain))
length(unique(paste(data$Batch, data$Strain)))
# Format data
data <- data %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line", "EvoTime"), sep = "_", 
                  remove = FALSE, convert = FALSE)

# Load other data
data2 <- read.csv("data_ploidy_20230211_combine.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::filter(!Strain %in% c("PA_4_t1000", "PA_0_t0", "PM_0_t0"))
# Check each strain comes from only one batch
length(unique(data2$Strain))
length(unique(paste(data2$Batch, data2$Strain)))

# Combine data
data <- rbind(data, data2)
remove(data2)
# Check each strain comes from only one batch
length(unique(data$Strain))
length(unique(paste(data$Batch, data$Strain)))

# Format data
data <- data %>%
  dplyr::arrange(factor(Condition, levels = conditions), Line, 
                 factor(EvoTime, levels = evotimes))
data$Condition <- factor(data$Condition, levels = conditions)
data$Line <- factor(data$Line, levels = lines)
data$EvoTime <- factor(data$EvoTime, levels = evotimes)
summary(data)

##### Plot ploidy #####

# Get G1 peak for each strain (G1 peak is usually the highest peak)
data_summary <- data %>%
  dplyr::group_by(Batch, Strain, Condition, Line, EvoTime) %>%
  dplyr::summarise(Count = n(),  
                   G1Peak = density(DNA_content)$x[which.max(density(DNA_content)$y)])
# Adjust for some strains by finding all peaks and choosing the right peak
# PA_2_t600
x <- data[data$Strain == "PA_2_t600", "DNA_content"]
d <- density(x)
d$x[c(F, diff(diff(d$y)>=0)<0)]  # 3.950631
data_summary[data_summary$Strain == "PA_2_t600", "G1Peak"] <- 3.950631
# PA_2_t1000
x <- data[data$Strain == "PA_2_t1000", "DNA_content"]
d <- density(x)
d$x[c(F, diff(diff(d$y)>=0)<0)]  # 4.7448181
data_summary[data_summary$Strain == "PA_2_t1000", "G1Peak"] <- 4.7448181

# Plot ploidy evolution
data_summary_t0 <- dplyr::filter(data_summary, EvoTime == "t0")
data_summary_plot <- dplyr::filter(data_summary, EvoTime != "t0")
for (line in c(lines, NA)) {  # final NA is for plotting gray color at the top of overlapping points at t0
  data_summary_plot <- rbind(data_summary_plot, 
                             dplyr::mutate(data_summary_t0, Line = line))
}
ggplot(data_summary_plot %>% 
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()), 
       aes(x = EvoTime2, y = G1Peak, color = Line)) +
  facet_wrap(~Condition, nrow = 1) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0, 5), breaks = 0:5, labels = c(0, paste0(1:5, "N")), 
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days of evolution", y = "DNA content", color = "Line") +
  ggplot_custom_theme4 +
  theme(panel.grid.major.y = element_line(size = 0.25, color = "gray"))
save_ggplot("ploidy_evo", width = 7, height = 4)
save_ggplot("ploidy_evo", width = 7, height = 4, mode = "paper")

# Save data
write.csv(data, file = "data_ploidy.csv", row.names = FALSE)
saveRDS(data, file = "data_ploidy.rds")
write.csv(data_summary, file = "data_ploidy_summary.csv", row.names = FALSE)
saveRDS(data_summary, file = "data_ploidy_summary.rds")
