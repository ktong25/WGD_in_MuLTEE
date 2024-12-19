library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20230920_agar_rev_combine")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
in_dir <- "data_final_csvs"
out_fig_path <- "."
control_strains <- c("2N", "4N")  # first 2N then 4N
conditions <- c("PM", "PA")
#conditions_color <- get_pal_colors("Set1")[c(2,1)]  # 3 # like in Bozdag2021, blue, red, green # "#377EB8" "#E41A1C" "#4DAF4A"
#names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
#conditions_label <- c("Mixotrophic", "Anaerobic")  #, "Aerobic")
#lines_shape <- c("circle", "triangle", "square", "diamond", "square cross")

# Load data
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
  dplyr::filter(!( (Batch == "20230713_1" & Strain == "8B_agarexpt70") | 
                   (Batch == "20230714_1" & Strain == "2151D_agarexp_t70") | 
                   (Batch == "20230915" & Strain == "GOB8B_D63")
                ))
summary(data)
# Check each strain comes from only one batch
length(unique(data$Strain))
length(unique(paste(data$Batch, data$Strain)))

# Change strain name
unique(data$Strain)
# "8A_agarexpt70" ~ "8A_agarexp_t70"
# "8C_agarexp_t70"
# "8D_agarexp_t70"
# "21A_agarexp_t70"
# "21B_agarexp_t70"
# "21C_agarexp_t70"
# "21D_agarexp_t70"
# "2153C_agarexp_t70"
# "2153D_agarexp_t70"
# "GOB2151D_D70" ~ "2151D_agarexp_t70"
# "GOB8B_D70" ~ "8B_agarexp_t70"
data <- data %>%
  dplyr::mutate(Strain = case_when(
    Strain == "8A_agarexpt70" ~ "8A_agarexp_t70", 
    Strain == "GOB2151D_D70" ~ "2151D_agarexp_t70", 
    Strain == "GOB8B_D70" ~ "8B_agarexp_t70",
    TRUE ~ Strain
  )) 
unique(data$Strain)
data$Strain <- data$Strain %>%
  gsub("_agarexp_t70", "", ., fixed = TRUE) %>%
  gsub("8", "PM_NA_t0_", ., fixed = TRUE) %>% 
  gsub("2151", "PA_3_t1000_", ., fixed = TRUE) %>% 
  gsub("2153", "PA_5_t1000_", ., fixed = TRUE) %>%
  gsub("21", "PA_NA_t0_", ., fixed = TRUE)
unique(data$Strain)
data <- data %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line", "EvoTime", "Replicate"), sep = "_", 
                  remove = FALSE, convert = FALSE)
data <- data %>%
  dplyr::mutate(Strain = Strain %>%
                  gsub("_A", "", ., fixed = TRUE) %>%
                  gsub("_B", "", ., fixed = TRUE) %>%
                  gsub("_C", "", ., fixed = TRUE) %>%
                  gsub("_D", "", ., fixed = TRUE)) %>%
  dplyr::arrange(factor(Condition, levels = conditions), Line, EvoTime, Replicate)
unique(data$Strain)

# Combine data
data1 <- data
data2 <- read.csv("20221122_combine_data_rev.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::filter(!( (Strain == "PA_3_t1000" & Replicate == "D") | 
                   (Strain == "PA_5_t1000" & Replicate %in% c("C", "D"))
  )) %>%
  dplyr::select(!PloidyReduction)
length(unique(data2$Strain))  #10
length(unique(paste(data2$Strain, data2$Replicate)))  #37
# t0
data_t0 <- data1 %>% dplyr::filter(EvoTime == "t0")
length(unique(data_t0$Strain))  #2
length(unique(paste(data_t0$Strain, data_t0$Replicate)))  #8
length(unique(paste(data_t0$Strain, data_t0$Replicate, data_t0$Batch)))  #8
data_t0 <- data_t0 %>% dplyr::arrange(factor(Condition, levels = conditions), Line, EvoTime, Replicate)
# t1000
data_t1000 <- rbind(data1 %>% dplyr::filter(EvoTime == "t1000"), data2)
length(unique(data_t1000$Strain))  #10
length(unique(paste(data_t1000$Strain, data_t1000$Replicate)))  #40
length(unique(paste(data_t1000$Strain, data_t1000$Replicate, data_t1000$Batch)))  #40
data_t1000 <- data_t1000 %>% dplyr::arrange(factor(Condition, levels = conditions), Line, EvoTime, Replicate)

# Save data
write.csv(data_t0, file = "data_t0.csv", row.names = FALSE)
write.csv(data_t1000, file = "data_t1000.csv", row.names = FALSE)
write.csv(rbind(data_t0, data_t1000), file = "data_t0_t1000.csv", row.names = FALSE)

