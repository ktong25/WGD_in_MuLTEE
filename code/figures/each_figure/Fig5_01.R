library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/chromosome_cnv/20231230_donutspread")  ###
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
in_cov_path <- "03_coverage"
in_cov_suffix <- ".coverage"
out_fig_path <- "."
bin_length <- 1000
conditions <- c("PA")
lines <- as.character(1:5)
evotimes <- c("t600", "t1000")
reps <- paste0("rep", 1:3)
reps_label <- gsub("rep", "", reps, fixed = TRUE)
colonymorphs <- c("D", "S") 
chrs_color <- hue_pal()(16)
chrs_color <- chrs_color[c(seq(1,16,3), seq(2,16,3), seq(3,16,3))]
#show_col(chrs_color, labels = TRUE)
chr_copy_nums_color <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(9)[2:9]  #PuOr
names(chr_copy_nums_color) <- 1:8
#show_col(chr_copy_nums_color, labels = TRUE)

# Load sample metadata
metadata <- read.csv(file.path(in_cov_path, "metadata.csv"), 
                     row.names = NULL, colClasses = "character")
metadata$Ploidy <- as.numeric(metadata$Ploidy)
metadata_vars <- colnames(metadata)  # first column is "Sample_ID", last column is "Ploidy"
#ref_sample <- "PM_t0"

# Load and process data
load_one_coverage <- function(file) {
  # Import data
  df <- read.table(file, sep = "\t", row.names = NULL, 
                   header = FALSE, col.names = c("Chr", "Base", "Coverage"), 
                   stringsAsFactors = FALSE)
  # Extract sample ID
  filename <- basename(file)
  sample_id <- substr(filename, 1, stri_length(filename) - stri_length(in_cov_suffix))
  # Pre-process data
  df <- df %>%
    # Remove MT rows
    dplyr::filter(Chr != "ref|NC_001224|") %>%
    # Convert chromosome ref # to chromosome #
    dplyr::mutate(Chr = Chr %>%
                    substr(., stri_length("ref|NC_00")+1, stri_length(.)-1) %>% 
                    as.integer() - 1132, 
                  Chr = factor(Chr, levels = 1:16, labels = as.character(as.roman(1:16))) )
  # Get per-bin coverage (average per-base coverage in each non-overlapping sliding window) for each chromosome
  df_bin <- df %>%
    dplyr::group_by(Chr) %>%
    dplyr::mutate(Bin = Base %/% bin_length) %>%  # 0-based
    dplyr::group_by(Chr, Bin) %>%
    dplyr::summarise(Bin_coverage = sum(Coverage) / bin_length) %>%  # not mean(Coverage) because some bases may be absent
    as.data.frame()
  # Normalize by whole-genome coverage (median of per-bin coverages)
  df_bin <- df_bin %>%
    dplyr::mutate(Bin_coverage_norm = Bin_coverage / median(Bin_coverage))
  # Add metadata
  cur_vars <- colnames(df_bin)
  for (v in metadata_vars) {
    df_bin <- df_bin %>%
      dplyr::mutate("{v}" := metadata[metadata$Sample_ID == sample_id, v])
  }
  id_vars <- metadata_vars[!metadata_vars %in% c("Sample_ID", "Ploidy")]
  df_bin <- df_bin %>% 
    tidyr::unite(., Sample, all_of(id_vars), sep = "_", remove = FALSE, na.rm = TRUE) %>%
    dplyr::select(all_of(c("Sample_ID", "Sample", id_vars, "Ploidy", cur_vars)))
  return(df_bin)
}
data <- ldply(.data = list.files(path = in_cov_path, 
                                 pattern = paste0("*", in_cov_suffix), full.names = TRUE),
              .fun = load_one_coverage)
summary(data)
# The above step is slow, thus save data now in case you need to start from scratch later
write.csv(data, file = "data_load.csv", row.names = FALSE)
# Here is how you re-load the data (note the data type of some columns may be changed)
# data <- read.csv("data_load.csv", row.names = NULL, stringsAsFactors = FALSE)
# data$Chr <- factor(data$Chr, levels = as.character(as.roman(1:16)))
# summary(data)

# Remove and rename some samples
# Sample column and Strain column are equivalent here, Sample was used for chromosome_cnv scripts, Strain was used for donut-spread cluster/cell/ploidy scripts
data <- data %>%
  dplyr::mutate(Strain = Sample) %>%
  # Remove some samples
  dplyr::filter(!(
    Strain %in% c(paste0("PA_1_t600_", c("D_rep1", "S_rep1", "D_rep4", "S_rep4")), 
                  paste0("PA_4_t600_", c("D", "S"), "_rep3"))
  )) %>%
  # Rename samples
  # PA1 t600 rep2/3 -> rep1/2
  # PA4 t600 rep4 -> rep3
  dplyr::mutate(Strain = case_when(Strain == "PA_1_t600_D_rep2" ~ "PA_1_t600_D_rep1",
                                   Strain == "PA_1_t600_S_rep2" ~ "PA_1_t600_S_rep1",
                                   Strain == "PA_1_t600_D_rep3" ~ "PA_1_t600_D_rep2",
                                   Strain == "PA_1_t600_S_rep3" ~ "PA_1_t600_S_rep2",
                                   Strain == "PA_4_t600_D_rep4" ~ "PA_4_t600_D_rep3",
                                   Strain == "PA_4_t600_S_rep4" ~ "PA_4_t600_S_rep3",
                                   TRUE ~ Strain)) %>%
  dplyr::mutate(Sample = Strain)

# Data formatting
# Format data
# Note: now Sample_ID is based on input file names, Sample and Strain columns are correct, re-generate other sample metadata variables
data <- data %>%
  dplyr::select(Sample_ID, Sample, Strain, Ploidy, Chr, Bin, Bin_coverage, Bin_coverage_norm) %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line", "EvoTime", "ColonyMorph", "Replicate"), sep = "_", 
                  remove = FALSE, convert = FALSE) %>%
  tidyr::unite(col = "CondLine", Condition, Line, sep = "", remove = FALSE) %>%
  tidyr::unite(col = "StrainBackground", Condition, Line, EvoTime, sep = "_", remove = FALSE) %>%
  tidyr::unite(col = "StrainBackgroundRep", StrainBackground, Replicate, sep = "_", remove = FALSE) %>%
  tidyr::unite(col = "ColonyRep", ColonyMorph, Replicate, sep = "", remove = FALSE) %>%
  dplyr::mutate(ColonyRep = gsub("rep", "", ColonyRep, fixed = TRUE)) %>%
  dplyr::relocate(ColonyRep, .before = Condition) %>%
  dplyr::mutate(StrainBackground = StrainBackground %>%
                  gsub("_t", " t", ., fixed = TRUE) %>%
                  gsub("_", "", ., fixed = TRUE)) %>%
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph)
samples <- unique(data$Sample)
strains <- unique(data$Strain)
#samples <- c(ref_sample, samples[samples != ref_sample])
samples_label <- samples %>%  # PA_1_t600_D_rep1 -> PA1 t600 D1
  gsub("PA_", "PA", ., fixed = TRUE) %>%
  gsub("_t", " t", ., fixed = TRUE) %>%
  gsub("_rep", "", ., fixed = TRUE) %>%
  gsub("_", " ", ., fixed = TRUE)
condlines <- unique(data$CondLine)
strainbackgrounds <- unique(data$StrainBackground)
strainbackgroundreps <- unique(data$StrainBackgroundRep)
colonyreps <- unique(data$ColonyRep)
data <- data %>%
  dplyr::mutate(Sample = factor(Sample, levels = samples), 
                Strain = factor(Strain, levels = strains), 
                CondLine = factor(CondLine, levels = condlines), 
                StrainBackground = factor(StrainBackground, levels = strainbackgrounds),
                StrainBackgroundRep = factor(StrainBackgroundRep, levels = strainbackgroundreps), 
                ColonyRep = factor(ColonyRep, levels = colonyreps), 
                Line = factor(Line, levels = lines), 
                EvoTime = factor(EvoTime, levels = evotimes), 
                Replicate = factor(Replicate, levels = reps), 
                ColonyMorph = factor(ColonyMorph, levels = colonymorphs)
  )

# Manually adjust Ploidy value
# metadata.csv is only a preliminary guess
# Here aim to get most chr copy numbers in later bar plot to be around integers and match experimentally-measured ploidy level
# Run after data formatting
# Check if Ploidy_round later is 5 or 3 or any other unusual value
data <- data %>%
  dplyr::mutate(Ploidy = case_when(
    Sample == "PA_1_t600_S_rep2" ~ Ploidy + 0.05, #
    Sample == "PA_1_t1000_D_rep2" ~ Ploidy + 0.05, #
    Sample == "PA_1_t1000_S_rep2" ~ Ploidy + 0.1, #
    Sample == "PA_1_t1000_D_rep3" ~ Ploidy + 0.05, #
    Sample == "PA_1_t1000_S_rep3" ~ Ploidy + 0.1, #
    Sample == "PA_2_t600_D_rep1" ~ Ploidy + 0.1, #
    Sample == "PA_2_t600_D_rep2" ~ Ploidy + 0.1, #
    Sample == "PA_2_t600_D_rep3" ~ Ploidy + 0.1, #
    Sample == "PA_2_t1000_D_rep1" ~ Ploidy + 0.05, #
    Sample == "PA_2_t1000_D_rep2" ~ Ploidy + 0.05, #
    Sample == "PA_2_t1000_D_rep3" ~ Ploidy + 0.05, #
    Sample == "PA_3_t1000_S_rep2" ~ Ploidy - 0.05, #
    Sample == "PA_4_t600_D_rep1" ~ Ploidy + 0.1, #
    Sample == "PA_4_t600_D_rep2" ~ Ploidy + 0.15, #
    Sample == "PA_4_t600_S_rep2" ~ Ploidy + 0.2, #
    Sample == "PA_4_t600_D_rep3" ~ Ploidy + 0.1, #
    Sample == "PA_4_t600_S_rep3" ~ Ploidy + 0.21, #
    Sample == "PA_4_t1000_D_rep1" ~ Ploidy + 0.1, #
    Sample == "PA_4_t1000_S_rep1" ~ Ploidy + 0.1, #
    Sample == "PA_4_t1000_D_rep2" ~ Ploidy + 0.1, #
    Sample == "PA_4_t1000_S_rep2" ~ Ploidy + 0.1, #
    Sample == "PA_4_t1000_D_rep3" ~ Ploidy + 0.1, #
    Sample == "PA_5_t600_D_rep1" ~ Ploidy - 0.1, #
    Sample == "PA_5_t600_S_rep1" ~ Ploidy - 0.1, #
    Sample == "PA_5_t600_D_rep2" ~ Ploidy - 0.1, #
    Sample == "PA_5_t600_S_rep2" ~ Ploidy - 0.1, #
    Sample == "PA_5_t600_D_rep3" ~ Ploidy - 0.05, #
    Sample == "PA_5_t600_S_rep3" ~ Ploidy - 0.1, #
    Sample == "PA_5_t1000_D_rep2" ~ Ploidy - 0.05, #
    Sample == "PA_5_t1000_S_rep2" ~ Ploidy - 0.15, #
    TRUE ~ Ploidy
  ))

# Check
length(unique(data$Strain))
unique(data$Strain)
summary(data)
id_vars <- c("Sample_ID", "Sample", 
             "Strain", "CondLine", "StrainBackground", "StrainBackgroundRep", "ColonyRep", 
             "Condition", "Line", "EvoTime", "Replicate", "ColonyMorph")
id_vars_2 <- c(id_vars, "Ploidy")

# Check each sample: coverage, bin number
data %>% 
  dplyr::group_by(Sample) %>%
  dplyr::summarise(Median_coverage = median(Bin_coverage), 
                   Bin_count = n()) %>%
  #dplyr::arrange(Median_coverage) %>%
  #dplyr::arrange(Bin_count) %>%
  View()
## Median coverage: 25-85
## Bin number: two 12079, others around 12000

# # Get fold change relative to reference sample
# ref <- data %>%
#   dplyr::filter(Sample == ref_sample) %>%
#   dplyr::select(Chr, Bin, Bin_coverage, Bin_coverage_norm) %>%
#   dplyr::rename(Bin_coverage_ref = Bin_coverage, 
#                 Bin_coverage_norm_ref = Bin_coverage_norm)
# data_fc <- data.frame()
# for (sample in samples) {  # including reference sample
#   data_fc_sample <- 
#     dplyr::inner_join(x = dplyr::filter(data, Sample == sample), 
#                       y = ref, 
#                       by = c("Chr", "Bin")) %>%
#     dplyr::mutate(Bin_coverage_fc = Bin_coverage / Bin_coverage_ref, 
#                   Bin_coverage_norm_fc = Bin_coverage_norm / Bin_coverage_norm_ref)
#   data_fc <- bind_rows(data_fc, data_fc_sample)
# }
# data <- data_fc
# rm(data_fc, data_fc_sample, ref)

# Function for plotting
plot_bin_cov <- function(data, yvar) {
  ggplot(data, aes(x = Bin, y = .data[[yvar]])) +
    facet_grid(Sample~Chr, scales = "free", space = "free_x") +
    geom_point(size = 0.2, color = "skyblue4") + 
    ggplot_custom_theme +
    theme(panel.spacing = unit(0, "lines"), 
          strip.text.y = element_text(angle = 0), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = rel(0.5)))
}

# Plot bin coverage
plot_bin_cov(data, "Bin_coverage")
## A few bins with exccessive coverage in Chr XII (rDNA locus 451-468 1kb bins)
save_ggplot("bin_coverage", width = 15, height = length(samples) * 0.5)

# Plot bin coverage norm after removing rDNA locus (Chr XII 451-468 1kb bins)
plot_bin_cov(data %>% 
               dplyr::filter(!( Chr == "XII" & 451 <= Bin & Bin <= 468 )) %>%
               dplyr::filter(Bin_coverage_norm < 2.5),  # block out outliers for showing plot more clearly
             "Bin_coverage_norm")
save_ggplot("bin_coverage_norm_remove_rDNA", width = 15, height = length(samples) * 0.5)
## Check if smiley coverage

# # Plot FC of normalized bin coverage
# plot_bin_cov(data, "Bin_coverage_norm_fc") 
# plot_bin_cov(data, "Bin_coverage_norm_fc") +
#   scale_y_continuous(trans = "log2")  # to see more clearly (esp. bins with coverage FC < 1)
# ## Centered around 1 (due to normalization and FC), some outlier bins
# save_ggplot("bin_coverage_norm_fc_log2", width = 15, height = length(samples) * 0.5)

# Note: in our Y55 strains, there is a coverage valley in 14-23 1kb bins, suggesting loss in this region

# No filtering is required for this version of script
data_filter <- data

# Use sample ploidy to estimate bin copy number
data_filter <- data_filter %>%
  dplyr::mutate(Bin_copy_num = Bin_coverage_norm * Ploidy)

# Estimate chromosome copy number (median per-bin coverage of each chromosome)
data_chr <- data_filter %>%
  dplyr::group_by_at(c(id_vars_2, "Chr")) %>%
  dplyr::summarise(Chr_copy_num = median(Bin_copy_num), 
                   Chr_copy_num_round = round(Chr_copy_num)) %>%  # rounded to nearest integer
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph)

# Plot bin copy number
# Bin copy numbers should lie around integer lines
bin_copy_num_cutoff <- 12  # above which show triangle
data_bin_plot <- data_filter %>%
  dplyr::mutate(Sample = factor(Sample, levels = samples, labels = samples_label)) %>%
  dplyr::mutate(Type = if_else(Bin_copy_num <= bin_copy_num_cutoff, "Below", "Above"), 
                Type = factor(Type, levels = c("Below", "Above")), 
                Bin_copy_num = if_else(Bin_copy_num <= bin_copy_num_cutoff, Bin_copy_num, bin_copy_num_cutoff))
data_bin_plot_ploidy <- data_bin_plot %>%
  dplyr::ungroup() %>%
  dplyr::select(all_of(c("Sample", metadata_vars))) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Ploidy_round = round(Ploidy))  # rounded to "baseline" ploidy
ggplot(data_bin_plot, 
       aes(x = Bin, y = Bin_copy_num, color = Chr, shape = Type)) +
  facet_grid(Sample~Chr, scales = "free", space = "free_x") +
  geom_point(size = 0.2) +
  scale_color_manual(values = chrs_color) +
  scale_shape_manual(values = c("Below" = "circle", "Above" = "triangle down open")) + 
  scale_y_continuous(breaks = seq(2,8,2), labels = seq(2,8,2)) +
  geom_hline(yintercept = seq(2,8,2), size = 0.25, color = "gray60") +
  geom_hline(yintercept = seq(1,7,2), size = 0.25, color = "gray85") +
  geom_hline(data = data_bin_plot_ploidy, 
             mapping = aes(yintercept = Ploidy_round),  
             size = 1, color = "red2") +
  labs(x = "Chromosomal coordinates", y = "Copy number") + 
  ggplot_custom_theme3 +
  theme(panel.spacing = unit(0, "lines"), 
        strip.text.y = element_text(angle = 0), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = rel(0.75)),
        legend.position = "none", 
        panel.border = element_rect(size = 0.5))
save_ggplot("copy_num_bin", width = 15, height = length(samples) * 1, limitsize = FALSE)
save_ggplot("copy_num_bin", width = 15, height = length(samples) * 1, limitsize = FALSE, mode = "paper")

# Plot chromosome copy number (bar plot)
ggplot(data_chr,
       aes(x = Chr, y = Chr_copy_num)) +
  facet_wrap(~Sample, ncol = 5) +
  geom_col(size = 0.9, fill = "skyblue4", color = NA) +
  scale_y_continuous(breaks = seq(2,8,2), labels = seq(2,8,2), limits = c(0,9)) +
  labs(x = "Chromosome", y = "Copy number") +
  theme(strip.text = element_text(size = rel(1.25)),
        axis.title = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = rel(0.75), angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = rel(1)))
save_ggplot("copy_num_chr_barplot", width = 2 * 5, height = 2 * 10)

# Manual correction of chromosome copy number rounded
# Compare chromosome copy number heatmap below and bin copy number plot above
data_chr$Chr_copy_num_round_correct <- data_chr$Chr_copy_num_round
# data_chr <- dplyr::mutate(data_chr, Chr_copy_num_round_correct = case_when(
#   Sample %in% c("PM_1_t600", "PM_2_t600", "PM_3_t600", "PM_4_t600", "PM_5_t600") & Chr == "XII" ~ 4,
#   TRUE ~ Chr_copy_num_round
# ))

# Plot chromosome copy number rounded corrected (heatmap)
# See the other script

# Save data
write.csv(data_filter, "data_bin.csv", row.names = FALSE)
write.csv(data_chr, "data_chr.csv", row.names = FALSE)
write.csv(data_bin_plot, "data_bin_plot.csv", row.names = FALSE)
write.csv(data_bin_plot_ploidy, "data_ploidy.csv", row.names = FALSE)
saveRDS(data_filter, file = "data_bin.rds")
saveRDS(data_chr, file = "data_chr.rds")
saveRDS(data_bin_plot, file = "data_bin_plot.rds")
saveRDS(data_bin_plot_ploidy, file = "data_ploidy.rds")
save.image(file = "plot.RData")

