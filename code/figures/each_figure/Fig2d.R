library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/chromosome_cnv/20231214_summary_MuLTEE")  ###
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
in_cov_path <- "03_coverage"
in_cov_suffix <- ".coverage"
out_fig_path <- "."
bin_length <- 1000
conditions <- c("PM", "PA")  #, "PO")
conditions_color <- get_pal_colors("Set1")[c(2,1)]  #[c(2,1,3)]  # like in Bozdag2021, blue, red, green # "#377EB8" "#E41A1C" "#4DAF4A"
names(conditions_color) <- conditions
#show_col(conditions_color, labels = TRUE)
lines <- as.character(1:5)
evotimes <- c("t0", "t200", "t400", "t600", "t1000")
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

# Data formatting
data <- data %>%
  dplyr::arrange(factor(Condition, levels = conditions), 
                 factor(EvoTime, levels = evotimes),
                 factor(Line, levels = lines))
samples <- unique(data$Sample)
#samples <- c(ref_sample, samples[samples != ref_sample])
samples_label <- samples %>% 
  gsub("_t", " t", ., fixed = TRUE) %>%
  gsub("_", "", ., fixed = TRUE)
data <- data %>%
  dplyr::arrange(factor(Sample, levels = samples), 
                 factor(Condition, levels = conditions), 
                 factor(EvoTime, levels = evotimes), 
                 factor(Line, levels = lines))
data$Sample <- factor(data$Sample, levels = samples)
data$Condition <- factor(data$Condition, levels = conditions)
data$Line <- factor(data$Line, levels = lines)
data$EvoTime <- factor(data$EvoTime, levels = evotimes)
summary(data)
id_vars <- c("Sample_ID", "Sample", "Condition", "Line", "EvoTime")
id_vars_2 <- c(id_vars, "Ploidy")

# Check each sample: coverage, bin number
data %>% 
  dplyr::group_by(Sample) %>%
  dplyr::summarise(Median_coverage = median(Bin_coverage), 
                   Bin_count = n()) %>%
  #dplyr::arrange(Median_coverage) %>%
  #dplyr::arrange(Bin_count) %>%
  View()
## Median coverage: 30-200
## Bin number: some 12079, others around 12000

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
                   Chr_copy_num_round = round(Chr_copy_num))  # rounded to nearest integer

# Plot bin copy number
# Bin copy numbers should lie around integer lines
bin_copy_num_cutoff <- 12  # above which show triangle
data_bin_plot <- data_filter %>%
  #dplyr::filter(Sample != ref_sample) %>%
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
save_ggplot("copy_num_bin", width = 15, height = length(samples) * 1)
save_ggplot("copy_num_bin", width = 15, height = length(samples) * 1, mode = "paper")

# Plot chromosome copy number (bar plot)
data_chr_barplot <- data_chr #%>% 
  #dplyr::filter(Sample != ref_sample) %>%
ggplot(data_chr_barplot,
       aes(x = Chr, y = Chr_copy_num)) +
  facet_wrap(~Sample, ncol = 7) +
  geom_col(size = 0.9, fill = "skyblue4", color = NA) +
  scale_y_continuous(breaks = seq(2,8,2), labels = seq(2,8,2), limits = c(0,9)) +
  labs(x = "Chromosome", y = "Copy number") +
  theme(strip.text = element_text(size = rel(1.25)),
        axis.title = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = rel(0.75), angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = rel(1)))
save_ggplot("copy_num_chr_barplot", width = 2 * 7, height = 2 * 6)

# Manual correction of chromosome copy number rounded
# Compare chromosome copy number heatmap below and bin copy number plot above
data_chr$Chr_copy_num_round_correct <- data_chr$Chr_copy_num_round
# data_chr <- dplyr::mutate(data_chr, Chr_copy_num_round_correct = case_when(
#   Sample %in% c("PM_1_t600", "PM_2_t600", "PM_3_t600", "PM_4_t600", "PM_5_t600") & Chr == "XII" ~ 4,
#   TRUE ~ Chr_copy_num_round
# ))

# Plot chromosome copy number rounded corrected (heatmap)
data_chr_heatmap <- data_chr %>%
  #dplyr::filter(Sample != ref_sample) %>%
  dplyr::mutate(Sample = factor(Sample, levels = samples, labels = samples_label))

# # Plot by sample
# ggplot(data_chr_heatmap, 
#        aes(x = Chr, y = Sample, fill = factor(as.integer(Chr_copy_num_round_correct)))) +
#   geom_tile(color = "black") + 
#   scale_fill_manual(name = "Copy\nnumber", 
#                     values = chr_copy_nums_color, 
#                     breaks = min(data_chr_heatmap$Chr_copy_num_round_correct):max(data_chr_heatmap$Chr_copy_num_round_correct)) +
#   scale_y_discrete(limits = rev) +
#   coord_fixed(ratio = 4) +
#   labs(x = "Chromosome", y = NULL) +
#   ggplot_custom_theme3 +
#   theme(panel.border = element_blank(),
#         axis.ticks = element_blank(), 
#         legend.title = element_text(size = rel(1.5)), 
#         legend.text = element_text(size = rel(1.5)), 
#         axis.text.x = element_text(size = rel(0.4), angle = 90, hjust = 1, vjust = 0.5))
# save_ggplot("copy_number_chr_heatmap", 
#             width = 6, height = 10)

# Plot by evotimes
ggplot(data_chr_heatmap,
       aes(x = Chr, y = Line, fill = factor(as.integer(Chr_copy_num_round_correct)))) +
  facet_grid(EvoTime ~ Condition, scales = "free_y", space = "free_y") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "Copy\nnumber",
                    values = chr_copy_nums_color,
                    breaks = min(data_chr_heatmap$Chr_copy_num_round_correct):max(data_chr_heatmap$Chr_copy_num_round_correct)) +
  scale_y_discrete(limits = rev, breaks = lines) +
  #coord_fixed(ratio = 4) +
  labs(x = "Chromosome", y = "Line") +
  ggplot_custom_theme3 +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = rel(1.25)),
        legend.text = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = rel(0.6), angle = 90, hjust = 1, vjust = 0.5, 
                                   margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(hjust = 0.46))
save_ggplot("copy_number_chr_heatmap", 
            width = length(conditions) * 2.5 + 0.5, 
            height = length(evotimes) * length(lines) * 1/5 + 0.5)
save_ggplot("copy_number_chr_heatmap", 
            width = length(conditions) * 2.5 + 0.5, 
            height = length(evotimes) * length(lines) * 1/5 + 0.5, 
            mode = "paper")

# Plot by line
lines_label <- c("", lines)
names(lines_label) <- c("0", lines)
ggplot(data_chr_heatmap %>% 
         dplyr::mutate(Line = as.character(Line),
                       Line = if_else(is.na(Line), "0", Line)),
       aes(x = Chr, y = EvoTime, fill = factor(as.integer(Chr_copy_num_round_correct)))) +
  facet_grid(Line ~ Condition, scales = "free_y", space = "free_y", 
             labeller = labeller(Line = lines_label), switch = "y") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "Copy\nnumber",
                    values = chr_copy_nums_color,
                    breaks = min(data_chr_heatmap$Chr_copy_num_round_correct):max(data_chr_heatmap$Chr_copy_num_round_correct)) +
  scale_y_discrete(limits = rev, breaks = evotimes, position = "right") +
  #coord_fixed(ratio = 4) +
  labs(x = "Chromosome", y = NULL) +
  ggplot_custom_theme3 +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = rel(1.25)),
        legend.text = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = rel(0.6), angle = 90, hjust = 1, vjust = 0.5, 
                                   margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(hjust = 0.46))
save_ggplot("copy_number_chr_heatmap_by_line", 
            width = length(conditions) * 2.5 + 0.5, 
            height = length(evotimes) * length(lines) * 1/5 + 0.5)
save_ggplot("copy_number_chr_heatmap_by_line", 
            width = length(conditions) * 2.5 + 0.5, 
            height = length(evotimes) * length(lines) * 1/5 + 0.5, 
            mode = "paper")

# Save data
write.csv(data_filter, "data_bin.csv", row.names = FALSE)
write.csv(data_chr, "data_chr.csv", row.names = FALSE)
saveRDS(data_filter, file = "data_bin.rds")
saveRDS(data_chr, file = "data_chr.rds")

######################

# ggplot(data_chr_heatmap,
#        aes(x = Chr, y = Line, fill = factor(as.integer(Chr_copy_num_round_correct)))) +
#   facet_grid(Condition ~ EvoTime) +  # , scales = "free_y", space = "free_y"
#   geom_tile(color = "black") +
#   scale_fill_manual(name = "Copy\nnumber",
#                     values = chr_copy_nums_color,
#                     breaks = min(data_chr_heatmap$Chr_copy_num_round_correct):max(data_chr_heatmap$Chr_copy_num_round_correct)) +
#   scale_y_discrete(limits = rev) +
#   coord_fixed(ratio = 4) +
#   labs(x = "Chromosome", y = "Line") +
#   ggplot_custom_theme3 +
#   theme(panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title = element_text(size = rel(1.5)),
#         legend.text = element_text(size = rel(1.5)),
#         axis.text.x = element_text(size = rel(0.4), angle = 90, hjust = 1, vjust = 0.5, 
#                                    margin = margin(t = 0, r = 0, b = 0, l = 0)))
# save_ggplot("copy_number_chr_heatmap", 
#             width = length(evotimes) * 1.5 + 1, 
#             height = length(conditions) * length(lines) * 1/3 + 0.5)
# save_ggplot("copy_number_chr_heatmap", 
#             width = length(evotimes) * 1.5 + 1, 
#             height = length(conditions) * length(lines) * 1/3 + 0.5, 
#             mode = "paper")
