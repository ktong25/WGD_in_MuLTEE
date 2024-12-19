library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(stringr)
library(glue)
library(RColorBrewer)
library(scales)
library(ggh4x)
library(ggpubr)
library(ggforce)
library(ggbreak)

setwd("~/Documents/projects/R_data_analysis/mutations/20231230_donutspread")  ###
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")
load("plot.RData")  # load the entire workspace saved after plot.R run
rm(data_chr_diff, data_ds_diff_dchr)

# Parameters
in_chr_cnv_folder <- "~/Documents/projects/R_data_analysis/chromosome_cnv/20231230_donutspread"
read_chr_cnv_data <- function(rds_filename) {
  readRDS(file.path(in_chr_cnv_folder, rds_filename)) 
}
variant_impacts <- c("High", "Moderate", "Low", "Modifier")
variant_impacts_color <- viridis_pal(option = "mako", begin = 0.2, end = 1, direction = -1)(9)[seq(8,2,-2)]   # colorRampPalette(c("#1B7837", "white"))(25)[seq(1,19,6)]
names(variant_impacts_color) <- variant_impacts
#show_col(variant_impacts_color)
variant_sign2s <- c("Maintain", "Gain", "Loss")
dCCN_sign2s <- c("Maintain", "Gain", "Loss")
defreq_modes <- c("Gain", "Loss", "LOH", "Increase", "Decrease")
defreq_mode2s <- c(defreq_modes, "Same")
dACNs_color <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(9)[c(1:3,5,7:9)]  #(9)[2:8] # [c(11 - rev(c(2,4,6,8,9)), 11, 11 + c(2,4,6,8,9))]
dACNs_color[4] <- "gray80"
names(dACNs_color) <- seq(-3,3)
#show_col(dACNs_color)
np_combine_levels <- c("Gain", "Random mutation", "Loss", "LOH", "Increase", "Decrease", "Same", "Donut")
strainbackgrounds_color <- c("#FBAEAD", "#E9474A", "#B8D8E9", "#4D93C3", "#5CB356", 
                             "#D5C1DE", "#8764AE", "#FDCC8C", "#FF9933")  #get_pal_colors("Paired")[c(5,6,1,2,4,9,10,7,8)] # t600 light, t1000 dark
names(strainbackgrounds_color) <- strainbackgrounds_label
#show_col(strainbackgrounds_color)

##### Prepare data 

# Exclude two 3N strains for all subsequent analyses
# Include PA4 for all subsequent analyses

### Check if estimation of allele copy number is appropriate (ignore segmental aneuploidy for now)
# Exclude 3N strains
# Add chromosome information (if you need to account for segmental aneuploidy, this is where you need to adjust)
length(unique(data$StrainBackgroundRep))
data_chr <- data %>%
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::left_join(read_chr_cnv_data("data_ds.rds") %>% dplyr::select(Strain, Chr, CCN), 
                   by = c("Strain", "Chr")) #%>% 
  #dplyr::filter(Alt_freq > 1.5/CCN)
length(unique(data_chr$StrainBackgroundRep))
# Plot mutation allele frequency and compare with chromosome copy number
p <- ggplot(data_chr, aes(x = factor(CCN), y = Alt_freq))
ccns <- sort(unique(data_chr$CCN))
for (i in seq_along(ccns)) {
  xpos <- i
  ccn <- ccns[i]
  for (n in 0:ccn) {
    p <- p + geom_segment(x = xpos-0.475, xend = xpos+0.475, y = n/ccn, yend = n/ccn, 
                          color = "gray", linewidth = 0.25, linetype = "dashed")
  }
}
p <- p + 
  geom_violin(scale = "width", trim = TRUE, size = 0.5, color = "black", fill = "white", width = 0.7) + 
  geom_violin(scale = "width", trim = TRUE, size = 0.5, color = "black", fill = "skyblue3", alpha = 0.2, width = 0.7) + 
  geom_sina(scale = "width", seed = 42, size = 0.1,  color = "skyblue3", maxwidth = 0.7) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c(0L, 0.25, 0.5, 0.75, 1L), 
                     expand = expansion(mult = 0.025)) +
  labs(x = "Chromosome copy number", y = "Allele frequency") +
  ggplot_custom_theme4
p
save_ggplot("allele_frequency_CCN", width = 6, height = 4)
save_ggplot("allele_frequency_CCN", width = 6, height = 4, mode = "paper")

### Compare mutations between donut and spread
# Here still include 3N strains, exclude later
data_diff <- data %>%
  dplyr::select(-Sample_ID, -Strain, -ColonyRep, -Manual_add, 
                -Ref_Alt_depth, -Ref_depth, -Alt_depth, -Total_depth) %>%
  tidyr::spread(key = "ColonyMorph", value = "Alt_freq", fill = 0) %>%  # if gain/loss then assign as 0
  dplyr::rename(D_freq = D, S_freq = S) %>%
  dplyr::mutate(DS_dfreq = S_freq - D_freq, 
                DS_dmut_sign = ifelse(D_freq != 0 & S_freq != 0, "Maintain", 
                                      ifelse(D_freq == 0, "Gain", "Loss")) %>% factor(levels = variant_sign2s))
summary(data_diff)  # 1107 rows, 14 gain, 34 loss
sum(data_diff$D_freq == 0 & data_diff$S_freq == 0)  # 0
length(unique(data_diff$StrainBackgroundRep))
# Exclude 3N strains
# Add chromosome information
data_diff_chr <- data_diff %>%
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::left_join(read_chr_cnv_data("data_ds_diff.rds") %>%
                     dplyr::select(StrainBackgroundRep, Chr, D_CCN, S_CCN, DS_dCCN), 
                   by = c("StrainBackgroundRep", "Chr")) %>%
  dplyr::mutate(DS_dCCN_sign = case_when(DS_dCCN > 0 ~ "Gain",
                                         DS_dCCN < 0 ~ "Loss", 
                                         DS_dCCN == 0 ~ "Maintain") %>% 
                  factor(levels = dCCN_sign2s))
summary(data_diff_chr)
length(unique(data_diff_chr$StrainBackgroundRep))

### Estimate allele copy number (and estimated allele frequency and its change) based on allele frequency and chromosome copy number
# Find the closest allele copy number that fits allele frequency and chromosome copy number
# Exclude 3N strains
estimate_ACN <- function(freq, CCN, return_NA_if_tie = TRUE, other_ACN = NULL) {
  if (freq == 0) {  # if mutation loss
    return(0)
  } else if (freq == 1) {  # if mutation LOH
    return(CCN)
  }
  possible_ACN <- 1:(CCN-1)
  differences <- abs(freq - (possible_ACN / CCN))
  min_indices <- which(differences == min(differences))
  if (length(min_indices) > 1) {  # If there is more than one index (tie)
    if (return_NA_if_tie) {
      return(NA_integer_)
    } else {  # return the allele copy number with less change from the other allele copy number
      tie_ACNs <- possible_ACN[min_indices]
      differences <- abs(tie_ACNs - other_ACN)
      min_index <- which(differences == min(differences))
      return(tie_ACNs[min_index])
    }
  } else {
    return(possible_ACN[min_indices])
  }
}
data_diff_chr_acn <- data_diff_chr %>%
  dplyr::rowwise() %>%
  dplyr::mutate(D_ACN = estimate_ACN(D_freq, D_CCN), 
                S_ACN = estimate_ACN(S_freq, S_CCN), 
                D_ACN_tie = ifelse(is.na(D_ACN), TRUE, FALSE), 
                S_ACN_tie = ifelse(is.na(S_ACN), TRUE, FALSE)
                )
summary(data_diff_chr_acn) # 5 D_ACN ties, 2 S_ACN ties; 985 maintain, 13 gain, 20 loss, 1018 DS mutation pairs (rows)
sum(data_diff_chr_acn$D_freq != 0) + sum(data_diff_chr_acn$S_freq != 0) # 2003 mutations, i.e. 1018 * 2 - 13 gain - 20 loss
length(unique(data_diff_chr_acn$StrainBackgroundRep))
# Check
# data_diff_chr_acn %>%
#   dplyr::filter(D_ACN_tie | S_ACN_tie) %>%
#   View()
# Correct for ties
data_diff_chr_acn <- data_diff_chr_acn %>%
  dplyr::rowwise() %>%
  dplyr::mutate(D_ACN = ifelse(D_ACN_tie, 
                               estimate_ACN(D_freq, D_CCN, return_NA_if_tie = FALSE, other_ACN = S_ACN), 
                               D_ACN), 
                S_ACN = ifelse(S_ACN_tie, 
                               estimate_ACN(S_freq, S_CCN, return_NA_if_tie = FALSE, other_ACN = D_ACN), 
                               S_ACN)
  )
# Check
# data_diff_chr_acn %>%
#   dplyr::filter(D_ACN_tie | S_ACN_tie) %>%
#   View()
# data_diff_chr_acn %>%
#   dplyr::select(DS_dmut_sign, DS_dCCN_sign, D_freq, S_freq, D_CCN, S_CCN, D_ACN, S_ACN) %>%
#   View()
# table(paste(data_diff_chr_acn$D_ACN, data_diff_chr_acn$S_ACN)) %>%
#   data.frame() %>%
#   View()
# data_diff_chr_acn %>%
#   dplyr::filter(abs(D_ACN-S_ACN) > 1) %>%
#   dplyr::select(DS_dmut_sign, DS_dCCN_sign, D_freq, S_freq, D_CCN, S_CCN, D_ACN, S_ACN) %>%
#   View()  # 6 cases, all abs(D_ACN-S_ACN) == 2, some are close to ties, can be manually corrected if needed
# Calculate allele copy number change, estimated allele frequency (based on allele copy number and chromosome copy number) and its change
data_diff_chr_acn <- data_diff_chr_acn %>%
  dplyr::mutate(DS_dACN = S_ACN - D_ACN, 
                DS_dACN_abs = abs(DS_dACN), 
                D_efreq = D_ACN / D_CCN,
                D_efreq_text = paste0(D_ACN, "/", D_CCN),
                S_efreq = S_ACN / S_CCN, 
                S_efreq_text = paste0(S_ACN, "/", S_CCN),
                DS_defreq = S_efreq - D_efreq, 
                DS_defreq_abs = abs(DS_defreq))
# Check
# table(paste(data_diff_chr_acn$D_efreq_text, data_diff_chr_acn$S_efreq_text, data_diff_chr_acn$DS_defreq)) %>%
#   data.frame() %>%
#   View()
# table(paste(data_diff_chr_acn$DS_dACN, data_diff_chr_acn$DS_defreq)) %>%
#   data.frame() %>%
#   View()
# Histogram of estimated allele frequency change 
# ggplot(data_diff_chr_acn %>% dplyr::filter(!DS_defreq_mode %in% c("Gain", "Loss", "LOH")), 
#        aes(x = DS_defreq)) +
#   geom_histogram()  #binwidth = 0.5, fill = "blue", color = "black"

### Classify mutation changes into five types: gain/loss, LOH (<1 to 1), allele frequency increase/decrease
data_diff_chr_acn <- data_diff_chr_acn %>%
  dplyr::mutate(DS_defreq_mode = ifelse(D_efreq == 0 | S_efreq == 0, 
                                        ifelse(D_efreq == 0, "Gain", "Loss"), 
                                        ifelse(D_efreq < 1 & S_efreq == 1, "LOH", 
                                               ifelse(D_efreq == S_efreq, "Same", 
                                                      ifelse(D_efreq < S_efreq, "Increase", "Decrease"))
                                               )) %>% factor(levels = defreq_mode2s))
summary(data_diff_chr_acn)
table(data_diff_chr_acn$DS_defreq_mode)
# Gain     Loss      LOH Increase Decrease     Same 
# 13       20        4      123       86      772 

##### Plot

### Check allele copy number changes (almost all are by 1 or 0 copy)
# Correlate mutation AF changes (gain/loss/LOH/increase/decrease) with ACN changes (color by ACN changes)
# table(paste(data_diff_chr_acn$DS_defreq_mode, data_diff_chr_acn$DS_dACN)) %>%
#   data.frame() %>%
#   View()
# data_diff_chr_acn %>%
#   dplyr::filter(DS_defreq_mode == "Same" & DS_dACN != 0) %>%
#   dplyr::select(DS_dmut_sign, DS_dCCN_sign, D_freq, S_freq, D_CCN, S_CCN, D_ACN, S_ACN, D_efreq, S_efreq, DS_defreq) %>%
#   View()  # All 5 cases: AF=100% in donut and spread, with chromosome/allele copy number +1
dACNs_range <- rev(min(data_diff_chr_acn$DS_dACN):max(data_diff_chr_acn$DS_dACN))
ggplot(data_diff_chr_acn, aes(x = DS_defreq_mode, fill = factor(DS_dACN, levels = dACNs_range))) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = dACNs_color,
                    labels = dACNs_range, 
                    drop = FALSE) +
  scale_y_break(breaks = c(130, 745), expand = expansion(mult = c(0.01,0.05))) + 
  scale_y_continuous(limits = c(0,775), breaks = c(seq(0,125,25), seq(750,775,25))) +
  labs(x = "Change of corrected allele frequency", y = "# of mutations", fill = "∆Allele\ncopy number") +
  ggplot_custom_theme4 + 
  theme(axis.ticks.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.line.y.right = element_blank(), 
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("barplot_mut-mode_allele-copy-num-change", width = 4.5, height = 5)
save_ggplot("barplot_mut-mode_allele-copy-num-change", width = 4.5, height = 5, mode = "paper")
# Correlate mutation AF changes with ACN changes and CCN changes (color by ACN changes)
ggplot(data_diff_chr_acn %>%
         dplyr::mutate(DS_dCCN_sign = factor(DS_dCCN_sign, levels = dCCN_sign2s[c(2,1,3)])), 
       aes(x = DS_dCCN_sign, fill = factor(DS_dACN, levels = dACNs_range))) +
  facet_wrap(~DS_defreq_mode, ncol = 3, scales = "free_y") +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = dACNs_color,
                    labels = dACNs_range, 
                    drop = FALSE) +
  scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0.02,0.05))) +
  labs(x = "Change of chromosome copy number", y = "# of mutations", fill = "∆Allele\ncopy number") +
  ggplot_custom_theme4 +
  theme(legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)), 
        axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("barplot_mut-mode_chr-sign_allele-copy-num-change", width = 6.5, height = 5)
save_ggplot("barplot_mut-mode_chr-sign_allele-copy-num-change", width = 6.5, height = 5, mode = "paper")

### Plot for each replicate
# Number and impact of mutation gain/loss/LOH/increase/decrease/same
ggplot(data_diff_chr_acn %>%
         dplyr::mutate(StrainBackgroundRep = factor(StrainBackgroundRep, 
                                                    levels = strainbackgroundreps_no3N, strainbackgroundreps_label_no3N)), 
       aes(x = StrainBackgroundRep, fill = Variant_impact)) +
  facet_wrap(~DS_defreq_mode, ncol = 1, scales = "free_y") +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = variant_impacts_color) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0,0.05))) +
  labs(x = NULL, y = "# of mutations", fill = "Mutation Impact") +
  ggplot_custom_theme4 + 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)), 
        axis.text.x = element_text(size = rel(0.7), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("barplot_strainrep_mut-mode_impact", width = 6, height = 10)
save_ggplot("barplot_strainrep_mut-mode_impact", width = 6, height = 10, mode = "paper")
# Allele frequency change and impact of mutation increase/decrease
rects <- data.frame(xstart = seq(0.5, 0.5 + 2 * (round(length(strainbackgroundreps_no3N) / 2) - 1), 2),  # 0.5-1.5 gray, 2.5-3.5 gray
                    xend = seq(1.5, 1.5 + 2 * (round(length(strainbackgroundreps_no3N) / 2) - 1), 2))
set.seed(0)
ggplot(data_diff_chr_acn %>%
         dplyr::mutate(StrainBackgroundRep = factor(StrainBackgroundRep, 
                                                    levels = strainbackgroundreps_no3N, strainbackgroundreps_label_no3N)) %>%
         dplyr::filter(DS_defreq_mode %in% c("Increase", "Decrease")), 
       aes(x = StrainBackgroundRep, y = DS_defreq_abs, color = Variant_impact)) +
  facet_grid(Variant_impact ~ DS_defreq_mode) + 
  geom_rect(data = rects, 
            mapping = aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, 
            fill = "gray", alpha = 0.4) +
  geom_jitter(width = 0.4, alpha = 0.7, size = 1) +
  scale_color_manual(values = variant_impacts_color) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, NA), n.breaks = 4, 
                     expand = expansion(mult = c(0,0.05))) +
  labs(x = NULL, y = "Absolute value of change of corrected allele frequency", color = "Mutation Impact") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggplot_custom_theme4 + 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)), 
        axis.text.x = element_text(size = rel(0.7), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("mutation_strainrep_mut-increase-decrease_impact_AFchange", width = 12, height = 10)
save_ggplot("mutation_strainrep_mut-increase-decrease_impact_AFchange", width = 12, height = 10, mode = "paper")

### Donut mutation impact distribution
# Exclude 3N strains
# Calculate mutation impact number and percentage
data_donut_all_np <- data %>%
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::filter(ColonyMorph == "D") %>%
  dplyr::mutate(StrainBackgroundRep = factor(StrainBackgroundRep, levels = strainbackgroundreps_no3N)) %>%
  dplyr::group_by(StrainBackgroundRep, Variant_impact, .drop = FALSE) %>%
  dplyr::summarise(Mut_num = n()) %>%
  dplyr::mutate(Mut_num_total = sum(Mut_num), 
                Mut_perc = Mut_num / Mut_num_total * 100)
# Plot
# xlabel (show mutation number): e.g., PA1 t600 (D1) (n=39)
xlabels <- data_donut_all_np %>%
  dplyr::mutate(StrainBackgroundRep = factor(StrainBackgroundRep, 
                                             levels = strainbackgroundreps_no3N, 
                                             labels = strainbackgroundreps_label_no3N %>% 
                                               gsub("(", "(D", ., fixed = TRUE)), 
                StrainBackgroundRep = paste0(StrainBackgroundRep, " (n=", Mut_num_total, ")")) %>%
  dplyr::pull(StrainBackgroundRep) %>%
  unique()
ggplot(data_donut_all_np %>%
         dplyr::mutate(StrainBackgroundRep = factor(StrainBackgroundRep, 
                                                    levels = strainbackgroundreps_no3N, 
                                                    labels = xlabels)),
       aes(x = StrainBackgroundRep, y = Mut_perc, fill = Variant_impact)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = variant_impacts_color) +
  scale_y_continuous(limits = c(0,100.1), expand = expansion(mult = 0)) +
  labs(x = NULL, y = "% of mutations", fill = "Mutation impact") +
  ggplot_custom_theme4 + 
  theme(legend.position = "right", 
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)), 
        axis.text.x = element_text(size = rel(0.65), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("barplot_donut_mut_impact_perc", width = 7, height = 4)
save_ggplot("barplot_donut_mut_impact_perc", width = 7, height = 4, mode = "paper")

### Whole dataset: impact percentages of mutation gain/loss/LOH/increase/decrease
# Calculate numbers and percentages of mutations
data_diff_chr_acn_np <- data_diff_chr_acn %>%
  dplyr::group_by(DS_defreq_mode, Variant_impact, .drop = FALSE) %>%
  dplyr::summarise(Mut_num = n()) %>%
  dplyr::group_by(DS_defreq_mode) %>%
  dplyr::mutate(Mut_num_total = sum(Mut_num), 
                Mut_perc = Mut_num / Mut_num_total * 100)
summary(data_diff_chr_acn_np)
# All donut strains (exclude 3N strains)
data_donut_np <- data %>% ###
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::filter(ColonyMorph == "D") %>%
  dplyr::group_by(Variant_impact) %>%
  dplyr::summarise(Mut_num = n()) %>%
  dplyr::mutate(Mut_num_total = sum(Mut_num), 
                Mut_perc = Mut_num / Mut_num_total * 100) %>%
  dplyr::mutate(DS_defreq_mode = "Donut") %>%
  dplyr::relocate(DS_defreq_mode, .before = Variant_impact)
# Probabilities of nonsense/missense/synonymous/intergenic mutations, after introducing a single-nt substitution into yeast genome
# Note: here low-impact includes synonymous and intronic mutations
# Use S288c genome build R64-3-1_20210421, see random_mutation_prob folder
cds_length <- 9076860
cds_genomic_length <- 9189215
intron_length <- cds_genomic_length - cds_length
genome_length <- 12157105
cds_ratio <- cds_length / genome_length
intron_ratio <- intron_length / genome_length
nonsense_prob <- 0.0536
nonsyn_prob <- 0.7394
syn_prob <- 0.2070
df_impact_probs <- data.frame(
  DS_defreq_mode = "Random mutation",
  Variant_impact = variant_impacts, 
  Mut_perc = c(cds_ratio * c(nonsense_prob, nonsyn_prob, syn_prob) + c(0,0,intron_ratio), 1-cds_ratio-intron_ratio) * 100 
)
sum(df_impact_probs$Mut_perc)
# Combine
data_np_combine <- dplyr::bind_rows(data_diff_chr_acn_np %>% dplyr::mutate(DS_defreq_mode = as.character(DS_defreq_mode)), 
                                    data_donut_np, df_impact_probs) %>%  ###
  dplyr::mutate(DS_defreq_mode = factor(DS_defreq_mode, levels = np_combine_levels), 
                Variant_impact = factor(Variant_impact, levels = variant_impacts)) %>%
  dplyr::arrange(factor(DS_defreq_mode, levels = np_combine_levels))
summary(data_np_combine)
# Compare probability distributions
# Chi-squared test for given probabilities, goodness-of-fit test; Chi-squared approximation may be incorrect
# Compare Gain with Random mutation
x <- data_np_combine %>% dplyr::filter(DS_defreq_mode == "Gain") %>% dplyr::pull(Mut_num)
y <- data_np_combine %>% dplyr::filter(DS_defreq_mode == "Random mutation") %>% dplyr::pull(Mut_perc)
result <- chisq.test(x, p = y / sum(y))
print(result)
df_chi_test1 <- tibble::tribble(
  ~group1, ~group2, ~statistic, ~df, ~p, ~Variant_impact,
  "Gain", "Random mutation", result$statistic, result$parameter, result$p.value, "High",
)
# Compare Loss/LOH/Increase/Decrease/Same with Donut
df_chi_test2 <- tibble::tribble(
  ~group1, ~group2, ~statistic, ~df, ~p, ~Variant_impact,
)
for (m in c("Loss", "LOH", "Increase", "Decrease", "Same")) {
  x <- data_np_combine %>% dplyr::filter(DS_defreq_mode == m) %>% dplyr::pull(Mut_num)
  y <- data_np_combine %>% dplyr::filter(DS_defreq_mode == "Donut") %>% dplyr::pull(Mut_num)
  result <- chisq.test(x, p = y / sum(y))
  print(result)
  df_chi_test2 <- df_chi_test2 %>% 
    add_row(group1 = m, group2 = "Donut", 
            statistic = result$statistic, df = result$parameter, p = result$p.value, 
            Variant_impact = "High")
}
# Compare Donut and Random mutation
# x <- data_np_combine %>% dplyr::filter(DS_defreq_mode == "Donut") %>% dplyr::pull(Mut_num)
# y <- data_np_combine %>% dplyr::filter(DS_defreq_mode == "Random mutation") %>% dplyr::pull(Mut_perc)
# result <- chisq.test(x, p = y / sum(y))
# print(result)
# Plot
xlabels <- unique(paste0(data_np_combine$DS_defreq_mode, " (n=", data_np_combine$Mut_num_total, ")"))
xlabels[xlabels=="Random mutation (n=NA)"] <- "Random mutation"
ggplot(data_np_combine, aes(x = DS_defreq_mode, y = Mut_perc, fill = Variant_impact)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = variant_impacts_color) +
  scale_x_discrete(labels = xlabels, expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(breaks = seq(0,100,20), expand = expansion(mult = c(0,0.1))) +
  labs(x = NULL, y = "% of mutations", fill = "Mutation impact") +
  stat_pvalue_manual(df_chi_test1 %>% dplyr::mutate(p = round(p, 3)), 
                     label = "p = {p}", 
                     y.position = 105, step.increase = 0.14, 
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.3) +
  stat_pvalue_manual(df_chi_test2 %>% dplyr::mutate(p = round(p, 3)), 
                     label = "p = {p}", 
                     y.position = 105, step.increase = 0.14, 
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.3) +
  annotate("segment", x = 0.1, xend = 0.1, y = 0, yend = 101, linewidth = 1) +
  ggplot_custom_theme4 + 
  theme(legend.position = "right", 
        legend.justification = c(1, 0.2),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.line.y = element_blank(), 
        axis.title.y.left = element_text(hjust = 0.15)
  )
save_ggplot("barplot_mut-mode_impact-perc", width = 6, height = 6)
save_ggplot("barplot_mut-mode_impact-perc", width = 6, height = 6, mode = "paper")

### Whole dataset: for mutation increase/decrease, test if mutations with higher impacts tend to have larger AF changes
# Use Kruskal-Wallis rank sum test
# Also plot: for mutation increase/decrease, AF changes in each mutation impact
df_kw_test <- tibble::tribble(
    ~DS_defreq_mode, ~.y., ~.x., ~statistic, ~df, ~p
  )
for (m in c("Increase", "Decrease")) {
  result <- kruskal.test(DS_defreq_abs ~ Variant_impact, 
                         data = data_diff_chr_acn %>% dplyr::filter(DS_defreq_mode == m))
  print(result)  # e.g., Kruskal-Wallis chi-squared = 5.9229, df = 3, p-value = 0.1154
  df_kw_test <- df_kw_test %>% 
    add_row(DS_defreq_mode = m, .y. = "DS_defreq_abs", .x. = "Variant_impact", 
            statistic = result$statistic, df = result$parameter, p = result$p.value)
}
# ggpubr::compare_means(DS_defreq_abs ~ Variant_impact,  
#                       data = data_diff_chr_acn %>% dplyr::filter(DS_defreq_mode == "Increase"), 
#                       method = "kruskal.test")  # do not report statistic and df
set.seed(42)
ggplot(data_diff_chr_acn %>% dplyr::filter(DS_defreq_mode %in% c("Increase", "Decrease")),
       aes(x = Variant_impact, y = DS_defreq_abs, color = Variant_impact)) +
  facet_wrap(~DS_defreq_mode, nrow = 1) +
  geom_jitter(width = 0.3, alpha = 0.7, size = 1) +
  stat_summary(geom = "errorbar", fun = mean, mapping = aes(ymin = after_stat(y), ymax = after_stat(y)),
               width = 0.4, color = "gray70", size = 1) +
  stat_summary(geom = "errorbar", fun.data = "mean_sdl", fun.args = list(mult = 1),
               width = 0.2, color = "gray70", size = 0.75) +
  scale_color_manual(values = variant_impacts_color) + 
  scale_x_discrete(expand = expansion(mult = 0.2)) + 
  scale_y_continuous(limits = c(0, 0.37), expand = expansion(mult = c(0, 0.05))) + 
  labs(x = "Mutation impact", y = "Absolute value of change of\ncorrected allele frequency") +
  ggpubr::stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 5, hjust = 0.1) +
  ggplot_custom_theme4 +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("mutation_mut-increase-decrease_impact_defreq", width = 6, height = 5)
save_ggplot("mutation_mut-increase-decrease_impact_defreq", width = 6, height = 5, mode = "paper")

# Check
table(data_diff_chr_acn %>% dplyr::filter(DS_defreq_mode != "Same") %>% dplyr::pull(DS_dCCN_sign))
# Maintain     Gain     Loss 
# 110       35      101

### For each DS pair: compare high/moderate-impact percentage of AF increase/decrease to donut background
# Not do this for mutation gain/loss/LOH because there are very few mutations per DS pair
data_diff_chr_acn_all_np <- data_diff_chr_acn %>%
  dplyr::mutate(StrainBackgroundRep = factor(StrainBackgroundRep, levels = strainbackgroundreps_no3N)) %>%
  dplyr::group_by(StrainBackgroundRep, DS_defreq_mode, Variant_impact, .drop = FALSE) %>%
  dplyr::summarise(Mut_num = n()) %>%
  dplyr::group_by(StrainBackgroundRep, DS_defreq_mode) %>%
  dplyr::mutate(Mut_num_total = sum(Mut_num), 
                Mut_perc = Mut_num / Mut_num_total * 100)
data_diff_chr_acn_all_np_incdec <- data_diff_chr_acn_all_np %>%  # increase and decrease only
  dplyr::filter(DS_defreq_mode %in% c("Increase", "Decrease")) %>%
  dplyr::filter(Variant_impact %in% c("High", "Moderate")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(StrainBackgroundRep, DS_defreq_mode) %>%
  dplyr::summarize(Mut_num_total = mean(Mut_num_total), 
                   Mut_perc_hm = sum(Mut_perc))  # % of high/moderate-impact mutations
data_donut_all_np_incdec <- data_donut_all_np %>%
  dplyr::filter(Variant_impact %in% c("High", "Moderate")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(StrainBackgroundRep) %>%
  dplyr::summarize(Mut_num = sum(Mut_num), 
                   Mut_num_total = mean(Mut_num_total), 
                   Mut_perc_hm = sum(Mut_perc))
data_all_np_combine_incdec <- dplyr::left_join(
  x = data_diff_chr_acn_all_np_incdec, 
  y = data_donut_all_np_incdec %>% dplyr::rename(Mut_num_donut = Mut_num, 
                                                 Mut_num_total_donut = Mut_num_total,
                                                 Mut_perc_hm_donut = Mut_perc_hm), 
  by = "StrainBackgroundRep"
  ) %>%
  dplyr::mutate(Mut_perc_hm_diff = Mut_perc_hm - Mut_perc_hm_donut) %>%
  tidyr::separate(col = "StrainBackgroundRep", into = c("Condition", "Line", "EvoTime", "Replicate"), sep = "_", 
                  remove = FALSE, convert = FALSE) %>%
  tidyr::unite(col = "StrainBackground", Condition, Line, EvoTime, sep = "_", remove = FALSE) %>%
  dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds))
# Perform statistical tests (one-sample t-test vs. 0, one-sided greater)
data_all_np_combine_incdec_stat <- data_all_np_combine_incdec %>% 
  dplyr::group_by(DS_defreq_mode) %>%
  rstatix::t_test(Mut_perc_hm_diff ~ 1, mu = 0, alternative = "greater", detailed = TRUE)
# Simulate for % of high/moderate-impact mutations in AF increase/decrease
simulate_hm_perc <- function(N, A, M) {
  # N: Total number of mutations
  # A: Number of A mutations
  # M: Number of mutations to select
  objects <- c(rep(1, A), rep(0, N - A))
  selected_objects <- sample(objects, M, replace = FALSE)
  count_A <- sum(selected_objects)
  perc_A <- count_A / M * 100
  return(perc_A)
}
set.seed(1)
data_all_np_combine_incdec_sim <- data_all_np_combine_incdec %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Mut_perc_hm_sim = simulate_hm_perc(Mut_num_total_donut, Mut_num_donut, Mut_num_total), 
                Mut_perc_hm_sim_diff = Mut_perc_hm_sim - Mut_perc_hm_donut)
# Plot with simulated data
data_all_np_combine_incdec_sim_plot <- data_all_np_combine_incdec_sim %>%
  dplyr::select(StrainBackgroundRep, StrainBackground, DS_defreq_mode, Mut_perc_hm_diff, Mut_perc_hm_sim_diff) %>%
  tidyr::pivot_longer(c("Mut_perc_hm_diff", "Mut_perc_hm_sim_diff"), 
                      names_to = "Dataset", values_to = "Mut_perc_hm_diff") %>%
  dplyr::mutate(DS_defreq_mode_label = ifelse(Dataset == "Mut_perc_hm_sim_diff", 
                                              paste0(DS_defreq_mode, "\n(simulated)"), 
                                              as.character(DS_defreq_mode)) %>% 
                  factor(levels = c("Increase", "Increase\n(simulated)", 
                                    "Decrease", "Decrease\n(simulated)")))
set.seed(0)
ggplot(data_all_np_combine_incdec_sim_plot %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, labels = strainbackgrounds_label)),
       aes(x = DS_defreq_mode_label, y = Mut_perc_hm_diff)) +
  geom_jitter(aes(color = StrainBackground), width = 0.2, alpha = 0.5, size = 2.5) +
  stat_summary(geom = "errorbar", fun = mean, mapping = aes(ymin = after_stat(y), ymax = after_stat(y)),
               width = 0.4, color = "gray70", size = 1) +
  stat_summary(geom = "errorbar", fun.data = "mean_sdl", fun.args = list(mult = 1),
               width = 0.2, color = "gray70", size = 0.75) +
  #geom_boxplot(width = 0.5, size = 0.75, outlier.shape = NA, fill = NA, color = "gray60") +
  scale_color_manual(values = strainbackgrounds_color) + 
  labs(x = NULL, y = "% of high/moderate-impact mutations\n(relative to donut background)", color = NULL) +
  ggplot_custom_theme4
save_ggplot("mutation-highmoderate-perc_incdec-vs-donut_sim", width = 7, height = 5)
save_ggplot("mutation-highmoderate-perc_incdec-vs-donut_sim", width = 7, height = 5, mode = "paper")

##### Save data

save.image("analysis.RData")
for (df in c("data", "data_chr", 
             "data_diff", "data_diff_chr", "data_diff_chr_acn", 
             "data_diff_chr_acn_np", "data_donut_np", "data_np_combine", 
             "df_chi_test1", "df_chi_test2", "df_kw_test", 
             "data_donut_all_np", "data_donut_all_np_incdec", 
             "data_diff_chr_acn_all_np", "data_diff_chr_acn_all_np_incdec", 
             "data_all_np_combine_incdec_sim", "data_all_np_combine_incdec_sim_plot", "data_all_np_combine_incdec_stat"
             )) {
  write.csv(get(df), glue("{df}.csv"), row.names = FALSE)
  saveRDS(get(df), file = glue("{df}.rds"))
}

##### Brief checks

# Check mutations with 100% allele frequency in donuts
# data_diff_chr_acn %>%
#   dplyr::filter(D_efreq == 1) %>%
#   # dplyr::select(StrainBackgroundRep, Mut_ID, Variant_type, gene_name) %>%  # some genes appear in only one replicate
#   dplyr::select(StrainBackground, Mut_ID, Variant_type, gene_name, gene_id) %>%
#   dplyr::distinct() %>%
#   dplyr::left_join(y = read.csv("gene_annotations.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
#                      dplyr::select(gene_id, description, gene_biotype),  # gene_name,
#                    by = "gene_id") %>%
#   dplyr::arrange(StrainBackground, gene_name) %>%
#   dplyr::select(-Mut_ID, -gene_id, -gene_biotype) %>%
#   View()

# Check the number of point mutation AF changes with/without chromosome copy number changes
# data_diff_chr_acn %>%
#   dplyr::filter(!DS_defreq_mode %in% c("Same", "Gain")) %>%
#   dplyr::pull(DS_dCCN_sign) %>%
#   table()
# Maintain     Gain     Loss 
# 98       35      100
