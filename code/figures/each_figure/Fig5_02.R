library(dplyr)
library(ggplot2)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)
library(tidyr)
library(purrr)
library(cowplot)

setwd("~/Documents/projects/R_data_analysis/chromosome_cnv/20231230_donutspread")  ###
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")
load("plot.RData")  # load the entire workspace saved after plot.R run

# Parameters
strainbackgroundreps_3N <- c("PA_4_t600_rep1", "PA_4_t1000_rep3")
DS_dCCN_modes <- c("T", "A")  # change towards or away from 4
DS_dCCN_modeTs <- c("4", "n4")
#DS_dCCN_modeTs_color <- c("red", "black")
#names(DS_dCCN_modeTs_color) <- DS_dCCN_modeTs
chrs <- as.character(as.roman(1:16))

# Load and re-format data
# Data for donut/spread
data_ds <- data_chr %>%
  dplyr::ungroup() %>%
  dplyr::select(-Sample_ID, -Sample, -Chr_copy_num, -Chr_copy_num_round) %>%  # Sample and Strain columns are equivalent
  dplyr::rename(CCN = Chr_copy_num_round_correct) %>%
  dplyr::mutate(StrainBackground = paste(Condition, Line, EvoTime, sep = "_"))
strainbackgrounds <- unique(data_ds$StrainBackground)
data_ds$StrainBackground <- factor(data_ds$StrainBackground, levels = strainbackgrounds)
strainbackgrounds_label <- strainbackgrounds %>%
  gsub("_t", " t", ., fixed = TRUE) %>%
  gsub("_", "", ., fixed = TRUE)
summary(data_ds)
# Data for MuLTEE evolved isolates
data_evo <- readRDS("~/Documents/projects/R_data_analysis/chromosome_cnv/20231214_summary_MuLTEE/data_chr.rds") %>%
  dplyr::ungroup() %>%
  dplyr::select(-Sample_ID, -Chr_copy_num, -Chr_copy_num_round) %>%
  dplyr::rename(CCN = Chr_copy_num_round_correct, Strain = Sample) %>%
  dplyr::mutate(Strain = paste(Condition, Line, EvoTime, sep = "_"))
strains_evo <- unique(data_evo$Strain)
data_evo$Strain <- factor(data_evo$Strain, levels = strains_evo)
summary(data_evo)
# id_vars
id_vars <- c("Strain", "CondLine", "StrainBackground", "StrainBackgroundRep", "ColonyRep", 
             "Condition", "Line", "EvoTime", "Replicate", "ColonyMorph")
id_vars_diff <- id_vars[!id_vars %in% c("Strain", "ColonyRep", "ColonyMorph")]

# Parameters
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
strainbackgrounds_color <- c("#FBAEAD", "#E9474A", "#B8D8E9", "#4D93C3", "#5CB356", 
                             "#D5C1DE", "#8764AE", "#FDCC8C", "#FF9933")  #get_pal_colors("Paired")[c(5,6,1,2,4,9,10,7,8)] # t600 light, t1000 dark
names(strainbackgrounds_color) <- strainbackgrounds_label
#show_col(strainbackgrounds_color, labels = TRUE)
dCCN_signs <- c("Gain", "Loss")
dCCN_signs_color <- c("#FF9BB1", "#A9C7F0") # c("#F4B3C1", "#BFD1F0") # alpha = 0.7 # c("#EA738D", "#89ABE3")  # Bubblegum Pink, Sky Blue
names(dCCN_signs_color) <- dCCN_signs
#show_col(dCCN_signs_color, labels = TRUE)
dCCN_signs_evods <- c(paste0(c("Gain", "Loss"), "\n(evo)"), paste0(c("Gain", "Loss"), "\n(DS)"))
dCCN_signs_evods_color <- c(dCCN_signs_color, dCCN_signs_color) #c("#FFE7ED", "#E9F1FB", "#F394A8", "#A0BCE3") # c("#FFEBEF", "#EEF4FC", "#FF9BB1", "#A9C7F0")  # alpha = 0.2, alpha = 1
names(dCCN_signs_evods_color) <- dCCN_signs_evods
#show_col(dCCN_signs_evods_color, labels = TRUE)
DS_dCCN_mode2s <- c("Tn4", "T4", "An4", "A4")
DS_dCCN_mode2s_label <- c("Towards 4", "To 4", "Away from 4", "From 4")
DS_dCCN_mode2s_color <- c("#E7D1CE", "#E3A29A", "#C0D6CA", "#78ACA8")
names(DS_dCCN_mode2s_color) <- DS_dCCN_mode2s_label
#show_col(DS_dCCN_mode2s_color, labels = TRUE)

# Parameters
colonymorphs_shape2 <- c("circle", "bullet")
names(colonymorphs_shape2) <- colonymorphs

# Calculate the karyotype difference between each donut/spread and its corresponding evolved isolate (strain background)
# When plotting heatmap, label only the donut non-zero difference, label 1 as "+1" and -1 as "-1" etc
data_ds <- purrr::map_dfr(strains, function (strain) {
  df <- data_ds %>% dplyr::filter(Strain == strain)
  strainbackground <- unique(df$StrainBackground) %>% as.character()
  df$StrainBackground_CCN <- data_evo %>%
    dplyr::filter(Strain == strainbackground) %>%
    dplyr::pull(CCN)
  df$StrainBackground_dCCN <- df$CCN - df$StrainBackground_CCN
  return(df)
}) %>%
  dplyr::mutate(
    StrainBackground_dCCN_label = ifelse(ColonyMorph == "D" & StrainBackground_dCCN != 0, 
                                             ifelse(StrainBackground_dCCN > 0, 
                                                    paste0("+", StrainBackground_dCCN), 
                                                    as.character(StrainBackground_dCCN)), 
                                             NA)
  )

# Custom theme for chr copy number heatmaps
CCN_heatmap_theme <- ggplot_custom_theme3 +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        strip.placement = "outside", 
        strip.text.x = element_blank(), 
        strip.text.y = element_text(size = rel(0.8)), 
        axis.title.x = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(0.7), angle = 90, hjust = 1, vjust = 0.5, 
                                   margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.text.y = element_text(size = rel(0.7), hjust = 0), 
        legend.position = "bottom", 
        legend.title = element_text(size = rel(1.4)),
        legend.text = element_text(size = rel(0.9)),
        legend.key.size = unit(12, "pt"),
        legend.margin = margin(t = 0, b = 0), 
        title = element_text(size = rel(0.65)))

# Plot heatmap
# Show donut-background difference as text in heatmap grid
CCNs_color <- chr_copy_nums_color
CCNs_color[4] <- "gray95"
CCNs_range <- min(data_ds$CCN):max(data_ds$CCN)
p1 <- ggplot(data_ds %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, 
                                                 labels = gsub(" ", "\n", strainbackgrounds_label, fixed = TRUE)), 
                       CCN = factor(CCN, levels = CCNs_range)),
       aes(x = Chr, y = ColonyRep, fill = CCN)) +
  facet_grid(StrainBackground ~ Condition, 
             scales = "free_y", space = "free_y", switch = "y") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "Copy\nnumber\n",
                    values = CCNs_color,
                    labels = paste0(CCNs_range, "\n"), 
                    guide = guide_legend(direction = "horizontal", nrow = 1, 
                                         title.position = "left", label.position = "bottom",
                                         label.vjust = 0.5), 
                    drop = FALSE) +
  scale_y_discrete(limits = rev, breaks = colonyreps, position = "left") +
  labs(x = "Chromosome", y = NULL, title = "D & S") +
  #coord_fixed(ratio = 1) +
  CCN_heatmap_theme #+
  #geom_text(aes(label = StrainBackground_dCCN_label), color = "black", size = 3.5)
save_ggplot("copy_number_chr_heatmap", width = 3.5, height = 10)
save_ggplot("copy_number_chr_heatmap", width = 3.5, height = 10, mode = "paper")

# Calculate the karyotype difference between donut and spread
data_ds_diff <- data_ds %>%
  dplyr::select(all_of(c(id_vars, "Chr", "CCN"))) %>%
  dplyr::select(-Strain, -ColonyRep) %>%
  tidyr::spread(key = "ColonyMorph", value = "CCN") %>%
  dplyr::rename(D_CCN = D, S_CCN = S) %>%
  dplyr::mutate(DS_dCCN = S_CCN - D_CCN) %>%
  # Distinguish CCN change towards 4 ("T") and away from 4 ("A")
  # Distinguish CCN change to 4 i.e. restore 4 ("T4")
  # Remove 3N strains
  dplyr::mutate(DS_dCCN_mode = ifelse(D_CCN == S_CCN | StrainBackgroundRep %in% strainbackgroundreps_3N, # base ploidy is not 4N, thus not valid to say "change towards/from 4"
                                          NA,  
                                          ifelse(abs(S_CCN-4) < abs(D_CCN-4), "T", "A")
                ) %>% factor(., levels = DS_dCCN_modes), # CCN changes towards/from 4, note this does not apply to e.g. 5->3 but we do not have such cases (all change by 1, only two by 2 but 8->6)
                DS_dCCN_mode_label = ifelse(StrainBackgroundRep %in% strainbackgroundreps_3N, 
                                                "-", 
                                                ifelse(DS_dCCN_mode == "T", "T", NA)), 
                DS_dCCN_modeT = ifelse(DS_dCCN_mode == "T",
                                           ifelse(S_CCN == 4, "4", "n4"),  # CCN changes to (restored to) 4, otherwise non-4
                                           NA) %>% factor(., levels = DS_dCCN_modeTs), 
                DS_dCCN_3N_label = ifelse(StrainBackgroundRep %in% strainbackgroundreps_3N, 
                                            "-", 
                                            NA), 
                )

# Plot heatmap
DS_dCCNs_color <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(21)[c(11 - rev(c(2,4,6,8,9)), 11, 11 + c(2,4,6,8,9))]
DS_dCCNs_color[6] <- "gray95"
names(DS_dCCNs_color) <- seq(-5,5)
#show_col(DS_dCCNs_color, labels = TRUE)
DS_dCCNs_range <- min(data_ds_diff$DS_dCCN):max(data_ds_diff$DS_dCCN)
p2 <- ggplot(data_ds_diff %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, 
                                                 labels = gsub(" ", "\n", strainbackgrounds_label, fixed = TRUE)), 
                       Replicate2 = factor(Replicate, levels = reps, labels = reps_label), 
                       DS_dCCN = factor(DS_dCCN, levels = DS_dCCNs_range)),
       aes(x = Chr, y = Replicate2, fill = DS_dCCN)) +
  facet_grid(StrainBackground ~ Condition, 
             scales = "free_y", space = "free_y", switch = "y") +
  geom_tile(color = "black") +
  scale_fill_manual(name = "∆Copy\nnumber\n",
                    values = DS_dCCNs_color,
                    labels = paste0(DS_dCCNs_range, "\n"), 
                    guide = guide_legend(direction = "horizontal", nrow = 1, 
                                         title.position = "left", label.position = "bottom",
                                         label.vjust = 0.5), 
                    drop = FALSE) +
  scale_y_discrete(limits = rev, breaks = reps_label, position = "left") +
  labs(x = "Chromosome", y = NULL, title = "S - D") +
  #coord_fixed(ratio = 1) +
  CCN_heatmap_theme +
  geom_text(aes(label = DS_dCCN_3N_label), color = "black", size = 3.5)
  #geom_text(aes(label = DS_dCCN_mode_label, color = DS_dCCN_modeT), size = 3.5) +
  #scale_color_manual(values = DS_dCCN_modeTs_color) +
  #guides(color = "none")
save_ggplot("copy_number_chr_heatmap_diff", width = 3.5, height = 10)
save_ggplot("copy_number_chr_heatmap_diff", width = 3.5, height = 10, mode = "paper")

# Combine
p2_simp <- p2 + theme(strip.text.y = element_blank())
p3 <- plot_grid(p1, p2_simp, rel_widths = c(1,0.75))
save_ggplot("copy_number_chr_heatmap_combine", width = 3.5*1.75, height = 10)
save_ggplot("copy_number_chr_heatmap_combine", width = 3.5*1.75, height = 10, mode = "paper")

# Plot CCN changes in evolved macroscopic isolates and CCN changes in D-S pairs
# Evolved macroscopic isolates: show Gain/Loss as rectangle background
# D-S pairs: show Gain/Loss as positive/negative bars
# Color: Gain/Loss
data_macro_chr <- read.csv("data_macro_chr.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::arrange(factor(StrainBackground, levels = strainbackgrounds), Chr)
data_ds_diff_dchr <- data_ds_diff %>%
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::filter(DS_dCCN != 0) %>%
  dplyr::mutate(DS_dCCN_sign = ifelse(DS_dCCN > 0, "Gain", "Loss") %>% factor(levels = dCCN_signs), 
                DS_dCCN_mode2 = ifelse(DS_dCCN_mode == "T",  # only T or A in our dataset
                                       ifelse(S_CCN == 4, "T4", "Tn4"), # to 4, towards 4 but not to 4
                                       ifelse(D_CCN == 4, "A4", "An4")  # from 4, away from 4 but not from 4
                ) %>% factor(levels = DS_dCCN_mode2s)
  )
summary(data_ds_diff_dchr)  # 53 chr changes
data_ds_diff_dchr_summary_sb <- data_ds_diff_dchr %>%  # strain background
  dplyr::group_by_at(c("CondLine", "StrainBackground", "Condition", "Line", "EvoTime", 
                       "Chr", "DS_dCCN_sign")) %>%
  dplyr::summarize(Count = n())
data_ds_diff_dchr_summary_sb_blank <- data_ds %>%
  dplyr::select(StrainBackground, Condition, Replicate) %>%
  dplyr::distinct() %>%
  dplyr::group_by(StrainBackground, Condition) %>%
  dplyr::summarize(Count = n()) %>%
  dplyr::mutate(Chr = "I", 
                DS_dCCN_sign = "Gain")
p4 <- ggplot(data_ds_diff_dchr_summary_sb %>%
               dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, 
                                                       labels = gsub(" ", "\n", strainbackgrounds_label, fixed = TRUE)), 
                             DS_dCCN_sign = factor(DS_dCCN_sign, levels = dCCN_signs, labels = dCCN_signs_evods[3:4])), 
             aes(x = Chr, y = Count, fill = DS_dCCN_sign)) +
  facet_grid(StrainBackground ~ Condition, 
             scales = "free_y", space = "free_y", switch = "y") +
  geom_segment(data = data_macro_chr %>%
                 dplyr::inner_join(x = ., y = data_ds_diff_dchr_summary_sb %>% 
                                     dplyr::mutate(Chr = factor(Chr, levels = chrs, labels = 1:16) %>% as.integer()) %>% 
                                     dplyr::ungroup() %>%
                                     dplyr::select(StrainBackground, Chr), 
                                   by = c("StrainBackground", "Chr")) %>%
                 dplyr::left_join(x = ., y = data_ds_diff_dchr_summary_sb_blank %>% dplyr::select(StrainBackground, Count), 
                                  by = "StrainBackground") %>%
                 dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, 
                                                         labels = gsub(" ", "\n", strainbackgrounds_label, fixed = TRUE))), 
               mapping = aes(x = Chr, y = 0, xend = Chr, yend = Count + 0.5), inherit.aes = FALSE,
               linewidth = 0.5, linetype = "42", color = "grey") +
  geom_point(data = data_macro_chr %>%
               dplyr::left_join(x = ., y = data_ds_diff_dchr_summary_sb_blank %>% dplyr::select(StrainBackground, Count), 
                                by = "StrainBackground") %>%
               dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, 
                                                       labels = gsub(" ", "\n", strainbackgrounds_label, fixed = TRUE)), 
                             Macro_dCCN_sign = factor(Macro_dCCN_sign, levels = dCCN_signs, labels = dCCN_signs_evods[1:2])),
             mapping = aes(x = Chr, y = Count + 0.5, fill = Macro_dCCN_sign), inherit.aes = FALSE,
             shape = "diamond filled", stroke = 0, color = "white", size = 3.2) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) + ###
  scale_fill_manual(name = "Sign(∆Copy\nnumber)\n",
                    values = dCCN_signs_evods_color,
                    breaks = dCCN_signs_evods[c(3,4,1,2)],
                    #labels = dCCN_signs_evods,
                    guide = guide_legend(direction = "horizontal", nrow = 1, 
                                         title.position = "left", label.position = "bottom",
                                         label.vjust = 0.5), 
                    drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = 0:3, labels = 0:3, 
                     #expand = expansion(add = c(0,0.4)), 
                     position = "left") +
  labs(x = "Chromosome", y = "# of D-S pairs", title = "Summary") + 
  ggplot_custom_theme4 +
  CCN_heatmap_theme +
  theme(axis.line = element_blank(), 
        axis.ticks.y = element_line(color = "black", linewidth = 0.25)) + 
  geom_hline(yintercept = 0, linewidth = 0.25) + 
  geom_segment(data = data_ds_diff_dchr_summary_sb_blank %>%
                 dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, 
                                                         labels = gsub(" ", "\n", strainbackgrounds_label, fixed = TRUE))), 
               mapping = aes(x = 0, y = 0, xend = 0, yend = Count), linewidth = 0.25) +
  coord_cartesian(clip = "off")
save_ggplot("copy_number_chr_heatmap_summary", width = 3.5*1.1, height = 10)
save_ggplot("copy_number_chr_heatmap_summary", width = 3.5*1.1, height = 10, mode = "paper")

# Combine
p4_simp <- p4 + theme(strip.text.y = element_blank())
p5 <- plot_grid(p1, p2_simp, p4_simp, rel_widths = c(1,0.75, 0.88), nrow = 1)
save_ggplot("copy_number_chr_heatmap_combine_summary", width = 3.5*2.63, height = 10)
save_ggplot("copy_number_chr_heatmap_combine_summary", width = 3.5*2.63, height = 10, mode = "paper")

##### Draw stacked bar plots

# Exclude two 3N strains for all subsequent analyses
# Include PA4 for all subsequent analyses

# Plot histogram (stacked bar plot) of number of changed chromosomes per replicate
# Color by strain background
data_ds_diff_summary_rep <- data_ds_diff %>%
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::group_by_at(id_vars_diff) %>%
  dplyr::summarize(dChr_num = sum(DS_dCCN != 0))
summary(data_ds_diff_summary_rep)
dChr_nums <- unique(data_ds_diff_summary_rep$dChr_num) %>% sort()
ggplot(data_ds_diff_summary_rep %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, labels = strainbackgrounds_label), 
                       dChr_num = factor(dChr_num, levels = dChr_nums)), 
       aes(x = dChr_num, fill = StrainBackground)) + 
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = strainbackgrounds_color) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = seq(0,10,2), labels = seq(0,10,2), 
                     expand = expansion(mult = c(0,0.05))) +
  labs(x = "# of changed chromosomes", y = "# of D-S pairs", fill = NULL) +
  ggplot_custom_theme4 + 
  theme(panel.grid.major.y = element_line(size = 0.2, color = "gray80"), 
        legend.position = c(0.85,0.775), 
        legend.text = element_text(size = rel(0.75)), 
        legend.key.size = unit(12, "pt"))
save_ggplot("barplot_dChr_num", width = 4, height = 4)
save_ggplot("barplot_dChr_num", width = 4, height = 4, mode = "paper")

# Plot stacked bar plot of the number of times each chromosome changes and its direction
# Color by gain/loss
ggplot(data_ds_diff_dchr, aes(x = Chr, fill = DS_dCCN_sign)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = dCCN_signs_color, labels = paste0(dCCN_signs, " (DS)")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = seq(0,12,2), labels = seq(0,12,2), 
                     expand = expansion(mult = c(0,0.05))) +
  labs(x = "Chromosome", y = "# of D-S pairs", fill = NULL) + 
  ggplot_custom_theme4 + 
  theme(panel.grid.major.y = element_line(size = 0.2, color = "gray80"), 
        legend.position = c(0.135,0.9), 
        legend.text = element_text(size = rel(1)), 
        axis.text.x = element_text(size = rel(0.8)) #, angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("barplot_chr_sign", width = 6, height = 3)
save_ggplot("barplot_chr_sign", width = 6, height = 3, mode = "paper")
# Color by strain background
ggplot(data_ds_diff_dchr %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, labels = strainbackgrounds_label)), 
       aes(x = Chr, fill = StrainBackground)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = strainbackgrounds_color, 
                    guide = guide_legend(direction = "vertical", ncol = 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = seq(0,12,2), labels = seq(0,12,2), 
                     expand = expansion(mult = c(0,0.05))) +
  labs(x = "Chromosome", y = "# of D-S pairs", fill = NULL) + 
  ggplot_custom_theme4 + 
  theme(panel.grid.major.y = element_line(size = 0.2, color = "gray80"), 
        legend.position = c(0.205,0.825), 
        legend.text = element_text(size = rel(0.75)), 
        axis.text.x = element_text(size = rel(0.8)), 
        legend.key.size = unit(12, "pt"), 
        #legend.margin = margin(t = 1, b = 1, l = 1, r = 1)
  )
save_ggplot("barplot_chr_strainbkg", width = 6, height = 3)
save_ggplot("barplot_chr_strainbkg", width = 6, height = 3, mode = "paper")
# # Color by dCCN modes (towards or away from 4)
# ggplot(data_ds_diff_dchr %>%
#          dplyr::mutate(DS_dCCN_mode2 = factor(DS_dCCN_mode2, levels = DS_dCCN_mode2s, labels = DS_dCCN_mode2s_label)), 
#        aes(x = Chr, fill = DS_dCCN_mode2)) +
#   geom_bar(stat = "count", position = "stack") +
#   scale_fill_manual(values = DS_dCCN_mode2s_color) +
#   scale_x_discrete(drop = FALSE) +
#   scale_y_continuous(breaks = seq(0,12,2), labels = seq(0,12,2), 
#                      expand = expansion(mult = c(0,0.05))) +
#   labs(x = "Chromosome", y = "# of D-S pairs", fill = NULL) + 
#   ggplot_custom_theme4 + 
#   theme(panel.grid.major.y = element_line(size = 0.2, color = "gray80"), 
#         legend.position = c(0.15,0.8), 
#         legend.text = element_text(size = rel(1)), 
#         axis.text.x = element_text(size = rel(0.8)))
# save_ggplot("barplot_chr_mode", width = 6, height = 3)
# save_ggplot("barplot_chr_mode", width = 6, height = 3, mode = "paper")

# Plot stacked bar plot of the number of times each chromosome changes and its direction (macro vs micro isolates)
# Color by gain/loss
ggplot(data_macro_chr %>% 
         dplyr::mutate(Chr = factor(Chr, levels = 1:16, labels = chrs)), 
       aes(x = Chr, fill = Macro_dCCN_sign)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = dCCN_signs_color, labels = paste0(dCCN_signs, " (evo)")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0,9), breaks = seq(0,12,2), labels = seq(0,12,2), 
                     expand = expansion(mult = c(0,0.05))) +
  labs(x = "Chromosome", y = "# of macro isolates", fill = NULL) + 
  ggplot_custom_theme4 + 
  theme(panel.grid.major.y = element_line(size = 0.2, color = "gray80"), 
        legend.position = c(0.135,0.9), 
        legend.text = element_text(size = rel(1)), 
        axis.text.x = element_text(size = rel(0.8)) #, angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0))
  )
save_ggplot("barplot_macro_chr_sign", width = 6, height = 3)
save_ggplot("barplot_macro_chr_sign", width = 6, height = 3, mode = "paper")
# Color by strain background
ggplot(data_macro_chr %>%
         dplyr::mutate(Chr = factor(Chr, levels = 1:16, labels = chrs)) %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, labels = strainbackgrounds_label)), 
       aes(x = Chr, fill = StrainBackground)) +
  geom_bar(stat = "count", position = "stack") +
  scale_fill_manual(values = strainbackgrounds_color, 
                    guide = guide_legend(direction = "vertical", ncol = 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0,9), breaks = seq(0,12,2), labels = seq(0,12,2), 
                     expand = expansion(mult = c(0,0.05))) +
  labs(x = "Chromosome", y = "# of macro isolates", fill = NULL) + 
  ggplot_custom_theme4 + 
  theme(panel.grid.major.y = element_line(size = 0.2, color = "gray80"), 
        legend.position = c(0.205,0.825), 
        legend.text = element_text(size = rel(0.75)), 
        axis.text.x = element_text(size = rel(0.8)), 
        legend.key.size = unit(12, "pt"), 
        #legend.margin = margin(t = 1, b = 1, l = 1, r = 1)
  )
save_ggplot("barplot_macro_chr_strainbkg", width = 6, height = 3)
save_ggplot("barplot_macro_chr_strainbkg", width = 6, height = 3, mode = "paper")

##### PCA

# Combine karyotype data of donut-spread (data_ds) and evolved isolates (data_evo)
# Exclude two 3N strains for all subsequent analyses
data_ds_evo <- dplyr::bind_rows(
  data_ds %>%
    dplyr::select(-Ploidy, -StrainBackground_CCN, -StrainBackground_dCCN, -StrainBackground_dCCN_label, 
                  -CondLine, -ColonyRep),
  data_evo %>% 
    dplyr::select(-Ploidy) %>%
    dplyr::filter(Strain %in% strainbackgrounds) %>%
    dplyr::mutate(StrainBackground = Strain, 
                  Replicate = "rep0", 
                  StrainBackgroundRep = paste(StrainBackground, Replicate, sep = "_"), 
                  ColonyMorph = "D")) %>%
  dplyr::filter(!StrainBackgroundRep %in% strainbackgroundreps_3N) %>%
  dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds))
data_ds_evo$Strain <- factor(data_ds_evo$Strain, 
                             levels = as.character(unique(data_ds_evo$Strain)))
data_ds_evo$StrainBackgroundRep <- factor(data_ds_evo$StrainBackgroundRep, 
                                          levels = as.character(unique(data_ds_evo$StrainBackgroundRep)))
data_ds_evo$EvoTime <- factor(data_ds_evo$EvoTime, 
                              levels = as.character(unique(data_ds_evo$EvoTime)))
data_ds_evo$Replicate <- factor(data_ds_evo$Replicate, 
                                levels = as.character(unique(data_ds_evo$Replicate)))
data_ds_evo$ColonyMorph <- factor(data_ds_evo$ColonyMorph, 
                                  levels = as.character(unique(data_ds_evo$ColonyMorph)))
summary(data_ds_evo)

# Extract PCA's input data
data_pca <- data_ds_evo %>%
  dplyr::select(Strain, Chr, CCN) %>%
  tidyr::spread(key = "Chr", value = "CCN")
write.csv(data_pca[,-1], file = "pca_data.csv", row.names = FALSE)

# Pause running this R script
stop("Pause and now run PCA in pca.py Python script")
# Run PCA in Python script then come back here and proceed running
# Also copy here the % of variance explained by PC1/2

# Perform PCA in Python using sklearn
# Scaling the data is not necessary here as copy numbers of different chromosomes are on the same scale, and the 2D PCA plot looks more informative without scaling
pca_result <- read.csv("pca_result.csv", row.names = NULL)
pca_result$Strain <- as.character(data_pca$Strain)
# Combine with strain information
data_pca_result <- data_ds_evo %>%
  dplyr::select(-Chr, -CCN) %>%
  dplyr::distinct() %>%
  dplyr::left_join(pca_result, by = "Strain")
# % of variances explained by each PC
percentVar <- c(43.36, 22.78) %>% round()  # PC1 and PC2

# # Perform PCA in R
# pca_result <- prcomp(data_pca[,-1], center = TRUE, scale. = TRUE, rank. = 2)
# # % of variances explained by each PC
# summary(pca_result) 
# percentVar <- summary(pca_result)$importance[2,1:2]
# # Coordinates of each individual
# pca_result_df <- data.frame(pca_result$x[,1:2])
# pca_result_df$Strain <- as.character(unique(data_pca$Strain))
# # Combine with strain information
# data_pca_result <- data_ds_evo %>%
#   dplyr::select(-Chr, -CCN) %>%
#   dplyr::distinct() %>%
#   dplyr::left_join(pca_result_df, by = "Strain")

# Plot PCA
ggplot(data_pca_result %>%
         dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, labels = strainbackgrounds_label)), 
       aes(x = PC1, y = PC2, color = StrainBackground)) +
  # Draw thick lines between donuts of the same strain background
  geom_line(data = data_pca_result %>% 
              dplyr::filter(ColonyMorph == "D") %>%
              dplyr::mutate(StrainBackground = factor(StrainBackground, levels = strainbackgrounds, labels = strainbackgrounds_label)), 
            mapping = aes(group = StrainBackground),
            size = 1.5, alpha = 0.7) +
  # Draw thin lines between each donut-spread pair
  geom_line(mapping = aes(group = StrainBackgroundRep),
            size = 0.8, alpha = 0.3, color = "gray60") +
  # Draw data points of each donut and spread
  geom_point(mapping = aes(shape = ColonyMorph), 
             size = 3, alpha = 1, color = "white") + 
  geom_point(mapping = aes(shape = ColonyMorph), 
             size = 3, alpha = 0.7) + 
  scale_color_manual(values = strainbackgrounds_color, 
                     breaks = strainbackgrounds_label, labels = strainbackgrounds_label) +
  scale_shape_manual(values = colonymorphs_shape2) + 
  labs(#title = "Karyotype fitness landscape", 
       x = glue("PC1 ({percentVar[1]}% variance)"),
       y = glue("PC2 ({percentVar[2]}% variance)"), 
       color = NULL, shape = NULL) +
  ggplot_custom_theme4 +
  theme(legend.margin = margin(t = 0, b = 0), 
        title = element_text(size = rel(0.8)))
save_ggplot("pca", width = 6, height = 4)
save_ggplot("pca", width = 6, height = 4, mode = "paper")

##### Save data

save.image("analysis.RData")
for (df in c("data_ds", "data_evo", "data_ds_diff", "data_ds_diff_dchr", 
             "data_ds_diff_summary_rep", "data_ds_diff_dchr_summary_sb", 
             "data_ds_evo", "data_pca", "data_pca_result")) {
  write.csv(get(df), glue("{df}.csv"), row.names = FALSE)
  saveRDS(get(df), file = glue("{df}.rds"))
}

