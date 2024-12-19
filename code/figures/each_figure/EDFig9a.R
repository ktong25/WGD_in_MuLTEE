rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")
load("/Users/kaitong/Documents/projects/R_data_analysis/mutations/20231230_donutspread/analysis.RData")

lwpt <- 0.352/0.75  # line width: size 1 = 0.75mm, print 1pt = 0.352mm
adjust_theme <- function(p) {
  p +
    theme(text = element_text(family = "Helvetica"), # Arial requires additional code to work in Mac
          plot.title = element_text(size = 7, hjust = 0.5), 
          axis.title = element_text(size = 7), 
          axis.text = element_text(size = 6), 
          axis.text.x = element_text(color = "black"), 
          axis.text.y = element_text(color = "black"), 
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6), 
          strip.text = element_text(size = 7),
          axis.line = element_line(linewidth = lwpt * 0.5), 
          axis.ticks = element_line(linewidth = lwpt * 0.5), 
          axis.ticks.length = unit(0.75, units = "mm"), 
          legend.spacing = unit(0, 'points'),
          legend.box.spacing = unit(0, units = "points"),
    )
}

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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        legend.key.size = unit(8, "points"),
  ) %>% adjust_theme()
#save_ggplot("barplot_mut-mode_allele-copy-num-change", width = 60, height = 60, units = "mm")
save_ggplot("barplot_mut-mode_allele-copy-num-change", width = 60, height = 60, units = "mm", mode = "paper")

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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        legend.key.size = unit(8, "points"),
  ) %>% adjust_theme() + 
  theme(legend.position = "none") ###
save_ggplot("barplot_mut-mode_allele-copy-num-change_nolegend", width = 50, height = 60, units = "mm", mode = "paper")

