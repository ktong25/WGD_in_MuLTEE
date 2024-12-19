library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(rstatix)

setwd("~/Documents/projects/R_data_analysis/competition_assay/20220909_2N4N_corrected_remove_PO_polish_20230928")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
conditions_label <- c("Mixotrophic", "Anaerobic")
selections <- c("S", "NS")
selections_label <- c("w/", "w/o")
selections_label_color <- c("#B3197F", "#5B83BA")  # purple, blue
names(selections_label_color) <- selections_label
#show_col(selections_label_color, labels = TRUE)

# Load data
dfit <- read.csv("comp_dfit.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::mutate(Condition = factor(Condition, levels = conditions, labels = conditions_label), 
                Selection = factor(Selection, levels = selections, labels = selections_label))
summary(dfit)

# Perform statistical tests
dfit_stat <- dfit %>% 
  dplyr::group_by(Condition) %>%
  rstatix::t_test(SR ~ Selection, paired = TRUE, detailed = TRUE) %>%  # default var.equal = FALSE (Welch's t-test)
  #adjust_pvalue(method = "bonferroni") %>%
  #add_significance("p.adj") %>%
  add_xy_position(x = "Condition", fun = "max", dodge = 0.8)

# Plot selection rate
ggplot(dfit, aes(x = Condition, y = SR)) +
  stat_summary(mapping = aes(fill = Selection), 
               geom = "bar", fun = mean, color = NA, alpha = 0.3, width = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_point(mapping = aes(color = Selection), 
             alpha = 0.6, size = 2.5, 
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8, seed = 2)) +
  stat_summary(mapping = aes(color = Selection), 
               geom = "errorbar", fun.data = mean_se, width = 0.4, linewidth = 0.75,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = selections_label_color) +
  scale_color_manual(values = selections_label_color) +
  stat_pvalue_manual(dfit_stat, label = "p = {p}",  # Selection should not be in ggplot(aes()) but in separate geom_XXX(aes())
                     bracket.size = 0.5, tip.length = 0.03, label.size = 5, vjust = -0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0.05,0.15)), n.breaks = 5) +
  labs(x = NULL, y = "Selection rate (per day)", title = "4N vs. 2N", 
       fill = "Settling selection", color = "Settling selection") +
  ggplot_custom_theme4 +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(1.25)))
save_ggplot("selection_rate", width = 4, height = 5)
save_ggplot("selection_rate", width = 4, height = 5, mode = "paper")

# Perform statistical tests (one-sample t-test vs. 0)
dfit_stat_vs_0_alt_greater <- dfit %>% 
  dplyr::group_by(Condition, Selection) %>%
  rstatix::t_test(SR ~ 1, mu = 0, alternative = "greater", detailed = TRUE)
dfit_stat_vs_0_alt_less <- dfit %>% 
  dplyr::group_by(Condition, Selection) %>%
  rstatix::t_test(SR ~ 1, mu = 0, alternative = "less", detailed = TRUE)

# Save data
write.csv(dfit, file = "dfit.csv", row.names = FALSE)
write.csv(dfit_stat %>% apply(2, as.character), file = "dfit_stat.csv", row.names = FALSE)
write.csv(dfit_stat_vs_0_alt_greater %>% apply(2, as.character), file = "dfit_stat_vs_0_alt_greater.csv", row.names = FALSE)
write.csv(dfit_stat_vs_0_alt_less %>% apply(2, as.character), file = "dfit_stat_vs_0_alt_less.csv", row.names = FALSE)
saveRDS(dfit, file = "dfit.rds")
saveRDS(dfit_stat, file = "dfit_stat.rds")
saveRDS(dfit_stat_vs_0_alt_greater, file = "dfit_stat_vs_0_alt_greater.rds")
saveRDS(dfit_stat_vs_0_alt_less, file = "dfit_stat_vs_0_alt_less.rds")
