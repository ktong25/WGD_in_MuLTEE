# Run after Fig2d.R

# Split by condition
for (condition in conditions) {

  ggplot(data_bin_plot %>% dplyr::filter(Condition == condition),  ###
         aes(x = Bin, y = Bin_copy_num, color = Chr, shape = Type)) +
    facet_grid(Sample~Chr, scales = "free", space = "free_x") +
    geom_point(size = 0.2) +
    scale_color_manual(values = chrs_color) +
    scale_shape_manual(values = c("Below" = "circle", "Above" = "triangle down open")) + 
    scale_y_continuous(breaks = seq(2,8,2), labels = seq(2,8,2)) +
    geom_hline(yintercept = seq(2,8,2), size = 0.25, color = "gray60") +
    geom_hline(yintercept = seq(1,7,2), size = 0.25, color = "gray85") +
    geom_hline(data = data_bin_plot_ploidy %>% dplyr::filter(Condition == condition),  ###
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
  save_ggplot(glue("copy_num_bin_{condition}"), width = 15, height = length(samples) / 2 * 1)  ###
  save_ggplot(glue("copy_num_bin_{condition}"), width = 15, height = length(samples) /2 * 1, mode = "paper")  ###

}

# Split by condition and group by line

for (condition in conditions) {
  
  ggplot(data_bin_plot %>% dplyr::filter(Condition == condition) %>%  ###
           dplyr::mutate(Line = as.character(Line)) %>%
           dplyr::mutate(Line = if_else(is.na(Line), "0", Line)) %>%
           dplyr::arrange(factor(Condition, levels = conditions), 
                          factor(Line, levels = c("0", lines)), 
                          factor(EvoTime, levels = evotimes)) %>%  ### this puts t0 in the front
           dplyr::mutate(Sample = factor(Sample, levels = unique(Sample))), 
         aes(x = Bin, y = Bin_copy_num, color = Chr, shape = Type)) +
    facet_grid(Sample~Chr, scales = "free", space = "free_x") +
    geom_point(size = 0.2) +
    scale_color_manual(values = chrs_color) +
    scale_shape_manual(values = c("Below" = "circle", "Above" = "triangle down open")) + 
    scale_y_continuous(breaks = seq(2,8,2), labels = seq(2,8,2)) +
    geom_hline(yintercept = seq(2,8,2), size = 0.25, color = "gray60") +
    geom_hline(yintercept = seq(1,7,2), size = 0.25, color = "gray85") +
    geom_hline(data = data_bin_plot_ploidy %>% dplyr::filter(Condition == condition),  ###
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
  save_ggplot(glue("copy_num_bin_{condition}_by_line"), width = 15, height = length(samples) / 2 * 1)  ###
  save_ggplot(glue("copy_num_bin_{condition}_by_line"), width = 15, height = length(samples) /2 * 1, mode = "paper")  ###
  
}
  

