ggplot_custom_theme <- theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = rel(2), hjust = 0.5), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)), 
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.25)), 
        strip.text = element_text(size = rel(1.5)))

ggplot_custom_theme2 <- theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = rel(2), hjust = 0.5), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)), 
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.25)), 
        strip.text = element_text(size = rel(1.5)), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

ggplot_custom_theme3 <- theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = rel(2), hjust = 0.5), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)), 
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.25)), 
        strip.text = element_text(size = rel(1.5)), 
        strip.background = element_blank())

ggplot_custom_theme4 <- theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = rel(2), hjust = 0.5), 
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)), 
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.25)), 
        strip.text = element_text(size = rel(1.5)), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(), 
        strip.background = element_blank())

get_pal_colors <- function(pal) {
  brewer.pal(n = brewer.pal.info[pal, "maxcolors"], name = pal)
}

save_ggplot <- function(filename, 
                        plot = last_plot(), 
                        path = get("out_fig_path", 
                                   envir = knitr::knit_global()), 
                        format = "png",  # "png", "pdf", "rds"
                        mode = "slide",  # "slide", "paper" 
                        ...) {
  if (mode == "paper") {
    path <- file.path(path, "paper")
    if (!file.exists(path)) {dir.create(path)}
    format <- c("pdf", "rds")
  }
  if ("rds" %in% format) {
    saveRDS(plot, file = file.path(path, paste0(filename, ".rds")))
  }
  for (device in intersect(c("png", "pdf"), format)) {
    ggsave(paste(filename, device, sep = "."), plot = plot, path = path, ...)
  }
}

# save_ggplot <- function(filename, 
#                         plot = last_plot(), 
#                         path = get("out_fig_path", 
#                                    envir = knitr::knit_global()), 
#                         format = c("png", "pdf"),  # "png", "pdf", "rds"
#                         ...) {
#   if ("rds" %in% format) {
#     saveRDS(plot, file = file.path(path, paste0(filename, ".rds")))
#   }
#   for (device in intersect(c("png", "pdf"), format)) {
#     ggsave(paste(filename, device, sep = "."), plot = plot, path = path, ...)
#   }
# }