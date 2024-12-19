library(plyr)  # load before dplyr
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)
library(RColorBrewer)
library(scales)
library(mixtools)
library(effsize)

setwd("~/Documents/projects/R_data_analysis/ploidy_measurement/20241021_MA_Sayantan")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup2.R")

# Parameters
in_dir <- "data_final_csvs"
out_fig_path <- "."
control_strains <- c("2N", "4N")  # first 2N then 4N
conditions <- c("PM", "PA")
lines <- as.character(1:5)
evotimes <- c("t0", "t1000", "4N")
evotimes_label <- c("t0 (2N)", "t1000 (4N)", "Engineered 4N")
replicates <- c("anc", "A", "B", "C")
line_evotimes <- c("t0", lines, "4N")
line_evotimes_label <- c("t0", lines, "Engineered 4N")
ploidy_changes <- c("up", "small change", "down")
ploidy_changes_color <- c("#F8766D", "gray", "#00BFC4") #c("#1B9E77", "lightgreen", "#D95F02")
names(ploidy_changes_color) <- ploidy_changes
#show_col(ploidy_changes_color)
ploidy_changes_color2 <- c("#C65E57", "gray10", "#00999D") # c("#F8766D", "gray10", "#00BFC4")  #c("#1B9E77", "lightgreen", "#D95F02")
names(ploidy_changes_color2) <- ploidy_changes
#show_col(ploidy_changes_color2)

# Load data
# PA_t0_rep2 is abnormal thus absent here
# Input file name: e.g. "<batch>.csv"
load_one_csv <- function(file) {
  df <- read.csv(file, row.names = NULL, stringsAsFactors = FALSE) %>%
    dplyr::filter(!Strain %in% control_strains) %>%
    dplyr::select(Strain, RawIntDen_bgs, DNA_content) %>%
    dplyr::mutate(Batch = substr(basename(file), 1, stri_length(basename(file))-stri_length(".csv"))) %>%
    dplyr::select(Batch, Strain, RawIntDen_bgs, DNA_content)
  return(df)
}
data <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
              .fun = load_one_csv)
# Filter data
repeats_strains <- data %>%
  dplyr::filter(Batch == "Repeats") %>%
  dplyr::pull(Strain) %>%
  unique()
repeats_strains_notused <- c("PA_5_t1000_D0", "PA_5_t1000_rep1", "PM_t0_rep2")
repeats_strains_used <- repeats_strains[!repeats_strains %in% repeats_strains_notused]
data <- data %>%
  dplyr::filter(!(Batch != "Repeats" & Strain %in% repeats_strains_used)) %>%
  dplyr::filter(!(Batch == "Repeats" & Strain %in% repeats_strains_notused)) %>%
  dplyr::filter(!(Batch != "Repeats_2" & Strain == "PM_t0_rep2"))
length(unique(data$Strain)) 
length(unique(paste(data$Batch, data$Strain)))
# Replace with PA5 t1000 data from earlier dataset
df <- read.csv("~/Documents/projects/R_data_analysis/ploidy_measurement/20230606_PMPA_evo_combine/data_ploidy.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::filter(Strain == "PA_5_t1000") %>%
  dplyr::select(Strain, RawIntDen_bgs, DNA_content) %>%
  dplyr::mutate(Batch = "20230606_PMPA_evo_combine") %>%
  dplyr::select(Batch, Strain, RawIntDen_bgs, DNA_content)
unique(df$Strain)
data <- rbind(data %>% dplyr::filter(Strain != "PA_5_t1000_D0"), 
              df %>% dplyr::mutate(Strain = "PA_5_t1000_D0"))
rm(df)
length(unique(data$Batch)) # 7
length(unique(data$Strain)) # 55
length(unique(paste(data$Batch, data$Strain))) # 55

# Format data
data <- data %>%
  dplyr::mutate(Strain = gsub("t0", "NA_t0", Strain, fixed = TRUE) %>%
                  gsub("D0", "anc", ., fixed = TRUE) %>%
                  gsub("4N_M", "PM_NA_4N", ., fixed = TRUE) %>%
                  gsub("4N_A", "PA_NA_4N", ., fixed = TRUE) %>%
                  gsub("rep1", "A", ., fixed = TRUE) %>%
                  gsub("rep2", "B", ., fixed = TRUE) %>%
                  gsub("rep3", "C", ., fixed = TRUE)
                ) %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line", "EvoTime", "Replicate"), sep = "_", 
                  remove = FALSE, convert = FALSE) %>%
  tidyr::unite(col = "Background", Condition, Line, EvoTime, sep = "_", remove = FALSE) %>%
  dplyr::relocate(Background, .before = "Condition") %>%
  dplyr::mutate(Line_EvoTime = ifelse(EvoTime == "t1000", as.character(Line), as.character(EvoTime))) %>%
  dplyr::relocate(Line_EvoTime, .after = "EvoTime") %>%
  dplyr::mutate(Condition = factor(Condition, levels = conditions), 
                Line = factor(Line, levels = lines), 
                EvoTime = factor(EvoTime, levels = evotimes), 
                Line_EvoTime = factor(Line_EvoTime, levels = line_evotimes),  
                Replicate = factor(Replicate, levels = replicates)) %>%
  dplyr::arrange(factor(Condition, levels = conditions),
                 factor(Line_EvoTime, levels = line_evotimes),
                 factor(Replicate, levels = replicates))
length(unique(data$Background)) # 14
length(unique(data$Line_EvoTime)) # 7
backgrounds <- unique(data$Background)
id_vars <- colnames(data)[1:(which(colnames(data) == "RawIntDen_bgs")-1)]

# Get G1 peak for each strain (G1 peak is usually the highest peak)
data_summary <- data %>%
  dplyr::group_by_at(id_vars) %>%
  dplyr::summarise(Count = n(),  
                   G1Peak = density(DNA_content)$x[which.max(density(DNA_content)$y)]) %>%
  dplyr::arrange(factor(Condition, levels = conditions),
                 factor(Line_EvoTime, levels = line_evotimes),
                 factor(Replicate, levels = replicates)) %>%
  dplyr::ungroup()
# Adjust for some strains by finding all peaks and choosing the right peak
# PM_t0_rep2
x <- data[data$Strain == "PM_NA_t0_B", "DNA_content"]
d <- density(x)
d$x[c(F, diff(diff(d$y)>=0)<0)]  # 2.00288598
data_summary[data_summary$Strain == "PM_NA_t0_B", "G1Peak"] <- 2.00288598
# Subtract by ancestor ploidy
data_summary <- dplyr::left_join(x = data_summary, 
                                 y = data_summary %>% 
                                   dplyr::filter(Replicate == "anc") %>%
                                   dplyr::select(Background, G1Peak) %>%
                                   dplyr::rename(G1Peak_anc = G1Peak), 
                                 by = "Background") %>%
  dplyr::mutate(G1Peak_diff = G1Peak - G1Peak_anc)

### Determine ploidy change

# Function to fit GMM and extract posterior probabilities
fit_and_extract_G1 <- function(data) {
  # Fit GMM
  set.seed(42)
  gmm_fit <- normalmixEM(data, k = 2) # set max ploidy  # if needed, set initial guesses for means and weights
  n <- ifelse(gmm_fit$mu[1] < gmm_fit$mu[2], 1, 2)  # component ID for G1
  # Extract and return parameters and data of G1 peak
  return(list(
    G1_mean = gmm_fit$mu[n], 
    G1_sd = gmm_fit$sigma[n], 
    G1_weight = gmm_fit$lambda[n],
    G1_data = data[gmm_fit$posterior[,n] > gmm_fit$posterior[,3-n]], # when n is 1, 3-n is 2; when n is 2, 3-n is 1
    G2_data = data[gmm_fit$posterior[,n] < gmm_fit$posterior[,3-n]]
  ))
}

# Iterate across each background then each replicate
data_summary_ploidy_diff <- map_dfr(
  .x = backgrounds, 
  .f = function(background) {
    replicates <- data %>% dplyr::filter(Background == background) %>% dplyr::pull(Replicate) %>% as.character() %>% unique()
    min_ploidy <- ifelse(endsWith(background, "t0"), 1, 2)
    max_ploidy <- ifelse(endsWith(background, "t0"), 5.5, 
                         ifelse(startsWith(background, "PM") | endsWith(background, "4N"), 10, 11.5))
    data2 <- data %>% dplyr::filter(Background == background & Replicate == "anc") %>% 
      dplyr::filter(DNA_content > min_ploidy & DNA_content < max_ploidy) %>% dplyr::pull(DNA_content)
    map_dfr(
      .x = replicates, 
      .f = function(replicate) {
        min_ploidy <- ifelse(background == "PA_5_t1000" & replicate %in% c("B", "C"), 1.5, min_ploidy) ###
        data1 <- data %>% dplyr::filter(Background == background & Replicate == replicate) %>% 
          dplyr::filter(DNA_content > min_ploidy & DNA_content < max_ploidy) %>% dplyr::pull(DNA_content)
        # Fit GMM to both datasets
        fit1 <- fit_and_extract_G1(data1)
        fit2 <- fit_and_extract_G1(data2)
        # Perform t-test 
        t_test_result <- t.test(fit1$G1_data, fit2$G1_data) #format.pval(t_test_result$p.value, digits = 3)
        # Calculate effect size (Cohen's d)
        effect_size <- cohen.d(fit1$G1_data, fit2$G1_data) #print(effect_size)
        # Report
        data.frame(Background = background, 
                   Replicate = replicate, 
                   G1_mean = fit1$G1_mean,
                   G1_mean_anc = fit2$G1_mean,
                   G1_sd = fit1$G1_sd, 
                   G1_sd_anc = fit2$G1_sd,
                   t_pval = t_test_result$p.value, 
                   cohend_d = effect_size$estimate, 
                   cohend_mag = effect_size$magnitude
                   )
    })
  })
data_summary <- dplyr::left_join(x = data_summary, y = data_summary_ploidy_diff, by = c("Background", "Replicate"))
data_summary <- data_summary %>%
  dplyr::mutate(Ploidy_change = ifelse(Replicate == "anc", NA, 
                                       ifelse(abs(cohend_d) < 0.8, "small change", 
                                              ifelse(G1_mean > G1_mean_anc, "up", "down")) ) %>% factor(levels = ploidy_changes))
# data_summary_ploidy_diff %>%
#   dplyr::filter(cohend_mag %in% c("medium", "large")) %>%
#   View()

# Add ploidy change levels to data
data2 <- data
data <- dplyr::left_join(x = data, 
                         y = data_summary %>% dplyr::select(Strain, Ploidy_change), 
                         by = "Strain")

# Plot distributions with ploidy change level assigned
ggplot(data %>% dplyr::mutate(Line_EvoTime = factor(Line_EvoTime, levels = line_evotimes, labels = line_evotimes_label)), 
       aes(x = DNA_content, y = Replicate, fill = Ploidy_change)) +
  facet_grid(Line_EvoTime ~ Condition) +  #labeller = labeller(Line = lines_label)
  geom_density_ridges(scale = 2, size = lwpt * 0.4, alpha = 0.2) +
  scale_x_continuous(limits = c(0, 11.5), 
                     breaks = seq(0, 11, 2),
                     #minor_breaks = seq(0, 11, 1), 
                     labels = c("0", paste0(seq(2, 11, 2), "N")), 
                     expand = expansion(mult = c(0,0.05))) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = ploidy_changes_color2) +
  labs(x = "DNA content", y = "Replicate population", fill = "Ploidy") +
  ggplot_custom_theme2 +
  theme(legend.position = "bottom", 
        legend.key.size = unit(8, "points"),
        panel.grid.major.x = element_line(linewidth = lwpt * 0.25, colour = "gray80"),
        #panel.grid.minor.x = element_line(linewidth = lwpt * 0.25, colour = "gray80"),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_text(vjust = 0, margin = margin(r = 2)), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(), 
        axis.title.y = element_text(margin = margin(r = 5)),
  )
save_ggplot("MA", width = 100, height = 120, units = "mm")
save_ggplot("MA", width = 100, height = 120, units = "mm", mode = "paper")

# Plot ploidy change with levels assigned
set.seed(42)
ggplot(data_summary %>% dplyr::filter(Replicate != "anc") %>%
         dplyr::mutate(EvoTime = factor(EvoTime, levels = evotimes, labels = evotimes_label)), 
       aes(x = Condition, y = G1Peak_diff, color = Ploidy_change)) +
  facet_wrap(~EvoTime, nrow = 1) +
  #geom_boxplot(width = 0.5, size = lwpt * 0.5, fill = NA, outlier.shape = NA, color = "skyblue3") +
  geom_hline(yintercept = 0, color = "gray", linewidth = lwpt * 0.5, linetype = "dashed") +
  geom_jitter(width = 0.2, alpha = 0.5, size = pspt * 2) +
  # stat_summary(geom = "errorbar", fun = mean, mapping = aes(ymin = after_stat(y), ymax = after_stat(y)),
  #              width = 0.4, size = lwpt * 1, color = "skyblue2") +
  # stat_summary(geom = "errorbar", fun.data = "mean_sdl", fun.args = list(mult = 1),
  #              width = 0.2, size = lwpt * 0.75, color = "skyblue2") +
  scale_color_manual(values = ploidy_changes_color) +
  scale_y_continuous(limits = c(-1,0.5), breaks = c(-1, -0.5, 0, 0.5), labels = c("-1N", "-0.5N", 0, "0.5N")) +
  labs(x = NULL, y = "Ploidy change", color = "Ploidy") +  ###
  ggplot_custom_theme2 +
  theme(legend.position = "bottom", 
        legend.key.size = unit(5, "points"))
save_ggplot("MA_diff", width = 75, height = 50, units = "mm")
set.seed(42)
save_ggplot("MA_diff", width = 75, height = 50, units = "mm", mode = "paper")

# Save data
for (df in c("data", "data_summary")) {
  write.csv(get(df), glue("{df}.csv"), row.names = FALSE)
  saveRDS(get(df), file = glue("{df}.rds"))
}

