library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

setwd("~/Documents/projects/R_data_analysis/chromosome_cnv/20231214_summary_MuLTEE")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
out_fig_path <- "."
conditions <- c("PM", "PA")
lines <- 1:5 %>% as.character()
lines_color <- get_pal_colors("Paired")[c(6,2,4,10,8)]  # consistent with Bozdag2023
names(lines_color) <- lines
#show_col(lines_color, labels = TRUE)
evotimes <- c("t0", "t200", "t400", "t600", "t1000")
evotimes2 <- gsub("t", "", evotimes, fixed = TRUE) %>% as.numeric()

# Load data
data_chr <- readRDS("data_chr.rds")
summary(data_chr)

# For each strain, calculate the 
data_chr$Chr_copy_num_diff <- data_chr$Chr_copy_num_round_correct - round(data_chr$Ploidy)
data_summary <- data_chr %>% 
  dplyr::group_by(Condition, Line, EvoTime) %>%
  dplyr::summarize(Num_chr_diff = sum(Chr_copy_num_diff != 0) %>% as.integer(), # number of chromosomes with copy number deviating from 4
                   Sum_chr_diff = sum(abs(Chr_copy_num_diff)) %>% as.integer(), # sum of abs(chromosome copy number - baseline copy number)
                   CV_chr = sd(Chr_copy_num_round_correct) / mean(Chr_copy_num_round_correct), # Coefficient of Variation
                   )  

# Plot
data_summary_t0 <- dplyr::filter(data_summary, EvoTime == "t0")
data_summary_plot <- dplyr::filter(data_summary, EvoTime != "t0")
for (line in c(lines, NA)) {  # final NA is for plotting gray color at the top of overlapping points at t0
  data_summary_plot <- rbind(data_summary_plot, 
                             dplyr::mutate(data_summary_t0, Line = line))
}
ggplot(data_summary_plot %>% 
         dplyr::filter(Condition == "PA") %>%
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()), 
       aes(x = EvoTime2, y = Num_chr_diff, color = Line)) +
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0,8,2), expand = expansion(mult = c(0.01, 0.05))) +  ###
  labs(title = "PA", x = "Days of evolution", y = "Number of chr with\ndeviated copy number", color = "Line") + 
  ggplot_custom_theme4
save_ggplot("num_chr_diff_evo", width = 6, height = 5)
save_ggplot("num_chr_diff_evo", width = 6, height = 5, mode = "paper")
ggplot(data_summary_plot %>% 
         dplyr::filter(Condition == "PA") %>%
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()), 
       aes(x = EvoTime2, y = Sum_chr_diff, color = Line)) +   ###
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))) +
  labs(title = "PA", x = "Days of evolution", y = "Sum of chr copy number deviation", color = "Line") +   ###
  ggplot_custom_theme4
save_ggplot("sum_chr_diff_evo", width = 6, height = 5)   ###
save_ggplot("sum_chr_diff_evo", width = 6, height = 5, mode = "paper")   ###
ggplot(data_summary_plot %>% 
         dplyr::filter(Condition == "PA") %>%
         dplyr::mutate(EvoTime2 = gsub("t", "", EvoTime, fixed = TRUE) %>% as.numeric()), 
       aes(x = EvoTime2, y = CV_chr, color = Line)) +   ###
  geom_line(size = 0.75, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = lines_color) +
  scale_x_continuous(breaks = evotimes2, expand = expansion(mult = c(0.05, 0.1))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))) +
  labs(title = "PA", x = "Days of evolution", y = "CV of chr copy number", color = "Line") +   ###
  ggplot_custom_theme4
save_ggplot("cv_chr_evo", width = 6, height = 5)   ###
save_ggplot("cv_chr_evo", width = 6, height = 5, mode = "paper")   ###

# Linear regression
data_summary$EvoTime2 <- factor(data_summary$EvoTime, 
                                levels = evotimes, 
                                labels = evotimes2) %>% as.numeric()
result <- lm(Num_chr_diff ~ EvoTime2, data = data_summary)
summary(result) %>% print()
# Residual standard error: 2.541 on 40 degrees of freedom
# Multiple R-squared:  0.1231,	Adjusted R-squared:  0.1012 
# F-statistic: 5.618 on 1 and 40 DF,  p-value: 0.02269
result <- lm(Sum_chr_diff ~ EvoTime2, data = data_summary)
summary(result) %>% print()
# Residual standard error: 3.588 on 40 degrees of freedom
# Multiple R-squared:  0.1335,	Adjusted R-squared:  0.1119 
# F-statistic: 6.164 on 1 and 40 DF,  p-value: 0.01733
result <- lm(CV_chr ~ EvoTime2, data = data_summary)
summary(result) %>% print()
# Residual standard error: 0.07993 on 40 degrees of freedom
# Multiple R-squared:  0.1284,	Adjusted R-squared:  0.1066 
# F-statistic: 5.894 on 1 and 40 DF,  p-value: 0.01979

# Save data
write.csv(data_summary, file = "data_chr_summary.csv", row.names = FALSE)