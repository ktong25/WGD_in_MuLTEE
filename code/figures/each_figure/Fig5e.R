# Run after running Fig4cdef.R

### Correlate convergent karyotype changes with cell volume/AR changes
# Convergent karyotype changes: 
# - donut/spread points (modes in colors, cases in shapes, not distinguish donut/spread because the direction always goes down except for one that is not involved here)
# - connecting line (modes in colors, cases in linetypes) 
# The others show grey points and gray lines, as the bottom layer in the figure
# Parameters
dKar_modes <- c("Same background",  
                "Same line", 
                "Different lines", 
                "Other")
dKar_ids <- c("1", "2", "3", "0")
dKar_modes_color <- c("#FF0000", "#33B050", "#00B0F0", "gray90")  # red, green, blue
names(dKar_modes_color) <- dKar_modes
#show_col(dKar_modes_color, labels = TRUE)
dKar_ids_shape <- c("circle", "square", "triangle", "circle small")
names(dKar_ids_shape) <- dKar_ids
dKar_ids_linetype <- c("solid", "32", "12", "solid")
names(dKar_ids_linetype) <- dKar_ids
# Annotate cases
data_summary_dKar <- data_summary %>%
  dplyr::mutate(
    dKar_mode = case_when(
      StrainBackgroundRep %in% c("PA_2_t600_rep1", "PA_2_t600_rep2") ~ dKar_modes[1], 
      StrainBackgroundRep %in% c("PA_2_t1000_rep2", "PA_2_t1000_rep3") ~ dKar_modes[1], 
      StrainBackgroundRep %in% c("PA_5_t600_rep2", "PA_5_t600_rep3") ~ dKar_modes[1], 
      StrainBackgroundRep %in% c("PA_1_t600_rep2", "PA_1_t1000_rep2") ~ dKar_modes[2], 
      StrainBackgroundRep %in% c("PA_1_t1000_rep1", "PA_4_t1000_rep1") ~ dKar_modes[3], 
      StrainBackgroundRep %in% c("PA_2_t1000_rep1", "PA_3_t1000_rep2") ~ dKar_modes[3], 
      TRUE ~ dKar_modes[4]
    ) %>% factor(levels = dKar_modes), 
    dKar_id = case_when(
      StrainBackgroundRep %in% c("PA_2_t600_rep1", "PA_2_t600_rep2") ~ dKar_ids[1], 
      StrainBackgroundRep %in% c("PA_2_t1000_rep2", "PA_2_t1000_rep3") ~ dKar_ids[2], 
      StrainBackgroundRep %in% c("PA_5_t600_rep2", "PA_5_t600_rep3") ~ dKar_ids[3], 
      StrainBackgroundRep %in% c("PA_1_t600_rep2", "PA_1_t1000_rep2") ~ dKar_ids[1], 
      StrainBackgroundRep %in% c("PA_1_t1000_rep1", "PA_4_t1000_rep1") ~ dKar_ids[1], 
      StrainBackgroundRep %in% c("PA_2_t1000_rep1", "PA_3_t1000_rep2") ~ dKar_ids[2], 
      TRUE ~ dKar_ids[4]
    ) %>% factor(levels = dKar_ids)
  )
# Plot
ggplot(data_summary_dKar, aes(x = Mean_volume, y = Mean_AR)) +
  # Plot "Other" data points at the bottom
  geom_line(data = data_summary_dKar %>% dplyr::filter(dKar_mode == "Other"), 
            mapping = aes(group = StrainBackgroundRep, color = dKar_mode, linetype = dKar_id),
            size = 0.6, alpha = 1) +
  geom_point(data = data_summary_dKar %>% dplyr::filter(dKar_mode == "Other"),
             mapping = aes(color = dKar_mode, shape = dKar_id), 
             size = 2, alpha = 1) +
  # Then plot the relevant data points
  geom_line(data = data_summary_dKar %>% dplyr::filter(dKar_mode != "Other"), 
            mapping = aes(group = StrainBackgroundRep, color = dKar_mode, linetype = dKar_id),
            size = 0.6, alpha = 1) +
  geom_point(data = data_summary_dKar %>% dplyr::filter(dKar_mode != "Other"),
             mapping = aes(color = dKar_mode, shape = dKar_id), 
             size = 2, alpha = 1) +
  scale_color_manual(values = dKar_modes_color, breaks = dKar_modes[1:3], labels = dKar_modes[1:3]) +
  scale_shape_manual(values = dKar_ids_shape, breaks = dKar_ids[1:3], labels = dKar_ids[1:3]) + 
  scale_linetype_manual(values = dKar_ids_linetype, breaks = dKar_ids[1:3], labels = dKar_ids[1:3]) + 
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5) +
  labs(x = expression(Cell~volume~(mu*m^3)), y = "Cell aspect ratio",
       color = "Convergent\nkaryotype changes", shape = "Case", linetype = "Case") +
  guides(color = guide_legend(order = 1)) +
  ggplot_custom_theme4 +
  theme(legend.position = "right", 
        legend.key.width = unit(30, "pt"))
save_ggplot("cell_volume_vs_AR_dKar", width = 6, height = 4)
save_ggplot("cell_volume_vs_AR_dKar", width = 6, height = 4, mode = "paper")
# Save data
write.csv(data_summary_dKar, file = "data_summary_dKar.csv", row.names = FALSE)
saveRDS(data_summary_dKar, file = "data_summary_dKar.rds")
