library(dplyr)

setwd("~/Documents/projects/R_data_analysis/papers/WGD")
rm(list=ls())

# Parameters
# Input
in_folder_global <- "~/Documents/projects/R_data_analysis"
read_data_rds <- function(rds_path) {
  readRDS(file.path(in_folder_global, rds_path)) 
}
read_data_csv <- function(csv_path) {
  read.csv(file.path(in_folder_global, csv_path), row.names = NULL, stringsAsFactors = FALSE)
}
# Output
out_source_data_path <- "data/source_data"
out_raw_data_path <- "data/raw_data"
if (!file.exists(out_source_data_path)) {dir.create(out_source_data_path, recursive = TRUE)}
if (!file.exists(out_raw_data_path)) {dir.create(out_raw_data_path, recursive = TRUE)}
write_source_data <- function(filename, df) {
  write.csv(df, file = file.path(out_source_data_path, paste0(filename, ".csv")), row.names = FALSE)
}
write_raw_data <- function(filename, df) {
  write.csv(df, file = file.path(out_raw_data_path, paste0(filename, ".csv")), row.names = FALSE)
}
sup_fig_prefix <- "ED"

# Write data
write_source_data(
  "Fig1cde_evo", 
  read_data_csv("cluster_cell_correlation/20230926/data_evo.csv") %>%
    dplyr::select(Condition, Line, EvoTime, 
                  Cluster.Weighted_mean_radius, Cell.Mean_volume, Cell.Mean_AR)
)
write_source_data(
  "Fig1cde_ploidy_Fig3f", 
  read_data_csv("cluster_cell_correlation/20230926/data_artificial.csv") %>%
    dplyr::select(Condition, Ploidy,  
                  Cluster.Weighted_mean_radius, Cell.Mean_volume, Cell.Mean_AR)
)
write_source_data(
  paste0("Fig1c_raw_", sup_fig_prefix, "Fig2ab"), 
  read_data_csv("cluster_size/20230926_all/data.csv")
)
write_source_data(
  paste0("Fig1de_raw_", sup_fig_prefix, "Fig2cd"), 
  read_data_csv("cell_size_aspect_ratio/20230926_all/data.csv") %>%
    dplyr::select(-ImageID, -Circ., -Solidity)
)
write_source_data(
  "Fig1fgh", 
  read_data_csv("snowflake_physics_modeling/20231211/data_summary.csv") %>%
    dplyr::select(-Count, -ends_with("_se"))
)
write_raw_data(
  "Fig1fgh", 
  read_data_csv("snowflake_physics_modeling/20231211/data.csv") %>%
    dplyr::select(-Cluster_volume) %>%
    dplyr::relocate(Cluster_radius, .before = Cluster_ncells)
)
write_source_data(
  "Fig2a", 
  read_data_csv("allele_frequency/20231223_MuLTEE_summary/data.csv")
)
write_source_data(
  paste0(sup_fig_prefix, "Fig3b"), 
  read_data_csv("ploidy_measurement/20220923_trial/data_controls.csv") %>%
    dplyr::select(Strain, RawIntDen, cluster_ID, cluster_cyto_BackInt, RawIntDen_bgs) %>%
    dplyr::mutate(Strain = case_when(Strain == "2Nm" ~ "2N", 
                                     Strain == "4Nm" ~ "4N"), 
                  Area_px = (RawIntDen - RawIntDen_bgs) / cluster_cyto_BackInt) %>% 
    dplyr::relocate(Area_px, .before = cluster_ID)
)
write_source_data(
  "Fig2b", 
  read_data_csv("ploidy_measurement/20230606_PMPA_evo_combine/data_ploidy_summary.csv") %>%
    dplyr::select(-Batch, -Strain, -Count)
)
write_raw_data(
  "Fig2b", 
  read_data_csv("ploidy_measurement/20230606_PMPA_evo_combine/data_ploidy.csv") %>%
    dplyr::select(-Batch, -Strain, -RawIntDen_bgs)
)
write_source_data(
  "Fig2c", 
  read_data_csv("ploidy_measurement/20230606_PMPA_early_combine/data.csv") %>%
    dplyr::select(-Batch, -Strain, -RawIntDen_bgs, -Ploidy)
)
write_source_data(
  "Fig2d", 
  read_data_csv("chromosome_cnv/20231214_summary_MuLTEE/data_chr.csv") %>%
    dplyr::select(-Sample_ID, -Sample, -Ploidy, -Chr_copy_num, -Chr_copy_num_round) %>%
    dplyr::rename(CCN = Chr_copy_num_round_correct)
)
write_source_data(
  paste0(sup_fig_prefix, "Fig4ab"), 
  read_data_csv("chromosome_cnv/20231214_summary_MuLTEE/data_bin.csv") %>%
    dplyr::select(-Sample_ID, -Sample)
)
write_source_data(
  "Fig3bce", 
  dplyr::full_join(
    x = read_data_csv("cluster_size/20230926_artificial/data_summary_rep.csv") %>%
      dplyr::filter(Time == "24h") %>%
      dplyr::select(Condition, Ploidy, Replicate, Weighted_mean_radius) %>%
      dplyr::rename(Cluster.Weighted_mean_radius = Weighted_mean_radius), 
    y = read_data_csv("cell_size_aspect_ratio/20230926_artificial/data_summary_rep.csv") %>%
      dplyr::select(Condition, Ploidy, Replicate, Mean_volume, Mean_AR) %>%
      dplyr::rename(Cell.Mean_volume = Mean_volume, Cell.Mean_AR = Mean_AR), 
    by = c("Condition", "Ploidy", "Replicate"))
)
write_source_data(
  paste0("Fig3bcf_raw_", sup_fig_prefix, "Fig5de"), 
  read_data_csv("cell_size_aspect_ratio/20230926_artificial/data.csv") %>%
    dplyr::select(-ImageID, -Circ., -Solidity)
)
write_source_data(
  paste0("Fig3e_raw_", sup_fig_prefix, "Fig5bc"), 
  read_data_csv("cluster_size/20230926_artificial/data.csv")
)
write_source_data(
  "Fig3g", 
  read_data_csv("competition_assay/20220909_2N4N_corrected_remove_PO/comp_dfit.csv") %>%
    dplyr::rename(Freq_4N_Day0 = StrainB_Perc_0, Freq_4N_Day3 = StrainB_Perc_3)
)
write_source_data(
  paste0(sup_fig_prefix, "Fig6b"), 
  read_data_csv("competition_assay/20220909_2N4N_corrected_remove_PO/cs_data.csv") %>%
    dplyr::filter(Replicate == "Rep1") %>%
    dplyr::select(Condition, Ploidy, Area, cell_top5_mean_area)
)
write_source_data(
  "Fig3h", 
  read_data_csv("ploidy_measurement/20230920_agar_rev_combine/data_plus_anc_t0.csv") %>%
    dplyr::select(-Batch, -Strain, -RawIntDen_bgs, -Ploidy)
)
write_source_data(
  "Fig3i", 
  read_data_csv("ploidy_measurement/20230920_agar_rev_combine/data_plus_anc_t1000.csv") %>%
    dplyr::select(-Batch, -Strain, -RawIntDen_bgs, -PloidyReduction)
)
write_source_data(
  paste0(sup_fig_prefix, "Fig7b"), 
  read_data_csv("cluster_size/20231201_agar/data.csv") %>%
    dplyr::select(-PloidyReduction)
)
write_source_data(
  paste0(sup_fig_prefix, "Fig7b_mean"), 
  read_data_csv("cluster_size/20231201_agar/data_summary.csv") %>%
    dplyr::select(-Count, -PloidyReduction)
)
write_source_data(
  "Fig3j", 
  read_data_csv("ploidy_measurement/20241021_MA_Sayantan/data_summary.csv") %>%
    dplyr::select(-Batch, -Strain, -Background, -Line_EvoTime)
)
write_source_data(
  paste0(sup_fig_prefix, "Fig8a"), 
  read_data_csv("ploidy_measurement/20241021_MA_Sayantan/data.csv") %>%
    dplyr::select(-Batch, -Strain, -Background, -Line_EvoTime, -RawIntDen_bgs, -Ploidy_change)
)
write_source_data(
  "Fig4bcdef_Fig5e", 
  dplyr::full_join(
    x = read_data_csv("cluster_size/20231227_donutspread/data_summary.csv") %>%
      dplyr::select(Condition, Line, EvoTime, Replicate, ColonyMorph, Weighted_mean_radius) %>%
      dplyr::rename(Cluster.Weighted_mean_radius = Weighted_mean_radius), 
    y = read_data_csv("cell_size_aspect_ratio/20231227_donutspread/data_summary.csv") %>%
      dplyr::select(Condition, Line, EvoTime, Replicate, ColonyMorph, Mean_volume, Mean_AR, 
                    Mean_volume.SignPadj, Mean_AR.SignPadj, Mean_volume_AR.Signs) %>%
      dplyr::rename(Cell.Mean_volume = Mean_volume, Cell.Mean_AR = Mean_AR, 
                    Cell.Mean_volume.SignPadj = Mean_volume.SignPadj, 
                    Cell.Mean_AR.SignPadj = Mean_AR.SignPadj, 
                    Cell.Mean_volume_AR.Signs = Mean_volume_AR.Signs), 
    by = c("Condition", "Line", "EvoTime", "Replicate", "ColonyMorph")) #%>% View()
)
write_raw_data(
  "Fig4b",
  read_data_csv("cluster_size/20231227_donutspread/data.csv") %>%
    dplyr::select(-Strain, -CondLine, -StrainBackgroundRep, -StrainBackground)
)
write_raw_data(
  "Fig4cdef_Fig5e",
  read_data_csv("cell_size_aspect_ratio/20231227_donutspread/data.csv") %>%
    dplyr::select(-File, -Strain, -CondLine, -StrainBackgroundRep, -StrainBackground, 
                  -ImageID, -Crop, -Circ., -Solidity, -Round)
)
write_raw_data(
  "Fig4f_dV_stat", 
  read_data_csv("cell_size_aspect_ratio/20231227_donutspread/data_volume_stat.csv")
)
write_raw_data(
  "Fig4f_dAR_stat", 
  read_data_csv("cell_size_aspect_ratio/20231227_donutspread/data_AR_stat.csv")
)
write_source_data(
  "Fig5a", 
  read_data_csv("chromosome_cnv/20231230_donutspread/data_ds.csv") %>%
    dplyr::select(Condition, Line, EvoTime, Replicate, ColonyMorph, 
                  Chr, CCN, StrainBackground, StrainBackground_CCN, StrainBackground_dCCN, StrainBackground_dCCN_label)
)
write_source_data(
  "Fig5b", 
  read_data_csv("chromosome_cnv/20231230_donutspread/data_ds_diff.csv") %>%
    dplyr::select(Condition, Line, EvoTime, Replicate,
                  Chr, D_CCN, S_CCN, DS_dCCN, DS_dCCN_3N_label)
)
write_source_data(
  "Fig5c_evo", 
  read_data_csv("chromosome_cnv/20231230_donutspread/data_macro_chr.csv") %>%
    dplyr::select(Condition, Line, EvoTime, 
                  Chr, Macro_dCCN_sign)
)
write_source_data(
  "Fig5c_DS", 
  read_data_csv("chromosome_cnv/20231230_donutspread/data_ds_diff_dchr_summary_sb.csv") %>%
    dplyr::select(Condition, Line, EvoTime, 
                  Chr, DS_dCCN_sign, Count)
)
write_raw_data(
  paste0(sup_fig_prefix, "Fig9_mut_include3Nstrains"), 
  read_data_csv("mutations/20231230_donutspread/data.csv") %>%
    dplyr::select(-Sample_ID, -Strain, -CondLine, -StrainBackground, -StrainBackgroundRep, -ColonyRep, 
                  -Mut_ID, -Ref_Alt_depth, -Manual_add)
)
write_raw_data(
  paste0(sup_fig_prefix, "Fig9ab"), 
  read_data_csv("mutations/20231230_donutspread/data_diff_chr_acn.csv") %>%
    dplyr::select(-CondLine, -StrainBackground, -StrainBackgroundRep, 
                  -Mut_ID, -DS_dmut_sign)
)
write_source_data(
  paste0(sup_fig_prefix, "Fig9a"), 
  read_data_rds("mutations/20231230_donutspread/data_diff_chr_acn.rds") %>%
    dplyr::group_by(DS_defreq_mode, DS_dACN_abs) %>%
    dplyr::summarize(Mut_num = n()) #%>% View()
)
write_source_data(
  paste0(sup_fig_prefix, "Fig9b"), 
  read_data_rds("mutations/20231230_donutspread/data_diff_chr_acn.rds") %>%
    dplyr::mutate(DS_dCCN_sign = factor(DS_dCCN_sign, levels = c("Gain", "Maintain", "Loss"))) %>%
    dplyr::group_by(DS_defreq_mode, DS_dCCN_sign, DS_dACN_abs) %>%
    dplyr::summarize(Mut_num = n()) #%>% View()
)
write_source_data(
  paste0(sup_fig_prefix, "Fig9c"), 
  read_data_csv("mutations/20231230_donutspread/data_diff_chr_acn_all_np.csv")
)
write_source_data(
  paste0(sup_fig_prefix, "Fig9d"), 
  read_data_csv("mutations/20231230_donutspread/data_donut_all_np.csv")
)
write_source_data(
  paste0(sup_fig_prefix, "Fig9e"), 
  read_data_csv("mutations/20231230_donutspread/data_np_combine.csv")
)
write_source_data(
  paste0(sup_fig_prefix, "Fig9f"), 
  read_data_csv("mutations/20231230_donutspread/data_all_np_combine_incdec_sim.csv") %>%
    dplyr::mutate(Mut_num_hm	= round(Mut_num_total	* Mut_perc_hm / 100)) %>%
    dplyr::relocate(Mut_num_hm, .before = Mut_num_total) %>%
    dplyr::rename(Mut_num_hm_donut = Mut_num_donut) #%>% View()
)
