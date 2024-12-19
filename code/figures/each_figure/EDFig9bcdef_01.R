library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(stringr)
library(glue)
library(RColorBrewer)
library(scales)
library(ggh4x)

setwd("~/Documents/projects/R_data_analysis/mutations/20231230_donutspread")
rm(list=ls())
source("~/Documents/projects/R_data_analysis/global_setup.R")

# Parameters
in_folder_path <- "snpEff_vcf"
in_file_suffix <- "_denovo_final_snpEff.vcf"
out_fig_path <- "."
conditions <- c("PA")
lines <- as.character(1:5)
evotimes <- c("t600", "t1000")
reps <- paste0("rep", 1:3)
reps_label <- gsub("rep", "", reps, fixed = TRUE)
reps_label2 <- gsub("rep", "Rep", reps, fixed = TRUE)
colonymorphs <- c("D", "S")
colonymorphs_shape <- c("circle", "bullet")
names(colonymorphs_shape) <- colonymorphs
chrs <- as.character(as.roman(1:16))
chrs_color <- hue_pal()(16)
chrs_color <- chrs_color[c(seq(1,16,3), seq(2,16,3), seq(3,16,3))]
names(chrs_color) <- chrs
#show_col(chrs_color, labels = TRUE)
chr_copy_nums_color <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(9)[2:9]  #PuOr
names(chr_copy_nums_color) <- 1:8
#show_col(chr_copy_nums_color, labels = TRUE)
strainbackgroundreps_3N <- c("PA_4_t600_rep1", "PA_4_t1000_rep3")

# Load data
data_ds_diff_dchr <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/chromosome_cnv/20231230_donutspread/data_ds_diff_dchr.rds")
data_chr_diff <- readRDS("/Users/kaitong/Documents/projects/R_data_analysis/chromosome_cnv/20231230_donutspread/data_ds_diff.rds")

##### Load and pre-process data

# Load mutation data
# Input file name: e.g., "t600_PA1_D1_denovo_final_snpEff.vcf"
load_one_vcf <- function(file) {
  # Read and process file
  df <- read.table(file, header = FALSE, sep = "\t", 
                   row.names = NULL, stringsAsFactors = FALSE) %>%  # this somehow ignores the first rows in VCF file
    dplyr::select(1,2,4,5,8,10) %>%
    dplyr::mutate(V1 = gsub("chr", "", V1, fixed = TRUE))
  colnames(df) <- c("Chr", "Pos", "Ref", "Alt", "Info", "Info2")
  ## Process Info column
  df$Info <- sapply(df$Info,
                    FUN = function(x) {strsplit(x, ",", fixed = TRUE)[[1]][1]},
                    simplify = TRUE, USE.NAMES = FALSE)
  df <- df %>%
    dplyr::mutate(Info = gsub("|", ",", Info, fixed = TRUE)) %>%  # somehow separating by "|" causes error thus first convert "|" to ","
    tidyr::separate(col = "Info", into = paste0("V", 1:16),
                    sep = ",", remove = TRUE, convert = FALSE) %>%
    dplyr::select(!all_of(paste0("V", c(1,6:9,16))))
  colnames(df) <- c(colnames(df)[1:4],
                    c("Variant_type", "Variant_impact", "gene_name", "gene_id", "HGVS.c", "HGVS.p",
                      "cDNA.pos/length", "CDS.pos/length", "AA.pos/length",
                      "Distance"),
                    "Info2")
  ## Process Info2 column
  ## Normally, it should be like "0/1:65,23:88:99:662,0,2316"
  ## Sometimes, it can be like "1/2:0,46,16:62:99:2604,549,380,1705,0,1582", due to multiple Alt nucleotides detected
  ## Handle such cases: extract first value as Ref_depth, sum up remaining values (no matter one or more) as Alt_depth
  ## Also extract Total_depth reported by VCF file to compare later
  df <- df %>%
    tidyr::separate(col = "Info2", into = c("X1", "Ref_Alt_depth", "Total_depth_vcf", "X4"),
                    sep = ":", remove = TRUE, convert = TRUE, extra = "merge") %>%
    dplyr::select(-X1, -X4) %>%
    dplyr::mutate(
      split_values = str_split(Ref_Alt_depth, ",", simplify = FALSE) %>% map(~as.numeric(.x)),
      Ref_depth = map_dbl(split_values, ~.x[1]),
      Alt_depth = map_dbl(split_values, ~sum(.x[-1], na.rm = TRUE))
    ) %>%
    dplyr::select(-split_values) %>%
    dplyr::mutate(Total_depth = Ref_depth + Alt_depth, 
                  Alt_freq = Alt_depth / Total_depth)
  ## Record current column names
  measure_vars <- colnames(df)
  # Extract metadata
  filename <- basename(file)
  sample_id <- substr(filename, 1, stri_length(filename) - stri_length(in_file_suffix))
  metadata <- strsplit(sample_id, split = "_", fixed = TRUE)[[1]]
  df$Sample_ID <- sample_id
  df$Condition <- substr(metadata[2], 1, 2)
  df$Line <- substr(metadata[2], 3, 3)
  df$EvoTime <- metadata[1]
  df$ColonyMorph <- substr(metadata[3], 1, 1)
  df$Replicate <- paste0("rep", substr(metadata[3], 2, 2))
  # Reorder column names
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}
data <- ldply(.data = list.files(path = in_folder_path, 
                                 pattern = paste0("*", in_file_suffix), full.names = TRUE),
              .fun = load_one_vcf)
# Check
# data %>%
#   dplyr::filter(Total_depth_vcf != Total_depth) %>%
#   View()
# Good
data <- data %>%
  dplyr::select(-Total_depth_vcf)

# Combine with mutations_AF_list.csv
# Sayantan generated a sheet that contains manually curated false-gain/loss mutations and their allele frequencies
# Combine it with the mutations in the VCF files to generate a final data table containing all mutations detected in donut/spread genomes
# Load data
mutations_AF_list <- read.csv("mutations_AF_list.csv", row.names = NULL, stringsAsFactors = FALSE) %>%
  dplyr::mutate(ColonyMorph = ifelse(gain.loss == "Lost", "S", "D")) %>%  # it is thought to be lost in spread, but it is actually still in spread
  dplyr::mutate(Sample_ID = paste(Time, Strain, paste0(ColonyMorph, Rep), sep = "_")) %>%
  dplyr::select(-Gene_name, -Gene_ID, -Variant_ID, -ID, -gain.loss) %>%
  dplyr::rename(EvoTime = Time, Replicate = Rep, Ref_Alt_depth = AD) %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line"), sep = 2, 
                  remove = TRUE, convert = FALSE) %>%
  tidyr::separate(col = "Mut_ID", into = c("Chr", "Pos"), sep = "_", 
                  remove = TRUE, convert = TRUE) %>%
  tidyr::separate(col = "Ref_Alt_depth", into = c("Ref_depth", "Alt_depth"), sep = ",", 
                  remove = FALSE, convert = TRUE) %>%
  dplyr::mutate(Replicate = paste0("rep", Replicate),
                Chr = gsub("chr", "", Chr, fixed = TRUE), 
                Total_depth = Ref_depth + Alt_depth, 
                Alt_freq = Alt_depth / Total_depth
                )
# Fill in mutation info from corresponding donut or spread in the same replicate
columns_to_fill <- setdiff(colnames(data), colnames(mutations_AF_list))
columns_for_join <- c("Condition", "Line", "EvoTime", "Replicate", "Chr", "Pos")
mutations_AF_list <- mutations_AF_list %>%
  dplyr::left_join(data %>% dplyr::select(all_of(c(columns_for_join, columns_to_fill))), 
                   by = columns_for_join)
# Combine with data
data <- dplyr::bind_rows(data, mutations_AF_list)
data$Manual_add <- c(rep(FALSE, nrow(data) - nrow(mutations_AF_list)), 
                     rep(TRUE, nrow(mutations_AF_list)))

# Remove and rename some samples
length(unique(data$Sample_ID))
data <- data %>%
  tidyr::unite(Strain, all_of(c("Condition", "Line", "EvoTime", "ColonyMorph", "Replicate")), 
               sep = "_", remove = FALSE, na.rm = TRUE) %>%
  # Remove some samples
  dplyr::filter(!(
    Strain %in% c(paste0("PA_1_t600_", c("D_rep1", "S_rep1", "D_rep4", "S_rep4")), 
                  paste0("PA_4_t600_", c("D", "S"), "_rep3"))
  )) %>%
  # Rename samples
  # PA1 t600 rep2/3 -> rep1/2
  # PA4 t600 rep4 -> rep3
  dplyr::mutate(Strain = case_when(Strain == "PA_1_t600_D_rep2" ~ "PA_1_t600_D_rep1",
                                   Strain == "PA_1_t600_S_rep2" ~ "PA_1_t600_S_rep1",
                                   Strain == "PA_1_t600_D_rep3" ~ "PA_1_t600_D_rep2",
                                   Strain == "PA_1_t600_S_rep3" ~ "PA_1_t600_S_rep2",
                                   Strain == "PA_4_t600_D_rep4" ~ "PA_4_t600_D_rep3",
                                   Strain == "PA_4_t600_S_rep4" ~ "PA_4_t600_S_rep3",
                                   TRUE ~ Strain))
length(unique(data$Sample_ID))

# Data formatting
# Note: now Sample_ID is based on input file names, Strain column is correct, re-generate other sample metadata variables
id_vars_core <- c("Condition", "Line", "EvoTime", "ColonyMorph", "Replicate")
data <- data %>%
  dplyr::select(!all_of(id_vars_core)) %>%
  tidyr::separate(col = "Strain", into = c("Condition", "Line", "EvoTime", "ColonyMorph", "Replicate"), sep = "_", 
                  remove = FALSE, convert = FALSE) %>%
  tidyr::unite(col = "CondLine", Condition, Line, sep = "", remove = FALSE) %>%
  tidyr::unite(col = "StrainBackground", Condition, Line, EvoTime, sep = "_", remove = FALSE) %>%
  tidyr::unite(col = "StrainBackgroundRep", StrainBackground, Replicate, sep = "_", remove = FALSE) %>%
  dplyr::relocate(StrainBackground, .before = StrainBackgroundRep) %>%
  tidyr::unite(col = "ColonyRep", ColonyMorph, Replicate, sep = "", remove = FALSE) %>%
  dplyr::mutate(ColonyRep = gsub("rep", "", ColonyRep, fixed = TRUE)) %>%
  dplyr::relocate(ColonyRep, .before = Condition) %>%
  dplyr::arrange(Line, factor(EvoTime, levels = evotimes), Replicate, ColonyMorph, 
                 factor(Chr, levels = chrs))
strains <- unique(data$Strain)
condlines <- unique(data$CondLine)
strainbackgrounds <- unique(data$StrainBackground)
strainbackgrounds_label <- strainbackgrounds %>%
  gsub("_t", " t", ., fixed = TRUE) %>%
  gsub("_", "", ., fixed = TRUE)
strainbackgroundreps <- unique(data$StrainBackgroundRep)
strainbackgroundreps_label <- strainbackgroundreps %>%
  gsub("PA_", "PA", ., fixed = TRUE) %>%
  gsub("_t", " t", ., fixed = TRUE) %>%
  gsub("_rep", " (", ., fixed = TRUE) %>%
  paste0(., ")")
strainbackgroundreps_no3N <- strainbackgroundreps[!strainbackgroundreps %in% strainbackgroundreps_3N]
strainbackgroundreps_label_no3N <- strainbackgroundreps_label[!strainbackgroundreps %in% strainbackgroundreps_3N]
colonyreps <- unique(data$ColonyRep)
data <- data %>%
  dplyr::mutate(Strain = factor(Strain, levels = strains), 
                CondLine = factor(CondLine, levels = condlines), 
                StrainBackground = factor(StrainBackground, levels = strainbackgrounds),
                StrainBackgroundRep = factor(StrainBackgroundRep, levels = strainbackgroundreps), 
                ColonyRep = factor(ColonyRep, levels = colonyreps), 
                Line = factor(Line, levels = lines), 
                EvoTime = factor(EvoTime, levels = evotimes), 
                Replicate = factor(Replicate, levels = reps), 
                ColonyMorph = factor(ColonyMorph, levels = colonymorphs), 
                Chr = factor(Chr, levels = chrs)
  )

# Check
length(unique(data$Strain))
unique(data$Strain)
summary(data)

# Group variant types into certain categories
table(data$Variant_type)
table(paste(data$Variant_type, data$Variant_impact)) #%>% View()
# conservative_inframe_deletion MODERATE 6 
# disruptive_inframe_deletion MODERATE 6 
# frameshift_variant&stop_lost&splice_region_variant HIGH 6 
# missense_variant MODERATE 1314
# splice_region_variant&intron_variant LOW 6 
# splice_region_variant&stop_retained_variant LOW 2 
# stop_gained HIGH 53 
# synonymous_variant LOW 342 
# upstream_gene_variant MODIFIER 431
# Summary: 
# HIGH: stop_gained, frameshift_variant&stop_lost&splice_region_variant
# MODERATE: missense_variant, conservative_inframe_deletion, disruptive_inframe_deletion
# LOW: synonymous_variant, splice_region_variant&intron_variant, splice_region_variant&stop_retained_variant
# MODIFIER: upstream_gene_variant
variant_types <- unique(data$Variant_type)
data <- data %>% 
  dplyr::mutate(Variant_impact = factor(Variant_impact, levels = toupper(variant_impacts), labels = variant_impacts))

# Create a column for mutation ID
data <- data %>%
  tidyr::unite(col = "Mut_ID", Chr, Pos, Ref, Alt, sep = "_", remove = FALSE)

##### Save data

save.image("plot.RData")

