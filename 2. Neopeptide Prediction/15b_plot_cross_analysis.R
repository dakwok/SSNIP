# Title: "Step 26: Pan-algorithm analysis"
# July 31, 2022 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Analyze the candidate lists from all three algorithm approaches:
#             - HLAthena
#             - MHCFlurry 2.0
#             - antigen.garnish

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggsci)
library(data.table)
library("Biostrings") # For opening FASTA


#  Load Directories -------------------------------------------------------
directory_15 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/15_cross_analyze_algorithms/out"

# Load Files --------------------------------------------------------------
# HLAthena
setwd(directory_15)
ha_10percentile = read_tsv("2023_0812_hlathena_01percentile_all_nmers.tsv", col_names = T)
mf_10percentile = read_tsv("2023_0812_mhcflurry_01percentile_all_nmers.tsv", col_names = T)
df_10percentile = read_tsv("2023_0812_cross_alg_01percentile_all_nmers.tsv", col_names = T)
df_all =          read_tsv("2023_0812_cross_alg_all_nmers.tsv", col_names = T)


# WT peptide sequences
# From Step 11: All confirmed aa sequences (WT and MUT)
setwd("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/11_aaseq_prediction/output")
filename_aa = "20201026_res_aa_prediction_confirmed.tsv"
aa_table = read_tsv(filename_aa, col_names = T)

# From External: FASTA File containing all proteins
setwd("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/external/fasta_uniprot")
filename_fasta = "UP000005640_9606.fasta"
fastaFile <- readAAStringSet(filename_fasta)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
fasta <- data.frame(seq_name, sequence)

###########################################################################
#  Step 2: Get 10th percentile candidates from HLAthena and MHCFlurry -----
###########################################################################
ha_10percentile[is.na(ha_10percentile)] <- ""
mf_10percentile[is.na(mf_10percentile)] <- ""

# Combine the two dataframes
ha_mf_combined = subset(full_join(ha_10percentile, mf_10percentile), select = c(allele, peptide, n_flank, c_flank, hlathena_presentation_score, mhcflurry_presentation_score))
ha_mf_combined = ha_mf_combined[!duplicated(ha_mf_combined[c(1,2,3,4)]),]

# Combine the two algorithm's TOP candidates
ha_mf_combined$shared = "Top 10%tile in One Alg. (HA or MF)"
ha_mf_combined$shared[is.na(ha_mf_combined$hlathena_presentation_score) == F & is.na(ha_mf_combined$mhcflurry_presentation_score) == F] = "Top 10%tile in Both Alg. (HA and MF)"

# Fill in the missing data
df_all[is.na(df_all)] <- ""

plot_ha_mf_combined = NULL
for (i in 1:nrow(ha_mf_combined)){
  print((i/nrow(ha_mf_combined)*100))
  
  row_i = ha_mf_combined %>% dplyr::slice(i)
  
  HLA     = row_i %>% dplyr::pull(allele)
  NMER    = row_i %>% dplyr::pull(peptide)
  N_FLANK = row_i %>% dplyr::pull(n_flank)
  C_FLANK = row_i %>% dplyr::pull(c_flank)
  
  HLATHENA  = row_i %>% dplyr::pull(hlathena_presentation_score)
  MHCFLURRY = row_i %>% dplyr::pull(mhcflurry_presentation_score)
  
  if (is.na(HLATHENA) == T){
    search_i = df_all %>% 
      dplyr::filter(allele == HLA & peptide == NMER & n_flank == N_FLANK & c_flank == C_FLANK) %>% 
      dplyr::pull(hlathena_presentation_score)
    row_i$hlathena_presentation_score = search_i
  }
  
  if (is.na(MHCFLURRY) == T){
    search_i = df_all %>% 
      dplyr::filter(allele == HLA & peptide == NMER & n_flank == N_FLANK & c_flank == C_FLANK) %>% 
      dplyr::pull(mhcflurry_presentation_score)
    row_i$mhcflurry_presentation_score = search_i
  }
  
  plot_ha_mf_combined = rbind(plot_ha_mf_combined, row_i)
}

###########################################################################
#  Step 4: Generate plots -------------------------------------------------
###########################################################################
# Plot the combined data sets and highlight the candidates that are shared between both candidates
ggplot(plot_ha_mf_combined, aes(hlathena_presentation_score, mhcflurry_presentation_score)) +
  geom_point(aes(color = shared), size = 3, alpha = 0.8) +
  theme_minimal() + 
  labs(title="") +
  #scale_color_manual(values = c("gray70", "#00AFBB", "#FC4E07", "#E7B800")) + 
  scale_color_manual(values = c("#489fb5", "gray70", "#ff5400", "#E7B800")) + 
  theme(plot.title = element_text(size = 16,  family="Helvetica", face = "bold", hjust = 0.5),
        text = element_text(size = 16,  family="Helvetica", face = "bold"),
        legend.position="bottom") +
  labs(x = "HLAthena Presentation Score", y = "MHCFlurry 2.0 Presentation Score") + 
  guides(color = guide_legend(nrow = length(unique(plot_ha_mf_combined$shared)))) +
  
  labs(colour = "") 
# Dimensions (2022_1009): 5.75 x 5.25

# Count the number of shared candidates
shared_2 = nrow(plot_ha_mf_combined %>% dplyr::filter(shared != "1 Hit (HA or MF)"))         # 636  # Updated 2023_0812 -> 4304

setwd(directory_15)
write_tsv(plot_ha_mf_combined,   "2023_0812_cross_analysis_summary_ha_mf.tsv", na = "NA", col_names = T, quote_escape = "double")

# test = read_tsv("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/12_nmer_generation/output_2022_0801/2022_0801_complete_list_all_nmers.tsv", col_names = T)
test2 = test[!duplicated(test[c(13,15,16)]),]  # 281,302

test3 = read_tsv("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/12_nmer_generation/output_2022_0801/2022_0801_hlathenalist_msic_allmers.tsv", col_names = T)
test4 = test3[!duplicated(test3[c(1,2,3)]),]  # 17,562
