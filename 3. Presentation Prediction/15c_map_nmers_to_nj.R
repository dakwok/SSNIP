# Step 15c: Map back to the original neojunctions
# March 07, 2021 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Map back the original amino acid sequences to identify "immunogenic neojunctions"
###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggsci)
library(data.table)

#  Load Directories -------------------------------------------------------
directory_11 = "{PATH_TO_STEP11_OUTPUT}"
directory_15 = "{PATH_TO_STEP15_OUTPUT}"

# Load Files --------------------------------------------------------------
setwd(directory_11)
filename_peptide = "2023_0710_res_aa_prediction_confirmed.tsv"
aa_conversion = read_tsv(filename_peptide, col_names = T)

setwd(directory_15)
df_combined = read_tsv("2023_0812_cross_analysis_summary_ha_mf.tsv", col_names = T)

df_top = df_combined %>% dplyr::filter(shared == "Top 10%tile in Both Alg. (HA and MF)")


###########################################################################
#  Step 1: Map back to the original AA ------------------------------------
###########################################################################
# Remove the STOP codons from the alt aa.seq in aa_conversion
aa_valid = NULL
for (i in 1:nrow(aa_conversion)){
  print(i)
  nmer_i = aa_conversion %>% dplyr::slice(i)
  peptide_full = nmer_i %>% dplyr::pull(aa.seq.alt)
  
  peptide_valid = gsub("\\*.*", "*", peptide_full)  # Remove everything after the first STOP codon
  print(peptide_valid)
  
  peptide_valid = gsub("[*]", "", peptide_valid) # Remove the STOP codon
  
  nmer_i = nmer_i %>% 
    mutate(alt.aa.valid = peptide_valid)
  
  aa_valid = rbind(aa_valid, nmer_i)
}

# Map the n-mers to their original AA sequence 
aa_shared = NULL
for (i in 1:nrow(df_top)){
  print(i)
  
  # Iterate through each peptide candidate
  PEPTIDE = df_top %>% dplyr::slice(i) %>% pull(peptide)
  N_FLANK = df_top %>% dplyr::slice(i) %>% pull(n_flank)
  C_FLANK = df_top %>% dplyr::slice(i) %>% pull(c_flank)
  ALLELE  = df_top %>% dplyr::slice(i) %>% pull(allele)
  SCORE_H = df_top %>% dplyr::slice(i) %>% pull(hlathena_presentation_score)
  SCORE_M = df_top %>% dplyr::slice(i) %>% pull(mhcflurry_presentation_score)
  SCORE_A = mean(c(SCORE_H, SCORE_M))
  SHARED  = df_top %>% dplyr::slice(i) %>% pull(shared)
  
  candidate_i = paste0(N_FLANK, PEPTIDE, C_FLANK)

  # Search for the candidate in the converted amino acid list
  peptide_i            = aa_valid[grep(candidate_i, aa_valid$alt.aa.valid),]
  peptide_i$nmer       = PEPTIDE
  peptide_i$n_flank    = N_FLANK
  peptide_i$c_flank    = C_FLANK
  peptide_i$hla_allele = ALLELE
  peptide_i$score_hlathena  = SCORE_H
  peptide_i$score_mhcflurry = SCORE_M
  peptide_i$score_average   = SCORE_A
  peptide_i$shared     = SHARED
  
  aa_shared = rbind(aa_shared, peptide_i)
  
}

nj_immuno = as.data.frame(t(table(aa_shared$junc.id)))[,c(2,3)]
nj_immuno = nj_immuno[order(nj_immuno$Freq, decreasing = TRUE),]

aa_map = subset(aa_shared, select = c(nmer, n_flank, c_flank, hla_allele, score_hlathena, score_mhcflurry, score_average, junc.id, symbol, ensg, fs, type, shared))

setwd(directory_15)
write_tsv(nj_immuno, "2023_0812_immunogenic_njs_n422.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(aa_map, "2023_0812_neoA_to_neoJ_map_ALL_summary.tsv", na = "NA", col_names = T, quote_escape = "double")
