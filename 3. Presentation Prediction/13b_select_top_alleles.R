# Step 13b: Select TOP HLA-allele presentation score for each n-mer
# July 31, 2022 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: For each of the n-mers, identify the best HLA-allele that it binds to the n-mer.

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)

#  Load Directories -------------------------------------------------------
directory_13 = "{PATH_TO_STEP13_OUTPUT}"

# Load Files --------------------------------------------------------------
setwd(directory_13)
nmer_08 = read_tsv("2023_0812_hlathena_08mer_edited.tsv", col_names = T)
nmer_09 = read_tsv("2023_0812_hlathena_09mer_edited.tsv", col_names = T)
nmer_10 = read_tsv("2023_0812_hlathena_10mer_edited.tsv", col_names = T)
nmer_11 = read_tsv("2023_0812_hlathena_11mer_edited.tsv", col_names = T)


###########################################################################
#  Step 1: Select highest HLA-allele and associated score -----------------
###########################################################################
nmer_08[is.na(nmer_08)] <- ""
nmer_09[is.na(nmer_09)] <- ""
nmer_10[is.na(nmer_10)] <- ""
nmer_11[is.na(nmer_11)] <- ""

for (i in 1:4){
  if (i == 1){nmer_i = nmer_08}
  if (i == 2){nmer_i = nmer_09}
  if (i == 3){nmer_i = nmer_10}
  if (i == 4){nmer_i = nmer_11}
  
  nmer_i_edit = NULL
  for (j in 1:(nrow(nmer_i)/5)){
    print(paste0(i, " ", j))
    nmer_j = nmer_i %>% dplyr::slice(j)
    
    PEPTIDE = nmer_j %>% dplyr::pull(peptide)
    N_FLANK = nmer_j %>% dplyr::pull(n_flank)
    C_FLANK = nmer_j %>% dplyr::pull(c_flank)
    
    decision_j = nmer_i %>% 
      dplyr::filter(peptide == PEPTIDE) %>% 
      dplyr::filter(n_flank == N_FLANK) %>% 
      dplyr::filter(c_flank == C_FLANK)
    
    max_i = decision_j[which.max(decision_j$hlathena_presentation_score),]
    
    nmer_i_edit = rbind(nmer_i_edit, max_i)
  }
  
  if (i == 1){nmer_08_edit = nmer_i_edit}
  if (i == 2){nmer_09_edit = nmer_i_edit}
  if (i == 3){nmer_10_edit = nmer_i_edit}
  if (i == 4){nmer_11_edit = nmer_i_edit}
}

# Export Files ------------------------------------------------------------
setwd(directory_13)
write_tsv(nmer_08_edit, "2023_0812_hlathena_08mer_selected_alleles.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nmer_09_edit, "2023_0812_hlathena_09mer_selected_alleles.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nmer_10_edit, "2023_0812_hlathena_10mer_selected_alleles.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nmer_11_edit, "2023_0812_hlathena_11mer_selected_alleles.tsv", na = "NA", col_names = T, quote_escape = "double")
