# Step 13a: Edit HLA Dataframes
# July 30, 2022 | Darwin Kwok M.S. | University of California San Francisco

# Modify the HLAthena output results to make an editted dataframe

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)

#  Load Directories -------------------------------------------------------
directory_hlathena_output = "{PATH_TO_HLATHENA_OUTPUT}"
directory_out             = "{PATH_TO_STEP13_OUTPUT}"


# Load Files --------------------------------------------------------------
setwd(directory_hlathena_output)
nmer_08 = read_tsv("2023_0812_08mers_hlathena_msic.tsv", col_names = T)
nmer_09 = read_tsv("2023_0812_09mers_hlathena_msic.tsv", col_names = T)
nmer_10 = read_tsv("2023_0812_10mers_hlathena_msic.tsv", col_names = T)
nmer_11 = read_tsv("2023_0812_11mers_hlathena_msic.tsv", col_names = T)



###########################################################################
#  Step 1: Cross Analyze Algorithms (HLAthena and MHCflurry 2.0) ----------
###########################################################################
# Edit the dataframe to match up with the MHCFlurry output as much as possible
# Remove the "-" in the upstream and downstream flanking aa sequences
nmer_08$ctex_up <- gsub("-", "", as.character(nmer_08$ctex_up))
nmer_08$ctex_dn <- gsub("-", "", as.character(nmer_08$ctex_dn))
nmer_09$ctex_up <- gsub("-", "", as.character(nmer_09$ctex_up))
nmer_09$ctex_dn <- gsub("-", "", as.character(nmer_09$ctex_dn))
nmer_10$ctex_up <- gsub("-", "", as.character(nmer_10$ctex_up))
nmer_10$ctex_dn <- gsub("-", "", as.character(nmer_10$ctex_dn))
nmer_11$ctex_up <- gsub("-", "", as.character(nmer_11$ctex_up))
nmer_11$ctex_dn <- gsub("-", "", as.character(nmer_11$ctex_dn))

# Isolate only the MSiC score column for each HLA allele and make a new column for HLA allele
# 1. Iterate through each n_mer (i)
for (i in 1:4){
  if (i == 1){nmer_i = nmer_08}
  if (i == 2){nmer_i = nmer_09}
  if (i == 3){nmer_i = nmer_10}
  if (i == 4){nmer_i = nmer_11}
  
  nmer_final = NULL
  # 2. Iterate through each HLA allele (j)
  for (j in 1:5){
    if (j == 1){hla_j = tibble("allele" = "HLA-A0101"); nmer_j = subset(nmer_i, select = c(n_mer, ctex_up, ctex_dn, MSiC_A0101))}
    if (j == 2){hla_j = tibble("allele" = "HLA-A0201"); nmer_j = subset(nmer_i, select = c(n_mer, ctex_up, ctex_dn, MSiC_A0201))}
    if (j == 3){hla_j = tibble("allele" = "HLA-A0301"); nmer_j = subset(nmer_i, select = c(n_mer, ctex_up, ctex_dn, MSiC_A0301))}
    if (j == 4){hla_j = tibble("allele" = "HLA-A1101"); nmer_j = subset(nmer_i, select = c(n_mer, ctex_up, ctex_dn, MSiC_A1101))}
    if (j == 5){hla_j = tibble("allele" = "HLA-A2402"); nmer_j = subset(nmer_i, select = c(n_mer, ctex_up, ctex_dn, MSiC_A2402))}
    
    nmer_append = cbind(hla_j, nmer_j)
    colnames(nmer_append) = c("allele", "peptide", "n_flank", "c_flank", "hlathena_presentation_score")
    
    nmer_final = rbind(nmer_final, nmer_append)
  }
  
  if (i == 1){nmer_08_edit = nmer_final}
  if (i == 2){nmer_09_edit = nmer_final}
  if (i == 3){nmer_10_edit = nmer_final}
  if (i == 4){nmer_11_edit = nmer_final}
}

setwd(directory_out)
write_tsv(nmer_08_edit, "2023_0812_hlathena_08mer_edited.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nmer_09_edit, "2023_0812_hlathena_09mer_edited.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nmer_10_edit, "2023_0812_hlathena_10mer_edited.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nmer_11_edit, "2023_0812_hlathena_11mer_edited.tsv", na = "NA", col_names = T, quote_escape = "double")

