# 14. Generate input CSV file for running MHCFlurry 2.0
# July 30, 2022 | Darwin Kwok M.S. | University of California San Francisco
#     The input CSV file is expected to contain columns 
#     “allele”, “peptide”, and, optionally, “n_flank”, and “c_flank”.

# To run on command line: $ mhcflurry-predict INPUT.csv –out RESULT.csv

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)

#  Load Directories -------------------------------------------------------
directory_12 = "{PATH_TO_STEP12_OUTPUT}"
directory_14 = "{PATH_TO_STEP14_OUTPUT}"

#  Load Files -------------------------------------------------------------
setwd(directory_12)
nmers_08 = read_tsv("2023_0812_hlathenalist_msic_08mers.tsv", col_names = T)
nmers_09 = read_tsv("2023_0812_hlathenalist_msic_09mers.tsv", col_names = T)
nmers_10 = read_tsv("2023_0812_hlathenalist_msic_10mers.tsv", col_names = T)
nmers_11 = read_tsv("2023_0812_hlathenalist_msic_11mers.tsv", col_names = T)


###########################################################################
#  Step 1: Edit dataframes ------------------------------------------------
###########################################################################
# Remove the TPM columns
nmers_08 = subset(nmers_08, select = -c(TPM))
nmers_09 = subset(nmers_09, select = -c(TPM))
nmers_10 = subset(nmers_10, select = -c(TPM))
nmers_11 = subset(nmers_11, select = -c(TPM))

# Change the column names to the suitable ones needed to run MHCFLurry 2.0
colnames(nmers_08) = c("peptide", "n_flank", "c_flank")
colnames(nmers_09) = c("peptide", "n_flank", "c_flank")
colnames(nmers_10) = c("peptide", "n_flank", "c_flank")
colnames(nmers_11) = c("peptide", "n_flank", "c_flank")

# Remove all "-" from columns
nmers_08$n_flank <- gsub("-", "", as.character(nmers_08$n_flank))
nmers_08$c_flank <- gsub("-", "", as.character(nmers_08$c_flank))
nmers_09$n_flank <- gsub("-", "", as.character(nmers_09$n_flank))
nmers_09$c_flank <- gsub("-", "", as.character(nmers_09$c_flank))
nmers_10$n_flank <- gsub("-", "", as.character(nmers_10$n_flank))
nmers_10$c_flank <- gsub("-", "", as.character(nmers_10$c_flank))
nmers_11$n_flank <- gsub("-", "", as.character(nmers_11$n_flank))
nmers_11$c_flank <- gsub("-", "", as.character(nmers_11$c_flank))

allele_a0101 = tibble("allele" = c("HLA-A0101"))
allele_a0201 = tibble("allele" = c("HLA-A0201"))
allele_a0301 = tibble("allele" = c("HLA-A0301"))
allele_a1101 = tibble("allele" = c("HLA-A1101"))
allele_a2402 = tibble("allele" = c("HLA-A2402"))

for (i in 1:4){
  if (i == 1){nmer_i = nmers_08}
  if (i == 2){nmer_i = nmers_09}
  if (i == 3){nmer_i = nmers_10}
  if (i == 4){nmer_i = nmers_11}
  
  nmer_i_a0101 = cbind(allele_a0101, nmer_i)
  nmer_i_a0201 = cbind(allele_a0201, nmer_i)
  nmer_i_a0301 = cbind(allele_a0301, nmer_i)
  nmer_i_a1101 = cbind(allele_a1101, nmer_i)
  nmer_i_a2402 = cbind(allele_a2402, nmer_i)
  
  n_mer_i_all = rbind(nmer_i_a0101, nmer_i_a0201, nmer_i_a0301, nmer_i_a1101, nmer_i_a2402)
  
  if (i == 1){nmers_08_final = n_mer_i_all}
  if (i == 2){nmers_09_final = n_mer_i_all}
  if (i == 3){nmers_10_final = n_mer_i_all}
  if (i == 4){nmers_11_final = n_mer_i_all}
}

setwd(directory_14)
write_csv(nmers_08_final, "2023_0812_08mer_mhcflurry_input.csv", na = "NA", col_names = T, quote_escape = "double")
write_csv(nmers_09_final, "2023_0812_09mer_mhcflurry_input.csv", na = "NA", col_names = T, quote_escape = "double")
write_csv(nmers_10_final, "2023_0812_10mer_mhcflurry_input.csv", na = "NA", col_names = T, quote_escape = "double")
write_csv(nmers_11_final, "2023_0812_11mer_mhcflurry_input.csv", na = "NA", col_names = T, quote_escape = "double")

