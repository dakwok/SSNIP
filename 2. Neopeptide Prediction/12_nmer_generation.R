# Title: "Step 12: Generate AA n-mers from Predicted AA Sequences"
# July 10, 2023 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: The previous step generated all of the predicted amino acid sequences that would be generated
#          from the neojunctions. This step will generate all of the predicted n-mers (from 8-mers to 10-mers),
#          and return only the n-mers that cannot be generated in the respective WT amino acid sequence.

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)

# Establish Directories ---------------------------------------------------

directory_00 = "{PATH_TO_META_FILES}"
directory_01 = "{PATH_TO_STEP01_OUTPUT}"
directory_11 = "{PATH_TO_STEP11_OUTPUT}"
directory_12 = "{PATH_TO_STEP12_OUTPUT}"

# Load Files --------------------------------------------------------------

# From External: TPM Ves
setwd(directory_00)
filename_tpm = "TcgaTargetGtex_rsem_isoform_log2tpm0001_TCGA_n656__protein_coding_137618tx.txt"
tpm = read_tsv(filename_tpm, na = c("", "NA"), col_types = cols(.default = col_double(), enst = col_character()))

# From Step 01: List of Patients
setwd(directory_01)
filename_tcga = "TCGA_Patient_List_Post_TumorPurity_Filter_20200604.tsv"
tcga = read_tsv(filename_tcga, na = c("", "NA"))
  
# From Step 11: List of Predicted AA Sequences
setwd(directory_11)
filename_aa = "2023_0710_res_aa_prediction_confirmed.tsv"
aa = read_tsv(filename_aa, na = c("", "NA"))
# 


###########################################################################
#  Step 1: Edit the TPM dataframe -----------------------------------------
###########################################################################

# The following steps are taken from Step 3 (03_tpm_filter_10.R)
# For greater detail on each step, refer to that step.

colnames(tpm) = substr(colnames(tpm), 1, 12)
tpm.edit = tpm %>% dplyr::filter(substr(enst, 1, 4) == "ENST")
tpm.edit = tpm.edit %>% dplyr::select(enst, colnames(tpm)[is.element(colnames(tpm), tcga$case)])
tpm.edit = tpm.edit %>% gather("case", "log2tpm", 2:ncol(tpm.edit))
tpm.edit = tpm.edit %>% mutate(tpm = round((2^as.numeric(log2tpm)) - 0.001, 4))
tpm.edit = tpm.edit %>% mutate(tpm = ifelse(is.na(tpm), 0, tpm))
tpm.edit = tpm.edit %>% mutate(tpm = ifelse(tpm < 0, 0, tpm))
tpm.edit = tpm.edit %>% dplyr::select( - log2tpm)
tpm.edit = tpm.edit %>% spread(case, tpm)
tpm.edit = tpm.edit %>%  mutate(enst.1 = substr(enst, 1, 15))

# Select out only the median TPM for all of the genes
tpm.edit = tpm.edit %>% 
  mutate(TPM = apply(tpm.edit[, -1], 1, median)) %>% 
  dplyr::select(enst.1, TPM)


###########################################################################
#  Step 2: Edit and Filter AA Sequences -----------------------------------
###########################################################################

# Remove NA sequences or any sequences that do not start with a start codon (M)

aa_valid = aa %>% 
  dplyr::select(junc.id, symbol, enst, enst.model, aa.change, aa.seq.wt, aa.seq.alt, ln.wt, ln.alt, ln.diff, sc, type) %>% 
  dplyr::filter(!is.na(aa.seq.wt) & !is.na(aa.seq.alt)) %>%              # Remove any sequences with NA in either the WT or Mut variants
  dplyr::filter(aa.seq.wt != aa.seq.alt) %>%                             # Ensure that the WT and Mut variants are not the same
  dplyr::filter(grepl("^M", aa.seq.wt) & grepl("^M", aa.seq.alt)) %>%    # Ensure that both AA sequences start with M
  mutate(aa.seq.alt = gsub("\\*.*", "*", aa.seq.alt))                    # Remove the first stop codon (*) and any amino acids following it


# aa %>% dim() %>% print()       --> 467  17
# aa_valid %>% dim() %>% print() --> 390  12

# Add the TPM values as the final column
aa_valid = aa_valid %>% 
  left_join(tpm.edit, by = c("enst" = "enst.1"))

# Turn all NA values under TPM into blanks
aa_valid$TPM[is.na(aa_valid$TPM)] = ""

# Remove duplicate aa.alts
aa_valid.1 = aa_valid[!duplicated(aa_valid$aa.seq.alt), ]
 
# ###########################################################################
# #  Step 3: Edit the AA sequences to only Include 8-11 mers ----------------
# ###########################################################################
# # Generate list of n-mers for MUTANT peptides (aa.seq.alt)
# setwd(directory_12)
# for (h in 8:11){
#   summary_nmer_h = NULL
#   for (i in 1:nrow(aa_valid.1)){
#     aa_i = aa_valid.1 %>% dplyr::slice(i)
#     aa_sequence = aa_valid.1 %>% dplyr::slice(i) %>% pull(aa.seq.alt)
#     
#     # Iterate through the entire length of the peptide
#     for (j in 1:nchar(aa_sequence)){
#       print(paste(h, round((i/nrow(aa_valid.1)*100), digits = 2), round((j/nchar(aa_sequence)*100), digits = 2)))
#       
#       # Generate the n-mer and their flanks
#       nmer_h = substr(aa_sequence, start = j, stop = j+h-1)
#       flank_10n = substr(aa_sequence, start = j-10, stop = j-1)
#       flank_10c = substr(aa_sequence, start = j+h, stop = j+h+10-1)
#       flank_30n = substr(aa_sequence, start = j-30, stop = j-1)
#       flank_30c = substr(aa_sequence, start = j+h, stop = j+h+30-1)
#       
#       append_row = aa_i %>% 
#         mutate(n_flank = flank_30n) %>%
#         mutate(n_mer = nmer_h) %>%
#         mutate(c_flank = flank_30c)
#       
#       if(nchar(nmer_h) == h){summary_nmer_h = rbind(summary_nmer_h, append_row)}
#     }
#   }
#   if (h == 8){nmer_8_mut = summary_nmer_h ; write_tsv(nmer_8_mut, "2023_0802_all_iterations_alt_list_08mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
#   if (h == 9){nmer_9_mut = summary_nmer_h ; write_tsv(nmer_9_mut, "2023_0802_all_iterations_alt_list_09mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
#   if (h == 10){nmer_10_mut = summary_nmer_h ; write_tsv(nmer_10_mut, "2023_0802_all_iterations_alt_list_10mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
#   if (h == 11){nmer_11_mut = summary_nmer_h ; write_tsv(nmer_11_mut, "2023_0802_all_iterations_alt_list_11mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
# }


# Generate list of n-mers for WILD-TYPE peptides (aa.seq.alt)
for (h in 8:11){
  summary_nmer_h = NULL
  for (i in 1:nrow(aa_valid.1)){
    aa_i = aa_valid.1 %>% dplyr::slice(i)
    aa_sequence = aa_valid.1 %>% dplyr::slice(i) %>% pull(aa.seq.wt)
    
    # Iterate through the entire length of the peptide
    for (j in 1:nchar(aa_sequence)){
      print(paste(h, round((i/nrow(aa_valid.1)*100), digits = 2), round((j/nchar(aa_sequence)*100), digits = 2)))
      
      # Generate the n-mer and their flanks
      nmer_h = substr(aa_sequence, start = j, stop = j+h-1)
      flank_10n = substr(aa_sequence, start = j-10, stop = j-1)
      flank_10c = substr(aa_sequence, start = j+h, stop = j+h+10-1)
      flank_30n = substr(aa_sequence, start = j-30, stop = j-1)
      flank_30c = substr(aa_sequence, start = j+h, stop = j+h+30-1)
      
      append_row = aa_i %>% 
        mutate(n_flank = flank_30n) %>%
        mutate(n_mer = nmer_h) %>%
        mutate(c_flank = flank_30c)
      
      if(nchar(nmer_h) == h){summary_nmer_h = rbind(summary_nmer_h, append_row)}
    }
  }
  if (h == 8){nmer_8_wt = summary_nmer_h ; write_tsv(nmer_8_wt,  "2023_0802_all_iterations_wt_list_08mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
  if (h == 9){nmer_9_wt = summary_nmer_h ; write_tsv(nmer_9_wt,  "2023_0802_all_iterations_wt_list_09mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
  if (h == 10){nmer_10_wt = summary_nmer_h ; write_tsv(nmer_10_wt,  "2023_0802_all_iterations_wt_list_10mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
  if (h == 11){nmer_11_wt = summary_nmer_h ; write_tsv(nmer_11_wt,  "2023_0802_all_iterations_wt_list_11mers.tsv", na = "NA", col_names = T, quote_escape = "double")}
}


###########################################################################
#  Step 4: Remove aa.seq.alt n-mers plus flanks found in aa.seq.wt --------
###########################################################################
setwd(directory_12)
nmer08_mut = read_tsv("2023_0802_all_iterations_alt_list_08mers.tsv", col_names = T)
nmer09_mut = read_tsv("2023_0802_all_iterations_alt_list_09mers.tsv", col_names = T)
nmer10_mut = read_tsv("2023_0802_all_iterations_alt_list_10mers.tsv", col_names = T)
nmer11_mut = read_tsv("2023_0802_all_iterations_alt_list_11mers.tsv", col_names = T)

nmer08_wt = read_tsv("2023_0802_all_iterations_wt_list_08mers.tsv", col_names = T)
nmer09_wt = read_tsv("2023_0802_all_iterations_wt_list_09mers.tsv", col_names = T)
nmer10_wt = read_tsv("2023_0802_all_iterations_wt_list_10mers.tsv", col_names = T)
nmer11_wt = read_tsv("2023_0802_all_iterations_wt_list_11mers.tsv", col_names = T)

nmer_08_neo <- anti_join(nmer08_mut, nmer08_wt, by = c("n_flank", "n_mer", "c_flank")) # 50563 neoantigens
nmer_09_neo <- anti_join(nmer09_mut, nmer09_wt, by = c("n_flank", "n_mer", "c_flank")) # 50801 neoantigens
nmer_10_neo <- anti_join(nmer10_mut, nmer10_wt, by = c("n_flank", "n_mer", "c_flank")) # 51039 neoantigens
nmer_11_neo <- anti_join(nmer11_mut, nmer11_wt, by = c("n_flank", "n_mer", "c_flank")) # 51275 neoantigens

nmer_all = rbind((rbind(rbind(nmer_08_neo, nmer_09_neo), nmer_10_neo)), nmer_11_neo)

###########################################################################
#  Step 5: Export Complete Files ------------------------------------------
###########################################################################
setwd(directory_12)
# Export the Complete n-mer files
filename_8mer  = "2023_0812_complete_list_08mers.tsv"
filename_9mer  = "2023_0812_complete_list_09mers.tsv"
filename_10mer = "2023_0812_complete_list_10mers.tsv"
filename_11mer = "2023_0812_complete_list_11mers.tsv"
filename_all   = "2023_0812_complete_list_all_mers.tsv"

write_tsv(nmer_08_neo, filename_8mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_09_neo, filename_9mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_10_neo, filename_10mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_11_neo, filename_11mer, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(nmer_all, filename_all, na = "NA", col_names = T, quote_escape = "double") ;


###########################################################################
#  Step 6: Filter for cancer-specific n-mers ------------------------------
###########################################################################
library("Biostrings") # For opening FASTA
setwd(directory_12)
nmer_8  = read_tsv(filename_8mer, na = c("", "NA"),)
nmer_9  = read_tsv(filename_9mer, na = c("", "NA"),)
nmer_10 = read_tsv(filename_10mer, na = c("", "NA"),)
nmer_11 = read_tsv(filename_11mer, na = c("", "NA"),)

# From External: FASTA File containing all proteins
setwd("{PATH_TO_FASTA_UNIPROT")
filename_fasta = "UP000005640_9606.fasta"
fastaFile <- readAAStringSet(filename_fasta)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
fasta <- data.frame(seq_name, sequence)

for (h in 4:4){
  if (h == 1){df_h = nmer_8}
  if (h == 2){df_h = nmer_9}
  if (h == 3){df_h = nmer_10}
  if (h == 4){df_h = nmer_11}
  nmers_valid = NULL
  for (i in 1:nrow(df_h)){
    print(h)
    print((i/nrow(df_h)*100))
    
    nmer_i = df_h[i,]
    PEPTIDE = as.character(nmer_i %>% dplyr::pull(n_mer))
    WT_PEPTIDE = fasta[grep(PEPTIDE, fasta$sequence),]
    
    if (nrow(WT_PEPTIDE) < 1){nmers_valid = rbind(nmers_valid, nmer_i)}
  }
  if (h == 1){df_valid_08 = nmers_valid} # 6808 --> 5315  # 2023_0812 = 16158
  if (h == 2){df_valid_09 = nmers_valid} # 6809 --> 5427  # 2023_0812 = 16704
  if (h == 3){df_valid_10 = nmers_valid} # 6641 --> 5518  # 2023_0812 = 17077
  if (h == 4){df_valid_11 = nmers_valid} # 6725 --> 5527  # 2023_0812 = 17425
}

setwd(directory_12)
write_tsv(df_valid_08, "2023_0812_cancer_specific_08mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(df_valid_09, "2023_0812_cancer_specific_09mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(df_valid_10, "2023_0812_cancer_specific_10mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(df_valid_11, "2023_0812_cancer_specific_11mers.tsv", na = "NA", col_names = T, quote_escape = "double") ;



# NOTE: After changing to readAA rather than readDNA, the total changed from 21489 --> 21787






setwd(directory_12)
df_valid_08 = read_tsv("2023_0812_cancer_specific_08mers.tsv", col_names = T); df_valid_08$c_flank <- ifelse(is.na(df_valid_08$c_flank), "", df_valid_08$c_flank); df_valid_08$n_flank <- ifelse(is.na(df_valid_08$n_flank), "", df_valid_08$n_flank); df_valid_08$n_flank <- gsub("\\*", "", df_valid_08$n_flank); df_valid_08$c_flank <- gsub("\\*", "", df_valid_08$c_flank); df_valid_08 <- df_valid_08[!grepl("\\*", df_valid_08$n_mer), ]
df_valid_09 = read_tsv("2023_0812_cancer_specific_09mers.tsv", col_names = T); df_valid_09$c_flank <- ifelse(is.na(df_valid_09$c_flank), "", df_valid_09$c_flank); df_valid_09$n_flank <- ifelse(is.na(df_valid_09$n_flank), "", df_valid_09$n_flank); df_valid_09$n_flank <- gsub("\\*", "", df_valid_09$n_flank); df_valid_09$c_flank <- gsub("\\*", "", df_valid_09$c_flank); df_valid_09 <- df_valid_09[!grepl("\\*", df_valid_09$n_mer), ]
df_valid_10 = read_tsv("2023_0812_cancer_specific_10mers.tsv", col_names = T); df_valid_10$c_flank <- ifelse(is.na(df_valid_10$c_flank), "", df_valid_10$c_flank); df_valid_10$n_flank <- ifelse(is.na(df_valid_10$n_flank), "", df_valid_10$n_flank); df_valid_10$n_flank <- gsub("\\*", "", df_valid_10$n_flank); df_valid_10$c_flank <- gsub("\\*", "", df_valid_10$c_flank); df_valid_10 <- df_valid_10[!grepl("\\*", df_valid_10$n_mer), ]
df_valid_11 = read_tsv("2023_0812_cancer_specific_11mers.tsv", col_names = T); df_valid_11$c_flank <- ifelse(is.na(df_valid_11$c_flank), "", df_valid_11$c_flank); df_valid_11$n_flank <- ifelse(is.na(df_valid_11$n_flank), "", df_valid_11$n_flank); df_valid_11$n_flank <- gsub("\\*", "", df_valid_11$n_flank); df_valid_11$c_flank <- gsub("\\*", "", df_valid_11$c_flank); df_valid_11 <- df_valid_11[!grepl("\\*", df_valid_11$n_mer), ]

# Function to pad strings with hyphens to reach a total length of 30 characters
pad_to_30_front <- function(s) {
  str_pad(s, width = 30, side = "left", pad = "-")
}

pad_to_30_back <- function(s) {
  str_pad(s, width = 30, side = "right", pad = "-")
}

df_valid_08$ctex_up <- sapply(df_valid_08$n_flank, pad_to_30_front); df_valid_08$ctex_dn <- sapply(df_valid_08$c_flank, pad_to_30_back)
df_valid_09$ctex_up <- sapply(df_valid_09$n_flank, pad_to_30_front); df_valid_09$ctex_dn <- sapply(df_valid_09$c_flank, pad_to_30_back)
df_valid_10$ctex_up <- sapply(df_valid_10$n_flank, pad_to_30_front); df_valid_10$ctex_dn <- sapply(df_valid_10$c_flank, pad_to_30_back)
df_valid_11$ctex_up <- sapply(df_valid_11$n_flank, pad_to_30_front); df_valid_11$ctex_dn <- sapply(df_valid_11$c_flank, pad_to_30_back)


# Generate the HLAthena Compatible n-mer files for each of the algorithms (MSi, MSiC, and MSiCE)
# Files for the MSiCE algorithm
MSiCE_8  = df_valid_08 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM") 
MSiCE_9  = df_valid_09 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM")
MSiCE_10 = df_valid_10 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM")
MSiCE_11 = df_valid_11 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn", "TPM")

# Files for the MSiC algorithm
MSiC_8  = df_valid_08 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")
MSiC_9  = df_valid_09 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")
MSiC_10 = df_valid_10 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")
MSiC_11 = df_valid_11 %>% dplyr::select("n_mer", "ctex_up", "ctex_dn") %>% 
  mutate(TPM = "")

# Files for the MSiC algorithm
MSi_8 = df_valid_08 %>% dplyr::select("n_mer") %>%
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "") 
MSi_9 = df_valid_09 %>% dplyr::select("n_mer") %>% 
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "")
MSi_10 = df_valid_10 %>% dplyr::select("n_mer") %>% 
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "")
MSi_11 = df_valid_11 %>% dplyr::select("n_mer") %>% 
  mutate(ctex_up = "") %>% 
  mutate(ctex_dn = "") %>% 
  mutate(TPM = "")

# Export the HLAthena Compatible n-mer files
filename_MSi_8  = "2023_0812_hlathenalist_msi_08mers.tsv"
filename_MSi_9  = "2023_0812_hlathenalist_msi_09mers.tsv"
filename_MSi_10 = "2023_0812_hlathenalist_msi_10mers.tsv"
filename_MSi_11 = "2023_0812_hlathenalist_msi_11mers.tsv"

filename_MSiC_8  = "2023_0812_hlathenalist_msic_08mers.tsv"
filename_MSiC_9  = "2023_0812_hlathenalist_msic_09mers.tsv"
filename_MSiC_10 = "2023_0812_hlathenalist_msic_10mers.tsv"
filename_MSiC_11 = "2023_0812_hlathenalist_msic_11mers.tsv"

filename_MSiCE_8  = "2023_0812_hlathenalist_msice_08mers.tsv"
filename_MSiCE_9  = "2023_0812_hlathenalist_msice_09mers.tsv"
filename_MSiCE_10 = "2023_0812_hlathenalist_msice_10mers.tsv"
filename_MSiCE_11 = "2023_0812_hlathenalist_msice_11mers.tsv"

write_tsv(MSi_8, filename_MSi_8, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSi_9, filename_MSi_9, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSi_10, filename_MSi_10, na = "NA", col_names = T, quote_escape = "double")
write_tsv(MSi_11, filename_MSi_11, na = "NA", col_names = T, quote_escape = "double")

write_tsv(MSiC_8, filename_MSiC_8, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiC_9, filename_MSiC_9, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiC_10, filename_MSiC_10, na = "NA", col_names = T, quote_escape = "double")
write_tsv(MSiC_11, filename_MSiC_11, na = "NA", col_names = T, quote_escape = "double")

write_tsv(MSiCE_8, filename_MSiCE_8, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiCE_9, filename_MSiCE_9, na = "NA", col_names = T, quote_escape = "double") 
write_tsv(MSiCE_10, filename_MSiCE_10, na = "NA", col_names = T, quote_escape = "double")
write_tsv(MSiCE_11, filename_MSiCE_11, na = "NA", col_names = T, quote_escape = "double")

