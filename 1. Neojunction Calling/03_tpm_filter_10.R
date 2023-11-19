# Title: "Step 3 - TPM Filter > 10"
# July 05, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Identify transcripts with TPM values > 10 within TCGA-data with tumor purity > 0.60
# Pipeline: 1. 

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

install.packages("readxl")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("tibble")
install.packages("ggsci")

library(readxl)
library(tidyverse)
library(ggsci)

#  Import Files -----------------------------------------------------------
directory_in  = "/Users/darwinkwok/Documents/UCSF\ 2020\ Okada-Costello\ Lab/Project\ A\ -\ IDH1\ Alternative\ Splicing-Associated\ Neojunctions/03_tpm_filter_10/input"
directory_out = "/Users/darwinkwok/Documents/UCSF\ 2020\ Okada-Costello\ Lab/Project\ A\ -\ IDH1\ Alternative\ Splicing-Associated\ Neojunctions/03_tpm_filter_10/output"
setwd(directory_in)

filename_tpm          = "TcgaTargetGtex_rsem_isoform_log2tpm0001_TCGA_n656__protein_coding_137618tx.txt"  # TPM (Isoform Expression) Table
filename_tcga.tp.0.60 = "TCGA_Patient_List_Post_TumorPurity_Filter_20200604"                              # TCGA samples with purity greater than 0.60 from Step 01
filename_gtf.coding   = "GTF_Protein_Coding_Genes_20200606"                                               # GTF (Protein Coding Genes) from Step 02

dataframe_tpm   = as_tibble(data.table::fread(filename_tpm))
dataframe_tcga  = as_tibble(data.table::fread(filename_tcga.tp.0.60))
dataframe_gtf   = as_tibble(data.table::fread(filename_gtf.coding))


# Check TPM File ----------------------------------------------------------
# Check the dimensions
dataframe_tpm %>% dim()
# [1] 137631    657

# Check to see if there are recurring samples (will end with "-02") and also if "enst" is the first and only column that is not a sample
# The function "!grepl("-01$")" combined with colnames() will retain any columns that does NOT END with "-01"
colnames(dataframe_tpm)[!grepl("-01$", colnames(dataframe_tpm))]
# [1] "enst"

# Check to see if there are any duplicated cases (same Sample ID in the first 12 characters)
# The function "substr" extracts a portion of the string
# The function "duplicated" identifies how many strings are the same in a list
(substr(colnames(dataframe_tpm), 1, 12))[duplicated(substr(colnames(dataframe_tpm), 1, 12))]
# character(0)


# Check TCGA File ---------------------------------------------------------
# Check the dimensions
dataframe_tcga %>% dim()
# [1] 429   8

# Check to see how many LGG and GBM cases are within the TCGA dataframe
dataframe_tcga %>% pull(study) %>% table()
# TCGA-GBM TCGA-LGG 
# 115      314 

# Check to see how many IDH subclasses are within the TCGA dataframe
dataframe_tcga %>% pull(cls) %>% table()
# IDH-A IDH-O IDHwt 
# 140   123   166 


# Check GTF File ----------------------------------------------------------
# Check the dimensions
dataframe_gtf %>% dim()
# [1] 145548     12


###########################################################################
#  Step 1: Filter for TCGA Samples and Convert log2tpm to Standard TPM ----
###########################################################################

# 1. Remove the "-01" from the back of every sample name by only retaining the first 12 symbols in each column name
# FUNCTION substr() extracts the specified characters from a string
# EXAMPLE  'TCGA-19-1787-01' --> 'TCGA-19-1787'
colnames(dataframe_tpm) = substr(colnames(dataframe_tpm), 1, 12)

# 2. Retain genes (rows) that have ENST ID's
#    Remove any rows that do not have "ENST" as the first 4 characters
# FUNCTION filter() retains all rows that satisfies the condition (to be retained, the row must produce a value of TRUE)
dataframe_tpm.edit = dataframe_tpm %>% 
  dplyr::filter(substr(enst, 1, 4) == "ENST")

# 3. Extract samples (columns) shared with TCGA samples
#    Remove any columns that are not found in the TCGA data frame
# FUNCTION select() selects variables in a data frame (e.g. a:f selects columns a through f), can predicate functions (e.g. is.element()) to select variable for a property
#          is.element() performs set union/intersection between two vectors
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  dplyr::select(enst, colnames(dataframe_tpm)[is.element(colnames(dataframe_tpm), dataframe_tcga$case)])

# 4.  Reformat the TPM data into a dataframe "list"
#     Pair ENST value with Sample # to their log2tpm value as each row in the new dataframe
# FUNCTION gather() takes multiple columns and collapses them into their individual rows
#          e.g.         Sample 1    Sample 2             ENST 1   Sample 1   1
#               ENST 1         1           2     -->     ENST 1   Sample 2   2
#               ENST 2         3           4             EnST 2   Sample 1   3
#                                                        ENST 2   Sample 2   4
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  gather("case", "log2tpm", 2:ncol(dataframe_tpm.edit))

# 5. Convert log2 TPM values into standard TPM values
#    Convert log2(tpm + 0.001) into tpm
# FUNCTION as.numeric() converts a variable into a numeric class
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  mutate(tpm = round((2^as.numeric(log2tpm)) - 0.001, 4))

# 6. Convert NA values to 0
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  mutate(tpm = ifelse(is.na(tpm), 0, tpm))

# 7. Convert negative values to 0
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  mutate(tpm = ifelse(tpm < 0, 0, tpm)) # fix negative values
  
# 8. Remove the log2tpm column.
#    Now that we have the tpm column, this column is no longer necessary
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  dplyr::select( - log2tpm)

# 9. Revert the table format to its original form
#    gather() originally collapsed the table into a list format
# FUNCTION spread() reverts gather() and transforms the dataframe to back into the original TPM format
dataframe_tpm.edit = dataframe_tpm.edit %>% 
  spread(case, tpm)

# 10. Check the dimensions 
dataframe_tpm.edit %>% dim()
# [1] 137618    425

# 11. Check to see if all transcripts in the TPM frame are protein-coding
#     There should be no missing or different genes between the TPM and GTF dataframe
# FUNCTION anti_join() returns all of the values from the Table 1 without a match in Table 2
dataframe_tpm.edit %>% 
  dplyr::select(1:5) %>%                                 # Select the first 5 samples to not overwhelm the system
  mutate(enst = substr(enst, 1, 15)) %>%                 # Extract the ENST ID by extracting the first 15 characters
  anti_join(dataframe_gtf, by = "enst")                  # Compare the ENST ID between the TPM and GTF dataframe
# A tibble: 0 x 5


###########################################################################
#  Step 2: Identify and Retain Expressed Genes ----------------------------
###########################################################################

dataframe_tpm.median.final = list()
transcript_tpm.median.passed = list()

# 1. Generate a list of cases and their associated samples
#    The list will contain 6 groups, and each group will contain a list of their associated glioma classifications
for (i in 1:6) {
  print(i) ;
  cases = list(
    dataframe_tcga$case,
    dataframe_tcga$case[dataframe_tcga$study == "TCGA-GBM"],
    dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDHwt"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDH-A"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDH-O"]
  )[[i]]
  print(cases)
  
  # 2. Retrieve the columns (and enst) from the TPM dataframe that corresponds to the cases within each iteration
  dataframe_tpm.i = dataframe_tpm.edit %>% 
    dplyr::select(enst, colnames(dataframe_tpm.edit[is.element(colnames(dataframe_tpm.edit), cases)]))
  
  # 3. Generate a new dataframe that only contains the ENST ID and the median TPM value
  dataframe_tpm.median = dataframe_tpm.i %>% 
    mutate(median = apply(dataframe_tpm.i[, -1], 1, median)) %>% 
    dplyr::select(enst, median)
  
  # 4. Retain transcripts that have a median tpm ≥ 10 within the subgroup
  list_tpm.median.passed = dataframe_tpm.median %>% 
    dplyr::filter(median >= 10) %>% 
    pull(enst)
  
  transcript_tpm.median.passed[[i]] = list_tpm.median.passed
  
  colnames(dataframe_tpm.median)[2] = c("all", "TCGA-GBM", "TCGA-LGG", "IDHwt", "IDH-A", "IDH-O")[i]
  dataframe_tpm.median.final[[i]] = dataframe_tpm.median
}


###########################################################################
#  Step 3: Generate the Final TPM Table Passing TPM ≥ 10 ------------------
###########################################################################
# Merge the medians into one data frame from the different disease groups

for (i in 1:6) {
  # 1. Generate a new data frame by which everything is merged into
  #    This new data frame is generated when i = 1
  if (i ==1) {
    transcript_uniq = transcript_tpm.median.passed[[i]]
    dataframe_tpm.merged = dataframe_tpm.median.final[[i]]
  }
  else {
    transcript_uniq = unique(c(transcript_uniq, transcript_tpm.median.passed[[i]]))
    dataframe_tpm.merged = dataframe_tpm.merged %>% 
      full_join(dataframe_tpm.median.final[[i]], by = "enst")
  }
}

# 2. Only retain transcripts that have a TPM ≥ 10 in at least 1 disease category in dataframe_tpm.edit
dataframe_tpm.final = dataframe_tpm.edit %>% 
  dplyr::filter(is.element(enst, transcript_uniq))


###########################################################################
#  Step 4: Generate the TPM Summary Tables (Raw and Filtered) -------------
###########################################################################

# Generate the raw summary table of the median TPM in each disease group (merged)
dataframe_tpm.summary.raw = dataframe_gtf %>% 
  dplyr::select(enst, symbol, ensg) %>% 
  right_join(dataframe_tpm.merged %>% 
               mutate(enst = substr(enst, 1, 15)), by = "enst") %>% 
  arrange(desc(all))

# Generate the summary table of the median TPM with the filter of TPM > 10 applied
dataframe_tpm.summary.passed = dataframe_tpm.summary.raw %>% 
  dplyr::filter(is.element(enst, substr(c(transcript_uniq), 1, 15)))


###########################################################################
#  Step 5: Generate the Final GTF Table Passing TPM ≥ 10 ------------------
###########################################################################

dataframe_gtf.passed = dataframe_gtf %>%
  semi_join(dataframe_tpm.summary.passed, by = "enst") %>% 
  dplyr::select(chr = X1, source = X2, start = X4, end = X5, strand = X7, enst, symbol, ensg)


###########################################################################
#  Step 6: Output Data ----------------------------------------------------
###########################################################################

setwd(directory_out)

# TPM Dataframe (Passed TPM ≥ 10 Filter)
filename_tpm.final = "TPM_Filter10_ProteinCodingTx_n11617_20200717"
write_tsv(dataframe_tpm.final, filename_tpm.final, na = "NA", col_names = T, quote_escape = "double")

# TPM Summary Table (Raw)
filename_tpm.summary.raw = "TPM_Summary_Table_Complete_20200717"
write_tsv(dataframe_tpm.summary.raw, filename_tpm.summary.raw, na = "NA", col_names = T, quote_escape = "double")

# TPM Summary Table (Passed TPM ≥ 10 Filter)
filename_tpm.summary.passed = "TPM_Summary_Table_FilterTPM10_20200717"
write_tsv(dataframe_tpm.summary.passed, filename_tpm.summary.passed, na = "NA", col_names = T, quote_escape = "double")

# GTF Dataframe (Passed TPM ≥ 10 Filter)
filename_gtf.passed = "GTF_ProteinCoding_FilterTPM10_20200717"
write_tsv(dataframe_gtf.passed, filename_gtf.passed, na = "NA", col_names = T, quote_escape = "double")

