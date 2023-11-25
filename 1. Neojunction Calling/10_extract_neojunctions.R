# Title: "Step 10 Extract Neojunctions"
# July 08, 2023 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Extract Neojunctions by filtering with GTEx junctions that have a PSR 
#          of < 0.01 (list derived from Step 09)

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))
library(dplyr)
library(tidyverse)
library(readxl)
library(ggsci)


# Establish Directories ---------------------------------------------------

directory_external   = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/tnejo/altspl/00_external"
directory_out.step07 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/07_count_depth_freq_judge_psr/output"
directory_out.step09 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/09_count_depth_freq_judge_psr_gtex/out"
directory_out.step10 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/10_extract_neojunctions/out"


# Load Files --------------------------------------------------------------

# From Step 7: PSR Tables
setwd(directory_out.step07)
filename_psr.tcga   = "P1_PSR_20200818.tsv"
dataframe_psr.tcga  = read_tsv(filename_psr.tcga, na = c("", "NA"))
# junc.id                    junc.id.overlap           
# <chr>                      <chr>                     
# 1 chr1:-:101342417-101343031 chr1:-:101342420-101343031
# 2 chr1:-:101467100-101479265 chr1:-:101467143-101479265
# 3 chr1:-:109608854-109618542 chr1:-:109608854-109617608
# 4 chr1:-:109822778-109823398 chr1:-:109822778-109823017
# 5 chr1:-:109823111-109823391 chr1:-:109823111-109823398
# 6 chr1:-:109823630-109823729 chr1:-:109823668-109823758

# From Step 9: PSR Tables for GTEx
setwd(directory_out.step09)
filename_psr.gtex   = "2023_0704_psr_retained_and_passed_gtex.tsv"
dataframe_psr.gtex  = read_tsv(filename_psr.gtex, na = c("", "NA"))
# Dimensions
# dataframe_psr.p1       = 6943 x 8 (confirmed) 
# dataframe_psr.idh      = 7104 x 8 (confirmed)
# dataframe_psr.gtex.p1  = 241 x 5 (confirmed) 
# dataframe_psr.gtex.idh = 267 x 5 (confirmed)


# From Step 7: Count, Depth, Freq, Judge Tables for General (1) and IDH1 Subclasses (2)
setwd(directory_out.step07)
filename_count  = "P1_Count_Table_Retained_and_Passed_Junctions_20200818.tsv"
filename_depth  = "P1_Depth_Table_Retained_and_Passed_Junctions_20200818.tsv"
filename_freq   = "P1_Freq_Table_Retained_and_Passed_Junctions_20200818.tsv"
filename_judge  = "P1_Judgement_Table_Retained_and_Passed_Junctions_20200818.tsv"
dataframe_count = read_tsv(filename_count, na = c("", "NA"))
dataframe_depth = read_tsv(filename_depth, na = c("", "NA"))
dataframe_freq  = read_tsv(filename_freq, na = c("", "NA"))
dataframe_judge = read_tsv(filename_judge, na = c("", "NA"))

###########################################################################
#  Step 1: Filter for Neojunctions: PSR -----------------------------------
###########################################################################

# Step 09 generated a list of neojunctions with a GTEx PSR that is less than 0.01
# Use that list to filter for the neojunctions from the TCGA

dataframe_psr.neo = right_join(
  dataframe_psr.tcga %>% dplyr::rename(max.tcga = max),
  dataframe_psr.gtex %>% dplyr::rename(min.gtex = min),
  by = "junc.id")



setwd(directory_out.step10)
filename_psr.p1  = "2023_0708_psr_neojunctions.tsv"
write_tsv(dataframe_psr.neo, filename_psr.p1, na = "NA", col_names = T, quote_escape = "double")



###########################################################################
#  Step 2: Filter for Neojunctions: Count, Depth, Freq, Judge -------------
###########################################################################

# From the above filtered neojunctions, this step simply filters out the Count, Depth, Frequency,
# and Judgement Tables generated in Step 7

list_tables = list(dataframe_count, dataframe_depth, dataframe_freq, dataframe_judge)

list_filenames = list("2023_0708_count_neojunctions.tsv", "2023_0708_depth_neojunctions.tsv", "2023_0708_freq_neojunctions.tsv", "2023_0708_judge_neojunctions.tsv")
i = 1
for (i in 1:length(list_tables)){
  print(i)
  dataframe.i = list_tables[[i]] %>% 
    semi_join(dataframe_psr.neo, by = "junc.id")
  
  # Output the Count, Depth, Freq, and Judge tables
  setwd(directory_out.step10)
  write_tsv(dataframe.i, list_filenames[[i]], na = "NA", col_names = T, quote_escape = "double")
}
















