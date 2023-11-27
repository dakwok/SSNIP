# Step 15a: Cross Analyze Algorithms
# July 31, 2022 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Generate a dataframe cross-analyzing the MHCFlurry 2.0 and HLAthena scores

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))
library(readxl)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(ggsci)
library(data.table)

# Font Import
font_import()
loadfonts(device = "pdf")

#  Load Directories -------------------------------------------------------
directory_13 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/13_hlathena_analysis/out_processed_files"
directory_14 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out"
directory_15 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/15_cross_analyze_algorithms/out"

# Load Files --------------------------------------------------------------
setwd(directory_13)
hla_08 = read_tsv("2023_0812_hlathena_08mer_edited.tsv", col_names = T)
hla_09 = read_tsv("2023_0812_hlathena_09mer_edited.tsv", col_names = T)
hla_10 = read_tsv("2023_0812_hlathena_10mer_edited.tsv", col_names = T)
hla_11 = read_tsv("2023_0812_hlathena_11mer_edited.tsv", col_names = T)

setwd(directory_14)
mhc_08 = read_csv("2023_0812_08mers_flank_mhcflurry.csv", col_names = T)
mhc_09 = read_csv("2023_0812_09mers_flank_mhcflurry.csv", col_names = T)
mhc_10 = read_csv("2023_0812_10mers_flank_mhcflurry.csv", col_names = T)
mhc_11 = read_csv("2023_0812_11mers_flank_mhcflurry.csv", col_names = T)


###########################################################################
#  Step 1: Generate dataframe combining HLAthena and MHCFlurry results ----
###########################################################################
hla_all = rbind(hla_08, hla_09, hla_10, hla_11)
mhc_all = rbind(mhc_08, mhc_09, mhc_10, mhc_11)

# Remove all NA's
hla_all[is.na(hla_all)] <- ""
mhc_all[is.na(mhc_all)] <- ""

# Generate a column for n_mer ID by combining n_mer, upstream, downstream sequences
combined_all <- dplyr::left_join(hla_all, mhc_all, by=c("allele", "peptide", "n_flank", "c_flank"))
combined_all = combined_all[!duplicated(combined_all[c(1,2,3,4)]),]

combined_all_1 = subset(combined_all, select = c("allele", "peptide", "n_flank", "c_flank", "hlathena_presentation_score", "mhcflurry_presentation_score"))
combined_all_1$combined_score = rowMeans(combined_all[,c('hlathena_presentation_score', 'mhcflurry_presentation_score')], na.rm=TRUE)

combined_all_1 = combined_all_1[order(combined_all_1$combined_score, decreasing = TRUE),]


###########################################################################
#  Step 2: Generate the top 10th percentile dataframe ---------------------
###########################################################################
# Calculate the 10th percentile number
#total_nmers = nrow(hla_08)/5 + nrow(hla_09)/5 + nrow(hla_10)/5 + nrow(hla_11)/5
#percentile_10 = as.integer(total_nmers/10)

# EDIT 2023_0428: The 10th percentile should be based on the total number of nmer:HLA combinations, not nmers themself
total_nmers = nrow(hla_08) + nrow(hla_09) + nrow(hla_10) + nrow(hla_11)
percentile_10 = as.integer(total_nmers/100)

# Order by score
hla_all = hla_all[order(hla_all$hlathena_presentation_score, decreasing = TRUE),]
mhc_all = mhc_all[order(mhc_all$mhcflurry_presentation_score, decreasing = TRUE),]

# Generate a new dataframe for the 10th percentile
hla_10percentile = hla_all[1:percentile_10,]
mhc_10percentile = mhc_all[1:percentile_10,]

cross_10percentile <- dplyr::semi_join(hla_10percentile, mhc_10percentile, by=c("allele", "peptide", "n_flank", "c_flank"))
cross_10percentile <- dplyr::inner_join(combined_all_1, cross_10percentile, by=c("allele", "peptide", "n_flank", "c_flank"))
cross_10percentile$hlathena_presentation_score.y <- NULL

cross_10percentile_unique <- cross_10percentile[!duplicated(cross_10percentile[, c("peptide", "n_flank", "c_flank")]), ]


###########################################################################
#  Step 3: Output files ---------------------------------------------------
###########################################################################
setwd(directory_15)
write_tsv(hla_10percentile,   "2023_0812_hlathena_01percentile_all_nmers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(mhc_10percentile,   "2023_0812_mhcflurry_01percentile_all_nmers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(cross_10percentile, "2023_0812_cross_alg_01percentile_all_nmers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(cross_10percentile_unique, "2023_0812_cross_alg_01percentile__unique_all_nmers.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(combined_all_1,     "2023_0812_cross_alg_all_nmers.tsv", na = "NA", col_names = T, quote_escape = "double")

