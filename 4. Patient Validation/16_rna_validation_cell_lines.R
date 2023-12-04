# Step 27: RNA Validation of NJ via cell lines
# February 13, 2022 | Darwin Kwok M.S. | University of California San Francisco

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggsci)
library(data.table)

#  Load Directories -------------------------------------------------------
directory_15  = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/15_cross_analyze_algorithms/out_2022_0802"
directory_out = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/27_rna_validation_cell_lines/out"

# Load Files --------------------------------------------------------------
setwd(directory_15)
candidates = read_tsv("2022_0815_cross_analysis_summary_ha_mf_ag.tsv", col_names = T)


###########################################################################
#  Step 2: Identified shared HLA-A0201 Candidates -------------------------
###########################################################################
candidates_2alg_hlaa0201 = candidates %>% 
  dplyr::filter(shared != "1 Hit (HA or MF)") %>% 
  dplyr::filter(allele == "HLA-A0201")

setwd(directory_out)
filename_hlaa0201 = "2022_0815_shared_candidates_hlathena_mhcflurry_HLA_A0201.tsv"
write_tsv(candidates_2alg_hlaa0201, filename_hlaa0201, na = "NA", col_names = T, quote_escape = "double")



###########################################################################
#  Step 3: Identify Candidates that are found in Cell Lines ---------------
###########################################################################
# Cell lines to look at RNAseq with: SF10417 and SF10602 (both HLA*A0201)

# Load Neojunction Count Table for both cell lines
directory_00 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/external/rnaseq_sj.out.tab_20201120"
directory_11 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/11_aaseq_prediction/output"
filename_sj10417 = "SF10417-IDH1KO_RNAseq_hg19_SJ.out.tab"
filename_sj10602 = "SF10602p_RNAseq_hg19_SJ.out.tab"
filename_peptide = "20201026_res_aa_prediction_confirmed.tsv"

setwd(directory_00)
sj_10417 = read_tsv(filename_sj10417, col_names = F)
sj_10602 = read_tsv(filename_sj10602, col_names = F)

setwd(directory_11)
aa_conversion = read_tsv(filename_peptide, col_names = T)

colnames(sj_10417) = c("chr", "int.start", "int.end", "strand", "int.motif", "annot", "n.uniq.map", "n.mult.map", "max.spl.overhang")
colnames(sj_10602) = c("chr", "int.start", "int.end", "strand", "int.motif", "annot", "n.uniq.map", "n.mult.map", "max.spl.overhang")

# Put all of the SJ.out.tab files into one compact list to run
list_sj = list(sj_10417, sj_10602)

for (i in 1:length(list_sj)){
  sj_i = list_sj[[i]]
  sj_i = sj_i %>% 
    # Remove splicing junctions in each file that are not found in the raw reference splicing junction file: dataframe_sj
    # Make a new entry within list_sj for each iteration of files
    mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
    dplyr::filter(strand !=	"undefined") %>%
    dplyr::filter(chr != "M") %>% 
    mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
    dplyr::select(chr, int.start, int.end, strand, n.uniq.map, junc.id)
    #anti_join(dataframe_sj, by = "junc.id") %>%
    # Extract junctions with 10 or more supportive read counts
    #dplyr::filter(n.uniq.map >= 10) %>%
    #dplyr::select(chr, int.start, int.end, strand, junc.id)
  
  list_sj[[i]] = sj_i
}

# Reassign the filtered SJ.out.tab dataframes to their original variables
sj_10417_edit = list_sj[[1]]    # Passed = 29,653 junctions (Nonannotated, Supportive read count > 10)
sj_10602_edit = list_sj[[2]]  # Passed = 20,255 junctions (Nonannotated, Supportive read count > 10)


###########################################################################
# Identify putative RNA-level expression in cell lines --------------------
###########################################################################
# Identify which of the 249 NJs are expressed on an RNA-level in the cell lines
# Get the 249 NJs from the Judgement Table
directory_nj = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/10_extract_neojunctions/output"
filename_nj  = "PSR_Neojunctions_20201020.tsv"
setwd(directory_nj)
nj_249 = subset(read_tsv(filename_nj, col_names = T), select = "junc.id")

# Identify the number of unique neojunction hits with at least 1 HIT
nj_sf10417_00 = semi_join(sj_10417_edit, nj_249, by = "junc.id") %>% 
  dplyr::filter(n.uniq.map > 0)
nj_sf10602_00 = semi_join(sj_10602_edit, nj_249, by = "junc.id") %>% 
  dplyr::filter(n.uniq.map > 0)

# Identify the number of unique neojunction hits with at least 10 HITS
nj_sf10417_10 = semi_join(sj_10417_edit, nj_249, by = "junc.id") %>% 
  dplyr::filter(n.uniq.map > 9)
nj_sf10602_10 = semi_join(sj_10602_edit, nj_249, by = "junc.id") %>% 
  dplyr::filter(n.uniq.map > 9)

# Identify the neojunctions COMMON between cell lines
nj_shared_00 = subset(semi_join(nj_sf10417_00, nj_sf10602_00, by = "junc.id"), select = "junc.id")
nj_shared_10 = subset(semi_join(nj_sf10417_10, nj_sf10602_10, by = "junc.id"), select = "junc.id")
# nj_shared_00 = 31
# nj_shared_10 = 14


# Identify the neojunctions found in EITHER cell line
nj_either_00 = subset(full_join(nj_sf10417_00, nj_sf10602_00, by = "junc.id"), select = "junc.id")
nj_either_10 = subset(full_join(nj_sf10417_10, nj_sf10602_10, by = "junc.id"), select = "junc.id")
# nj_either_00 = 82
# nj_either_10 = 45

# Total the neojunctions found in both cell lines

setwd(directory_out)
write_tsv(nj_sf10417_00, "2022_0401_nj_in_sf10417_atleast01_hit.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_sf10602_00, "2022_0401_nj_in_sf10602_atleast01_hit.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_sf10417_10, "2022_0401_nj_in_sf10417_atleast10_hits.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_sf10602_10, "2022_0401_nj_in_sf10602_atleast10_hits.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_shared_00,  "2022_0401_shared_nj_in_sf10417_sf10602_01_hit.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_shared_10,  "2022_0401_shared_nj_in_sf10417_sf10602_10_hits.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_either_00,  "2022_0404_either_nj_in_sf10417_sf10602_01_hit.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_either_10,  "2022_0404_either_nj_in_sf10417_sf10602_10_hits.tsv", na = "NA", col_names = T, quote_escape = "double")


###########################################################################
# Validate expression with detected CPTAC NJs -----------------------------
###########################################################################
directory_ms = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/19_cptac_output/out_rscript_dk_detection"
filename_ms  = "2021_1216_detected_neojunction_derived_peptides.tsv"
setwd(directory_ms)
nj_cptac = subset(read_tsv(filename_ms, col_names = T), select = "Proteins")

nj_cptac_list = tibble("junc.id" = character())
# Split the NJs apart and generate a new table 
for (i in 1:nrow(nj_cptac)){
  peptide_i = as.character(nj_cptac %>% dplyr::slice(i))
  list_nj = strsplit(peptide_i, ";")[[1]]
  for (j in 1:length(list_nj)){
    nj_j = list_nj[j]
    id_j = strsplit(nj_j, "_")[[1]][2]
    append_j = tibble("junc.id" = id_j)
    nj_cptac_list = rbind(nj_cptac_list, append_j)
  }
}

# Filter out any duplicate NJs
nj_cptac_uniq = nj_cptac_list[!duplicated(nj_cptac_list[,c('junc.id')]),]

nj_shared_rna_pep_01 = semi_join(nj_shared_00, nj_cptac_uniq) # 5 NJs
nj_shared_rna_pep_10 = semi_join(nj_shared_10, nj_cptac_uniq) # 3 NJs
nj_either_rna_pep_00 = semi_join(nj_either_00, nj_cptac_uniq) # 7 NJs
nj_either_rna_pep_10 = semi_join(nj_either_10, nj_cptac_uniq) # 6 NJs

setwd(directory_out)
write_tsv(nj_shared_rna_pep_01, "2022_0404_nj_in_rnaseq_bothcelllines_and_cptac_01hit.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_shared_rna_pep_10, "2022_0404_nj_in_rnaseq_bothcelllines_and_cptac_10hits.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_shared_rna_pep_10, "2022_0404_nj_in_rnaseq_eithercelllines_and_cptac_00hits.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(nj_shared_rna_pep_10, "2022_0404_nj_in_rnaseq_eithercelllines_and_cptac_10hits.tsv", na = "NA", col_names = T, quote_escape = "double")


setwd("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/19_cptac_output/out_rscript_dk_detection")
write_tsv(nj_cptac_uniq, "2022_0401_unique_njs_cptac.tsv", na = "NA", col_names = T, quote_escape = "double")







# -------------------------------------------------------------------------
# Remove the STOP codons from the alt aa.seq in aa_conversion
# -------------------------------------------------------------------------

aa_valid = NULL
for (i in 1:nrow(aa_conversion)){
  print(i)
  nmer_i = aa_conversion %>% dplyr::slice(i)
  peptide_full = nmer_i %>% pull(aa.seq.alt)
  
  peptide_valid = gsub("\\*.*", "*", peptide_full)  # Remove everything after the first STOP codon
  print(peptide_valid)
  
  peptide_valid = gsub("[*]", "", peptide_valid) # Remove the STOP codon
  
  nmer_i = nmer_i %>% 
    mutate(alt.aa.valid = peptide_valid)
  
  aa_valid = rbind(aa_valid, nmer_i)
}


# -------------------------------------------------------------------------
# Map the n-mers to their original AA sequence 
# -------------------------------------------------------------------------

aa_shared_0201 = NULL
for (i in 1:nrow(candidates_2alg_hlaa0201)){
  print(i)
  
  # Iterate through each peptide candidate
  candidate_i = candidates_2alg_hlaa0201 %>% dplyr::slice(i)
  candidate_i = candidate_i %>% pull(peptide)
  
  # Search for the candidate in the converted amino acid list
  peptide_i = aa_valid[grep(candidate_i, aa_valid$alt.aa.valid),]
  peptide_i$peptide = candidate_i
  
  aa_shared_0201 = rbind(aa_shared_0201, peptide_i)
}


candidates_10602 = inner_join(aa_shared_0201, sj_10602_edit)
candidates_10602_complete = candidates_10602[!duplicated(candidates_10602$junc.id),]
candidates_10602_summary = subset(candidates_10602_complete, select = c(junc.id, symbol, ensg, alt.aa.valid, peptide, n.uniq.map))
# SF10602 Candidates = 12


candidates_10417 = inner_join(aa_shared_0201, sj_10417_edit)
candidates_10417_complete = candidates_10417[!duplicated(candidates_10417$junc.id),]
candidates_10417_summary = subset(candidates_10417_complete, select = c(junc.id, symbol, ensg, alt.aa.valid, peptide, n.uniq.map))
# SF10417 Candidates = 14

shared_sfline_candidates = inner_join(candidates_10602_summary, candidates_10417_summary, by = "junc.id")
shared_sfline_candidates_full = full_join(candidates_10602_summary, candidates_10417_summary, by = "junc.id")
# Shared Candidates = 7
# All found in cell lines = 19

# Change colnames to match
colnames(candidates_10602_summary)[5] = "n_mer"
colnames(candidates_10417_summary)[5] = "n_mer"


# -------------------------------------------------------------------------
# Reassign the HLAthena and MHCFlurry Scores to the Candidates
# -------------------------------------------------------------------------
# Reassign HLAthena scores
candidates_10602_summary = inner_join(candidates_10602_summary, ha_1000, by = "n_mer")
candidates_10602_summary = inner_join(candidates_10602_summary, mf_1000, by = "n_mer")
candidates_10602_summary = candidates_10602_summary[!duplicated(candidates_10602_summary$junc.id),]

# Reassign MHCflurry scores
candidates_10417_summary = inner_join(candidates_10417_summary, ha_1000, by = "n_mer")
candidates_10417_summary = inner_join(candidates_10417_summary, mf_1000, by = "n_mer")
candidates_10417_summary = candidates_10417_summary[!duplicated(candidates_10417_summary$junc.id),]

# Combine both cell lines into a COMBINED file

colnames(candidates_10417_summary)[6] = "SF10417_count"
colnames(candidates_10602_summary)[6] = "SF10602_count"
candidates_all_summary = full_join(candidates_10417_summary, candidates_10602_summary)
candidates_all_edit = candidates_all_summary[, c(1,2,3,4,5,6,11,7,10,8,9)]
candidates_all_edit[is.na(candidates_all_edit)] <- 0     # Change all NAs into 0 
candidates_all_final = NULL
for (i in 1:nrow(candidates_all_edit)){
  print(i)
  row_i = candidates_all_edit %>% dplyr::slice(i)
  sj_i = row_i %>% pull(junc.id)
  
  # Look for duplicates
  sj_filter = candidates_all_edit %>% dplyr::filter(junc.id == sj_i)
  
  # If there are duplicates, then combine them into one row
  if (nrow(sj_filter > 1)){
    row_i[1,6] = sum(sj_filter$SF10417_count)
    row_i[1,7] = sum(sj_filter$SF10602_count)
    candidates_all_final = rbind(candidates_all_final, row_i)
  }
  if (nrow(sj_filter == 1)){candidates_all_final = rbind(candidates_all_final, row_i)}
}

candidates_all_final = candidates_all_final[!duplicated(candidates_all_final$junc.id),]


setwd(directory_out)
filename_sf10417 = "2022_0213_SF10417_Neojunction_Candidates.tsv"
filename_sf10602 = "2022_0213_SF10602_Neojunction_Candidates.tsv"
filename_combine = "2022_0213_Combined_SF10417_SF10602_Neojunction_Candidates.tsv"
write_tsv(candidates_10417_summary, filename_sf10417, na = "NA", col_names = T, quote_escape = "double")
write_tsv(candidates_10602_summary, filename_sf10602, na = "NA", col_names = T, quote_escape = "double")
write_tsv(candidates_all_final, filename_combine, na = "NA", col_names = T, quote_escape = "double")



