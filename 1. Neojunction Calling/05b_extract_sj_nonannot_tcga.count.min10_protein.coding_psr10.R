# Title: "Step 5b - Further Filter SJ for Unique Protein-Coding Transcripts with GTF and a PSR > 10"
# July 31, 2020 | Darwin Kwok M.S. | University of California San Francisco
# Last updated: 10-23-2020

# Part A - Purpose: Filter splicing junctions for:
#             1. Non-annotated
#             2. Read count > 10
# Part B - Purpose: Filter splicing junctions for:
#             3. Protein Coding
#             4. Uniquely Expressed
#             5. PSR > 10 in at least one Disease Group

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE)) ;

library(tidyverse) ;
library(readxl) ;
library(ggsci) ;

directory_05a_in  = "PATH_TO_INPUT"
directory_05a_out = "PATH_TO_OUTPUT"
directory_in_GBM  = "PATH_TO_GBM_INPUT" ; # sj.out.tab of TCGA-GBM
directory_in_LGG  = "PATH_TO_LGG_INPUT"; # sj.out.tab of TCG-LGG


setwd(directory_05a_out) ; 
filename_dataframe_sj.combined = "SJ_List_NonAnnotated_Candidates_20200627.tsv" ; 
dataframe_sj.combined = read_tsv(filename_dataframe_sj.combined, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))

setwd(directory_05a_in) ; 
filename_gtf = "GTF_ProteinCoding_FilterTPM10_20200717" ; 
dataframe_gtf = read_tsv(filename_gtf, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character())) ;

filename_tcga = "TCGA_Patient_List_Post_TumorPurity_Filter_20200604" ; # TCGA high-tumor purity, n = 429
dataframe_tcga = read_tsv(filename_tcga, na = c("", "NA"), col_names = T) ;



###########################################################################
#  Step 1: Filter for SJ in Protein Coding Transcripts --------------------
###########################################################################
# The previous list of splicing junctions generated (in Step 2) are filtered for the following:
#         1. Are non-annotated
#         2. Have a read count >= 10
# The goal of this step is to further filter for splicing junctions for:
#         3. Within protein-coding transcript regions

start.time = proc.time() ;

for(i in 1:nrow(dataframe_sj.combined)){
  # Progress Bar
  print((i/nrow(dataframe_sj.combined))*100)
  
  # Iterate through each row to work with them individually
  dataframe_sj.combined.i = dataframe_sj.combined %>%
    dplyr::slice(i) ;
  
  CHR    = dataframe_sj.combined.i %>% pull(chr) ;
  STRAND = dataframe_sj.combined.i %>% pull(strand) ;
  START  = dataframe_sj.combined.i %>% pull(int.start) ;
  END    = dataframe_sj.combined.i %>% pull(int.end) ;
  
  # 1. Generate a new variable RETAIN that is 1 for retain or 0 for do NOT retain
  # This is calculated by determining whether junction(i) is within any transcript
  RETAIN = dataframe_gtf %>%
    dplyr::filter(chr == CHR & strand == STRAND) %>%
    dplyr::filter(start < START & END < end) %>%
    head(1) %>%
    nrow() ;
  
  # 2. Create a dataframe that adds a splicing junction if it is found in the protein-coding
  # GTF list (RETAIN == 1)
  dataframe_sj.combined.i = dataframe_sj.combined.i %>%
    mutate(retain = RETAIN) ;
  
  if(i == 1){
    dataframe_sj.combined.edited = dataframe_sj.combined.i ;
  }else{
    dataframe_sj.combined.edited = dataframe_sj.combined.edited %>%
      bind_rows(dataframe_sj.combined.i) ;
  }
}

RUNTIME = proc.time() - start.time ; # sec
RUNTIME %>% print() ;

#     user   system  elapsed 
# 2130.026   35.707 2115.852

# 3. Retain splicing junctions with RETAIN == 1
dataframe_sj.combined.retain = dataframe_sj.combined.edited %>%
  dplyr::filter(retain == 1) %>%
  dplyr::select(junc.id, chr, strand, int.start, int.end) ;

# List of 79,953 splicing junctions



###########################################################################
#  Step 2: Output SJ Data (Filtered for TPM and PC) -----------------------
###########################################################################

directory_out = "PATH_TO_OUTPUT"
setwd(directory_out)

filename_output.sj.filtered = "SJ_List_NonAnnotated_Candidates_&_Protein_Coding_20200630"
write_tsv(dataframe_sj.combined.retain, filename_output.sj.filtered, na = "NA", col_names = T, quote_escape = "double")



###########################################################################
#  Step 3: Generate SJ Count Table from sj.out.tab with the SJ List -------
###########################################################################
# From the list of splicing junctions generated in Step 4, this step generates the "count" 
# by overlaying the list of splicing junctions generated in Step 4 (79,953 junctions) with the 
# TCGA dataframes for LGG and GBMs

setwd(directory_out)
dataframe_sj.combined.retain =  read_tsv("SJ_List_NonAnnotated_Candidates_&_Protein_Coding_20200630.tsv", na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))


# 1. Initiate the count dataframe by taking only the junction coordinates 
dataframe_sj.count = dataframe_sj.combined.retain %>% 
  dplyr::select(junc.id) ; 


start.time = proc.time() ; 

# 2. Filter for the count data using the junction IDs and cross analyzing it with the TCGA count tables
for(i in 1:2){
  setwd(c(directory_in_GBM, directory_in_LGG)[i]) ; # Set i for the 1. GBM directory and 2. LGG director
  list_files = list.files()[is.element(substr(list.files(), 1, 12), dataframe_tcga$case)] ; 
  
  for(j in 1:length(list_files)){
    in.f = list_files[j] ; 
    dataframe_sj.combined.j = read_tsv(in.f, na = c("", "NA"), col_names = F, col_type = cols(X1 = col_character()))
    
    # 3. Retain junctions that are uniquely expressed (e.g. splicing junction only occurs in one event
    #    throughout the whole genome)
    dataframe_sj.count = dataframe_sj.count %>% 
      left_join(dataframe_sj.combined.j %>% 
                  dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
                  mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
                  dplyr::filter(strand !=	"undefined") %>%
                  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
                  dplyr::select(junc.id, n.uniq.map), by = "junc.id") %>% 
      mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map)); 
    
    colnames(dataframe_sj.count)[ncol(dataframe_sj.count)] = substr(list_files[j], 1, 12) ; 
    
    round((ncol(dataframe_sj.count) - 1)/429, 4) %>% print() ; # to see the progress
  }
}

RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print() ; 

# check -------------------------------------------------------------------

dataframe_sj.count %>% dim() %>% print() ; 


dataframe_sj.count %>% dplyr::select(1:5) %>% head(10) %>% print() ; 



###########################################################################
#  Step 6: Filter for Junctions with Read Count > 10 per Disease Group ----
###########################################################################
# From the previous step, I will filter for junctions that have a minimum Read Count of 10 in 
# at least of 10% of each disease groups

# 1. Manipulate the dataframe with FUNCTION gather() to have all the counts for each junction+sample in one column
dataframe_sj.count.warp = dataframe_sj.count %>% 
  gather("sample.id", "count", 2:ncol(dataframe_sj.count)) %>% 
  # Create a new column "judge" for that reads "1" whenever count > 10, and "0" if count < 10
  mutate(judge = ifelse(count < 10, 0, 1)) %>% 
  # Remove the count column
  dplyr::select( - count) 


# 2. Filter for functions with a minimum read count >= 10 in at least 10% within each disease group
for(i in 1:6){
  print(i) ;
  CASES = list(
    dataframe_tcga$case, 
    dataframe_tcga$case[dataframe_tcga$study == "TCGA-GBM"], 
    dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG"], 
    dataframe_tcga$case[dataframe_tcga$cls == "IDHwt"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDH-A"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDH-O"]
  )[[i]] ;
  
  # 3. Warp dataframe to be based on judgement in each case
  dataframe_sj.count.i = dataframe_sj.count.warp %>% 
    dplyr::filter(is.element(sample.id, CASES)) %>% 
    spread(sample.id, judge) 
  
  # 4. Generate a new table for positive sample rate (based on judgement) for each disease group
  dataframe_sj.psr.i = dataframe_sj.count.i %>% 
    mutate(n.pos = apply(dataframe_sj.count.i %>% dplyr::select(- junc.id), 1, sum)) %>% 
    mutate(pos.sample.rate = round(n.pos / length(CASES), 3)) %>% 
    dplyr::select(junc.id, pos.sample.rate) ; 
  
  colnames(dataframe_sj.psr.i)[2] = paste0(c("all", "TCGA-GBM", "TCGA-LGG", "IDHwt", "IDH-A", "IDH-O")[i], "_n", length(CASES)) ; 
  
  # 5. Combine the PSR columns for each disease group to the final dataframe_sj.psr
  if(i == 1){
    dataframe_sj.psr = dataframe_sj.psr.i ;
  }else{
    dataframe_sj.psr = dataframe_sj.psr %>% 
      bind_cols(dataframe_sj.psr.i %>% dplyr::select(- junc.id)) ; 
  }
}

# 6. Determine the junctions to retain by determining whether the maximum PSR in any disease group is greater than 10%
# A. Generate the Junction PSR List
dataframe_sj.psr = dataframe_sj.psr %>% 
  mutate(max = apply(dataframe_sj.psr[, -1], 1, max)) %>% 
  mutate(retain = ifelse(max >= 0.1, 1, 0))

dataframe_sj.psr.retain = dataframe_sj.psr %>% 
  dplyr::filter(retain == 1) %>% 
  dplyr::select( - retain)

# B. Generate the Junction Count List
dataframe_sj.count.retain = dataframe_sj.count %>% 
  semi_join(dataframe_sj.psr.retain, by = "junc.id") ; 


###########################################################################
#  Step 5: Output Data ----------------------------------------------------
###########################################################################

setwd(directory_out) ; 
filename_sj.retained.psr = "SJ_PSR_NonAnnotated_Candidates_&_Protein_Coding_20200630" 
write_tsv(dataframe_sj.psr.retain, filename_sj.retained.psr, na = "NA", col_names = T, quote_escape = "double") ;

filename_sj.retained.count = "SJ_Count_NonAnnotated_Candidates_&_Protein_Coding_20200630" 
write_tsv(dataframe_sj.count.retain, filename_sj.retained.count, na = "NA", col_names = T, quote_escape = "double") ;


# si ----------------------------------------------------------------------

sessionInfo() ;

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# 
# Matrix products: default
# BLAS:   /home/shared/cbc/software_cbc/R/R-3.6.3/lib64/R/lib/libRblas.so
# LAPACK: /home/shared/cbc/software_cbc/R/R-3.6.3/lib64/R/lib/libRlapack.so
# 
# locale:
# [1] C
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#  [1] ggsci_2.9       readxl_1.3.1    forcats_0.5.0   stringr_1.4.0
#  [5] dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.3
#  [9] tibble_3.0.1    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.4.6     cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3
#  [5] dbplyr_1.4.3     tools_3.6.3      jsonlite_1.6.1   lubridate_1.7.8
#  [9] lifecycle_0.2.0  nlme_3.1-144     gtable_0.3.0     lattice_0.20-38
# [13] pkgconfig_2.0.3  rlang_0.4.6      reprex_0.3.0     cli_2.0.2
# [17] rstudioapi_0.11  DBI_1.1.0        haven_2.2.0      withr_2.2.0
# [21] xml2_1.3.2       httr_1.4.1       fs_1.4.1         generics_0.0.2
# [25] vctrs_0.3.0      hms_0.5.3        grid_3.6.3       tidyselect_1.1.0
# [29] glue_1.4.0       R6_2.4.1         fansi_0.4.1      modelr_0.1.7
# [33] magrittr_1.5     backports_1.1.6  scales_1.1.1     ellipsis_0.3.0
# [37] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 utf8_1.1.4
# [41] stringi_1.4.6    munsell_0.5.0    broom_0.5.6      crayon_1.3.4


# end ---------------------------------------------------------------------

