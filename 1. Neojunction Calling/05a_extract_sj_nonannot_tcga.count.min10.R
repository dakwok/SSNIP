# Title: "Step 5 - Select for Non-annotated Junctions with a Read Count > 10"
# July 27, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Part A - Purpose: Filter splicing junctions for two things:
#             1. Non-annotated
#             2. Read count > 10
# Pipeline: 1. Filter the folder directory for TCGA patients/samples in the TCGA list generated in Step 1 (Tumor Purity > 0.60)
#           2. From each file in the directory, select for junctions that have a read count > 10
#           3. Combine the all junctions from each file into a comprehensive list of Junctions with Read Count > 10

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggsci)

directory_in = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/05_extract_non_annot_sj_min_count_10_in_tcga_n429/input"
setwd(directory_in)

# Directory for the TCGA GBM SJ RNA-seq data
# Contains RNA seq for every TCGA GBM sample in individual files
directory_in_GBM = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/05_extract_non_annot_sj_min_count_10_in_tcga_n429/input/tcga_gbm_sj_rnaseq"

# Directory for the TCGA LGG SJ RNA-seq data
# Contains RNA seq for every TCGA LGG sample in individual files
directory_in_LGG = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/05_extract_non_annot_sj_min_count_10_in_tcga_n429/input/tcga_lgg_sj_rnaseq"

#  Load Files --------------------------------------------------------------

# Splicing Junction Coordinates (Raw)
filename_sj  = "sjdbList.fromGTF.out.tab"
dataframe_sj = read_tsv(filename_sj, na = c("", "NA"), col_names = F, col_types = cols(X1 = col_character()))

# Patient List (From Step #1) Samples with a Tumor Purity > 60%
filename_tcga = "TCGA_Patient_List_Post_TumorPurity_Filter_20200604"
dataframe_tcga = read_tsv(filename_tcga, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))

# GTF File - (From Step #3) Protein-Coding Transcripts with Median TPM > 10
filename_gtf  = "GTF_ProteinCoding_FilterTPM10_20200717"
dataframe_gtf = read_tsv(filename_gtf, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))

dataframe_tcga %>% head(3) %>% print() ;


###########################################################################
#  Step 1: Rename the Columns of the SJ Coordinates Table -----------------
###########################################################################

dataframe_sj = dataframe_sj %>% 
  dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4) %>% 
  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) ; 


###########################################################################
#  Step 2: Generate the List of Non-Annotated SJ --------------------------
###########################################################################
# Goal: This step generates a list of splicing junctions that:
#         1. Are non-annotated
#         2. Have a read count >= 10

start.time = proc.time()

# Create an empty list to hold all of the "non-annotated" junction counts
list_sj = list()

for (i in 1:2) {
  # 0. Work through both the TCGA GBM and LGG folders 
  # FUNCTION: list.file() creates a character vector of the names of files in the current directory
  setwd(c(directory_in_GBM, directory_in_LGG)[i])
  
  # 1. Select for TCGA samples in each directory that are found in the tumor purity patient list (generated from Step #1)
  list_rnaseq.files = list.files()[is.element(substr(list.files(), 1, 12), dataframe_tcga$case)]
  
  for (j in 1:length(list_rnaseq.files)) {
    filename_j = list_rnaseq.files[j]
    dataframe_j = read_tsv(filename_j, na = c("", "NA"), col_names = F, col_type = cols(X1 = col_character()))
    
    # 2. Remove splicing junctions in each file that are not found
    #    in the raw reference splicing junction file: dataframe_sj
    # Make a new entry within list_sj for each iteration of files
    list_sj[[(length(list_sj) + 1)]] = dataframe_j %>%          # Starts with list_sj[[1]] and onwards
      dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
      mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
      dplyr::filter(strand !=	"undefined") %>%
      dplyr::filter(chr != "M") %>% 
      mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
      anti_join(dataframe_sj, by = "junc.id")
    
    # Extract junctions with 10 or more supportive read counts
    dataframe_j.read10 = list_sj[[length(list_sj)]] %>%
      dplyr::filter(n.uniq.map >= 10) %>%
      dplyr::select(chr, int.start, int.end, strand, junc.id)
    
    # Combine the junctions found in the current file with the junctions found in previous iterations
    if (i == 1 & j == 1) {
      dataframe_combined = dataframe_j.read10}
    else {
      dataframe_combined = dataframe_combined %>% 
        bind_rows(dataframe_j.read10) %>% 
        distinct(junc.id, .keep_all = TRUE)}
    
    # Progress Bar
    print(i)
    print((j/length(list_rnaseq.files))*100)
    }
  }
}

# 3. Edit and finalize the combined SJ dataframe
dataframe_combined = dataframe_combined %>% 
  arrange(chr, int.start, int.end)

RUNTIME = proc.time() - start.time 
RUNTIME %>% print()
#     user   system  elapsed 
#  476.438    7.180 2499.749

# List of 158,388 splicing junctions


###########################################################################
#  Step 3: Output Data ----------------------------------------------------
###########################################################################

directory_output = "/Users/darwinkwok/Documents/UCSF\ 2020\ Okada-Costello\ Lab/Project\ A\ -\ IDH1\ Alternative\ Splicing-Associated\ Neojunctions/05_extract_non_annot_sj_min_count_10_in_tcga_n429/output"
setwd(directory_output)

filename_output = "SJ_List_NonAnnotated_Candidates_20200627"
write_tsv(dataframe_combined, filename_output, na = "NA", col_names = T, quote_escape = "double")


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# Part B - Purpose: Filter splicing junctions for two things
# Pipeline: 1. 

#  Load Files --------------------------------------------------------------

directory_output = "/Users/darwinkwok/Documents/UCSF\ 2020\ Okada-Costello\ Lab/Project\ A\ -\ IDH1\ Alternative\ Splicing-Associated\ Neojunctions/05_extract_non_annot_sj_min_count_10_in_tcga_n429/output"
setwd(directory_output)

filename_sj.combined = "SJ_List_NonAnnotated_Candidates_20200627"
dataframe_sj.combined = read_tsv(filename_sj.combined, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))


###########################################################################
#  Step 4: Filter for Splicing Junctions in Protein Coding Transcripts ----
###########################################################################
# The previous list of splicing junctions generated (in Step 2) are filtered for the following:
#         1. Are non-annotated
#         2. Have a read count >= 10
# The goal of this step is to further filter for splicing junctions for:
#         3. Within protein-coding transcript regions

start.time = proc.time()

for (i in 1:nrow(dataframe_sj.combined)) {
  
  # Progress Bar
  print((i/nrow(dataframe_sj.combined))*100)
  
  # Iterate through each row to work with them individually
  dataframe_sj.combined.i = dataframe_sj.combined %>% 
    dplyr::slice(i)
  
  CHR    = dataframe_sj.combined.i %>% pull(chr)
  STRAND = dataframe_sj.combined.i %>% pull(strand)
  START  = dataframe_sj.combined.i %>% pull(int.start)
  END    = dataframe_sj.combined.i %>% pull(int.end)
  
  # Generate a new variable RETAIN that is 1 for retain or 0 for do NOT retain
  # This is calculated by determining whether junction(i) is within any transcript
  RETAIN = dataframe_gtf %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>%
    dplyr::filter(start < START & END < end) %>%
    head(1) %>%
    nrow() ;
  
  # Create a dataframe that adds a splicing junction if it is found in the protein-coding
  # GTF list (RETAIN == 1)
  if(i == 1) {
    dataframe_sj.filtered = dataframe_sj.combined.i}
  else{
    if(RETAIN == 1) {
      dataframe_sj.filtered = dataframe_sj.filtered %>% 
        bind_rows(dataframe_sj.combined.i)
    }
  }
}

RUNTIME = proc.time() - start.time ; # sec
RUNTIME %>% print() ;

#     user   system  elapsed 
# 2130.026   35.707 2115.852

# List of 79,953 splicing junctions


###########################################################################
#  Step 5: Filter for Splicing Junctions in Protein Coding Transcripts ----
###########################################################################
# From the list of splicing junctions generated in Step 4, this step generates the "count" 
# by overlaying the list of splicing junctions generated in Step 4 (79,953 junctions) with the 
# TCGA dataframes for LGG and GBMs

# 1. Initiate the count dataframe by taking only the junction coordinates 
dataframe_sj.count = dataframe_sj.filtered %>% 
  dplyr::select(junc.id)

# 2. Filter for the count data using the junction IDs and cross analyzing it with the TCGA count tables
start.time = proc.time() ; 

for(j in 1:2){
  setwd(c(directory_in_GBM, directory_in_LGG)[j]) ; # gbm, lgg ; 
  lf = list.files()[is.element(substr(list.files(), 1, 12), dataframe_tcga$case)] ; 
  
  for(i in 1:length(lf)){
    in.f = lf[i] ; 
    sj.i = read_tsv(in.f, na = c("", "NA"), col_names = F, 
                    col_type = cols(X1 = col_character())
    ) %>% 
      dplyr::filter(X1 == CHR) ;
    
    # 3. Retain junctions that are uniquely expressed (e.g. splicing junction only occurs in one event
    #    throughout the whole genome)
    dataframe_sj.count = dataframe_sj.count %>% 
      left_join(sj.i %>% 
                  dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
                  mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
                  dplyr::filter(strand !=	"undefined") %>%
                  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
                  dplyr::select(junc.id, n.uniq.map), by = "junc.id") %>% 
      mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map)); 
    
    colnames(dataframe_sj.count)[ncol(dataframe_sj.count)] = substr(lf[i], 1, 12) ; 
    
    round((ncol(dataframe_sj.count) - 1)/429, 4) %>% print() ; # to see the progress
  }
}

RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print() ; 

#    user  system elapsed 
# 107.855   6.372 118.933


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
for (i in 1:6) {
  print(i)
  CASES = list(
    dataframe_tcga$case,
    dataframe_tcga$case[dataframe_tcga$study == "TCGA-GBM"],
    dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDHwt"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDH-A"],
    dataframe_tcga$case[dataframe_tcga$cls == "IDH-O"])[[i]]
  
  # 3. Warp dataframe to be based on judgement in each case
  dataframe_sj.count.case = dataframe_sj.count.warp %>% 
    dplyr::filter(is.element(sample.id, CASES)) %>% 
    spread(sample.id, judge)
  
  # 4. Generate a new table for positive sample rate (based on judgement) for each disease group
  dataframe_sj.psr.i = dataframe_sj.count.case %>% 
    mutate(sum_positive.judgement = apply(dataframe_sj.count.case %>% 
                                            dplyr::select(- junc.id), 1, sum)) %>% 
    mutate(positive.sample.rate = round(sum_positive.judgement/length(CASES), 3)) %>% 
    dplyr::select(junc.id, positive.sample.rate)
  
  colnames(dataframe_sj.psr.i)[2] = paste0(c("all", "TCGA-GBM", "TCGA-LGG", "IDHwt", "IDH-A", "IDH-O")[i], "_n", length(CASES))
  
  # 5. Combine the PSR columns for each disease group to the final dataframe_sj.psr
  if (i == 1) {
    dataframe_sj.psr = dataframe_sj.psr.i}
  else {dataframe_sj.psr = dataframe_sj.psr %>% 
    bind_cols(dataframe_sj.psr.i %>% 
                dplyr::select(- junc.id))}
}

# 6. Determine the junctions to retain by determining whether the maximum PSR in any disease group is greater than 10%
dataframe_sj.psr = dataframe_sj.psr %>% 
  mutate(max = apply(dataframe_sj.psr[, -1], 1, max)) %>%    # Make a new column "max" that is the greatest PSR value from all disease groups for the junction of interest
  mutate(retain = ifelse(max >= 0.1, 1, 0))                  # Make a new column "retain" that is "1" if the maximum PSR is greater than 10% and "0" if less


dataframe_sj.psr.retained = dataframe_sj.psr %>% 
  dplyr::filter(retain == 1) %>% 
  dplyr::select( - retain)

dataframe_sj.psr %>% pull(retain) %>% table() %>% print()
###########################################################################
#  Step 5: Output Data ----------------------------------------------------
###########################################################################

directory_output = "/Users/darwinkwok/Documents/UCSF\ 2020\ Okada-Costello\ Lab/Project\ A\ -\ IDH1\ Alternative\ Splicing-Associated\ Neojunctions/05_extract_non_annot_sj_min_count_10_in_tcga_n429/output"
setwd(directory_output)

filename_output.sj.filtered = "SJ_List_NonAnnotated_Candidates_&_Protein_Coding_20200630"
write_tsv(dataframe_sj.filtered, filename_output.sj.filtered, na = "NA", col_names = T, quote_escape = "double")

dataframe_sj.filtered = read_tsv(filename_output.sj.filtered, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))



