# Title: "Step 9: Count, Depth, Frequency, and Judgement Tables for GTEx"
# August 26, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Note: SAME script as Step 7
# Purpose: Generate and prepare the following tables for the GTEx data:
#            1. Count table
#            1. Depth table
#            1. Frequency table
#            1. Judgement table

###########################################################################
#  Step 00: Connect to the TIPCC ------------------------------------------
###########################################################################
# This code will need to run through 9166 GTEx datasets which are over 19 Tb in size.
# Therefore, the TIPCC is necessary to efficiently process everything 

#ssh dkwok@cclc01.som.ucsf.edu
# <pw required>
#ssh n6
#module load CBC r
#R

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(tidyverse)
library(readxl)
library(ggsci)

# Establish Directories ---------------------------------------------------

directory_external   = "{PATH_TO_GTEX_ANNOTATION}"
directory_07  = "{PATH_TO_STEP07_OUTPUT}"
directory_08  = "{PATH_TO_STEP08_OUTPUT}"
directory_out = "{PATH_TO_STEP09_OUTPUT}"

# Load Files --------------------------------------------------------------
# From External: GTEx Information Dataframe
setwd(directory_external)
filename_gtex.info = "07102020_gtex_id_annotation_9166.tsv"
dataframe_gtex.info = read_tsv(filename_gtex.info, na = c("", "NA"))

# head(dataframe_gtex.info) %>% print()
# sample.id  sample.name          body.site      cls      
# <chr>      <chr>                <chr>          <chr>    
# 1 SRR1404062 GTEX-12696-1726-SM-… Stomach        Stomach  
# 2 SRR1404084 GTEX-139T8-0326-SM-… Nerve - Tibial Nerve    
# 3 SRR1404105 GTEX-1399R-0626-SM-… Artery - Aorta Blood Ve…
# 4 SRR1404123 GTEX-Y9LG-0526-SM-4… Lung           Lung     
# 5 SRR1404147 GTEX-YEC3-1426-SM-5… Stomach        Stomach  
# 6 SRR1404168 GTEX-YFC4-3026-SM-5… Brain - Cereb… Brain  

# From Step 7: Overlap Tables for General (1) and IDH1 Subclasses (2)
#setwd(directory_out.step07)
#filename_overlap.table.1 = "P1_Overlap_Table_Retained_and_Passed_Junctions_20200818"
#dataframe_overlap.table.1 = read_tsv(filename_overlap.table.1, na = c("", "NA"))

# NJs that passed the PSR > 0.1 filter in TCGA
setwd(directory_07)
filename_nj_passed_tcga_psr = "P1_PSR_20200818.tsv"
dataframe_passed_tcga_junc = read_tsv(filename_nj_passed_tcga_psr , col_names = T)

# From Step 8: GTEx Overlap Table
setwd(directory_08)
filename_overlap = "2023_0623_gtex_sj_overlap_table.tsv"
dataframe_overlap.table.1 = read_tsv(filename_overlap, na = c("", "NA"))
dataframe_overlap.table.1 = subset(inner_join(dataframe_overlap.table.1, dataframe_passed_tcga_junc), select = c("junc.id", "junc.id.overlap"))

# From Step 8: GTEx Count Tables for General (1) and IDH1 Subclasses (2)
setwd(directory_08)
filename_gtex.count.1 = "2023_0623_gtex_sj_count_table_sj.out.tab_retained.tsv"
dataframe_gtex.count.1 = read_tsv(filename_gtex.count.1, na = c("", "NA"))
# rows(dataframe_gtex.count.1) == 12369
# rows(dataframe_gtex.count.2) == 12617



###########################################################################
#  Step 2: Obtain the Count, Depth, Frequency, and Judgement Values -------
###########################################################################
# Make lists for dataframes for each of the project goals: General (1) and IDH1 Subclasses (2)
list_overlap.table = list(dataframe_overlap.table.1)
list_sj.count      = list(dataframe_gtex.count.1)

list_count = list()
list_depth = list()
list_freq  = list()
list_judge = list()
 
for (h in 1:1){
  dataframe_overlap.table.h = list_overlap.table[[h]]
  dataframe_sj.count.h      = list_sj.count[[h]]
  for (i in 1:nrow(dataframe_overlap.table.h)){
    # Progress Bar
    print(h)
    print(i/nrow(dataframe_overlap.table.h)*100)
    
    # Select out the junction of interest (JUNC.ID) and the overlapping junctions (JUNC.ID.OVERLAP)
    JUNC.ID = dataframe_overlap.table.h %>% 
      dplyr::slice(i) %>% 
      pull(junc.id)
    JUNC.ID.OVERLAP = dataframe_overlap.table.h %>% 
      dplyr::slice(i) %>% 
      pull(junc.id.overlap)
    
    
    # If there are NO overlapping junctions in the current junction of interest, then 
    if (is.na(JUNC.ID.OVERLAP)){
      res.i = dataframe_sj.count.h %>%                            # Within the GTEx Count Table:
        dplyr::filter(junc.id == JUNC.ID) %>%                     # Filter for only the current junction of interest in the iteration
        gather(case, value, 2:ncol(dataframe_sj.count.h)) %>%     # Warp dataframe to 3 columns: 1) junc.id, 2) case, 3) count
        spread(junc.id, value) %>%                                # Remove the "junc.id" column
        dplyr::rename(count = JUNC.ID) %>%                        # Rename the count column as "count"
        mutate(count = ifelse(is.na(count), 0, count)) %>%        # Within the "count" column, convert all NA values to 0
        mutate(depth = count) %>%                                 # Make "depth" column. Since depth is the sum of the counts of the junction and its overlapping junction, in the case in which there are no overlapping junctions (JUNC.ID.O == NA), then Depth = count + 0
        mutate(freq = round(count / depth, 4)) %>%                # Make "freq" column. Frequency is defined as the ratio of the junction count and depth (Frequency = count/(count + count.o))
        mutate(freq = ifelse(is.na(freq), 0, freq)) %>%           # Within the "freq" column, convert all NA values to 0
        mutate(judge = ifelse(count >= 2 & depth >= 10 
                              & freq >= 0.01, 1, 0)) %>%          # Make "judge" column. Judge for GTEx is counted as having count >= 2 & depth >= 10 & freq >= 0.01
        mutate(junc.id = JUNC.ID) %>% 
        dplyr::select(junc.id, case, count, depth, freq, judge)}
   
     # If there is an overlapping junction, bind the count data for both the junction of interest (JUNC.ID) and (JUNC.ID.OVERLAP)
    if (!is.na(JUNC.ID.OVERLAP)){
      res.i = bind_rows(dataframe_sj.count.h %>% dplyr::filter(junc.id == JUNC.ID),              # Filter for the counts of the junction of interest 
                        dataframe_sj.count.h %>% dplyr::filter(junc.id == JUNC.ID.OVERLAP)) %>%  # Filter for the counts of the overlapping junctions, bind all the rows.
        gather(case, value, 2:ncol(dataframe_sj.count.h)) %>%                 # Warp dataframe to 3 columns: 1) junc.id, 2) case, 3) count
        spread(junc.id, value) %>%                                # Remove the "junc.id" column
        dplyr::rename(count = JUNC.ID, 
                      count.o = JUNC.ID.OVERLAP) %>%              # Rename the junction count column as "count" and overlapping junction count column as "count.o"
        mutate(count = ifelse(is.na(count), 0, count)) %>%        # Within the "count" column, convert all NA values to 0
        mutate(count.o = ifelse(is.na(count.o), 0, count.o)) %>%  # Within the "count.o" column, convert all NA values to 0
        mutate(depth = count + count.o) %>%                       # Make "depth" column. Depth is defined as the Sum of the counts of both junctions (Depth = count + count.o)
        mutate(freq = round(count / depth, 4)) %>%                # 
        mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
        mutate(judge = ifelse(count >= 2 & depth >= 10 & freq >= 0.01, 1, 0)) %>% 
        mutate(junc.id = JUNC.ID) %>% 
        dplyr::select(junc.id, case, count, depth, freq, judge)}
    count.i = res.i %>% 
      dplyr::select(junc.id, case, count) %>% 
      spread(case, count)
    depth.i = res.i %>% 
      dplyr::select(junc.id, case, depth) %>% 
      spread(case, depth)
    freq.i = res.i %>% 
      dplyr::select(junc.id, case, freq) %>% 
      spread(case, freq)
    judge.i = res.i %>% 
      dplyr::select(junc.id, case, judge) %>% 
      spread(case, judge)
    
    # Generate the Count, Depth, Freq, and Judge Tables
    if(i == 1){
      count.h = count.i; 
      depth.h = depth.i; 
      freq.h  = freq.i; 
      judge.h = judge.i}
    # Each iteration (junction) is a new row, and therefore bind the rows
    if(i != 1){
      count.h = bind_rows(count.h, count.i); 
      depth.h = bind_rows(depth.h, depth.i);
      freq.h  = bind_rows(freq.h, freq.i);
      judge.h = bind_rows(judge.h, judge.i);}
  }
  list_count[[h]] = count.h
  list_depth[[h]] = depth.h
  list_freq[[h]]  = freq.h
  list_judge[[h]] = judge.h
}

setwd(directory_out)
filename_count   = "2023_0704_count_ALL_gtex.tsv"
filename_depth   = "2023_0704_depth_ALL_gtex.tsv"
filename_freq    = "2023_0704_freq_ALL_gtex.tsv"
filename_judge   = "2023_0704_judge_ALL_gtex.tsv"
write_tsv(list_count[[1]], filename_count, na = "NA", col_names = T, quote_escape = "double")
write_tsv(list_depth[[1]], filename_depth, na = "NA", col_names = T, quote_escape = "double")
write_tsv(list_freq[[1]], filename_freq, na = "NA", col_names = T, quote_escape = "double")
write_tsv(list_judge[[1]], filename_judge, na = "NA", col_names = T, quote_escape = "double")


###########################################################################
#  Step 3: Prepare 3 Subgroups --------------------------------------------
###########################################################################
setwd(directory_out)
count = read_tsv("2023_0704_count_ALL_gtex.tsv", col_names = T)
freq  = read_tsv("2023_0704_freq_ALL_gtex.tsv", col_names = T)
depth = read_tsv("2023_0704_depth_ALL_gtex.tsv", col_names = T)
judge = read_tsv("2023_0704_judge_ALL_gtex.tsv", col_names = T)

CASES = list(dataframe_gtex.info$sample.id, 
             dataframe_gtex.info$sample.id[dataframe_gtex.info$cls != "Testis"], 
             dataframe_gtex.info$sample.id[dataframe_gtex.info$cls == "Brain"])

names(CASES) = c("GTEx_all", "GTEx_all_but_ts", "Brain")

list_psr  = list()
list_pass = list()

colnames(judge)[-1] <- sub(".*/(SRR[0-9]+).*", "\\1", colnames(judge)[-1])

# Positive Sample Rate
for (i in 1:3) {
  print(i)
  judge.i = judge %>% 
    dplyr::select(junc.id, colnames(judge)[is.element(colnames(judge), CASES[[i]])])
  
  judge.i = judge.i %>% 
    mutate(pos = apply(judge.i[, -1], 1, sum)) %>% 
    mutate(psr = round(pos / length(CASES[[i]]), 4)) %>% 
    dplyr::select(junc.id, psr)
  
  colnames(judge.i)[2] = paste0(names(CASES)[[i]], "_n", length(CASES[[i]]))
  
  if(i == 1){
    judge.2 = judge.i}
  if(i != 1){
    judge.2 = judge.2 %>% 
      full_join(judge.i, by = "junc.id")}
}

judge.2 = judge.2 %>% 
  mutate(min = apply(judge.2[, 2:3], 1, min)) # retain if the psr doesn't exceed 0.01 (psr < 0.01) in GTEx all OR GTEx-all but Testis. not take into acount Brain here.

judge.pass = judge.2 %>% 
  dplyr::filter(min < 0.01) ;

h = 1

list_psr[[h]]  = judge.2
list_pass[[h]] = judge.pass



count.pass = count %>% semi_join(list_pass[[1]], by = "junc.id")
depth.pass = depth %>% semi_join(list_pass[[1]], by = "junc.id")
freq.pass  = freq  %>% semi_join(list_pass[[1]], by = "junc.id")
judge.pass = judge %>% semi_join(list_pass[[1]], by = "junc.id")
overlap.table.pass = dataframe_overlap.table.1 %>% semi_join(list_pass[[1]], by = "junc.id") 



setwd(directory_out)
filename_count   = "2023_0704_count_retained_and_passed_gtex.tsv"
filename_depth   = "2023_0704_depth_retained_and_passed_gtex.tsv"
filename_freq    = "2023_0704_freq_retained_and_passed_gtex.tsv"
filename_judge   = "2023_0704_judge_retained_and_passed_gtex.tsv"
filename_overlap = "2023_0704_overlaptable_retained_and_passed_gtex.tsv"
write_tsv(count.pass, filename_count, na = "NA", col_names = T, quote_escape = "double")
write_tsv(depth.pass, filename_depth, na = "NA", col_names = T, quote_escape = "double")
write_tsv(freq.pass, filename_freq, na = "NA", col_names = T, quote_escape = "double")
write_tsv(judge.pass, filename_judge, na = "NA", col_names = T, quote_escape = "double")
write_tsv(overlap.table.pass, filename_overlap, na = "NA", col_names = T, quote_escape = "double")

# Positive Sample Rate (PSR)
filename_psr     = "2023_0704_psr_ALL_gtex.tsv"
write_tsv(list_psr[[1]], filename_psr, na = "NA", col_names = T, quote_escape = "double")

filename_pass     = "2023_0704_psr_retained_and_passed_gtex.tsv"
write_tsv(list_pass[[1]], filename_pass, na = "NA", col_names = T, quote_escape = "double")

