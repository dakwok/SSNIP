# Title: "Step 7: Count, Depth, Frequency, and Judgement Tables for Retained Splicing Junctions"
# August 7, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Generate and prepare the following tables for the retained SPLICING JUNCTIONS data:
#            1. Count table
#            1. Depth table
#            1. Frequency table
#            1. Judgement table

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE)) ;

library(tidyverse) ;
library(readxl) ;
library(ggsci) ;
library(data.table)

# Establish Directories ---------------------------------------------------

directory_out.step01 = "PATH_TO_STEP01_OUTPUT"
directory_out.step06 = "PATH_TO_STEP06_OUTPUT"

directory_out = "PATH_TO_STEP07_OUTPUT"

# Load Files --------------------------------------------------------------

# From Step 1: TCGA Patient/Sample List (Tumor Purity > 0.60)
setwd(directory_out.step01)
filename_tcga = "TCGA_Patient_List_Post_TumorPurity_Filter_20200604.tsv"
dataframe_tcga = read_tsv(filename_tcga, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 

# From Step 6: SJ Overlap Table
setwd(directory_out.step06)
filename_overlap.table = "SJ_Overlap_Table_20200818.tsv"
dataframe_overlap.table = read_tsv(filename_overlap.table, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 

# From Step 6: SJ Count Table
filename_count.table = "SJ.out.tab_Retained_Junctions_Count_Table_20200818"
dataframe_count.table = read_tsv(filename_count.table, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 



###########################################################################
#  Step 1: Obtain the Count | Depth | Freq | Judge Tables -----------------
###########################################################################

start.time = proc.time()

for (i in 1:nrow(dataframe_overlap.table)){
  
  # Status Bar
  print(i/nrow(dataframe_overlap.table))
  
  # Define the junction ID and overlapping junction ID for every iteration
  JUNC.ID = dataframe_overlap.table %>% 
    dplyr::slice(i) %>% 
    pull(junc.id)
  
  JUNC.ID.OVERLAP = dataframe_overlap.table %>% 
    dplyr::slice(i) %>% 
    pull(junc.id.overlap)
  
  if (is.na(JUNC.ID.OVERLAP)){
    res.i = dataframe_count.table %>% 
      # Filter out for only the current junction of the iteration
      dplyr::filter(junc.id == JUNC.ID) %>% 
      
      # Warp dataframe to have three columns: 1) junc.id, 2) case, 3) count
      gather(case, value, 2:ncol(dataframe_count.table)) %>% 
      
      # Remove the junc.id column and only retain: 1) case, 2) count
      spread(junc.id, value) %>% 
      
      # Properly rename the count column as "count"
      dplyr::rename(count = JUNC.ID) %>% 
      
      # Remove any NA's by converting them into 0's
      mutate(count = ifelse(is.na(count), 0, count)) %>% 
      
      # Since depth is the sum of the counts of the junction and its overlapping junction, in the case in which there are
      # no overlapping junctions (JUNC.ID.O == NA), then Depth = count + 0
      mutate(depth = count) %>%
      
      # Frequency is defined as the ratio of the junction count and depth (Frequency = count/(count + count.o))
      mutate(freq = round(count / depth, 4)) %>% 
      mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
      
      mutate(judge = ifelse(count >= 10 & depth >= 20 & freq >= 0.01, 1, 0)) %>% 
      mutate(junc.id = JUNC.ID) %>% 
      dplyr::select(junc.id, case, count, depth, freq, judge)
  }
  else{
    # If there is an overlapping junction (e.g. if JUNC.ID.OVERLAP is not == NA),
    # bind the count data for both the junction of interest (JUNC.ID) and the overlapping junction (JUNCTION.ID.OVERLAP)
    res.i = bind_rows(
      dataframe_count.table %>% dplyr::filter(junc.id == JUNC.ID), 
      dataframe_count.table %>% dplyr::filter(junc.id == JUNC.ID.OVERLAP)) %>% 
      
      # Warp the dataframe to flip the columns and the row values: Junction IDs of the junction and its overlapping junction are columns and TCGA ID are rows
      # Columns should now be 1) TCGA Cases, 2) Junction ID Count (count), 3) Overlapping Junction ID Count (count.o)
      gather(case, value, 2:ncol(dataframe_count.table)) %>% 
      spread(junc.id, value) %>% 
      dplyr::rename(count = JUNC.ID, count.o = JUNC.ID.OVERLAP) %>%
      
      # Remove any NA values by converting any NAs to 0
      mutate(count = ifelse(is.na(count), 0, count)) %>% 
      mutate(count.o = ifelse(is.na(count.o), 0, count.o)) %>% 
      
      # Depth is defined as the Sum of the counts of both junctions (Depth = count + count.o)
      mutate(depth = count + count.o) %>% 
      
      # Frequency is defined as the ratio of the junction count and depth (Frequency = count/(count + count.o))
      mutate(freq = round(count / depth, 4)) %>% 
      mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
      
      # Judgement is passed if count is at least 10, depth is at least 20, and frequency is greater than 1%
      mutate(judge = ifelse(count >= 10 & depth >= 20 & freq >= 0.01, 1, 0)) %>% 
      
      # Generate a column for the junction ID (junc.ID)
      mutate(junc.id = JUNC.ID) %>% 
      
      dplyr::select(junc.id, case, count, depth, freq, judge)
  }
  
  # Warp the rows back so that the TCGA ID's are columns once again (429 columns) and the quality is the row
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
  
  # Combine the rows together into separate tables for Count, Depth, Frequency, and Judgement
  if(i == 1){
    count = count.i; 
    depth = depth.i; 
    freq = freq.i; 
    judge = judge.i
  }else{
    count = bind_rows(count, count.i); 
    depth = bind_rows(depth, depth.i);
    freq = bind_rows(freq, freq.i);
    judge = bind_rows(judge, judge.i);
  }
}

RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print() ;

# Export files without PSR filter
setwd(directory_out)
filename_raw_count   = "2021_0827_count_pass_filter_in_tcga_n429_without_PSR_0.1_filter.tsv"
filename_raw_depth   = "2021_0827_depth_pass_filter_in_tcga_n429_without_PSR_0.1_filter.tsv"
filename_raw_freq    = "2021_0827_freq_pass_filter_in_tcga_n429_without_PSR_0.1_filter.tsv"
filename_raw_judge   = "2021_0827_judge_pass_filter_in_tcga_n429_without_PSR_0.1_filter.tsv"
write_tsv(count, filename_raw_count, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(depth, filename_raw_depth, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(freq,  filename_raw_freq, na = "NA", col_names = T, quote_escape = "double") ;
write_tsv(judge, filename_raw_judge, na = "NA", col_names = T, quote_escape = "double") ;


###########################################################################
#  Step 2: Calculate Judgement and Filter out Junctions Based on PSR > 0.10
###########################################################################
# 2a. Prepare a list of all disease subtypes (General list)
cases = list(
  dataframe_tcga$case,
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-GBM"],
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG"], 
  dataframe_tcga$case[dataframe_tcga$cls == "IDHwt"],
  dataframe_tcga$case[dataframe_tcga$cls == "IDH-A"],
  dataframe_tcga$case[dataframe_tcga$cls == "IDH-O"]) 

names(cases) = c("all", "TCGA-GBM", "TCGA-LGG", "IDHwt", "IDH-A", "IDH-O")

# 2a. Prepare a list of all disease subtypes (IDH Subtypes list)
cases_idh = list(
  dataframe_tcga$case,
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG" & dataframe_tcga$cls == "IDH-A"],
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG" & dataframe_tcga$cls == "IDH-O"],
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-LGG" & dataframe_tcga$cls == "IDHwt"],
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-GBM" & dataframe_tcga$cls != "IDHwt"],
  dataframe_tcga$case[dataframe_tcga$study == "TCGA-GBM" & dataframe_tcga$cls == "IDHwt"]
)

names(cases_idh) = c("all", "IDHmut-A", "IDHmut-O", "IDHwt-A", "IDHmut-GBM", "IDHwt-GBM")

# Depending on which project is being studied:
# project[1] refers to the general disease types looked by Taka and Darwin
# project[2] refers to the IDH subtypes looked by Darwin
list_project = list(cases, cases_idh)
list_judge.all = list(judge.all_cases = NULL, judge.all_cases_idh = NULL)

# 2b. Calculate the positive sample rates for all of the disease subtypes
for (h in 1:length(list_project)){
  project = list_project[[h]]
  for(i in 1:length(project)){
    # Status bar
    print(h)
    print(i)
    
    # Select the disease type for each iteration (e.g. 1 = "all", 2 = "TCGA-GBM")
    # From the overall Judgement Table, retain TCGA samples that are covered by the disease group
    # Generate a new Judgement Table with only these TCGA samples
    CASES = project[[i]]
    judge.i = judge %>%
      dplyr::select(junc.id, colnames(judge)[is.element(colnames(judge), CASES)])
    
    # Calculate the positive sample rate (PSR)
    judge.i = judge.i %>% 
      # Calculate the number of positive samples per junction and make a column for it
      mutate(total_positive = apply(judge.i[, -1], 1, sum)) %>%
      # Calculate the PSR = (# of Positive Samples/Total # of Samples)
      mutate(psr = round(total_positive / length(CASES), 4)) %>% 
      # Remove everything except for the Junction ID and its associated PSR
      dplyr::select(junc.id, psr)
    
    # Change the 2nd Column (the name of disease group)
    if (h == 1){colnames(judge.i)[2] = paste0(names(cases)[[i]], "_n", length(CASES))}
    if (h == 2){colnames(judge.i)[2] = paste0(names(cases_idh)[[i]], "_n", length(CASES))
}
    if (i == 1){
      list_judge.all[[h]] = judge.i}
    else {
      list_judge.all[[h]] = list_judge.all[[h]] %>% 
        full_join(judge.i, by = "junc.id")}
  }
}

judge.all_cases = list_judge.all[[1]] %>% 
  mutate(max = apply(judge.all_cases[, -1], 1, max)) %>% 
  dplyr::filter(max >= 0.1)

# This added variable is the df that contains all of the NJs regardless of PSR > 0.1 (used to graph PSR dot plot of TCGA vs GTEx)
judge.all_cases_allnjs = list_judge.all[[1]] 
judge.all_cases_allnjs = judge.all_cases_allnjs %>% mutate(max = apply(judge.all_cases_allnjs[, -1], 1, max))
setwd(directory_out)
write_tsv(judge.all_cases_allnjs, "2023_0223_psr_all_njs_tcga_LGG_GBM.tsv", na = "NA", col_names = T, quote_escape = "double") ;


judge.noIDHmutGBM = list_judge.all[[2]] %>% 
  select(-c("IDHmut-GBM_n5"))

judge.all_cases.idh = list_judge.all[[2]] %>% 
  # Since IDHmut-GBM has only a sample size of 5, I can remove it from judgement consideration since if only one sample 
  # For now I will include it so that the average/median number of junctions can be compared with all groups
  # select(-c("IDHmut-GBM_n5")) %>%    
  mutate(max = apply(judge.all_cases.idh[, -1], 1, max)) %>% 
  dplyr::filter(max >= 0.1)


###########################################################################
#  Step 3: Generate the Final Count, Depth, Freq, Judge, and Overlap Tables
###########################################################################

# Generate the Tables for Project 1 (General Cases)
dataframe_count.pass.1 = count %>%
  semi_join(judge.all_cases, by = "junc.id")
dataframe_depth.pass.1 = depth %>% 
  semi_join(judge.all_cases, by = "junc.id")
dataframe_freq.pass.1 = freq %>%
  semi_join(judge.all_cases, by = "junc.id")
dataframe_judge.pass.1 = judge %>%
  semi_join(judge.all_cases, by = "junc.id")
dataframe_overlap.table.1 = dataframe_overlap.table %>% 
  semi_join(judge.all_cases, by = "junc.id")

# Generate the Tables for Project 2 (IDH1 Subtype Cases)
dataframe_count.pass.2 = count %>% 
  semi_join(judge.all_cases.idh, by = "junc.id")
dataframe_depth.pass.2 = depth %>%
  semi_join(judge.all_cases.idh, by = "junc.id")
dataframe_freq.pass.2 = freq %>%
  semi_join(judge.all_cases.idh, by = "junc.id")
dataframe_judge.pass.2 = judge %>%
  semi_join(judge.all_cases.idh, by = "junc.id")
dataframe_overlap.table.2 = dataframe_overlap.table %>% 
  semi_join(judge.all_cases.idh, by = "junc.id")

# Check for NA values - should all equate to NA's
dataframe_count.pass.1 %>% anyNA() %>% print()
dataframe_depth.pass.1 %>% anyNA() %>% print()
dataframe_freq.pass.1 %>% anyNA() %>% print()
dataframe_judge.pass.1 %>% anyNA() %>% print()

dataframe_count.pass.2 %>% anyNA() %>% print()
dataframe_depth.pass.2 %>% anyNA() %>% print()
dataframe_freq.pass.2 %>% anyNA() %>% print()
dataframe_judge.pass.2 %>% anyNA() %>% print()



###########################################################################
#  Step 4: Output Files ---------------------------------------------------
###########################################################################

setwd(directory_out)

# Export Project 1 --------------------------------------------------------
# Count Table
filename_count.1 = "P1_Count_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_count.pass.1, filename_count.1, na = "NA", col_names = T, quote_escape = "double") ;

# Depth Table
filename_depth.1 = "P1_Depth_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_depth.pass.1, filename_depth.1, na = "NA", col_names = T, quote_escape = "double") ;

# Frequency Table
filename_freq.1 = "P1_Freq_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_freq.pass.1, filename_freq.1, na = "NA", col_names = T, quote_escape = "double") ;

# Judgement Table (Same as Taka's judge.2.pass)
filename_judge.1 = "P1_Judgement_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_judge.pass.1, filename_judge.1, na = "NA", col_names = T, quote_escape = "double") ;

# Judgement Table (Raw without Filtering out PSR < 0.10)
filename_judge.all.1 = "P1_PSR_20200818"
write_tsv(judge.all_cases, filename_judge.all.1, na = "NA", col_names = T, quote_escape = "double") ;

# Overlap Table
filename_overlap.1 = "P1_Overlap_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_overlap.table.1, filename_overlap.1, na = "NA", col_names = T, quote_escape = "double") ;



# Export Project 2 --------------------------------------------------------
# Count Table
filename_count.2 = "IDH_Count_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_count.pass.2, filename_count.2, na = "NA", col_names = T, quote_escape = "double") ;

# Depth Table
filename_depth.2 = "IDH_Depth_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_depth.pass.2, filename_depth.2, na = "NA", col_names = T, quote_escape = "double") ;

# Frequency Table
filename_freq.2 = "IDH_Freq_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_freq.pass.2, filename_freq.2, na = "NA", col_names = T, quote_escape = "double") ;

# Judgement Table (Same as Taka's judge.2.pass)
filename_judge.2 = "IDH_Judgement_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_judge.pass.2, filename_judge.2, na = "NA", col_names = T, quote_escape = "double") ;

# Judgement Table (Raw without Filtering out PSR < 0.10)
filename_judge.all.2 = "IDH_PSR_20200818"
write_tsv(judge.all_cases.idh, filename_judge.all.2, na = "NA", col_names = T, quote_escape = "double") ;

# Overlap Table
filename_overlap.2 = "IDH_Overlap_Table_Retained_and_Passed_Junctions_20200818"
write_tsv(dataframe_overlap.table.2, filename_overlap.2, na = "NA", col_names = T, quote_escape = "double") ;
