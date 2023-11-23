# Title: "Step 6: SJ Overlap Table"
# July 31, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Generate and prepare an Overlap Table for each candidate splicing junction event
#          This step requires both the annotated (Step 4) and non-annotated (Step 5) SJs in order
#          to identify ALL overlapping junctions

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE)) ;

library(tidyverse) ;
library(readxl) ;
library(ggsci) ;
library(data.table)

# Establish Directories ---------------------------------------------------

directory_samples = "{PATH_TO_INPUT}"

directory_04 = "{PATH_TO_STEP04_OUTPUT}"
directory_05 = "{PATH_TO_STEP05_OUTPUT}"
directory_06 = "{PATH_TO_STEP06_OUTPUT}"

# Load Files --------------------------------------------------------------
# From Step 4: SJ ID's (Annotated)
setwd(directory_04)
filename_sj.annot = "SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_20200623.tsv"
dataframe_sj.annot = read_tsv(filename_sj.annot, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 

# From Step 5b: SJ PSR's (Non-Annotated with Read Count > 10, PSR > 10, and Protein Coding)
setwd(directory_05)
filename_sj.nonannot = "2022_0804_sj_list_non_annotated_protein_coding.tsv"
dataframe_sj.nonannot = read_tsv(filename_sj.nonannot, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 



###########################################################################
#  Step 1: Generate a Complete List of Junction ID's to Analyze -----------
###########################################################################

# 1. Take the Annotated SJ List and edit it such that it matches the Non-Annotated SJ List
#    with "Junc.ID", "Strand", "Int.Start", "Int.End"

dataframe_sj.nonannot = dataframe_sj.nonannot %>% 
  dplyr::select(junc.id) %>% 
  mutate(chr = (gsub ("chr", "", sapply(strsplit(junc.id, ":"), "[[", 1)))) %>% 
  mutate(strand = sapply(strsplit(junc.id, ":"), "[[", 2)) %>% 
  mutate(int.start = as.numeric(gsub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3))) + 1) %>% 
  mutate(int.end = as.numeric(gsub(".*-", "", sapply(strsplit(junc.id, ":"), "[[", 3)))) ;

# 2. Identify which "annotated" and "non-annotated" junction counts to collect
start.time = proc.time()

list_sj = list()
junc.ids = NULL

for(i in 1:nrow(dataframe_sj.nonannot)) {
  
  # Progress Bar
  print(i/nrow(dataframe_sj.nonannot)*100)
  
  # Slice out row(i) for every iteration in the for-loop
  sj.i = dataframe_sj.nonannot %>% dplyr::slice(i)
  
  # Extract the JUNC.ID, STRAND, START, and END for junction(i)
  CHR     = sj.i %>% pull(chr)
  JUNC.ID = sj.i %>% pull(junc.id) ; 
  STRAND  = sj.i %>% pull(strand) ; 
  START   = sj.i %>% pull(int.start) ; 
  END     = sj.i %>% pull(int.end) ; 
  
  # Find all of the overlapping junctions from the "annotated" SJ list
  JUNC.ID.OVERLAP = dataframe_sj.annot %>% 
    dplyr::filter(chr == CHR, strand == STRAND & int.start < END & START < int.end) %>% 
    pull(junc.id) ;
  
  # Create a list of overlapping junctions for each of the "non-annotated" junctions 
  if(length(JUNC.ID.OVERLAP) == 0){
    list_sj[[i]] = tibble(
      junc.id = JUNC.ID, 
      junc.id.overlap = NA)}
  else{
    list_sj[[i]] = tibble(
      junc.id = JUNC.ID, 
      junc.id.overlap = JUNC.ID.OVERLAP)}
  
  
  junc.ids.i = list_sj[[i]] %>% 
    # Gather and list all of the junctions (junction of interest including overlapping junctions)
    gather("label", "junc.id", 1:2) %>% 
    
    # If there are no overlapping junctions, remove all NA's
    dplyr::filter(!is.na(junc.id)) %>% 
    
    # Select only unique/distinct junc.id
    distinct(junc.id) %>% 
    pull(junc.id) ;
  
  junc.ids = unique(c(junc.ids, junc.ids.i)) ;
}

RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print() ;

# user    system   elapsed 
# 11262.674   548.539 11893.480 

head(junc.ids)


###########################################################################
#  Step 2: Analyze the Count Table with the Complete Junction ID List -----
###########################################################################

# 1. Generate an ordered list of junctions to bind to the sj.out.tab files
# 1a. Create the dataframe for SJ Counts with the first column being all of the Junction ID's
#     from the junc.id list made from Step 1
dataframe_sj.count = tibble(junc.id = junc.ids) %>% # for sorting
  # 1b. Make a new column that is the starting coordinate of the junction
  mutate(start = as.numeric(gsub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3)))) %>% 
  # 1c. Order the rows numerically by their starting coordinates
  arrange(start) %>% 
  # 1d. Select out the now ordered junction ID's
  dplyr::select(junc.id)

dataframe_sj.count %>% dim() %>% print()


# 2. Bind the ordered list of junctions to the sj.out.tab files for GBM and LGG
setwd(directory_samples)
list_rnaseq.files = list.files()[grepl("SJ.out.tab",list.files())] #2022_0627 -- total of 527 samples out of 535

start.time = proc.time()


for (j in 1:length(list_rnaseq.files)){
  # 2b. Work with each file in the directory iteratively
  filename_j = list_rnaseq.files[j]
  dataframe_sj.j = read_tsv(filename_j, na = c("", "NA"), col_names = F, col_type = cols(X1 = col_character()))
  
  append_i = dataframe_sj.j %>% 
    # Rename the columns of the Count Table
    dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
    # Change the "strand" values from "2" to "-", "1" to "+", and neither to "undefined"
    mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
    # Remove any junctions that have "strand" == undefined
    dplyr::filter(strand !=	"undefined") %>%
    # Generate a new column for "junc.id"
    mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
    dplyr::select(junc.id, n.uniq.map)
  
  dataframe_sj.count = left_join(dataframe_sj.count, append_i) %>% 
    mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map))
    
  colnames(dataframe_sj.count)[ncol(dataframe_sj.count)] = strsplit(list_rnaseq.files[j], "_SF.out.tab")[[1]][1]

  # Status Bar
  print(((ncol(dataframe_sj.count) - 1)/length(list_rnaseq.files))*100)
}


RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print()

#    user  system elapsed 
# 385.068   9.433 406.379 


# 3. Identify the "Most-Dominant" Junctions Among the Multiple Overlaps
start.time = proc.time()

dataframe_overlap.table = NULL

for(i in 1:length(list_sj)){
  # 3a. Use the first iteration to establish a new dataframe for dataframe_overlap.table
  if(nrow(list_sj[[i]]) == 1){
    dataframe_overlap.table = dataframe_overlap.table %>% 
      bind_rows(list_sj[[i]])}
  else{
    # sj_temp will be the current junction of interest with its overlapping junctions
    # sj_temp contains two columns: 1) "junc.id", 2) "junc.id.overlap"
    sj_temp = list_sj[[i]]  %>% 
      # Select out the "junc.id.overlap" column
      dplyr::select(junc.id.overlap) %>% 
      # 3b. Get the count data for the overlapping junctions by left_join() the overlapping junctions
      # with the count data for each overlapping junction from dataframe_sj.count
      left_join(dataframe_sj.count, by = c("junc.id.overlap" = "junc.id"))
    
    sj_temp = sj_temp %>% 
      # 3c. Generate a new column named "sum" that totals the count in each row 
      mutate(sum = apply(sj_temp %>% 
                           dplyr::select(-1), 1, sum)) %>% 
      # 3d. Remove the count data, and retain only the "junc.id.overlap" and "sum" columns
      dplyr::select(junc.id.overlap, sum) %>% 
      # 3e. Arrange the overlapping junctions in descending order and retain only the highest count junction
      arrange(desc(sum)) %>% 
      head(n = 1)
    
    # 3f. Add the overlapping junction to dataframe_overlap.table if it has not been added already
    dataframe_overlap.table = dataframe_overlap.table %>% 
      bind_rows(list_sj[[i]] %>% 
                  semi_join(sj_temp, by = "junc.id.overlap"))
    
    # Progress Report
    print((i/length(list_sj))*100)
  }
}

RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print() ; 


###########################################################################
#  Step 3: Edit the SJ Count Table ----------------------------------------
###########################################################################

# Retain the count data for the overlapping junctions
dataframe_sj.count.retain = dataframe_sj.count %>% 
  right_join(dataframe_overlap.table %>% 
               gather("label", "junc.id", 1:2) %>% 
               dplyr::filter(!is.na(junc.id)) %>% 
               distinct(junc.id),
             by = "junc.id")

dataframe_sj.count.retain %>% 
  distinct(junc.id, .keep_all = T) %>% 
  anyNA()



###########################################################################
#  Step 4: Export Files ---------------------------------------------------
###########################################################################

setwd(directory_06)
filename_output.overlap.table = "2022_0812_sj_overlap_table.tsv"
write_tsv(dataframe_overlap.table, filename_output.overlap.table, na = "NA", col_names = T, quote_escape = "double") ;


filename_output.retain = "2022_0812_sj_count_table_sj.out.tab_retained.tsv"
write_tsv(dataframe_sj.count.retain, filename_output.retain, na = "NA", col_names = T, quote_escape = "double") ;
