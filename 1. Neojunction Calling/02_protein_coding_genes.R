# Title: "Step 2 - Protein Coding Genes"
# July 05, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Identify and filter for protein-coding genes from GTF file.
# Pipeline: 1. Download and open the GTF file from the Box or online.
#           2. Filter for genes that are have a gene_biotype = "protein_coding"
#           3. Generate a list of protein-genes and their respective coordinates

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

install.packages("readxl")
install.packages("dplyr")
install.packages("fs")
install.packages("rlist")
install.packages("tidyverse")
install.packages("tibble")
install.packages("data.table", dependencies=TRUE)
install.packages("rapportools") # for is.empty

library(readxl)
library(tidyverse)
library(rlist)
library(data.table)
library(rapportools)

#  Download GTF File ------------------------------------------------------

directory_in = "{PATH_TO_INPUT}"
setwd(directory_in)

filename_gtf = "Homo_sapiens.GRCh37.87.chr.gtf"
data.table_gtf = as_tibble(data.table::fread(filename_gtf))    # 'data.table::fread()' reads .gtf files 


###########################################################################
#  Step 1: Rename Columns in the Original GTF File ------------------------
###########################################################################

# Rename the columns of the original GTF file to have columns 'X1-X9' (Match with Takahide's annotation)

names(data.table_gtf)[1] = "X1"
names(data.table_gtf)[2] = "X2"
names(data.table_gtf)[3] = "X3"
names(data.table_gtf)[4] = "X4"
names(data.table_gtf)[5] = "X5"
names(data.table_gtf)[6] = "X6"
names(data.table_gtf)[7] = "X7"
names(data.table_gtf)[8] = "X8"
names(data.table_gtf)[9] = "X9"


###########################################################################
#  Step 2: Extract the Protein Coding Genes from the Original GTF File ----
###########################################################################

data.table_gtf_protein_coding = data.table_gtf %>% 
  dplyr::filter(X1 != "MT") %>%                # Remove mitochondrial genes
  dplyr::filter(X3 == "transcript") %>%        # Select for transcript genes
  dplyr::filter(grepl("protein_coding", X9))   # Select for protein-coding genes

# Protein-coding genes = 145,548


###########################################################################
#  Step 3: Edit and Finalize GTF File -------------------------------------
###########################################################################

# Important information is in column X9 of the GTF file. This step extracts pertinent information,
# cleans them up, and puts them into their individual columns.

data.table_gtf_final = data.table_gtf_protein_coding %>% 
  # -----------------------------------------------------------------------
  # 1. Generate an Ensembl Transcript ID column ---------------------------
  # -----------------------------------------------------------------------
  mutate(enst = sapply(str_split(data.table_gtf_protein_coding$X9, "; "), "[[", 3)) %>%     # Create a new column "enst" that holds only the Transcript ID from Column X9
            # str_split() --> split the string into pieces based on a pattern ("; ")
            # saaply() --> retains the piece that corresponds to the index number (3)
  mutate(enst = gsub("transcript_id \"", "", enst)) %>%                                     # Remove "transcript_id" and retain only the ENST ID
            # gsub(x, y, source) --> replace x with y from the string source
  mutate(enst = gsub("\"", "", enst)) %>%                                                   # Remove the " from the end to finalize
  # -----------------------------------------------------------------------
  # 2. Generate a Gene Symbol column --------------------------------------
  # -----------------------------------------------------------------------
  mutate(symbol = sapply(str_split(data.table_gtf_protein_coding$X9, "; "), "[[", 5) ) %>%  # Create a new column "symbol" that holds the Gene Name from Column X9
  mutate(symbol = gsub("gene_name \"", "", symbol)) %>%                                     # Remove "gene_name" and retain only the Gene Symbol
  mutate(symbol = gsub("\"", "", symbol)) %>%                                               # Remove the " from the end to finalize 
  # -----------------------------------------------------------------------
  # 3. Generate an Ensembl Gene ID column ---------------------------------
  # -----------------------------------------------------------------------
  mutate(ensg = sapply(str_split(data.table_gtf_protein_coding$X9, "; "), "[[", 1) ) %>%    # Create a new column "ensg" that holds the Gene ID from Column X9
  mutate(ensg = gsub("gene_id \"", "", ensg)) %>%                                           # Remove "gene_id" and retain only the Gene ID 
  mutate(ensg = gsub("\"", "", ensg))                                                       # Remove the " from the end to finalize 


###########################################################################
#  Step 4: Output Data ----------------------------------------------------
###########################################################################

directory_output = "{PATH_TO_OUTPUT}"
setwd(directory_output)

filename_output = "GTF_Protein_Coding_Genes_20200606"
write.table(data.table_gtf_final, filename_output, col.names = T, row.names = F, sep = "\t", quote = F)
