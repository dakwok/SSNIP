# Title: "Step 1 - Tumor Purity"
# July 02, 2020 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Remove samples that have a tumor purity of lower than 60%.
# Pipeline: 1. Download and open the data summary table provided by Kahles 2018
#           2. Filter out patient/samples that have a tumor purity that is less than 0.60
#           3. Generate a new list of samples that are filtered out by this criteria

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

install.packages("readxl")
install.packages("dplyr")
install.packages("fs")
install.packages("rlist")
install.packages("tidyverse")
install.packages("tibble")
install.packages("data.table", dependencies=TRUE)

library(readxl)
library(tidyverse)
library(rlist)
library(data.table)

directory_in = "{PATH_TO_INPUT}"
setwd(directory_in)

filename_kahles  = "{PATH_TO_TUMOR_PURITY_FILE}"
dataframe_kahles = as_tibble(read_excel(filename_kahles))


###########################################################################
#  Step 1: Generate Glioma Specific Kahles 2018 Data Frame ----------------
###########################################################################

# Filter for RNA-seq samples specific to gliomas
# This data frame will contain 666 samples with columns for:
# Case, Study, Histology, Grade, IDH Status, 1p/19q Codeletion, IDH/codel Subtype, Original Subtype, Absolute Purity

dataframe_kahles_glioma = dataframe_kahles %>% 
  filter(RNAseq == "Yes") %>% 
  select(Case, Study, Histology, Grade, `IDH status`, `1p/19q codeletion`, `IDH/codel subtype`, `Original Subtype`, `ABSOLUTE purity`) %>% 
  mutate(Study = ifelse(Study == "Brain Lower Grade Glioma", "TCGA-LGG", "TCGA-GBM"))

# Unfortunately, certain patients in this data frame have an IDH/codel subtype of 'NA'
# Therefore, these samples must be assigned in order to preserve as many samples as possible
# The following are the 13 samples with NA values

dataframe_kahles_glioma %>% 
  filter( `IDH/codel subtype` == "NA" ) %>% 
  select(Case, Study, Histology, `IDH status`)

# As such, retain the 13 samples by manually checking the IDH status in RNAseq data using IGV 
dataframe_kahles_glioma = dataframe_kahles_glioma %>% 
  mutate(`IDH/codel subtype` = ifelse(Study == "TCGA-GBM" & `IDH/codel subtype` == "NA", 
                                      ifelse(Case == "TCGA-06-5417", "IDHmut-non-codel", "IDHwt"), 
                                      `IDH/codel subtype`)) %>% 
  filter(`IDH/codel subtype` != "NA")

table(dataframe_kahles_glioma$Study, dataframe_kahles_glioma$`IDH/codel subtype`)
#          IDHmut-codel IDHmut-non-codel IDHwt
# TCGA-GBM            0                9   145
# TCGA-LGG          168              247    94

table(dataframe_kahles_glioma$Study)
# TCGA-GBM TCGA-LGG 
#      154      509

table(dataframe_kahles_glioma$Study, dataframe_kahles_glioma$Grade)
#           G2  G3  G4  NA
# TCGA-GBM   0   0 153   1
# TCGA-LGG 214 238   0  57


###########################################################################
#  Step 2: Generate a Data Frame of Patients with Tumor Purity Parameter --
###########################################################################

dataframe_kahles_glioma_TP = dataframe_kahles_glioma %>% 
  mutate(cls = ifelse(`IDH/codel subtype` == "IDHmut-non-codel", "IDH-A", ifelse(`IDH/codel subtype` == "IDHmut-codel", "IDH-O", "IDHwt"))) %>% 
  rename(case = Case, study = Study, histology = Histology, grade = Grade, IDH = `IDH status`, codel = `1p/19q codeletion`, abs.purity = `ABSOLUTE purity`) %>% 
  dplyr::select(case, study, cls, grade, histology, IDH, codel, abs.purity) %>% 
  mutate(abs.purity = as.numeric(abs.purity))

dataframe_kahles_glioma_TP %>% 
  dplyr::filter(is.na(abs.purity)) %>% 
  nrow()

dim(dataframe_kahles_glioma_TP)
dataframe_kahles_glioma_TP %>% 
  head() %>% 
  print()
table(dataframe_kahles_glioma_TP$cls)
#  IDH-A IDH-O IDHwt 
#  256   168   239 

summary(dataframe_kahles_glioma_TP$abs.purity)


###########################################################################
#  Step 3: Generate List of Patients with Tumor Purity > 0.60 -------------
###########################################################################

threshold_TP = 0.60

dataframe_kahles_glioma_TP_0.60 = dataframe_kahles_glioma_TP %>% 
  dplyr::filter(!is.na(abs.purity)) %>% 
  dplyr::filter(abs.purity >= threshold_TP)

# check -------------------------------------------------------------------
dim(dataframe_kahles_glioma_TP_0.60)
table(dataframe_kahles_glioma_TP_0.60$cls)
# IDHwt IDH-A IDH-O
#   166   140   123
table(dataframe_kahles_glioma_TP_0.60$study, dataframe_kahles_glioma_TP_0.60$cls)
#          IDHwt IDH-A IDH-O
# TCGA-GBM   110     5     0
# TCGA-LGG    56   135   123


###########################################################################
#  Step 4: Output Data ----------------------------------------------------
###########################################################################

directory_output = "{PATH_TO_OUTPUT}"
setwd(directory_output)

filename_output = "TCGA_Patient_List_Post_TumorPurity_Filter_0.80_20201015"
write.table(dataframe_kahles_glioma_TP_0.60, filename_output, col.names = T, row.names = F, sep = "\t", quote = F)
