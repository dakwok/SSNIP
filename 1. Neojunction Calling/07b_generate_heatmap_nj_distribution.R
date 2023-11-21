# Title: "Step 7b: Heatmap for Read Frequency of TCGA samples"
# November 30, 2022 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Generate heatmaps for the read frequency table of TCGA

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE)) ;

library(tidyverse) ;
library(readxl) ;
library(ggsci) ;
library(data.table);
library("pheatmap")
library(RColorBrewer)
library(viridis)


# Establish Directories ---------------------------------------------------
directory_07   = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/07_count_depth_freq_judge_psr/output"
directory_10   = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/10_extract_neojunctions/output/"
directory_29 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/29_immunogenic_neojunctions/"
directory_meta = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/external/meta_tcga"

# Load Files --------------------------------------------------------------
# From Step 7: Read Frequency Table
setwd(directory_07)
filename_overlap.table = "P1_Freq_Table_Retained_and_Passed_Junctions_20200818.tsv"
freq = read_tsv(filename_overlap.table, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 

# From Step 7: Judge Frequency Table
filename_judge.table = "P1_Judgement_Table_Retained_and_Passed_Junctions_20200818.tsv"
judge = read_tsv(filename_judge.table, na = c("", "NA"), col_names = T, col_types = c(chr = col_character())) ; # non-annot pass-filter junc.ids 

# Final NJ list
setwd(directory_10)
df_nj = subset(read_tsv("Judgement_Neojunctions_20201020.tsv", col_names = T), select = c(junc.id))

# Immunogenic NJ list
setwd(directory_29)
df_immuno_nj = read_tsv("2022_1202_immunogenic_njs.txt", col_names = T)

# Annotation Metadata for TCGA and GTEx
setwd(directory_meta)
filename_tcga_meta = "07032020_tcga_glioma_annotation_n663.tsv"
meta_tcga = read_tsv(filename_tcga_meta, na = c("", "NA"), col_names = T)
#meta_tcga$case <- paste0(meta_tcga$case, "-01")   # Change the ID of the meta data to match the column names of the TPM dataframe

###########################################################################
#  Step 1: Generate Heatmap -----------------------------------------------
###########################################################################
# Heatmap meta data
meta_heatmap = subset(meta_tcga, select = c(case, cls))
meta_heatmap$cls <- gsub('-', '_', meta_heatmap$cls)
meta_heatmap2  = meta_heatmap
meta_heatmap   = subset(meta_heatmap2, select = -c(case))
meta_heatmap   = as.data.frame(meta_heatmap)
row.names(meta_heatmap) <- meta_heatmap2$case


# Filter for only the neojunctions in generating the heatmap
freq_nj = left_join(df_nj, freq) # 248/249 NJs were detected
judge_nj = left_join(df_nj, judge)
judge_nj_immuno = left_join(df_immuno_nj, judge)

# Convert NA to 0
freq_nj[is.na(freq_nj)] <- 0
judge_nj[is.na(judge_nj)] <- 0
judge_nj_immuno[is.na(judge_nj_immuno)] <- 0

# ***Important*** For the values with Judge = 0, turn to 0 in the Freq Table
for (i in 1:nrow(judge_nj)){
  for (j in 2:ncol(judge_nj)){
    print(paste(i, j))
    judge_sample = as.numeric(judge_nj[i,j])
    if (judge_sample == 0){freq_nj[i,j] = 0}
  }
}

heatmap_freq = subset(freq_nj, select = -c(junc.id))        # Remove the non-numerical columns from the second table to be bound
heatmap_rows = subset(freq_nj, select = c(junc.id))
row.names(heatmap_freq) <- heatmap_rows$junc.id

heatmap_judge = subset(judge_nj, select = -c(junc.id))        # Remove the non-numerical columns from the second table to be bound
heatmap_rows = subset(judge_nj, select = c(junc.id))
row.names(heatmap_judge) <- heatmap_rows$junc.id

heatmap_judge_immuno = subset(judge_nj_immuno, select = -c(junc.id))        # Remove the non-numerical columns from the second table to be bound
heatmap_rows_immuno = subset(judge_nj_immuno, select = c(junc.id))
row.names(heatmap_judge_immuno) <- heatmap_rows_immuno$junc.id

# Remove rows with standard deviation of zero
heatmap_corrected_freq = NULL
for (i in 1:nrow(heatmap_freq)){
  print(i)
  
  row_i = heatmap_freq %>% dplyr::slice(i)
  row_num = as.numeric(row_i)
  sd = sd(row_i)
  if (sd != 0){heatmap_corrected_freq = rbind(heatmap_corrected_freq, row_i)}
}

heatmap_corrected_judge = NULL
for (i in 1:nrow(heatmap_judge)){
  print(i)
  
  row_i = heatmap_judge %>% dplyr::slice(i)
  row_num = as.numeric(row_i)
  sd = sd(row_i)
  if (sd != 0){heatmap_corrected_judge = rbind(heatmap_corrected_judge, row_i)}
}

heatmap_corrected_judge_immuno = NULL
for (i in 1:nrow(heatmap_judge_immuno)){
  print(i)
  
  row_i = heatmap_judge_immuno %>% dplyr::slice(i)
  row_num = as.numeric(row_i)
  sd = sd(row_i)
  if (sd != 0){heatmap_corrected_judge_immuno = rbind(heatmap_corrected_judge_immuno, row_i)}
}

# Set the general annotation colors
annotation_colors_1 = list(cls = c(IDHwt = "#0b84a5", IDH_A = "#f6c85f", IDH_O = "#ca472f"))
metamatrix_heatmap <- pheatmap(mat = heatmap_corrected_judge_immuno,
         color = inferno(100),
         annotation_col = meta_heatmap,
         annotation_colors = annotation_colors_1,
         cluster_rows = T,
         scale = "none", 
         breaks=seq(0, 1, length.out = 100),
         border_color = NA,
         show_rownames = T,
         treeheight_row = 20,
         treeheight_col = 10,
         #cutree_rows = 3,
         #cutree_cols = 3,
         main = "Putative TNJs Passing Read Count, Depth, and Freq Criteria (TCGA GBM/LGG)", fontsize_col = 0.5, fontsize_row = 1)

# Export heatmap (dimension = 12 x 5 Landscape)


###########################################################################
#  Step 3: Generate list of NJ clusters -----------------------------------
###########################################################################
# Generate the list of NJs based on the 3 clusters that were generated
# The way we do this is by finding order by which the NJ junc.id is done from top to bottom in the heatmap
# and manually split based on the first and last NJ in each of the clusters

# 1. Generate the ordered numerical list
#    All junc.ids are numerically ordered, but in the heatmap they are obviously in different positions.
#    Use the following lines of code to determine what order the NJs are then generate a new list of ordered NJs.
list_order_num = metamatrix_heatmap[["tree_row"]][["order"]]
list_order_nj = c()
for (i in 1:length(list_order_num)){
  print(i)
  list_i = list_order_num[i]
  list_order_nj[i] = metamatrix_heatmap$tree_row[["labels"]][list_i]
}

# 2. Look at the heatmap and manually mark when clusters start and end
#    Cluster 1:  1 to  25   (IDH1wt-specific NJs)
#    Cluster 2: 26 to  92   (non-specific NJs)
#    Cluster 3: 93 to 249   (IDH1mut-specific NJs)

list_idh1wt_nj  = list_order_nj[1:25]
list_nonspec_nj = list_order_nj[26:92]
list_idh1mut_nj = list_order_nj[93:249]

df_idh1wt_nj  = tibble("junc.id" =  list_idh1wt_nj)
df_nonspec_nj = tibble("junc.id" =  list_nonspec_nj)
df_idh1mut_nj = tibble("junc.id" =  list_idh1mut_nj)

setwd(directory_07)
write_tsv(df_idh1wt_nj,  "2022_1201_nj_list_idh1wt_cluster.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(df_nonspec_nj, "2022_1201_nj_list_nonspec_cluster.tsv", na = "NA", col_names = T, quote_escape = "double")
write_tsv(df_idh1mut_nj, "2022_1201_nj_list_idh1mut_cluster.tsv", na = "NA", col_names = T, quote_escape = "double")


###########################################################################
#  Step 4: Calculate the  -----------------------------------------------
###########################################################################


