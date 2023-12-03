# Title: Step 28b: Visulalize expression level of neojunctions relative to normal junctions (GBM)

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))

library(readxl)
library(tidyverse)
library(ggsci)
library(data.table)

#  Load Directories -------------------------------------------------------
directory_00 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/external/rnaseq_sj.out.tab_20201120"
directory_fig = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/27_rna_validation_cell_lines/figures"
filename_sj10417 = "SF10417-IDH1KO_RNAseq_hg19_SJ.out.tab"
filename_sj10602 = "SF10602p_RNAseq_hg19_SJ.out.tab"

setwd(directory_00)
list_files = c(filename_sj10417, filename_sj10602)

count_table = tibble()
for (i in 1:length(list_files)){
  filename_i = list_files[i]
  if (grepl(".tab", filename_i) == T){
    print(i)
    
    # Get junction ID from sj.out.tab files
    mayo_i = read_tsv(filename_i, col_names = F)
    sj_i = mayo_i %>% 
      dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
      mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
      dplyr::filter(strand !=	"undefined") %>%
      mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
      dplyr::select(junc.id, n.uniq.map) %>% 
      mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map)); 
    
    colnames(sj_i)[2] = filename_i
    
    if (ncol(count_table) == 0){count_table = sj_i}
    if (ncol(count_table) != 0){count_table = full_join(count_table, sj_i)}
  }
}

count_table = count_table[!duplicated(count_table$junc.id),] # 533905
count_table[is.na(count_table)] <- 0  # Convert all NA values to 0

# log2 transform the data
log2_table = count_table
log2_table[,2:ncol(log2_table)] <- log2(log2_table[,2:ncol(log2_table)] + 0.01)

log2_values   = subset(log2_table, select = -c(junc.id))
log2_rownames = subset(log2_table, select = c(junc.id))
log2_values$avg_count = apply(log2_values,1,mean)
log2_mean = subset(log2_values, select = c(avg_count))
log2_summary = cbind(log2_rownames, log2_mean)

# Filter NJs in the dataset
nj = subset(read_tsv("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/28_rna_validation_cell_lines_gbm_mayo/out/2022_1020_detection_of_all_nj_candidates_in_mayo_gbm_lines_by_nj.tsv", col_names = T), select = c("junc.id"))
log2_nj = inner_join(log2_summary, nj, by = "junc.id") # 84 NJ found in the SF lines
log2_wt = anti_join(log2_summary, nj, by = "junc.id")  # 533821 normal junctions

# Box and whisker plot
log2_nj$detected_in = "mutant"
log2_wt$detected_in = "normal"

log2_all = rbind(log2_wt, log2_nj)

# Reorder the groups by changing the order of the factor log2_all$detected_in 

log2_all$detected_in <- factor(log2_all$detected_in , levels=c("normal", "mutant"))


my_colors = c("gray80", "#b5179e", "#ca472f")
ggplot(log2_all, aes(x=detected_in, y=avg_count, fill = detected_in)) + 
  theme_minimal() +
  scale_fill_manual(values = my_colors) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        axis.title.y = element_text(size=16, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, face = "bold")) +
  xlab("") +
  ylab("log2(junction reads)") +
  ggtitle("") +
  theme(plot.title = element_text(size = 16,  family="Helvetica", face = "bold", hjust = 0.5),  # Change the text of the title and adjust to center with hjust
        text = element_text(size = 15,  family="Helvetica"),                                    # Change the text of the legend and the axis  
        panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size = 0)) +               # Remove the background grids
  geom_boxplot(alpha = 1, width=0.4)
# Export size = 3.5 x 5.25

# Calculate wilcox test results between the normal junctions and neojunction peptides
list_normal = as.numeric(t((log2_wt$avg_count)))
list_tumor = as.numeric(t(as.list(log2_nj$avg_count)))

wilcox.test(list_normal, list_tumor) #  W = 148781, p-value = 0.3484    # As of 2022_1020
