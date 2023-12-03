directory_00 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/external/rnaseq_sj.out.tab_20201120"
directory_11 = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/11_aaseq_prediction/output"
filename_sj10417 = "SF10417-IDH1KO_RNAseq_hg19_SJ.out.tab"
filename_sj10602 = "SF10602p_RNAseq_hg19_SJ.out.tab"

setwd(directory_00)
sj_10417 = read_tsv(filename_sj10417, col_names = F)
sj_10602 = read_tsv(filename_sj10602, col_names = F)

# Convert the sj.out.tab files of the SF10417 and SF10602 into usable formats
sj_10417_edit = sj_10417 %>% 
  dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
  mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
  dplyr::filter(strand !=	"undefined") %>%
  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
  dplyr::select(junc.id, n.uniq.map) %>% 
  mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map)); 

sj_10602_edit = sj_10602 %>% 
  dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>% 
  mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>% 
  dplyr::filter(strand !=	"undefined") %>%
  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
  dplyr::select(junc.id, n.uniq.map) %>% 
  mutate(n.uniq.map = ifelse(is.na(n.uniq.map), 0, n.uniq.map)); 

# Generate tables for ALL 249 neojunctions
# Here we see how many putative neojunctions are detected in each of the SF cell lines
# Open the file with the 249 neojunctions
nj = subset(read_tsv("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/10_extract_neojunctions/output/PSR_Neojunctions_20201020.tsv", col_names = T), select = c("junc.id"))

hits_sf10417 = inner_join(sj_10417_edit, nj)
hits_sf10602 = inner_join(sj_10602_edit, nj)

hits_sf10417 = hits_sf10417 %>% filter(n.uniq.map > 0) # 54 NJs
hits_sf10602 = hits_sf10602 %>% filter(n.uniq.map > 0) # 56 NJs

colnames(hits_sf10417)[2] = "SF10417"
colnames(hits_sf10602)[2] = "SF10602"

complete_table = full_join(hits_sf10417, hits_sf10602)
complete_table[is.na(complete_table)] <- 0

directory_out = "/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/27_rna_validation_cell_lines/out"
setwd(directory_out)
write_tsv(complete_table, "2022_1020_count_of_all_nj_candidates_in_sf_lines.tsv", na = "NA", col_names = T, quote_escape = "double")

complete_table_id  = subset(complete_table, select = c(junc.id))
complete_table_num = subset(complete_table, select = -c(junc.id))
complete_table_num[complete_table_num > 0] <- 1
complete_table_num$sum = rowSums(complete_table_num)
complete_table_2 = cbind(complete_table_id, complete_table_num)

complete_table_3 = subset(complete_table_2, select = c("junc.id", "sum"))
setwd(directory_out)
write_tsv(complete_table_3, "2022_1020_detection_of_all_nj_candidates_in_sf_lines.tsv", na = "NA", col_names = T, quote_escape = "double")

