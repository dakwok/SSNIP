# Identify top candidates
# These are the candidates found in both RNA
rnaseq = read_tsv("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/27_rna_validation_cell_lines/out/2022_1020_detection_of_all_nj_candidates_in_sf_lines.tsv", col_names = T)

cptac = read_tsv("/Users/darwinkwok/Desktop/cclc01/okadalab/data1/dkwok/proj_01_altspl/19_cptac_output/out_rscript_dk_detection/2022_0401_unique_njs_cptac.tsv", col_names = T)

test = inner_join (rnaseq, cptac, by = "junc.id")
write_tsv(test, "2023_0426_validated_8_nj_candidates.tsv", na = "NA", col_names = T, quote_escape = "double") ;
