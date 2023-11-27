# Generate Figure 5i: Top validated NJs (n=8) and their corresponding neoantigens

# Purpose: Generate heatmaps for the read frequency table of TCGA

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################
rm(list = ls(all.names = TRUE)) ;

library(tidyverse) ;
library(readxl) ;
library(ggsci) ;
library(data.table);
library(RColorBrewer) ;
library(gridExtra)

# Establish Directories ---------------------------------------------------
directory_15 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/15_cross_analyze_algorithms/out"

# Load files --------------------------------------------------------------
setwd(directory_15)
df_all_map = read_tsv("2023_0812_neoA_to_neoJ_map_ALL_summary.tsv", col_names = T)            # 5211
df_na_nj_map = df_all_map %>% dplyr::filter(shared == "Top 10%tile in Both Alg. (HA and MF)") # 1187

###########################################################################
#  Step 1. Plot FS NJ's neoantigen scores for the top 92 neoantigens ------
###########################################################################
# Filter for the scores of frame-shifts (FS) and in-frame (IF)
df_fs = subset(df_na_nj_map %>% dplyr::filter(fs == "fs"), select = c("hla_allele", "score_average"))        # Total = 777
df_if = subset(df_na_nj_map %>% dplyr::filter(fs == "in-frame"), select = c("hla_allele", "score_average"))  # Total = 392

# Combine data frames
df_combined <- rbind(df_fs, df_if)
df_combined$type <- c(rep("Frame-shift (n=777)", nrow(df_fs)), rep("In-frame (n=392)", nrow(df_if)))

# Define color palette for HLA types
hla_colors <- c("HLA-A0101" = "#390099", "HLA-A0201" = "#9e0059", "HLA-A0301" = "#ff0054", "HLA-A1101" = "#ff5400", "HLA-A2402" = "#ffbd00")

ggplot(df_combined, aes(x = type, y = score_average, color = hla_allele)) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = hla_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 12),
        legend.title = element_text(family = "Helvetica", size = 12, face = "bold")) +
  labs(x = "Type of neoantigens", y = "Average presentation score", title = "")

# combine the two data frames
df_all <- rbind(transform(df_fs, type = "Frame-shift"),
                transform(df_if, type = "In-frame"))

# plot grouped boxplot
ggplot(df_all, aes(x = hla_allele, y = score_average, fill = type)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = "HLA Allele", y = "Immunogenicity Score") +
  scale_fill_manual(values = c("Frame-shift" = "#006e90", "In-frame" = "#f18f01")) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 10),
        legend.title = element_text(family = "Helvetica", size = 10, face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")


###########################################################################
#  Step 2. Plot FS NJ's neoantigen scores for the ALL neoantigens ---------
###########################################################################
# Since there are so many values, we might have to do a histogram
# Split into 5 different histograms based on HLA allele

hist = subset(df_all_map, select = c("hla_allele", "score_average", "fs"))
hist_0101 = hist %>% dplyr::filter(hla_allele == "HLA-A0101")
hist_0201 = hist %>% dplyr::filter(hla_allele == "HLA-A0201")
hist_0301 = hist %>% dplyr::filter(hla_allele == "HLA-A0301")
hist_1101 = hist %>% dplyr::filter(hla_allele == "HLA-A1101")
hist_2402 = hist %>% dplyr::filter(hla_allele == "HLA-A2402")

# Log2 transform the scores
hist_0101$score_log2 <- log2(hist_0101$score_average + 0.001)
hist_0201$score_log2 <- log2(hist_0201$score_average + 0.001)
hist_0301$score_log2 <- log2(hist_0301$score_average + 0.001)
hist_1101$score_log2 <- log2(hist_1101$score_average + 0.001)
hist_2402$score_log2 <- log2(hist_2402$score_average + 0.001)

# plot histogram with density curve
gg_0101 = ggplot(hist_0101, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*01:01)", y = "Density") +
  scale_fill_manual(values = c("#006e90", "#f18f01")) +
  xlim(-1, 0) +
  theme_bw()
gg_0201 = ggplot(hist_0201, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*02:01)", y = "Density") +
  scale_fill_manual(values = c("#006e90", "#f18f01")) +
  xlim(-1, 0) +
  theme_bw()
gg_0301 = ggplot(hist_0301, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*03:01)", y = "Density") +
  scale_fill_manual(values = c("#006e90", "#f18f01")) +
  xlim(-1, 0) +
  theme_bw()
gg_1101 = ggplot(hist_1101, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*11:01)", y = "Density") +
  scale_fill_manual(values = c("#006e90", "#f18f01")) +
  xlim(-1, 0) +
  theme_bw()
gg_2402 = ggplot(hist_2402, aes(x = score_log2, fill = fs)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*24:02)", y = "Density") +
  scale_fill_manual(values = c("#006e90", "#f18f01")) +
  xlim(-1, 0) +
  theme_bw()

grid.arrange(gg_0101, gg_0201, gg_0301, gg_1101, gg_2402, ncol = 1)


###########################################################################
#  Step 3. Plot TYPE NJ's neoantigen scores for the top 92 neoantigens ----
###########################################################################
# Filter for the scores of different splice types
df_a3loss = subset(df_na_nj_map %>% dplyr::filter(type == "A3.loss"), select = c("hla_allele", "score_average"))            # Total = 98  # 27830 (all)
df_a3gain = subset(df_na_nj_map %>% dplyr::filter(type == "A3.gain"), select = c("hla_allele", "score_average"))            # Total = 56  #  6730 (all) 
df_a5loss = subset(df_na_nj_map %>% dplyr::filter(type == "A5.loss"), select = c("hla_allele", "score_average"))            # Total = 26  #  7730 (all)
df_a5gain = subset(df_na_nj_map %>% dplyr::filter(type == "A5.gain"), select = c("hla_allele", "score_average"))            # Total = 86  #  7210 (all) 
df_juncin = subset(df_na_nj_map %>% dplyr::filter(type == "JUNC.WITHIN.EXON"), select = c("hla_allele", "score_average"))   # Total = 17  #  9130 (all)
df_juncex = subset(df_na_nj_map %>% dplyr::filter(type == "JUNC.WITHIN.INTRON"), select = c("hla_allele", "score_average")) # Total = 14  #  3360 (all) 
df_exskip = subset(df_na_nj_map %>% dplyr::filter(type == "ES"), select = c("hla_allele", "score_average"))                 # Total = 1   #   165 (all)
df_others = subset(df_na_nj_map %>% dplyr::filter(type == "OTHERS"), select = c("hla_allele", "score_average"))             # Total = 26  #  5690 (all) 

# Combine data frames
df_combined <- rbind(df_a3loss, df_a3gain, df_a5loss, df_a5gain, df_juncin, df_juncex, df_exskip, df_others)
df_combined$type <- c(rep("A3 loss (n=466)", nrow(df_a3loss)), rep("A3 gain (n=241)", nrow(df_a3gain)),
                      rep("A5 loss (n=106)", nrow(df_a5loss)), rep("A5 gain (n=212)", nrow(df_a5gain)),
                      rep("JWE (n=22)", nrow(df_juncin)), rep("JWI (n=41)", nrow(df_juncex)),
                      rep("ES (n=27)", nrow(df_exskip)), rep("OTHERS (n=54)", nrow(df_others)))

# Define color palette for HLA types
hla_colors <- c("HLA-A0101" = "#390099", "HLA-A0201" = "#9e0059", "HLA-A0301" = "#ff0054", "HLA-A1101" = "#ff5400", "HLA-A2402" = "#ffbd00")

ggplot(df_combined, aes(x = type, y = score_average, color = hla_allele)) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = hla_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold", angle = 0, hjust = 0.5),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 10),
        legend.title = element_text(family = "Helvetica", size = 10, face = "bold")) +
  labs(x = "Type of neoantigens", y = "Average presentation score", title = "")


# combine the two data frames
df_all <- rbind(transform(df_a3loss, type = "A3.loss"),
                transform(df_a3gain, type = "A3.gain"),
                transform(df_a5loss, type = "A5.loss"),
                transform(df_a5gain, type = "A5.gain"),
                transform(df_juncin, type = "JUNC.WITHIN.EXON"),
                transform(df_juncex, type = "JUNC.WITHIN.INTRON"),
                transform(df_exskip, type = "ES"),
                transform(df_others, type = "OTHERS"))

# plot grouped boxplot
ggplot(df_all, aes(x = hla_allele, y = score_average, fill = type)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = "HLA Allele", y = "Immunogenicity Score") +
  scale_fill_manual(values = c("A3.gain" = "#f94144", "A3.loss" = "#f3722c",
                               "A5.gain" = "#f8961e", "A5.loss" = "#f9c74f",
                               "ES" = "#90be6d", "JUNC.WITHIN.EXON" = "#43aa8b",
                               "JUNC.WITHIN.INTRON" = "#4d908e", "OTHERS" = "#577590")) +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        legend.text = element_text(family = "Helvetica", size = 10),
        legend.title = element_text(family = "Helvetica", size = 10, face = "bold")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom")


###########################################################################
#  Step 2. Plot FS NJ's neoantigen scores for the ALL neoantigens ---------
###########################################################################
# Since there are so many values, we might have to do a histogram
# Split into 5 different histograms based on HLA allele

hist = subset(df_all_map, select = c("hla_allele", "score_average", "type"))
hist_0101 = hist %>% dplyr::filter(hla_allele == "HLA-A0101")
hist_0201 = hist %>% dplyr::filter(hla_allele == "HLA-A0201")
hist_0301 = hist %>% dplyr::filter(hla_allele == "HLA-A0301")
hist_1101 = hist %>% dplyr::filter(hla_allele == "HLA-A1101")
hist_2402 = hist %>% dplyr::filter(hla_allele == "HLA-A2402")

# Log2 transform the scores
hist_0101$score_log2 <- log2(hist_0101$score_average + 0.001)
hist_0201$score_log2 <- log2(hist_0201$score_average + 0.001)
hist_0301$score_log2 <- log2(hist_0301$score_average + 0.001)
hist_1101$score_log2 <- log2(hist_1101$score_average + 0.001)
hist_2402$score_log2 <- log2(hist_2402$score_average + 0.001)

# plot histogram with density curve
gg_0101 = ggplot(hist_0101, aes(x = score_log2, fill = type)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*01:01)", y = "Density") +
  scale_fill_manual(values = c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590")) +
  theme_bw()
gg_0201 = ggplot(hist_0201, aes(x = score_log2, fill = type)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*02:01)", y = "Density") +
  scale_fill_manual(values = c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590")) +
  theme_bw()
gg_0301 = ggplot(hist_0301, aes(x = score_log2, fill = type)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*03:01)", y = "Density") +
  scale_fill_manual(values = c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590")) +
  theme_bw()
gg_1101 = ggplot(hist_1101, aes(x = score_log2, fill = type)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*11:01)", y = "Density") +
  scale_fill_manual(values = c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590")) +
  theme_bw()
gg_2402 = ggplot(hist_2402, aes(x = score_log2, fill = type)) +
  geom_density(alpha = 0.5) +
  labs(x = "log2(Average presentation score) (HLA-A*24:02)", y = "Density") +
  scale_fill_manual(values = c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590")) +
  theme_bw()

grid.arrange(gg_0101, gg_0201, gg_0301, gg_1101, gg_2402, ncol = 1)

