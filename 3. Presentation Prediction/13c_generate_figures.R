# Step 13c: Make histogram 
# July 31, 2022 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Generate a histogram illustrating the HLA-presentation scores for HLAthena

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

rm(list = ls(all.names = TRUE))
library(readxl)
library(tidyverse)
library(ggplot2)
library(extrafont)

# Font Import
font_import()
loadfonts(device = "pdf")

#  Load Directories -------------------------------------------------------
directory_13      = "{PATH_TO_STEP13_OUTPUT}"
directory_figures = "{PATH_TO_STEP13_FIGURES}"

# Load Files --------------------------------------------------------------
setwd(directory_13)
nmer_08 = read_tsv("2023_0812_hlathena_08mer_selected_alleles.tsv", col_names = T)
nmer_09 = read_tsv("2023_0812_hlathena_09mer_selected_alleles.tsv", col_names = T)
nmer_10 = read_tsv("2023_0812_hlathena_10mer_selected_alleles.tsv", col_names = T)
nmer_11 = read_tsv("2023_0812_hlathena_11mer_selected_alleles.tsv", col_names = T)

###########################################################################
#  Step 3: Remove duplicate n-mers ----------------------------------------
###########################################################################
nmer_08_unique = nmer_08[!duplicated(nmer_08[c(1,2,3,4)]),]
nmer_09_unique = nmer_09[!duplicated(nmer_09[c(1,2,3,4)]),]
nmer_10_unique = nmer_10[!duplicated(nmer_10[c(1,2,3,4)]),]
nmer_11_unique = nmer_11[!duplicated(nmer_11[c(1,2,3,4)]),]

nmer_all_unique = rbind(nmer_08_unique, nmer_09_unique, nmer_10_unique, nmer_11_unique)



###########################################################################
#  Step 2: Make stacked histogram with all n-mers and associated scores ---
###########################################################################
for (i in 1:5){
  if (i == 1){plot_i = nmer_08_unique;  title_i = "08-mers (n=12479)";  filename_i = "2023_0812_histogram_all_08mer_n13320.pdf"}
  if (i == 2){plot_i = nmer_09_unique;  title_i = "09-mers (n=12794)";  filename_i = "2023_0812_histogram_all_09mer_n13659.pdf"}
  if (i == 3){plot_i = nmer_10_unique;  title_i = "10-mers (n=13063)"; filename_i  = "2023_0812_histogram_all_10mer_n13945.pdf"}
  if (i == 4){plot_i = nmer_11_unique;  title_i = "11-mers (n=13085)"; filename_i  = "2023_0812_histogram_all_11mer_n14223.pdf"}
  if (i == 5){plot_i = nmer_all_unique; title_i = "All n-mers (n=51412)"; filename_i = "2023_0812_histogram_all_all_nmers_n55147.pdf"}
  
  ggplot(plot_i, aes(x = hlathena_presentation_score, fill = allele)) +
    geom_histogram(bins = 50) + 
    theme_minimal() + 
    scale_fill_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
    scale_color_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
    xlab("HLAthena Presentation Score") +
    ylab("Count") +
    ggtitle(title_i) +
    theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),  # Change the text of the title and adjust to center with hjust
          text = element_text(size = 20,  family="Helvetica"),                                    # Change the text of the legend and the axis  
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +               # Remove the background grids
    labs(fill = "HLA-A allele")
  
  setwd(directory_figures)
  ggsave(filename_i, limitsize = F, width = 8, height = 5)
}

###########################################################################
#  Step 3: Top 10 percentile -----------------------------------------------
###########################################################################
# 1. Generate a list of the top 01 percentile for each n_mer
nmer_08_ordered =  nmer_08_unique[order(nmer_08_unique$hlathena_presentation_score, decreasing = TRUE),]
nmer_09_ordered =  nmer_09_unique[order(nmer_09_unique$hlathena_presentation_score, decreasing = TRUE),]
nmer_10_ordered =  nmer_10_unique[order(nmer_10_unique$hlathena_presentation_score, decreasing = TRUE),]
nmer_11_ordered =  nmer_11_unique[order(nmer_11_unique$hlathena_presentation_score, decreasing = TRUE),]
nmer_all_ordered = nmer_all_unique[order(nmer_all_unique$hlathena_presentation_score, decreasing = TRUE),]

for (i in 1:5){
  if (i == 1){nmer_i = nmer_08_ordered}
  if (i == 2){nmer_i = nmer_09_ordered}
  if (i == 3){nmer_i = nmer_10_ordered}
  if (i == 4){nmer_i = nmer_11_ordered}
  if (i == 5){nmer_i = nmer_all_ordered}
  
  percentile_10 = as.integer(nrow(nmer_i)/10)
  
  if (i == 1){nmer_10percentile_08 = nmer_i[1:percentile_10,]}
  if (i == 2){nmer_10percentile_09 = nmer_i[1:percentile_10,]}
  if (i == 3){nmer_10percentile_10 = nmer_i[1:percentile_10,]}
  if (i == 4){nmer_10percentile_11 = nmer_i[1:percentile_10,]}
  if (i == 5){nmer_10percentile_all = nmer_i[1:percentile_10,]}
}

# 2. Generate a pie chart for the number of n_mers representing each allele in the top 10 percentile
pie_08 = as.data.frame(t(table(nmer_10percentile_08$allele)))[,c(2,3)]
pie_09 = as.data.frame(t(table(nmer_10percentile_09$allele)))[,c(2,3)]
pie_10 = as.data.frame(t(table(nmer_10percentile_10$allele)))[,c(2,3)]
pie_11 = as.data.frame(t(table(nmer_10percentile_11$allele)))[,c(2,3)]

# Generate a custom table for the all n-mers
nmer_10percentile_all_edit = NULL
for (i in 1:nrow(nmer_10percentile_all)){
  row_i = nmer_10percentile_all %>% dplyr::slice(i)
  peptide_i =  row_i %>% pull(peptide)
  LEN = nchar(peptide_i)
  row_i = row_i %>% mutate(len = LEN)
  nmer_10percentile_all_edit = rbind(nmer_10percentile_all_edit, row_i)
}
pie_all_nmer = as.data.frame(t(table(nmer_10percentile_all_edit$len)))[,c(2,3)]
pie_all_hla = as.data.frame(t(table(nmer_10percentile_all_edit$allele)))[,c(2,3)]


for (i in 1:6){
  if (i == 1){pie_i = pie_08; title_i = "8-mers (Top 10 Percentile)";  filename_i = "2023_0812_pie_chart_08mer_top10percentile.pdf"}
  if (i == 2){pie_i = pie_09; title_i = "9-mers (Top 10 Percentile)";  filename_i = "2023_0812_pie_chart_09mer_top10percentile.pdf"}
  if (i == 3){pie_i = pie_10; title_i = "10-mers (Top 10 Percentile)"; filename_i = "2023_0812_pie_chart_10mer_top10percentile.pdf"}
  if (i == 4){pie_i = pie_11; title_i = "11-mers (Top 10 Percentile)"; filename_i = "2023_0812_pie_chart_11mer_top10percentile.pdf"}
  if (i == 5){pie_i = pie_all_hla; title_i =  "All n-mers (Top 10 Percentile)"; filename_i =  "2023_0812_pie_chart_all_nmers_top10percentile_hla_distribution.pdf"}
  if (i == 6){pie_i = pie_all_nmer; title_i = "All n-mers (Top 10 Percentile)"; filename_i =  "2023_0812_pie_chart_all_nmers_top10percentile_len_distribution.pdf"}
  
  if (i < 6){
    ggplot(pie_i, aes(x = "", y = Freq, fill = Var2)) +
      scale_fill_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + # remove background, grid, numeric labels
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica"),) + 
      labs(fill = "HLA-A allele")
    
    setwd(directory_figures)
    ggsave(filename_i, limitsize = F, width = 8, height = 5)
  }
  
  if (i == 6){
    ggplot(pie_i, aes(x = "", y = Freq, fill = Var2)) +
      scale_fill_manual(values=c("#e6c229", "#f17105", "#d11149", "#6610f2")) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + # remove background, grid, numeric labels
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica"),) + 
      labs(fill = "n-mer Length")
    
    setwd(directory_figures)
    ggsave(filename_i, limitsize = F, width = 8, height = 5)
  }
}

nmer_10percentile_all_edit$len<-as.character(nmer_10percentile_all_edit$len)
nmer_10percentile_all_edit_len =nmer_10percentile_all_edit
nmer_10percentile_all_edit_len$len[nmer_10percentile_all_edit_len$len == "8"] <- "08"
nmer_10percentile_all_edit_len$len[nmer_10percentile_all_edit_len$len == "9"] <- "09"
nmer_10percentile_all_edit_len = nmer_10percentile_all_edit_len[order(nmer_10percentile_all_edit_len$len),]

# 3. Generate a histogram chart for the number of n_mers representing each allele in the top 10 percentile
for (i in 1:6){
  if (i == 1){plot_i = nmer_10percentile_08; title_i = "8-mers (Top 10 Percentile)";  filename_i = "2023_0812_histogram_top10percentile_08mer.pdf"}
  if (i == 2){plot_i = nmer_10percentile_09; title_i = "9-mers (Top 10 Percentile)";  filename_i = "2023_0812_histogram_top10percentile_09mer.pdf"}
  if (i == 3){plot_i = nmer_10percentile_10; title_i = "10-mers (Top 10 Percentile)"; filename_i = "2023_0812_histogram_top10percentile_10mer.pdf"}
  if (i == 4){plot_i = nmer_10percentile_11; title_i = "11-mers (Top 10 Percentile)"; filename_i = "2023_0812_histogram_top10percentile_11mer.pdf"}
  if (i == 5){plot_i = nmer_10percentile_all_edit; title_i = "All n-mers (Top 10 Percentile)"; filename_i = "2023_0812_histogram_top10percentile_all_nmers_hla_distribution.pdf"}
  if (i == 6){plot_i = nmer_10percentile_all_edit_len; title_i = "All n-mers (Top 10 Percentile)"; filename_i = "2023_0812_histogram_top10percentile_all_nmers_len_distribution.pdf"}
  
  if (i < 6){
    ggplot(plot_i, aes(x = hlathena_presentation_score, fill = allele)) +
      geom_histogram(bins = 30) + 
      theme_minimal() + 
      scale_fill_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
      scale_color_manual(values=c("#2c9061", "#147ab0", "#cf0e41", "#3d387e", "#dc7320")) + 
      xlab("HLAthena Presentation Score") +
      ylab("Count") +
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      labs(fill = "HLA-A allele")
    
    setwd(directory_figures)
    ggsave(filename_i, limitsize = F, width = 8, height = 5)
  }
  
  if (i == 6){
    ggplot(plot_i, aes(x = hlathena_presentation_score, fill = len)) +
      geom_histogram(bins = 30) + 
      theme_minimal() + 
      scale_fill_manual(values=c("#e6c229", "#f17105", "#d11149", "#6610f2")) + 
      scale_color_manual(values=c("#e6c229", "#f17105", "#d11149", "#6610f2")) + 
      xlab("HLAthena Presentation Score") +
      ylab("Count") +
      ggtitle(title_i) +
      theme(plot.title = element_text(size = 20,  family="Helvetica", face = "bold", hjust = 0.5),
            text = element_text(size = 20,  family="Helvetica"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      labs(fill = "n-mer Length")
    
    setwd(directory_figures)
    ggsave(filename_i, limitsize = F, width = 8, height = 5)
  }
}

# Export files
setwd(directory_13)
