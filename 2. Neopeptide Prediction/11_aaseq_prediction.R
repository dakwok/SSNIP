# Title: "Step 11: Amino-Acid Sequence Prediction"
# July 10, 2023 | Darwin Kwok M.S. | University of California San Francisco

# Purpose: Predict the "template" or "isoform" of each candidate junction derived from Step 10
#          This is based on the following two information/assumptions:
#              1. Genomic coordinate
#              2. Premise that one or two ends of each junction will overlap with that of the canonical junction
#                 In the case of Exon Skipping, both edges will end at either two of the canonical junctions.
#                 In the case of A3 and A5, either of the edges is expected to be located at the canonical junction.

###########################################################################
#  Step 0: Load Packages and Data -----------------------------------------
###########################################################################

# Install Packages if necessary 
# XML seems to only be able to be downloaded using this method
install.packages("XML", type = "binary")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("GenomicRanges")
BiocManager::install("AnnotationHub")
BiocManager::install("Rsamtools")

rm(list = ls(all.names = TRUE))

library(BiocManager)
library(tidyverse) ;
library(XML)
library(ensembldb) ;
library(EnsDb.Hsapiens.v75) ;
library(GenomicRanges) ;
library(AnnotationHub) ;
library(Rsamtools) ;

# library(BSgenome.Hsapiens.UCSC.hg19) ;


# Establish Directories ---------------------------------------------------

directory_00 = "{PATH_TO_EXTERNAL}"
directory_03 = "{PATH_TO_STEP_03}"
directory_04 = "{PATH_"
directory_10 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/10_extract_neojunctions/out"
directory_11 = "/Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/11_aaseq_prediction/out"


# Load Files --------------------------------------------------------------

# From Step 10: PSR Tables for Neojunctions Filtered and Generated from Step 10
setwd(directory_10)
filename_psr.neo     = "2023_0708_psr_neojunctions.tsv"
dataframe_psr.neo     = read_tsv(filename_psr.neo, na = c("", "NA"))      # 241 x 12

  # Extract only the Junction ID's from the neojunction list
  neo_ids = dataframe_psr.neo %>% 
    dplyr::select(junc.id)

  
# From Step 03: GTF Table of Junctions Filtered for 1) Protein-Coding and 2) TPM > 10
setwd(directory_03) ; 
filename_gtf = "GTF_ProteinCoding_FilterTPM10_20200717.tsv"
dataframe_gtf = read_tsv(filename_gtf, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))


# From Step 04: Annotated SJ Junction List of the TCGA Junction List Filtered by Step 03 
setwd(directory_04) ; 
filename_sj = "SJ_List_Filtered_by_GTF_ProteinCoding_ExpressedTranscripts_20200623.tsv" # "annotated" junction list related to the isoforms of the current interest only. 
dataframe_sj = read_tsv(filename_sj, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character()))


# From Step 00: Original Annotated Junctions List
setwd(directory_00) ; 
filename_sj.orig = "sjdbList.fromGTF.out.tab"
dataframe_sj.orig = read_tsv(filename_sj.orig, na = c("", "NA"), col_names = F, col_types = cols(X1 = col_character()))


# Check -------------------------------------------------------------------

neo_ids %>% dim() %>% print() ; # "junc.id" column only
# [1] 249   1
dataframe_gtf %>% dim() %>% print() ; # gtf.prot.code.tx.pass
# [1] 11617     8
dataframe_sj %>% dim() %>% print() ; # sj.ref
# [1] 112110      5
dataframe_sj.orig %>% dim() %>% print() ; # sj.ref.orig
# [1] 344400      4



###########################################################################
#  Step 1: Edit the Loaded Dataframes -------------------------------------
###########################################################################

# Add informational columns for the neojunction IDs (obtained from Step 10)
neo_ids.2 = neo_ids %>% 
  mutate(chr = gsub("chr", "", sapply(strsplit(junc.id, ":"), "[[", 1))) %>%                        # Make a new column "chr" for the chromosome number
  mutate(strand = sapply(strsplit(junc.id, ":"), "[[", 2)) %>%                                      # Make a new column "strand" for the strand (- or +)
  mutate(int.start = as.numeric(gsub("-.*", "", sapply(strsplit(junc.id, ":"), "[[", 3))) + 1) %>%  # Make a new column "int.start" for the starting coordinate of the transcript
  mutate(int.end = as.numeric(gsub(".*-", "", sapply(strsplit(junc.id, ":"), "[[", 3))))            # Make a new column "int.end" for the ending coordinate for the transcript

# head(neo_ids.2)
# junc.id                    chr   strand int.start   int.end
# <chr>                      <chr> <chr>      <dbl>     <dbl>
#   1 chr1:-:1235406-1235888     1     -        1235407   1235888
# 2 chr1:-:153507303-153507677 1     -      153507304 153507677
# 3 chr1:-:156383289-156383875 1     -      156383290 156383875
# 4 chr1:-:156384545-156396233 1     -      156384546 156396233
# 5 chr1:-:156385262-156391346 1     -      156385263 156391346
# 6 chr1:-:156396661-156399167 1     -      156396662 156399167


# Edit the column names of the complete junction list (obtained from Step 00)
dataframe_sj.orig = dataframe_sj.orig %>% 
  mutate(junc.id = paste0("chr", X1, ":", X4, ":", (X2 - 1), "-", X3)) %>%       # Make a new column "junc.id"
  dplyr::select(junc.id, chr = X1, strand = X4, int.start = X2, int.end = X3)    # Rename columns

# head(dataframe_sj.orig)
# <chr>              <chr> <chr>      <dbl>   <dbl>
#   1 chr1:+:12057-12178 1     +          12058   12178
# 2 chr1:+:12227-12594 1     +          12228   12594
# 3 chr1:+:12227-12612 1     +          12228   12612
# 4 chr1:+:12697-12974 1     +          12698   12974
# 5 chr1:+:12721-13220 1     +          12722   13220
# 6 chr1:+:12721-13224 1     +          12722   13224


# refer to gtf and dataframe_sj (prot-code expr-tx only and original) ------------------------------------------------------------

# Iterate through each neojunction from the Step 10 list
for(i in 1:nrow(neo_ids.2)){
  print(i/nrow(neo_ids.2)*100)   # Progress Bar
  
  # 1. Slice out row(i) from the list of neoantigens (Step 10)
  neo_ids.i = neo_ids.2 %>%                 
    dplyr::slice(i) ;       # if i = 1, print -->  # junc.id                chr   strand int.start int.end
                                                   # 1 chr1:-:1235406-1235888 1     -        1235407 1235888
  
  # 2. Assign variables for JUNC.ID, CHR, STRAND, START, END
  JUNC.ID = neo_ids.i %>% pull(junc.id)        # print --> [1] "chr1:-:1235406-1235888"
  CHR = neo_ids.i %>% pull(chr)                # print --> [1] "1" 
  STRAND = neo_ids.i %>% pull(strand)          # print --> [1] "-"
  INT.START = neo_ids.i %>% pull(int.start)    # print --> [1] 1235407
  INT.END = neo_ids.i %>% pull(int.end)        # print --> [1] 1235888
  
  # 3a. Identify whether the neojunction has a canonical junction from the annotated SJ's filtered for protein coding, TPM > 10 (Step 04)
  sj.hit = dataframe_sj %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>%                # Filter for the CHR and STRAND of slice(i)
    dplyr::filter(int.start == INT.START | int.end == INT.END) %>%  # Filter for the START and END of slice(i)
    pull(junc.id) ;                                                 # print --> [1] "chr1:-:1235582-1235888"
  
  # 3b. In the case where the neojunction does not hit any canonical junction, see if they fall in any of the non-protein coding transcripts, TPM < 10 
  # for the situation where the edge(s) of the junction hits those of out-of-inteneo_idst isoforms. 
  sj.orig.hit = dataframe_sj.orig %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>% 
    dplyr::filter(int.start == INT.START | int.end == INT.END) %>% 
    pull(junc.id) ;  
  
  # 4. Based on #3a/b, label the neojunctions as either:
  #    "no.hit", if the junction has neither ends that hit a protein-coding canonical junction nor ends that hit junctions that hit outside the scope of our study
  #    "out.of.scope", if the junction does not hit a protein-coding canonical junction, but has both ends hitting a junction that is non-protein-coding and/or has a TPM < 10
  sj.hit = if(length(sj.hit) == 0 & length(sj.orig.hit) == 0){
    "no.hit"
  }else if(length(sj.hit) == 0 & length(sj.orig.hit) > 0){ 
    "out.of.scope"
  }else{                                                
    sj.hit
  } ;
  
  # 5. Using the GTF file, identify the gene by which the neojunction's ends lie between
  gtf.i = dataframe_gtf %>% 
    dplyr::filter(chr == CHR & strand == STRAND) %>% 
    dplyr::filter(start < INT.START & INT.END < end) ; 
  
  # 6. Add columns for the Canonical ID, Gene Symbol, ENSG, and ENST 
  neo_ids.i = neo_ids.i %>% 
    mutate(canonical = paste(sj.hit, collapse = ";")) %>% # annotated junc.id
    mutate(symbol = paste(gtf.i %>% pull(symbol) %>% unique(), collapse = ";")) %>% # symbol
    mutate(ensg = paste(gtf.i %>% pull(ensg) %>% unique(), collapse = ";")) %>% # ensg
    mutate(enst = paste(gtf.i %>% pull(enst) %>% unique(), collapse = ";")); # enst(s)
  
  # 7. Concatenate rows together for every iteration
  if(i == 1){
    neo_ids.3 = neo_ids.i ;
  }else{
    neo_ids.3 = neo_ids.3 %>% 
      bind_rows(neo_ids.i) ;
  }
}

neo_ids.3 %>% head() %>% print() ;
# junc.id       chr   strand int.start  int.end canonical                                         symbol ensg    enst                                         
# 1 chr1:-:12354… 1     -        1235407   1.24e6 chr1:-:1235582-1235888                            ACAP3  ENSG00… ENST00000354700                              
# 2 chr1:-:15350… 1     -      153507304   1.54e8 chr1:-:153507306-153507677                        S100A6 ENSG00… ENST00000496817;ENST00000368719              
# 3 chr1:-:15638… 1     -      156383290   1.56e8 chr1:-:156377767-156383875;chr1:-:156383289-1563… C1orf… ENSG00… ENST00000608007;ENST00000400991;ENST00000368…
# 4 chr1:-:15638… 1     -      156384546   1.56e8 chr1:-:156384545-156385188;chr1:-:156384545-1563… C1orf… ENSG00… ENST00000400991;ENST00000471156;ENST00000497…
# 5 chr1:-:15638… 1     -      156385263   1.56e8 chr1:-:156374393-156391346;chr1:-:156384545-1563… C1orf… ENSG00… ENST00000400991;ENST00000471156;ENST00000497…
# 6 chr1:-:15639… 1     -      156396662   1.56e8 chr1:-:156396664-156399167                        C1orf… ENSG00… ENST00000400991;ENST00000497824;ENST00000310…

neo_ids.3 %>% nrow() %>% print() ;
# [1] 241

neo_ids.3 %>% dplyr::filter(canonical == "no.hit") %>% nrow() %>% print() ;
# [1] 70
neo_ids.3 %>% dplyr::filter(canonical == "out.of.scope") %>% nrow() %>% print() ;
# [1] 0

neo_ids.3 %>% dplyr::filter(grepl(";", symbol)) %>% nrow() %>% print() ; # those matched with multiple genes
# [1] 3
neo_ids.3 %>% dplyr::filter(grepl(";", symbol)) %>% dplyr::select(junc.id, symbol, enst) %>% print() ;
# # A tibble: 3 x 3
#   junc.id                symbol             enst
#   <chr>                  <chr>              <chr>
# 1 chr5:+:140855713-1408~ PCDHGB2;PCDHGB7;P~ ENST00000522605;ENST00000398594;ENS~
# 2 chr12:-:10586473-1058~ NKG2-E;KLRC2       ENST00000539033;ENST00000381902
# 3 chr18:-:47014927-4701~ RPL17-C18orf32;RP~ ENST00000584895;ENST00000579248;ENS~


# Export the list of neojunctions with their associated symbol, ENSG, ENST
setwd(directory_11) ;
filename_output = "2023_0710_junction_list_passfilter_symbol_ensg_enst.tsv"
write_tsv(neo_ids.3, filename_output, na = "NA", col_names = T, quote_escape = "double") ;





###########################################################################
#  Step 2: Prepare for Amino Acid Prediction ------------------------------
###########################################################################

# Define the ensembl database
database_ensembl = EnsDb.Hsapiens.v75 ; # ensembl database

# Check the organism of study
organism(database_ensembl) ; # check
# [1] "Homo sapiens"

supportedFilters(database_ensembl) ; # check
#                       filter                 field
# 1               EntrezFilter                entrez
# 2              ExonEndFilter              exon_end
# 3               ExonIdFilter               exon_id
# 4             ExonRankFilter             exon_rank
# 5            ExonStartFilter            exon_start
# 6              GRangesFilter                  <NA>
# 7          GeneBiotypeFilter          gene_biotype
# 8              GeneEndFilter              gene_end
# 9               GeneIdFilter               gene_id
# 10            GeneNameFilter             gene_name
# 11           GeneStartFilter            gene_start
# 12            GenenameFilter              genename
# 13           ProtDomIdFilter           prot_dom_id
# 14     ProteinDomainIdFilter     protein_domain_id
# 15 ProteinDomainSourceFilter protein_domain_source
# 16           ProteinIdFilter            protein_id
# 17             SeqNameFilter              seq_name
# 18           SeqStrandFilter            seq_strand
# 19              SymbolFilter                symbol
# 20           TxBiotypeFilter            tx_biotype
# 21               TxEndFilter                tx_end
# 22                TxIdFilter                 tx_id
# 23              TxNameFilter               tx_name
# 24             TxStartFilter              tx_start
# 25           UniprotDbFilter            uniprot_db
# 26             UniprotFilter               uniprot
# 27  UniprotMappingTypeFilter  uniprot_mapping_type

# bsg = BSgenome.Hsapiens.UCSC.hg19 ;
# unique(genome(bsg)) ; # check
# # [1] "hg19"


# Retrive the genomic nucleotide sequences
a.h = AnnotationHub() ;
query(a.h, c("Homo Sapiens", "TwoBit", "hg19")) ;
dna = a.h[["AH13964"]] ; # human data
head(seqlevels(dna)) ; # check 


# loop --------------------------------------------------------------------

setwd(directory_11) ;
filename_output = "2023_0804_log_sink_aaseq_prediction.txt" ; 
# Sink exports the Console output as an output file
sink(filename_output) ;

start.time = proc.time() ; 

neo_ids.4 = NULL ; 
for(i in 1:nrow(neo_ids.3)){
  
  # 0. Progress bar
  print((i / nrow(neo_ids.3))*100) ; 
  
  # 1. Slice out row(i) from the list of neoantigens (Step 10)
  neo_ids.i = neo_ids.3 %>% 
    dplyr::slice(i) ;   # print --> 1 chr1:-:1235406-1235888 1     -        1235407 1235888 chr1:-:1235582-1235888 ACAP3  ENSG00000131584 ENST00000354700
  
  JUNC.ID = neo_ids.i %>% pull(junc.id)                             # print --> [1] "chr1:-:1235406-1235888"
  ENSG = neo_ids.i %>% pull(ensg) %>% strsplit(";") %>% unlist()    # print --> [1] "ENSG00000131584"
  ENSTs = neo_ids.i %>% pull(enst) %>% strsplit(";") %>% unlist()   # print --> [1] "ENST00000354700" isoform(s) of current interest
  CHR = gsub("chr", "", sapply(str_split(JUNC.ID, ":"), "[[", 1))   # print --> [1] "1"
  STRAND =  sapply(str_split(JUNC.ID, ":"), "[[", 2)                # print --> 
  INT.START = as.numeric(gsub("-.*", "", sapply(str_split(JUNC.ID, ":"), "[[", 3))) + 1  # print --> 
  INT.END = as.numeric(gsub(".*-", "",sapply(str_split(JUNC.ID, ":"), "[[", 3))) ;       # print -->
  
  # load exonic info of all the protein-coding transcripts (isoforms) from the ENSG
  exons = exons(
    database_ensembl, 
    columns = c("tx_id", "exon_idx"), 
    filter = list(GeneIdFilter(ENSG), TxBiotypeFilter("protein_coding"))
  ) %>% 
    as_tibble() %>% 
    mutate(seqnames = paste0("chr", seqnames)) 
  
  exons = exons %>% 
    mutate(seq = as.vector(getSeq(dna, GRanges(exons) ) ) ) %>% # add nucleotide seq of each exon
    arrange(tx_id, exon_idx) %>% 
    dplyr::select(- exon_id, - gene_id, - tx_biotype)
  
  
  # check if 1 or 2 of the edge(s) of the junction will match the canonical junc --------
  
  ENSTs.test = exons %>% # tier 1 # retain ENSTs that fulfill both (highly-expneo_idssed / 1 or 2 edge(s) match)
    dplyr::filter(is.element(tx_id, ENSTs) & (end == INT.START - 1 | INT.END + 1 == start) ) %>% 
    pull(tx_id) %>% 
    unique() ; 
  
  if(length(ENSTs.test) == 0){ # tier 2 # retain ENSTs that fulfill either
    ENSTs.test = exons %>% 
      dplyr::filter(is.element(tx_id, ENSTs) | (end == INT.START - 1 | INT.END + 1 == start) ) %>% 
      pull(tx_id) %>% 
      unique() ; 
  }
  
  if(length(ENSTs.test) == 0){ # tier 3 # others
    ENSTs.test = exons %>% 
      pull(tx_id) %>% 
      unique() ; 
  }
  
  if(length(ENSTs.test) == 0){ # break in case of zero-hit
    print( "check required" ) ; next ;
  }
  
  # genomic range-based filter: further validation: exclude the isoforms of which genomic range don't include JUNC.ID
  ENSTs.test = transcripts(
    database_ensembl, 
    columns = c("tx_id"), 
    filter = list(TxIdFilter(ENSTs.test))
  ) %>% 
    as_tibble() %>% 
    dplyr::filter(start < INT.START & INT.END < end) %>% 
    pull(tx_id) ;
  
  
  # transcripts to be excluded (manually added) -----------------------------
  
  ENSTs.exclude = c(
    "ENST00000400991", # "no protein" in ensembl 
    "ENST00000368242", # "no protein" in ensembl
    "ENST00000404436", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000429422", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000423049", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000316851", # the beginning part of aaseq is a bit different from that in Ensembl. 
    "ENST00000084795", # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
    "ENST00000467825"  # "protein-coding" in ensembl, but has abnormal aaseq starting pattern (not with M). excluded.
  )
  
  ENSTs.test = ENSTs.test[!is.element(ENSTs.test, ENSTs.exclude)]; 
  
  
  # loop for each tx --------------------------------------------------------
  
  if(length(ENSTs.test) == 0){
    next ;
  }else{
    for(j in 1:length(ENSTs.test)){
      ENST.j = ENSTs.test[j] ;
      
      exons.j = exons %>% 
        dplyr::filter(tx_id == ENST.j) ; 
      
      
      # prep for AS type --------------------------------------------------------------------
      
      exons.j = exons.j %>% 
        arrange(start) ; 
      
      if(exons.j %>% nrow() == 1){
        exons.introns.j = exons.j %>% 
          mutate(lab = "exon") %>% 
          mutate(idx = "E", exon_idx) %>% 
          dplyr::select(start, end, lab, id = exon_idx, idx) ;
      }else{
        introns.j = exons.j %>% 
          dplyr::select(start, end, exon_idx) %>% 
          dplyr::slice(1:(nrow(exons.j) - 1)) %>% 
          mutate(int.start = end + 1) %>% 
          mutate(int.end = (exons.j$start)[2:nrow(exons.j)] - 1) %>% 
          mutate(int.idx = if(STRAND == "+"){1:(nrow(exons.j) - 1)}else{ (nrow(exons.j) - 1):1} ) %>% 
          dplyr::select(int.start, int.end, int.idx) ;            
        exons.introns.j = bind_rows(
          exons.j %>% 
            mutate(lab = "exon") %>% 
            mutate(idx = paste0("E", exon_idx)) %>% 
            dplyr::select(start, end, lab, id = exon_idx, idx), 
          introns.j %>% 
            mutate(lab = "intron") %>% 
            mutate(idx = paste0("I", int.idx)) %>% 
            dplyr::select(start = int.start, end = int.end, lab, id = int.idx, idx)
        )
      }
      
      
      # AS type -----------------------------------------------------------------
      
      # discriminate 9 different patterns
      
      # 1) ES; exon-skipping
      # - by confirming that both of "INT.START - 1" and "INT.END + 1" in the JUNC.ID match with "end" and "start" of the canonical exons, and that there are some skip(s) of "exon_idx".
      
      # 2) LEFT - GAIN -- A5 for plus-strand and A3 for minus-strand cases
      # - by confirming that "INT.END + 1" of the JUNC.ID match with "end" position of an exon, 
      # and that "INT.START - 1" is located in the intronic region between the exon and the next left exon. 
      
      # 3) LEFT - LOSS -- A5 for plus-strand and A3 for minus-strand cases
      # - by confirming that "INT.END + 1" of the JUNC.ID match with "end" position of an exon, 
      # and that "INT.START - 1" is not located in the intronic region mentioned above. 
      
      # 4) RIGHT - GAIN -- A3 for plus-strand and A5 for minus-strand cases
      # - by confirming that "INT.START - 1" of the JUNC.ID match with "start" position of an exon, 
      # and that "INT.END + 1" is located in the intronic region between the exon and the next right exon. 
      
      # 5) RIGHT - LOSS -- A3 for plus-strand and A5 for minus-strand cases
      # - by confirming that "INT.START - 1" of the JUNC.ID match with "start" position of an exon, 
      # and that "INT.END + 1" is not located in the intronic region mentioned above. 
      
      # 6) NEOJUNCTION WITHIN SINGLE EXONIC INTRON (NOVEL INTRON) -- cryptic junction within a single exon
      # - by confirming that both "INT.START - 1" and "INT.END + 1" are included within a single exon
      
      # 7) NEOJUNCTION IN INTRONIC REGION -- cryptic junction within a single intronic region
      # - by confirming that both "INT.START" and "INT.END" are included within a single intronic region
      
      # 8) THE OTHERS    
      
      # 1) ES
      
      if(
        exons.introns.j %>% 
        dplyr::filter(lab == "exon") %>% 
        dplyr::filter(end + 1 == INT.START | INT.END == start - 1) %>% 
        nrow() == 2
      ){
        if(
          exons.introns.j %>% 
          dplyr::filter(lab == "exon") %>% 
          dplyr::filter(INT.START < start & end < INT.END ) %>% 
          nrow() > 1 # no less than 1 exons are included between INT.START and INT.END
        ){
          TYPE = "ES" ; 
          IDX = exons.introns.j %>% 
            dplyr::filter(lab == "exon") %>% 
            dplyr::filter(end == INT.START - 1 | INT.END + 1 == start) %>% 
            pull(idx) ;
          NOTE = paste0("lt.", IDX[1], "-", "rt.", IDX[2]) ;
        }
        
      }else if( 
        
        # LEFT
        exons.introns.j %>% 
        dplyr::filter(lab == "exon") %>% 
        dplyr::filter(INT.END + 1 == start) %>% 
        nrow() == 1 # only the junction's "RIGHT-SIDE" matches "start" of any exons.
      ){
        
        RT.SIDE = exons.introns.j %>% 
          dplyr::filter(lab == "exon") %>% 
          dplyr::filter(INT.END + 1 == start) ;
        
        LT.SIDE = exons.introns.j %>% 
          dplyr::filter(start <= INT.START & INT.START <= end) ;
        
        if( # 2) LEFT.GAIN
          LT.SIDE %>% pull(lab) == "intron" & 
          ifelse(STRAND == "+", 
                 LT.SIDE %>% pull(id) == RT.SIDE %>% pull(id) - 1, 
                 LT.SIDE %>% pull(id) == RT.SIDE %>% pull(id) 
          ) 
        ){
          TYPE = ifelse(STRAND == "+", "A5.gain", "A3.gain") ;
        }else{ # 3) LEFT.LOSS
          TYPE = ifelse(STRAND == "+", "A5.loss", "A3.loss") ;
        }
        
        NOTE = paste0("lt.", (LT.SIDE %>% pull(idx) %>% tolower()), "-rt.", (RT.SIDE %>% pull(idx)) ) ;
        
      }else if(
        
        # RIGHT
        exons.introns.j %>% 
        dplyr::filter(lab == "exon") %>% 
        dplyr::filter(end == INT.START - 1) %>% 
        nrow() == 1 # only the junction's "LEFT-SIDE" matches "end" of any exons.
      ){
        
        LT.SIDE = exons.introns.j %>% 
          dplyr::filter(lab == "exon") %>% 
          dplyr::filter(end == INT.START - 1) ;
        
        RT.SIDE = exons.introns.j %>% 
          dplyr::filter(start <= INT.END & INT.END <= end) ;
        
        if( # 4) RIGHT.GAIN
          RT.SIDE %>% pull(lab) == "intron" & 
          ifelse(STRAND == "+", 
                 LT.SIDE %>% pull(id) == RT.SIDE %>% pull(id), 
                 LT.SIDE %>% pull(id) == RT.SIDE %>% pull(id) + 1 
          ) 
        ){
          TYPE = ifelse(STRAND == "+", "A3.gain", "A5.gain") ;
        }else{ # 5) RIGHT.LOSS
          TYPE = ifelse(STRAND == "+", "A3.loss", "A5.loss") ;
        }
        
        NOTE = paste0("lt.", (LT.SIDE %>% pull(idx)), "-rt.", (RT.SIDE %>% pull(idx) %>% tolower() ) ) ;
        
      }else if( 
        
        # 6) WITHIN SINGLE EXONIC REGION ;
        exons.introns.j %>% 
        dplyr::filter(lab == "exon") %>% 
        dplyr::filter( (start < INT.START & INT.START < end) & (start < INT.END & INT.END < end) ) %>% 
        nrow() == 1 # only the junction's "LEFT-SIDE" matches "end" of any exons.        
      ){
        TYPE = "JUNC.WITHIN.EXON" ; 
        NOTE = paste0("within.", exons.introns.j %>% dplyr::filter(lab == "exon") %>% dplyr::filter( (start < INT.START & INT.START < end) & (start < INT.END & INT.END < end) ) %>% pull(idx) %>% tolower() ) ; 
        
      }else if(
        
        # 7) WITHIN SINGLE INTRONIC REGION ;
        exons.introns.j %>% 
        dplyr::filter(lab == "intron") %>% 
        dplyr::filter( (start < INT.START & INT.START < end) & (start < INT.END & INT.END < end) ) %>% 
        nrow() == 1 # only the junction's "LEFT-SIDE" matches "end" of any exons.        
      ){
        TYPE = "JUNC.WITHIN.INTRON" ; 
        NOTE = paste0("within.", exons.introns.j %>% dplyr::filter(lab == "intron") %>% dplyr::filter( (start < INT.START & INT.START < end) & (start < INT.END & INT.END < end) ) %>% pull(idx) %>% tolower() ) ;   
        
      }else{
        # 8) THE OTHERS ;
        LT.SIDE = exons.introns.j %>% 
          dplyr::filter(start <= INT.START & INT.START <= end) ;
        
        RT.SIDE = exons.introns.j %>% 
          dplyr::filter(start <= INT.END & INT.END <= end) ;          
        
        TYPE = "OTHERS" ; 
        NOTE = paste0("lt.", (LT.SIDE %>% pull(idx) %>% tolower() ), "-rt.", (RT.SIDE %>% pull(idx) %>% tolower() ) ) ;   
      }
      
      
      # nucleotide sequence of “altered” tx -------------------------------------      
      
      # extract "involved" exons
      
      LT.SIDE = exons.introns.j %>% 
        arrange(start) %>% 
        mutate(nr = 1:nrow(exons.introns.j)) %>% 
        dplyr::filter(start <= INT.START & INT.START <= end) ; 
      
      RT.SIDE = exons.introns.j %>%   
        arrange(start) %>% 
        mutate(nr = 1:nrow(exons.introns.j)) %>% 
        dplyr::filter(start <= INT.END & INT.END <= end) ;       
      
      EXONS.NOT.INVOLVED = exons.introns.j %>% 
        arrange(start) %>% 
        mutate(nr = 1:nrow(exons.introns.j)) %>% 
        dplyr::filter(nr < (LT.SIDE %>% pull(nr)) | (RT.SIDE %>% pull(nr) < nr)) %>% 
        dplyr::filter(lab == "exon") ;
      
      # edit "involved" exons
      
      for(k in 1:2){ # 1: LT; 2: RT
        NEW.EXON = list(LT.SIDE, RT.SIDE)[[k]]  %>% 
          mutate(end = ifelse(k == 1, INT.START - 1, end)) %>% 
          mutate(start = ifelse(k == 2, INT.END + 1, start)) %>% 
          mutate(idx = id + 0.5) ;           
        
        if( NEW.EXON %>% pull(start) > NEW.EXON %>% pull(end) ){
          NEW.EXON = NULL ;
        }else{
          NEW.EXON = GRanges(
            seqnames = paste0("chr", CHR),
            strand = STRAND, 
            ranges = IRanges(
              NEW.EXON %>% pull(start), 
              NEW.EXON %>% pull(end)
            )
          ) %>% 
            as_tibble() %>% 
            mutate(tx_id = ENST.j) %>% 
            mutate(exon_idx = NEW.EXON %>% pull(idx)) %>% 
            mutate(seqnames = as.character(seqnames)) %>% 
            mutate(seq = as.vector(
              getSeq(
                dna, 
                GRanges(
                  seqnames = paste0("chr", CHR),
                  strand = STRAND, 
                  ranges = IRanges(
                    NEW.EXON %>% pull(start), 
                    NEW.EXON %>% pull(end)
                  )
                )
              )
            )
            ) ;
        }
        if(k == 1){
          LT.NEW.EXON = NEW.EXON ;
        }else{
          RT.NEW.EXON = NEW.EXON ;
        }
      }
      
      exons.new = exons.j %>% 
        semi_join(EXONS.NOT.INVOLVED, by = c("exon_idx" = "id")) %>% 
        bind_rows(LT.NEW.EXON) %>% 
        bind_rows(RT.NEW.EXON) ;
      
      if(STRAND == "+"){
        exons.new = exons.new %>% 
          arrange(start) ;
      }else{
        exons.new = exons.new %>% 
          arrange(desc(start)) ;
      } ;
      
      exons.j = exons.j %>% 
        arrange(exon_idx) ;
      
      
      # tx.seq ------------------------------------------------------------------
      
      # transcript dna seq - wildtype and altered
      tx.seq.wt = paste(exons.j$seq, collapse = "") ;
      tx.seq.alt = paste(exons.new$seq, collapse = "") ;
      
      
      # rm. 5'UTR / 3’UTR -----------------------------------------------------------
      
      # remove "untranslated regions" from the above dna seq
      utr5 = fiveUTRsByTranscript(database_ensembl, 
                                  columns = c("tx_id", "exon_idx"), 
                                  filter = list(TxIdFilter(ENST.j))
      ) ;
      if(is.null(utr5)){
        ln.utr5 = 0 ;
      }else{
        ln.utr5 = utr5 %>% 
          as_tibble() %>% 
          pull(width) %>% 
          sum() ; 
      }
      
      utr3 = threeUTRsByTranscript(database_ensembl, 
                                   columns = c("tx_id", "exon_idx"), 
                                   filter = list(TxIdFilter(ENST.j))
      ) ;
      if(is.null(utr3)){
        ln.utr3 = 0 ;
      }else{
        ln.utr3 = utr3 %>% 
          as_tibble() %>% 
          pull(width) %>% 
          sum() ; 
      }
      
      cds.seq.wt = substr(tx.seq.wt, (ln.utr5 + 1), (nchar(tx.seq.wt) - ln.utr3)) ;
      cds.seq.alt = substr(tx.seq.alt, (ln.utr5 + 1), (nchar(tx.seq.alt) - ln.utr3)) ;
      
      
      # translate ------------------------------------------------------------------
      
      aa.seq.wt = as.character( translate( DNAString(cds.seq.wt) ) ) ;
      aa.seq.alt = as.character( translate( DNAString(cds.seq.alt) ) ) ;
      
      
      # remove “*(stop-codon)” in the end ---------------------------------------
      
      aa.seq.wt = ifelse(str_sub(aa.seq.wt, -1, -1) == "*", str_sub(aa.seq.wt, 1, -2), aa.seq.wt) ;
      aa.seq.alt = ifelse(str_sub(aa.seq.alt, -1, -1) == "*", str_sub(aa.seq.alt, 1, -2), aa.seq.alt) ;
      
      
      # i = 224, SELENOW gene ---------------------------------------------------------------
      
      # "Selenocystein (Sec or U)" ref: https://en.wikipedia.org/wiki/Selenocysteine"
      aa.seq.wt = gsub("MALAVRVVYCGA\\*GYKSKYLQLKKKLE", "MALAVRVVYCGAUGYKSKYLQLKKKLE", aa.seq.wt) ;
      aa.seq.alt = gsub("MALAVRVVYCGA\\*GYKSKYLQLKKKLE", "MALAVRVVYCGAUGYKSKYLQLKKKLE", aa.seq.alt) ;
      
      # evaluate some elements --------------------------------------------------------------------
      
      # frame-shift or not 
      FS = ifelse(abs(nchar(cds.seq.wt) - nchar(cds.seq.alt) ) %% 3 == 0, "in-frame", "fs") ;
      # stop-codon creating or not
      SC = ifelse(grepl("\\*", aa.seq.alt), "sc", "no.sc") ;
      
      
      # check -------------------------------------------------------------------
      
      CHECK = if(aa.seq.wt != proteins(database_ensembl, 
                                       columns = c("protein_sequence"), 
                                       filter = TxIdFilter(ENST.j)) %>% 
                 as_tibble() %>% 
                 pull(protein_sequence)
      ){
        # print("check required") ; break ;
        "not.match"
      }else{
        "pass"
      } ;
      
      # check 
      paste0("i=", i, "; j=", j) %>% print(); 
      JUNC.ID %>% print() ;
      ENST.j %>% print() ;
      c(TYPE, NOTE, FS, SC, CHECK) %>% print() ;
      
      print("wt") ;
      exons.j %>% dplyr::select(start, end, strand, exon_idx) %>% anti_join(exons.new, by = "exon_idx") %>% print() ;
      print("alt") ;
      exons.new %>% dplyr::select(start, end, strand, exon_idx) %>% anti_join(exons.j, by = "exon_idx") %>% print();
      
      print("wt") ;      
      aa.seq.wt %>% print() ;
      print("alt") ;
      aa.seq.alt %>% print() ;
      
      
      # prep for output ------------------------------------------------------------
      
      neo_ids.4 = neo_ids.4 %>% 
        bind_rows(
          tibble(
            junc.id = JUNC.ID, 
            enst.model = ENST.j, 
            aa.change = ifelse(nchar(aa.seq.wt) > nchar(aa.seq.alt), "loss", "gain"), 
            aa.seq.wt = aa.seq.wt, 
            aa.seq.alt = aa.seq.alt, 
            cds.seq.wt = cds.seq.wt,
            cds.seq.alt = cds.seq.alt,
            ln.wt = nchar(aa.seq.wt), 
            ln.alt = nchar(aa.seq.alt), 
            ln.diff = ln.alt - ln.wt, 
            fs = FS, 
            sc = SC, 
            type = TYPE, 
            note = NOTE, 
            check = CHECK
          )
        )
    }
  }
}

RUNTIME = proc.time() - start.time ; # sec 
RUNTIME %>% print() ; 

# check -------------------------------------------------------------------

neo_ids.4 %>% dim() %>% print() ; 
neo_ids.4 %>% distinct(junc.id, .keep_all = T) %>% print() ;


# sink - out --------------------------------------------------------------

sink() ;


# edit --------------------------------------------------------------------

neo_ids.5 = neo_ids.3 %>% 
  dplyr::select(junc.id, canonical, symbol, ensg, enst) %>% 
  right_join(neo_ids.4, by = "junc.id") ; 

setwd(directory_11)
filename_out = "2023_0804_res_aa_prediction_confirmed.tsv"
write_tsv(neo_ids.5, filename_out, na = "NA", col_names = T, quote_escape = "double")
