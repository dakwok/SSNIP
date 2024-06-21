# SSNIP
The Spatial Splicing-derived Neoantigen Identifier Pipeline (SSNIP) allows for the precise characterization and spatial heterogeneity of neoantigens derived from cancer-specific splicing events (neojunctions). The code available in this repository are the neopeptide identifier stages: 1. Neojunction calling, 2. Neopeptide prediction, 3. Neopeptide presentation prediction, 4. Spatial heterogeneity characterization. The pipeline used in this GitHub repository is formatted for analysis of LGG and GBM; adjust the {PATH} directories and files to fit analyses to custom files or other cancer type analyses.

Sample input files that are not deposited in EGA can be found in _____ with the corresponding output files found in each step's parent folder.

## 1. Neojunction calling
The characterization of neojunctions is first performed by evaluating splice sites in tumor samples derived from TCGA and filtering against splice sites with putative read expression in a normal tissue repository, GTEx (n=9166).

### Step 01. Tumor purity
From meta data files, remove samples with a tumor purity of lower than 60%. Various TCGA publications will include meta data quantifying tumor purity. For our studies, we have used Kahles et al. 2018 (Cancer Cell).

### Step 02. Protein-coding genes
Identify and filter for protein-coding genes with a corresponding GTF file (Homo_sapiens.GRCh37.87.chr.gtf in our study). This will generate a dataframe of protein-encoding genes and their respectie coordinates.

### Step 03. TPM filter > 10
From tumor samples with a tumor purity > 0.60, identify transcripts with TPM values > 10.

### Step 04. Characterize annotated splicing junctions
Next characterize all splicing junctions that are identified at alignment. From STAR aligner's SJ.out.tab files, the splice junction site count data can be extracted. Based on Step 02 and Step 03, select for splice sites that are found in protein-coding genes with a TPM filter > 10.

### Step 05. Extract putative non-annotated splice junctions
Filter out annotated splicing junctions in the corresponding GTF file (sjdbList.fromGTF.out.tab) to identify non-annotated splicing junctions. Note that in this step, different GTF versions will affect the characterization of neojunctions that are identified. At the time of this code being used, our GRCh37.87 GTF sj.out.tab file version was GENCODE v33. 

### Step 06. Prepare splicing junction overlap table
