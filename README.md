# SSNIP
The Spatial Splicing-derived Neoantigen Identifier Pipeline (SSNIP) allows for the precise characterization and spatial heterogeneity of neoantigens derived from cancer-specific splicing events (neojunctions). The code available in this repository are the neopeptide identifier stages: 1. Neojunction calling, 2. Neopeptide prediction, 3. Neopeptide presentation prediction, 4. Spatial heterogeneity characterization. The pipeline used in this GitHub repository is formatted for analysis of LGG and GBM; adjust the {PATH} directories and files to fit analyses to custom files or other cancer type analyses.

Sample input files that are not deposited in EGA can be found in _____ with the corresponding output files found in each step's parent folder.

## 1. Neojunction calling
The characterization of neojunctions is first performed by evaluating splice sites in tumor samples derived from TCGA and filtering against splice sites with putative read expression in a normal tissue repository, GTEx (n=9166). 
