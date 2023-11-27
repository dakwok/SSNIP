# Title:  Step 14: MHCFlurry2.0 Analysis
# Author: Darwin Kwok M.S.
# Date:   July 30, 2022

# Purpose: MHCFlurry 2.0 allows us to assess protein processing likelihoods for the peptides we 
# predicted in the previous steps. This algorithm offers us two experimental predictors,
# both of which are trained on mass spec-identified MHC-ligands:
#    1. Antigen Processing - models MHC allele-independent effects such as proteosomal cleavage
#    2. Presentation - integrates processing predictions with binding affinity predictions to 
#     give a composite "presentation score".

# The following script is adapted and follows the tutorial given by MHCFlurry 2.0's page:
# https://openvax.github.io/mhcflurry/index.html

################################################################
# 0. Install and Activate conda (Terminal)
################################################################

# Downlod the Miniconda installer to the Home directory
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda.sh
bash ~/miniconda.sh -b                 # Install Miniconda quietly, accepting defaults, to the Home directory 
rm ~/miniconda.sh                         # Remove the Miniconda installer from the Home directory
source $HOME/miniconda3/bin/activate      # You will need to first activate conda by calling the activate command by its full system path to add it to your system's PATH environment variable



################################################################
# 0. Install MHCFlurry through conda environment (Terminal)
################################################################

# Install through conda to avoid tensorflow installation problems
conda create -q -n mhcflurry-env python=3.8 tensorflow
source activate mhcflurry-env

pip install mhcflurry        # Install the package
mhcflurry-downloads fetch    # Download their datasets and trained modules



################################################################
# 1. Download MHCFlurry models (Terminal)
################################################################

# Download pre-trained MHCFlurry models with the mhcflurry-downloads tool
# These models are distributed separately from the pip package

# models_class1_presentation is usually what we need as the presentation predictor includes
# a peptide/MHC I binding affinity predictor and an antigen processing predictor
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads path models_class1_presentation   # This gets the path to the downloaded files



################################################################
# 2. Generate Predictions (Terminal)
################################################################

# Generate the .cvs files of the prediction scores for all of the n-mers obtained from Step 11
# This method requires that you generate a CSV input file in order to run mhcflurry -- refer to the 14_mhcflurry2_input_df_generation.R script

# 8-mer
mhcflurry-predict /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_08mer_mhcflurry_input.csv --out /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_08mers_flank_mhcflurry.csv

# 9-mer
mhcflurry-predict /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_09mer_mhcflurry_input.csv --out /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_09mers_flank_mhcflurry.csv

# 10-mer
mhcflurry-predict /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_10mer_mhcflurry_input.csv --out /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_10mers_flank_mhcflurry.csv

# 11-mer
mhcflurry-predict /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_11mer_mhcflurry_input.csv --out /Users/darwinkwok/Desktop/cclc01/costellolab/data5/dkwok/proj_01_altspl/14_mhcflurry2_analysis/out/2023_0812_11mers_flank_mhcflurry.csv
