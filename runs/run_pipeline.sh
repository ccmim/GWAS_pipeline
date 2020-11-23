#!/bin/bash

# Move to the root directory of this repository
cd $(git rev-parse --show-toplevel)

RUNID=$1
LATENT=data/coma_output/${RUNID}/latent_space.csv
LATENT_ADJ=data/coma_output/${RUNID}/latent_space_adjusted.csv
GWAS_FP=${RUNID}/GWAS__{{phenotype}}_GBR

set -e
# Adjust for covariates 
#Rscript code/adjust_for_covariates.R --input_file $LATENT --output_file $LATENT_ADJ 

# Run GWAS 
#python main.py -c config_files/default_coma_config.yaml --phenotype_file $LATENT_ADJ --gwas_file $GWAS_FP

# Generate plots
Rscript analysis/gwas_analysis_new.R --gwas_folder $RUNID --cache_rds --qqplot_pooled
