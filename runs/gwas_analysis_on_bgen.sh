#!/bin/bash

cd $(git rev-parse --show-toplevel)

EXPERIMENT=$1

Rscript analysis/gwas_analysis.R --output_folder output/coma \
--gwas_folder ${EXPERIMENT}/BGEN \
--gwas_pattern "GWAS__{phenotype}__std_covariates__GBR__BGEN__qc" \
--phenotypes z0 z1 z2 z3 z4 z5 z6 z7
