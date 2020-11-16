#!/bin/bash

cd $(git rev-parse --show-toplevel)

for GWAS_FOLDER in `ls output/coma`; do
  Rscript analysis/gwas_analysis_new.R --gwas_folder $GWAS_FOLDER --cache_rds --qqplot_pooled
done
