#!/bin/bash

cd $(git rev-parse --show-toplevel)

PARTITIONS=("LV" "RV")
SCALING_=("scaled" "non_scaled")

for PARTITION in ${PARTITIONS[@]}; do
  for SCALING in ${SCALING_[@]}; do
    Rscript analysis/gwas_analysis_new.R --gwas_folder PCA__${PARTITION}__${SCALING}__5000_samples_adjusted --cache_rds --qqplot_pooled
  done
done
