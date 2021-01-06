#!/bin/bash

cd $(git rev-parse --show-toplevel)

for EXPERIMENT in `ls data/coma_output`; do
  Rscript code/adjust_for_covariates.R \
  -i data/coma_output/${EXPERIMENT}/latent_space.csv \
  -c data/covariates.csv \
  -o data/coma_output/${EXPERIMENT}/latent_space_adjusted.csv \
  --phenotypes_black_list ID subset
done
