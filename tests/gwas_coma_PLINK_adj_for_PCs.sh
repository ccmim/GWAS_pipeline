#!/bin/bash

# Move to the root directory of this repository
cd $(git rev-parse --show-toplevel)

EXPERIMENT="2020-09-11_02-13-41"

python main.py \
  --yaml_config_file config_files/config_coma.yaml \
  --covariates config_files/covariates/std_covariates_10PCs_and_LVEDV.yaml \
  --coma_experiment ${EXPERIMENT} \
  --gwas_file tests/examples/{experiment}/{suffix}/GWAS__{{phenotype}}__{suffix} \
  --suffix "test_{covariates_config}__{sample_white_lists}__{quality_control}" \
  --gwas_software plink \
  --chromosomes 1-22 
