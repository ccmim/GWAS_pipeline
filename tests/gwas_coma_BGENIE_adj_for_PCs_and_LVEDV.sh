#!/bin/bash

# Move to the root directory of this repository
cd $(git rev-parse --show-toplevel)

EXPERIMENTS=("2020-09-11_02-13-41" "2020-09-30_12-36-48")
EXPERIMENTS=("2020-09-11_02-13-41")

module load anaconda
source activate gwas

for EXPERIMENT in ${EXPERIMENTS[@]}; do
  python main.py \
    --yaml_config_file config_files/config_coma.yaml \
    --coma_experiment ${EXPERIMENT} \
    --phenotypes z5 \
    --gwas_file "output/coma/{experiment}/{suffix}/BGEN/GWAS__{{phenotype}}__{suffix}" \
    --covariates "config_files/covariates/std_covariates_10PCs_and_LVEDV.yaml" \
    --bgen_sample_file "~/data/ukbb/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3/ukb_imp_chr22_49824534-51243298__hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3.sample" \
    --phenotype_intermediate "data/coma_output/{experiment}/{suffix}/BGEN/latent_space.tsv" \
    --suffix "{covariates_config}__{sample_white_lists}__{quality_control}" \
    --gwas_software bgenie \
    --chromosomes 1-22 
done
