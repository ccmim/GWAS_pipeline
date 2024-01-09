#!/bin/bash

cd $(git rev-parse --show-toplevel)
MAF_THR=0.005
INFO_THR=0.3

for CHR in `seq 1 22`; do 
  SOURCE_FILE="$(pwd)/data/datasets/imputed/ukb_mfi_chr${CHR}_v3.txt"
  DEST_FILE="$(pwd)/data/transforms/snps_files/by_region/ukb_chr${CHR}__maf_gt_${MAF_THR}_info_gt_${INFO_THR}.txt"
  awk '$6 > '${MAF_THR}' && $8 > '${INFO_THR}' {print $2, $3}' ${SOURCE_FILE} > ${DEST_FILE};
done

