#!/bin/bash

for CHR in `seq 1 22`; do
  echo "Submitting job for chromosome ${CHR}..."
  qsub -v FILE=../ukb_imp_chr${CHR}_40k_maf_gt_1e-3.bgen -o logs/index_chr${CHR}.out -e logs/index_chr${CHR}.err -N index_chr${CHR} bgenix.sge
  sleep 0.01
done
