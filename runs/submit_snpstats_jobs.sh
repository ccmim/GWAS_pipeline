#!/bin/bash

for CHR in `seq 1 22`; do
  echo "Submitting job for chromosome ${CHR}..."
  qsub -v CHR=$CHR -o logs/chr${CHR}.out -e logs/chr${CHR}.err -N sumstats_chr${CHR} generate_sum_stats.sge 
  sleep 0.01
done
