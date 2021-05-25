#!/bin/bash

for CHR in `seq 1 22`; do
  echo "Submitting job for chromosome ${CHR}..."
  qsub -v CHR=$CHR -o logs/filter_chr${CHR}.out -e logs/filter_chr${CHR}.err -N filter_chr${CHR} filter_bgens.sge
  sleep 0.01
done
