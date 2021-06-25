for VAR in `seq 0 7`; do
  ./SPrediXcan.py \
  --model_db_path /home/home01/scrb/data/genetics/GTEx_v8/prediction_models/eqtl/mashr/mashr_Heart_Left_Ventricle.db \
  --covariance /home/home01/scrb/data/genetics/GTEx_v8/prediction_models/eqtl/mashr/mashr_Heart_Left_Ventricle.txt.gz \
  --gwas_folder /home/home01/scrb/src/GWAS_pipeline/output/coma/2020-09-11_02-13-41/BGEN \
  --gwas_file_pattern "GWAS__z${VAR}.*harmonized.txt.gz" \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --beta_column effect_size \
  --pvalue_column pvalue \
  --output_file results/z${VAR}.csv
done