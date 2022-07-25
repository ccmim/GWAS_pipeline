import os 

REGIONS_FILE = "data/ld_indep_regions/fourier_ls-all_EUR_hg19.bed"
INPUT_GENOTYPE_FILE_PATTERN = "data/datasets/imputed/ukb22828_c{chromosome}_b0_v3.bgen"
INDEX_FILE_PATTERN = "data/datasets/imputed/ukb_imp_chr{chromosome}_v3.bgen.bgi"
TEMP_GENOTYPE_FILE_PATTERN = "data/transforms/genotypes_by_region/ukb_chr{chromosome_adjusted}_{start_pos}-{end_pos}_all.bgen"
OUTPUT_GENOTYPE_FILE_PATTERN = "data/transforms/genotypes_by_region/ukb_chr{chromosome_adjusted}_{start_pos}-{end_pos}.bgen"
INCLUDE_SAMPLES = "data/ids_list/CMR_GBR.txt"
INCLUDE_SNP_FILES = [os.path.join("data/transforms/snps_files/by_chromosome", x) for x in ["ukb_chr{chromosome}__maf_gt_0.005_info_gt_0.3.txt"]]
REDUCED_SNP_FILE_PATTERN = "data/transforms/snps_files/by_region/ukb_chr{chromosome}_{start_pos}-{end_pos}__maf_gt_0.005_info_gt_0.3.txt"
SAMPLES_FILE = "data/datasets/imputed/ukb22828_c1_b0_v3_s487202.sample"

SNPS_STATS_FILE_PATTERN = OUTPUT_GENOTYPE_FILE_PATTERN[:-5] + "_snps_stats.txt"
SNPS_LIST_FILE_PATTERN = OUTPUT_GENOTYPE_FILE_PATTERN[:-5] + "_snps_passing_QC.txt"
FINAL_BGEN_FILE_PATTERN = OUTPUT_GENOTYPE_FILE_PATTERN[:-5] + "_snps_passing_QC.bgen"
