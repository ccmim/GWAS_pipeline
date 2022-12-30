library(dplyr)
library(ggplot2)
library(glue)

##########################################################################################

dataset_dir <- "data/datasets"
transforms_dir <- "data/transforms"

gbr_ids_file <- glue("{dataset_dir}/ids_list/british_fid_iid.txt")
cmr_gbr_ids_file <- glue("{dataset_dir}/ids_list/cmr_british_ids.txt")

bed_file_pattern <- "{dataset_dir}/calls/ukb22418_c{chromosome}_b0_v2.bed"
bim_file_pattern <- "{dataset_dir}/calls/ukb_snp_chr{chromosome}_v2.bim"
fam_file_pattern <- "{dataset_dir}/calls/ukb22418_c{chromosome}_b0_v2_s488170.fam"
long_range_LD_file <- glue("{dataset_dir}/long_range_LD_regions.bed")

bfile_filtered <- glue("{transforms_dir}/GenomicPCA/ukb_cal_all_chrs_v2_GBR_indiv_snps_for_PCA")
files_to_merge <- glue("{transforms_dir}/GenomicPCA/files_to_merge.txt")
snps_to_exclude_file <- glue("{transforms_dir}/GenomicPCA/exclude_snps_for_pca.txt")

out_bfile_pattern <- "{transforms_dir}/GenomicPCA/genotypes/ukb_cal_chr{chromosome}_v2_GBR_indiv"
geno_file_for_pca_chr1 <- glue("{transforms_dir}/GenomicPCA/genotypes/ukb_cal_chr1_v2_GBR_indiv")
geno_file_for_pca <- glue("{transforms_dir}/GenomicPCA/ukb_cal_all_chrs_v2_GBR_indiv")

maf_file_pattern <- paste0(out_bfile_pattern, ".frq")
missingness_file_pattern <- paste0(out_bfile_pattern, ".lmiss")
prunein_file_pattern <- paste0(out_bfile_pattern, ".prune.in")
hwe_file_pattern <- paste0(out_bfile_pattern, ".hwe")
prunein_snps_file <- glue("{transforms_dir}/snps_for_pca_prune.in")
genomic_pca_file <- glue("{transforms_dir}/genomicPCs_GBR.tsv")

##########################################################################################

extract_variants_in_range <- function(chromosome, start, end) {
  bim_file <- glue(bim_file_pattern)
  snp_df <- read.table(bim_file, header=FALSE)
  colnames(snp_df) <- c("chromosome", "rsid", "pos_cm", "position", "allele_1", "allele_2")
  snp_df <- snp_df %>% filter(position > start & position < end)
  snp_df$rsid
}

extract_variants_ambiguous_strand <- function(chromosome) {
  bim_file <- glue(bim_file_pattern)
  snp_df <- read.table(bim_file, header=FALSE)
  colnames(snp_df) <- c("chromosome", "rsid", "pos_cm", "position", "allele_1", "allele_2")
  snp_df <- snp_df %>% filter(allele_1 == "C" & allele_2 == "G" | allele_1 == "G" & allele_2 == "C" | allele_1 == "A" & allele_2 == "T" | allele_1 == "T" & allele_2 == "A")
  as.character(snp_df$rsid)
}

get_snps_below_maf_thr <- function(chromosome, maf_threshold=0.025) {
  maf_file <- glue(maf_file_pattern)
  print(maf_file)
  maf_df <- read.table(maf_file, header=TRUE)
  maf_df <- maf_df %>% filter(MAF < maf_threshold)
  as.character(maf_df$SNP)
}

get_snps_above_miss_thr <- function(chromosome, missingness_threshold=0.015) {
  missingness_file <- glue(missingness_file_pattern)
  missingness_df <- read.table(missingness_file, header=TRUE)
  missingness_df <- missingness_df %>% filter(F_MISS > missingness_threshold)
  as.character(missingness_df$SNP)
}

get_snps_below_hwe_p_thr <- function(chromosome, hwe_p_threshold=1e-5) {
  hwe_file <- glue(hwe_file_pattern)
  hwe_df <- read.table(hwe_file, header=TRUE)
  hwe_df <- hwe_df %>% filter(P < hwe_p_threshold)
  as.character(hwe_df$SNP)
}

get_snps_to_exclude <- function(ofile) {
  
  snps_below_maf_thr <- unlist(lapply(1:22, get_snps_below_maf_thr))
  snps_above_miss_thr <- unlist(lapply(1:22, get_snps_above_miss_thr))
  snps_below_hwe_p_thr <- unlist(lapply(1:22, get_snps_below_hwe_p_thr))

  snps_to_exclude <- unique(c(
    ambiguous_strand_snps, 
    snps_below_maf_thr, 
    snps_above_miss_thr, 
    snps_below_hwe_p_thr
  ))

  snps_exclude_df <- as.data.frame(snps_to_exclude)
  write.csv(snps_exclude_df, ofile, quote = FALSE, row.names = FALSE)
  snps_exclude_df

}

##########################################################################################

long_range_ld_df <- read.table(long_range_LD_file, sep = "\t", header=TRUE)
long_range_ld.lst <- long_range_ld_df %>% select(1:3) %>% split(., seq(nrow(.))) # dataframe to list by rows
snps_in_ld <- lapply(long_range_ld.lst, function(x) {x <- as.integer(x); extract_variants_in_range(x[1], x[2], x[3])})
snps_in_ld <- unlist(snps_in_ld)

ambiguous_strand_snps <- lapply(1:22, extract_variants_ambiguous_strand)
ambiguous_strand_snps <- unlist(ambiguous_strand_snps)

maf_command_pattern <- "plink --keep {gbr_ids_file} --freq {bedbimfam} --out {out_bfile}"
missing_command_pattern <- "plink --keep {gbr_ids_file} --missing {bedbimfam} --out {out_bfile}"
hwe_command_pattern <- "plink --keep {gbr_ids_file} --hardy {bedbimfam} --out {out_bfile}"

for (chromosome in 1:22) {
  bedfile <- glue(bed_file_pattern)
  bimfile <- glue(bim_file_pattern)
  famfile <- glue(fam_file_pattern)
  bedbimfam <- glue("--bed {bedfile} --bim {bimfile} --fam {famfile}")
  
  out_bfile <- glue(out_bfile_pattern)
  plink_command_maf <- glue(maf_command_pattern)
  #system(plink_command_maf)
  plink_command_missing <- glue(missing_command_pattern)
  #system(plink_command_missing)
  plink_command_hwe <- glue(hwe_command_pattern)
  #system(plink_command_hwe)
}

snps_exclude_df <- get_snps_to_exclude(snps_to_exclude_file)

plink_ld_prune_command <- "plink --exclude {snps_to_exclude_file} --indep-pairwise 100 10 0.1 {bedbimfam} --out {out_bfile}"

for(chromosome in 1:22) {
  bedfile <- glue(bed_file_pattern)
  bimfile <- glue(bim_file_pattern)
  famfile <- glue(fam_file_pattern)
  bedbimfam <- glue("--bed {bedfile} --bim {bimfile} --fam {famfile}")
  out_bfile <- glue(out_bfile_pattern)

  plink_ld_prune_command <- glue(plink_ld_prune_command)
  #system(plink_ld_prune_command)
}

prunein_snps <- unlist(lapply(1:22, function(chromosome) read.table(glue(prunein_file_pattern))[,1]))
prunein_snps_df <- as.data.frame(prunein_snps)
write.csv(prunein_snps_df, prunein_snps_file, quote = FALSE, row.names = FALSE)

if (!file.exists(files_to_merge)) {
  for (chromosome in 1:22) {
    file <- glue(out_bfile_pattern)
    system(glue("echo {file} >> {files_to_merge}"))
  }
}

merge_command <- "plink --bfile {geno_file_for_pca_chr1} --merge-list {files_to_merge} --make-bed --out {geno_file_for_pca}"

if (!file.exists(geno_file_for_pca)) {
  print(command)
  system(command)
}

extract_id <- function(x) strsplit(x, ":")[[1]][1]

if (!file.exists(genomic_pca_file)) {
  genomic_pca <- flashpcaR::flashpca(geno_file_for_pca, ndim=10)
  pca_proj <- genomic_pca$projection
  pca_proj_df <- as.data.frame(pca_proj)
  colnames(pca_proj_df) <- paste0("PC", 1:ncol(pca_proj_df))
  pca_proj_df$ID <- rownames(pca_proj_df)
  pca_proj_df <- cbind(pca_proj_df %>% select(ID), pca_proj_df %>% select(-ID)) # reorder columns
  write.table(pca_proj_df, file=genomic_pca_file , quote=FALSE, row.names=FALSE, sep = "\t")
} else {
  pca_proj_df <- read.table(genomic_pca_file, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
}

pca_proj_df$ID <- unname(sapply(pca_proj_df$ID, extract_id))
rownames(pca_proj_df) <- pca_proj_df$ID
