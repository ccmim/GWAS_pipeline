library(tidyverse)

SNPS_LOCATIONS = "~/01_repos/CardiacGWAS/results/snps_for_biomart__one_per_region.txt"
snps_df <- read.csv(SNPS_LOCATIONS)

gwas_info <- ieugwasr::gwasinfo()
gwas_info %>% filter(!grepl("ENSG", trait)) %>% View

cardiac_gwas_ids <- gwas_info %>% filter(grepl("[Hh]eart|[Cc]ardi|ECG", trait)) %>% .$id

# phewas <- ieugwasr::phewas(variants=snps_df$SNP[1])

# assocs <- lapply(snps_df$SNP, function(snp) ieugwasr::phewas(variants=snp))
# names(assocs) <- snps_df$SNP
# saveRDS(assocs, "phewas_associations_selected_snps.rds")
assocs <- readRDS(file = "phewas_associations_selected_snps.rds")


lapply(
  assocs, 
  function(df) tryCatch(expr = df %>% dplyr::filter(id %in% cardiac_gwas_ids), error=function(e) data.frame()
)