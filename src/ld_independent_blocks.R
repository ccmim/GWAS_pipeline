
setwd(system("git rev-parse --show-toplevel", intern = TRUE))

regions <- read.delim("data/ld_indep_regions/fourier_ls-all_EUR_hg19.bed", stringsAsFactors = F)
regions <- regions %>% group_by(chr) %>% mutate(id = paste0(row_number()))
regions$chr <-  sub("\\s+$", "", regions$chr)
regions <- regions %>% mutate(id=paste(chr, id, sep = "_")) %>% ungroup()

pheno <- "z0"
gwas <- read.delim(glue::glue("output/coma/2020-10-16_06-28-58/GWAS__{pheno}_adj_GBR.tsv"))
# gwas %>% group_by(CHR) %>% mutate(region=cut(BP, breaks=regions$start, labels=regions$id))

gwas_list <- list()
  
for(chr_ in 1:22) {
  reduced_regions <- regions %>% filter(chr==paste0("chr", chr_))
  reduced_gwas <- gwas %>% filter(CHR==chr_)
  reduced_gwas <- reduced_gwas %>% mutate(region=cut(BP, breaks=reduced_regions$start, labels=head(reduced_regions$id, -1)))
  print(reduced_gwas)
  gwas_list <- c(gwas_list, list(reduced_gwas))
}

best_p_per_region <- bind_rows(gwas_list) %>% filter(!is.na(region)) %>% group_by(region) %>% slice(which.min(P)) %>% ungroup()