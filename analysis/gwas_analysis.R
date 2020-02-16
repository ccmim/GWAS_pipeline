library(tidyverse)
library(qqman)

phenotypes <- read.table("../code/yaml_files/phenotypes.txt", stringsAsFactors = F)[,1]
phenotypes <- phenotypes[!grepl(pattern = "adj", x = phenotypes)]
# phenotypes <- phenotypes[grepl(pattern = "automatic", x = phenotypes)] %>% gsub(pattern = "_automatic", replacement = "", x = .)
methods <- c("automatic", "manual")

params <- list(
  gwas_results_dir="../output/inv_norm_only_white/",
  gwas_results_fp="gwas__white_indiv__invnorm__{pheno}.qassoc",
  gwas_results_adj_fp="gwas__white_indiv__invnorm__{pheno_adj}.qassoc",
  gwas_results_fp_rds="gwas__white_indiv__invnorm__{pheno}.rds",
  gwas_results_adj_fp_rds="gwas__white_indiv__invnorm__{pheno_adj}.rds",
  
  manhattan_dir="../figures/gwas/adj/inv-norm/",
  manhattan_fp="manhattan_{pheno}.png",
  manhattan_adj_fp="manhattan_{pheno_adj}.png"
  # ,
  # qqplot_dir="../figures/gwas/adj_vs_non-adj",
  # qq_comp_fp="qqplot__adj_vs_non-adj__{pheno}__inv-norm.png"
) 

# params <- list(
#   gwas_results_dir="../output/only_white/",
#   gwas_results_fp="gwas__white_indiv__{pheno}.qassoc",
#   gwas_results_adj_fp="gwas__white_indiv__{pheno_adj}.qassoc",
#   
#   manhattan_dir="../figures/gwas/only-white",
#   manhattan_fp="manhattan_{pheno}.png",
#   manhattan_adj_fp="manhattan_{pheno_adj}.png",
#   
#   qqplot_dir="../figures/gwas/adj_vs_non-adj",
#   qq_comp_fp="qqplot__adj_vs_non-adj__{pheno}__only_white.png"
# )

theme_bw(base_size = 30)

for (pheno in phenotypes) {
  
  print(pheno)
  pheno_adj <- paste(pheno, "adj", sep = "_")
  
  filename <- file.path(params$gwas_results_dir, glue::glue(params$gwas_results_fp))
  filename_rds <- file.path(params$gwas_results_dir, glue::glue(params$gwas_results_fp_rds))
  if (!file.exists(filename_rds)) {
    gwas_df <- read.table(file = filename, header=TRUE) %>% filter(!is.na(P))
    saveRDS(gwas_df, filename_rds)
  } else {
    readRDS(filename_rds)
  }
  
  filename <- file.path(params$gwas_results_dir, glue::glue(params$gwas_results_adj_fp))
  filename_rds <- file.path(params$gwas_results_dir, glue::glue(params$gwas_results_adj_fp_rds))
  if (!file.exists(filename_rds)) {
    gwas_adj_df <- read.table(file = filename, header=TRUE) %>% filter(!is.na(P))
    saveRDS(gwas_adj_df, filename_rds)
  } else {
    readRDS(filename_rds)
  }
  
  # df <- data.frame(P=gwas_df$P, P_adj=gwas_adj_df$P)
  
  png(file.path(params$manhattan_dir, glue::glue(params$manhattan_fp)), width = 2000, height = 1000)
  qqman::manhattan(gwas_df)
  dev.off()

  png(file.path(params$manhattan_dir, glue::glue(params$manhattan_adj_fp)), width = 2000, height = 1000)
  qqman::manhattan(gwas_adj_df)
  dev.off()
  
  # png(file.path(params$qqplot_dir, glue::glue(params$qq_comp_fp)), width = 1000, height = 1000)
  # pp <- ggplot(df %>% filter(!is.na(P) & !is.na(P_adj)), aes(-log10(sort(P)), -log10(sort(P_adj))))
  # pp <- pp + geom_point()
  # pp <- pp + theme_bw(base_size = 20)
  # pp <- pp + geom_point(size=0.5)
  # pp <- pp + geom_abline(slope = 1, intercept = 0)
  # pp <- pp + xlab("-log10(p-value): non-adjusted phenotype")
  # pp <- pp + ylab("-log10(p-value): adjusted phenotype")
  # print(pp)
  # dev.off()
}

# names(gwas_lst) <- phenotypes