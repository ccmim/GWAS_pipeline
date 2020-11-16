library(tidyverse)
library(qqman)
library(glue)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))

phenotypes <- c("EpicardialFatVolume", "EFV_to_LVEDV_plus_RVEDV_ratio", "EFV_to_LVESV_plus_RVESV_ratio", "EFV_to_mean_LVV_plus_RVV_ratio")
phenotypes <- c(phenotypes, "EFV_to_total_LVEDV_plus_RVEDV_ratio", "EFV_to_total_LVESV_plus_RVESV_ratio", "EFV_to_mean_total_LVV_plus_RVV_ratio")

phenotypes_ <- character()
for(pheno in phenotypes) {
  phenotypes_ <- c(phenotypes_, paste0(pheno,"_adj"))
}
phenotypes <- phenotypes_ # c(phenotypes, phenotypes_)

output_dir <- "output/EFT"
gwas_fp <- file.path(output_dir, "GWAS__{pheno}_GBR.tsv")
gwas_fp_rds <- file.path(output_dir, "GWAS__{pheno}_GBR.rds")
manhattan_fp <- file.path(output_dir, "GWAS__{pheno}__manhattan.png")
qqplot_fp <- file.path(output_dir, "GWAS__{pheno}__QQ-plot.png")

# params_df <- read.csv("~/data/coma/run_parameters.csv", header=TRUE) %>% filter(run_id %in% run_ids)
color1 <- "mediumblue"
color2 <- "lightblue2"

for (pheno in phenotypes) {
  
  if(file.exists(glue(gwas_fp_rds))) {
    gwas_f <- glue(gwas_fp_rds)
    gwas_df <- readRDS(gwas_f)
  } else {
    gwas_f <- glue(gwas_fp)
    gwas_df <- read_tsv(gwas_f, col_names = TRUE)
    saveRDS(gwas_df, glue(gwas_fp_rds))
  }
  
  print(glue("  {pheno}..."))
  
  gwas_df <- gwas_df %>% filter(!is.na(P) & P != 0)
  # pvals <- c(pvals, gwas_df$P)
  # nz <- params_df[params_df$run_id==run_id, "nz"]
  # chamber <- params_df[params_df$run_id==run_id, "partition"]
  
  plot_title <- glue("{pheno}")
  
  # MANHATTAN PLOT
  png(glue(manhattan_fp), res=100, width = 2000, height = 1000)
  pp <- qqman::manhattan(gwas_df, main=plot_title, col = c(color1, color2))
  print(pp)
  dev.off()
  
  # Q-Q PLOT
  png(glue(qqplot_fp), res=100, width = 1000, height = 1000)
  pp <- qqman::qq(gwas_df$P, main=plot_title, col = "blue4")
  print(pp)
  dev.off()
}