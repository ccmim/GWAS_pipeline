library(tidyverse)

concatenate_gwas <- function(experiment, z, results_dir, file_fp) {
  
  setwd(system("git rev-parse --show-toplevel", intern = TRUE))
  
  df_lst <- list()
  
  for (file in list.files(results_dir, full.names = TRUE)) {
      
    if (!grepl(pattern = "txt.gz", file))
      next
    
    df <- read.table(gzfile(file), header=TRUE, sep=" ")
    df_first <- df %>% select(chr, rsid, pos, af, a_0, a_1) %>% rename(CHR=chr, SNP=rsid, BP=pos, AF=af)
    
    df_ <- df %>% select(starts_with(z))
    
    beta_column <- paste0(z, "_beta")
    se_column <- paste0(z, "_se")
    t_column <- paste0(z, "_t")
    log10p_column <- paste0(z, ".log10p")
    
    col_names <- c(beta_column, se_column, t_column, log10p_column)
    
    df_ <- rename_(df_, .dots=setNames(col_names, c("BETA", "SE", "T", "LOG10P")))
    df_ <- mutate(df_, P=10^(-LOG10P))
    df_$LOG10P <- NULL
    df_ <- cbind(df_first, df_)
    
    df_lst <- c(df_lst, list(df_))
  }
  
  df <- bind_rows(df_lst)
  write_tsv(df, file.path(results_dir, paste0(glue::glue(file_fp), ".tsv")))
  write_rds(df, file.path(results_dir, paste0(glue::glue(file_fp), ".rds")))
}

library(argparse)

parser <- argparse::ArgumentParser()
parser$add_argument("--experiments", nargs="+")
parser$add_argument("--z", nargs="+")
parser$add_argument("--results_dir", default="output/coma/{experiment}/BGEN/")
parser$add_argument("--filename_pattern", default="GWAS__{z}__std_covariates__GBR__BGEN__qc")
args <- parser$parse_args()

for (experiment in args$experiments) {
  for (z in args$z) {
    print(glue::glue("Processing {experiment}/{z}..."))
    concatenate_gwas(
      experiment, z,
      results_dir = glue::glue(args$results_dir), 
      file_fp=glue::glue(args$filename_pattern)
    )
    print(glue::glue("Finished processing {experiment}/{z}..."))
  }
}
