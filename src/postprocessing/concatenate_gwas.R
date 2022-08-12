suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(argparse)
})

concatenate_gwas <- function(experiment, z, results_dir, output_file) {
  
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
  write_tsv(df, paste0(output_file, ".tsv"))
  write_rds(df, paste0(output_file, ".rds"))
}


parser <- argparse::ArgumentParser()
parser$add_argument("--experiment_ids", nargs="+")
parser$add_argument("--run_ids", nargs="+")
parser$add_argument("--phenotypes", nargs="+", default=NULL)
parser$add_argument("--input_results_dir")
parser$add_argument("--output_filename_pattern")
args <- parser$parse_args()

pending_jobs_file = file.path(dirname(args$input_results_dir), "pending_jobs")
while (file.info(pending_jobs_file)$size != 0)  {
  f <- file(pending_jobs_file, open="rb")
  nlines <- 0L
  while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
     nlines <- nlines + sum(chunk == as.raw(10L))
  }
  print(nlines)
  close(f)
  print(glue::glue("Waiting for {nlines} pending jobs..."))
  Sys.sleep(60)
}

for (experiment in args$experiment_ids) {
  for (run in args$run_ids) {
    for (phenotype in args$phenotypes) {
      logging::loginfo(glue("Processing {experiment}/{run}/{phenotype}..."))
      concatenate_gwas(
        experiment, phenotype,
        results_dir = glue(args$input_results_dir), 
        output_file = glue(args$output_filename_pattern)
      )
      logging::loginfo(glue("Finished processing {experiment}/{phenotype}..."))
  }
}
}
