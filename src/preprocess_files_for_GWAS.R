# install.packages("argparse")
library(argparse)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))

source("utils/preprocessing_for_gwas_helpers.R")

# This script loads a phenotype file, adjusts for covariates, 
# inverse-normalize the values and excludes samples

generate_adj_pheno <- function(
  pheno_file, 
  pheno_names, 
  exclude_columns, 
  samples_to_include, 
  samples_to_exclude, 
  ukb.sample=NULL,
  covariates_config, 
  gwas_software, 
  output_file=NULL, 
  overwrite_output_flag=FALSE, 
  fit_summaries_file=NULL) {
  
  if ( !is.null(output_file) && file.exists(output_file) ) {
    if (overwrite_output_flag) {
      logging::logwarn("Intermediate phenotype file, located at:\n\t{output_file}\nalready exists and will be overwritten." %>% glue)
    } else {
      logging::logwarn("Intermediate phenotype file, located at:\n\t{output_file}\nalready exists and won't be overwritten. If this isn't what you want, delete the file and run this R script again or pass the --overwrite flag." %>% glue)
      quit(save = "no", status = 0)
    }
  }
  
  logging::loginfo("Loading phenotype file {pheno_file}..." %>% glue)
  raw_pheno_df <- read_raw_pheno(pheno_file, pheno_names, exclude_columns)
  logging::loginfo("Excluding subjects...")
  raw_pheno_df <- raw_pheno_df %>% exclude_samples(samples_to_include, samples_to_exclude)
  # print(head(raw_pheno_df))
  
  logging::loginfo("Loading covariates file...")
  covariates_df <- generate_covariates_df(covariates_config)
  
  logging::loginfo("Generating covariate-adjusted phenotypes...")
  covariate_adj_fit = adj_by_covariates(raw_pheno_df, covariates_df)
  adj_pheno_df <- covariate_adj_fit$adj_pheno_df
  
  if (!is.null(fit_summaries_file)) {
      logging::loginfo(glue::glue("Saving a summary of the fitted model to file {fit_summaries_file}..."))
      fit_summaries_list <- covariate_adj_fit$fit_summaries
      saveRDS(fit_summaries_list, fit_summaries_file)
  }
  
  #create_adj_pheno_df(raw_pheno_df, covariates_df)
  
  adj_pheno_df <- format_df_for_tool(adj_pheno_df, gwas_software, ukb.sample)
  
  # Write output into file
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (tolower(gwas_software) == "plink") {
      logging::loginfo("Creating file of adjusted phenotypes (formatted for Plink)")
      readr::write_delim(adj_pheno_df, output_file, col_names = TRUE, delim = "\t", na = "NA")
    } else if (tolower(gwas_software) == "bgenie") {
      #TODO: support other NA strings
      logging::loginfo("Creating file of adjusted phenotypes (formatted for BGENIE)")
      readr::write_delim(adj_pheno_df, output_file, col_names = TRUE, delim = "\t", na = "-999")
    }
    logging::loginfo("File created successfully at:\n\t\t{output_file}" %>% glue)
  } else {
    # For testing it may be useful to return the dataframe without creating a file
    adj_pheno_df
  }
}



get_args <- function() {
  
  parser <- ArgumentParser()
  
  # Phenotype file
  parser$add_argument("-p", "--phenotype_file", required=TRUE, help="Path to the (original) phenotype file. This and all paths must be either absolute paths or paths relative to this repo's root directory")
  parser$add_argument("--phenotypes", default=NULL, nargs="+")
  parser$add_argument("--columns_to_exclude", default=NULL, nargs="+")
  
  # Samples
  parser$add_argument("--samples_to_include", nargs="+", default="data/datasets/ids_list/31838_unrelated_british_from_35610.txt")
  parser$add_argument("--samples_to_exclude", nargs="+", default=NULL)
  parser$add_argument("--bgen_sample_file", default=NULL)
  
  # Covariates
  parser$add_argument("-c", "--covariates_config_yaml", required=TRUE)
  # parser$add_argument("--impute_with_mean_for", nargs="+", default=c("X4079", "X4080"))
  
  # Output
  parser$add_argument("-o", "--output_file", required=TRUE)
  parser$add_argument("--gwas_software", default="plink", help="Which software to format the phenotype file for. Currently only Plink and BGENIE are supported.")
  parser$add_argument("--fit_summaries_file", default="data/fit_summaries.rds", help="An RDS file containing the summary of the fitted linear models.")
   
  parser$add_argument("--overwrite_output", default=FALSE, action="store_true", help="Flag indicating if this script should be re-run upon finding a previously generated file with the same name as --output_file.")
  
  args <- parser$parse_args()
  args$gwas_software <- tolower(args$gwas_software)
  
  check_args(args)
  args
  
}

check_args <- function(args) {
  if (args$gwas_software == "bgenie" && is.null(args$bgen_sample_file)) {
    logging::logerror("BGEN's sample file has not been provided and is required when running BGENIE. Aborting execution...")
    stop(2)
  }
}


main <- function(args) {
  generate_adj_pheno(
    args$phenotype_file, args$phenotypes, args$columns_to_exclude,
    args$samples_to_include, args$samples_to_exclude, args$bgen_sample_file,
    args$covariates_config_yaml, 
    args$gwas_software,
    args$output_file,
    args$overwrite_output,
    args$fit_summaries_file
  )
}

args <- get_args()
main(args)
