# install.packages("argparse")
library(argparse)
library(yaml)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))
source("analysis/helperFunctions.R")
source("data/preprocess_ukbb_data.R")

library(tidyverse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("-i", "--input_file", required=TRUE)
parser$add_argument("-o", "--output_file", required=TRUE)
parser$add_argument("-c", "--covariates_file")
parser$add_argument("--covariates", nargs="+", default=c("X50", "X4079", "X4080", "X21001", "X21003", "X31"))
parser$add_argument("--impute_with_mean_for", nargs="+", default=c("X4079", "X4080"))
parser$add_argument("--samples_white_list", nargs="+", default="data/ids_list/cmr_british_ids.txt")
parser$add_argument("--samples_black_list", default=NULL)
parser$add_argument("--phenotypes", default=NULL, nargs="+")
parser$add_argument("--phenotypes_black_list", default=NULL, nargs="+")
parser$add_argument("--gwas_software", default=NULL, nargs="+")
args <- parser$parse_args()

## AUXILIARY FUNCTIONS

get_adjusted_phenotype <- function(df, phenotype, covariates) {
  # Generate formula: phenotype ~ covariate_1 + covariate_2 + ...
  formula_as_text <- "{phenotype} ~ {paste(covariates, collapse = \" + \")}"
  fit <- lm(formula=as.formula(glue::glue(formula_as_text)), data=df, na.action = "na.exclude")
  resid(fit)
}

# Rank inverse-normalization
inverse_normalise <- function(x) 
  qnorm( (rank(x,na.last="keep")-0.5) / sum(!is.na(x)) )

## END AUXILIARY FUNCTIONS

## Load phenotype data
if (!file.exists(args$input_file)) {
  quit(status = 1)
} else {
  pheno_df <- read.table(args$input_file, sep = ",", header=TRUE)
  samples <- get_sample_list(args$samples_white_list, args$samples_black_list)
  pheno_df <- pheno_df %>% filter(ID %in% samples)
}

covariates <- read_yaml(args$covariates_file)
covariates_df <- data.frame()
covariate_names <- vector(length = 0)

for(covfile in names(covariates)) {
  df_ <- read.csv(covfile) %>% rename(ID=covariates[[covfile]][[1]][["id"]])
  df_$ID <- as.character(df_$ID)
  covariate_names <- c(covariate_names, unlist(covariates[[covfile]][2:length(covariates[[covfile]])]))
  # df_ <- select(df_, c("ID", covariate_names))
  
  if (nrow(covariates_df) == 0) {
    covariates_df <- df_  
  } else {
    covariates_df <- inner_join(covariates_df, df_, by="ID")
  }
}
  
covariates_df <- mean_across_visits(covariates_df, covariate_names, colnames(covariates_df)[startsWith(colnames(covariates_df), "X")])
covariates_df <- impute_na(covariates_df, args$impute_with_mean_for)
covariates_df <- na.omit(covariates_df)
covariates_df <- select(covariates_df, c("ID", covariate_names))

if (is.null(args$phenotypes)) {
  pheno_names <- colnames(pheno_df)
  pheno_names <- pheno_names[pheno_names != "ID"]
} else {
  pheno_names <- args$phenotypes
}

if (!is.null(args$phenotypes_black_list)) {
  pheno_names <- pheno_names[!pheno_names %in% args$phenotypes_black_list]
}

# Add covariates to the phenotype file to then perform linear regression
# ID | phenotype_1 | phenotype_2 | ... | covariate_1 | covariate_2 | ...
pheno_df$ID <- as.character(pheno_df$ID)
pheno_df <- pheno_df %>% left_join(covariates_df, by = "ID")

## Obtain covariate-adjusted and inverse-normalised phenotypes
for (i in seq_along(pheno_names)) {
  new_col <- get_adjusted_phenotype(pheno_df, pheno_names[i], covariate_names)
  new_col <- inverse_normalise(new_col)
  pheno_df[, pheno_names[i]] <- new_col
}

if (tolower(args$gwas_software) == "plink") {
  # For compatibility with PLINK
  pheno_df <- pheno_df %>% rename(IID=ID) %>% mutate(FID=IID)
  pheno_df <- pheno_df[,c("FID", "IID", pheno_names)]
} else if (tolower(args$gwas_software) == "bgenie") {
  pheno_df$ID <- NULL
  pheno_df <- pheno_df[samples, pheno_names]
}

# Write output into file
dir.create(dirname(args$output_file), recursive = TRUE)
write_delim(
  x = pheno_df,
  args$output_file,
  col_names = TRUE, delim = "\t", na = "NA"
)