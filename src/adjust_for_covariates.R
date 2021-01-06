#install.packages("argparse")
library(argparse)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))
source("analysis/helperFunctions.R")

library(tidyverse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("-i", "--input_file", required=TRUE)
parser$add_argument("-o", "--output_file", required=TRUE)
parser$add_argument("-c", "--covariates_file", default="data/covariates.csv")
parser$add_argument("--covariates", nargs="+", default=c("X50", "X4079", "X4080", "X21001", "X21003", "X31"))
parser$add_argument("--impute_with_mean_for", nargs="+", default=c("X4079", "X4080"))
parser$add_argument("--white_list", default=NULL)
parser$add_argument("--phenotypes", default=NULL, nargs="+")
parser$add_argument("--phenotypes_black_list", default=NULL, nargs="+")
parser$add_argument("--keep_non_cov_adj", action="store_true", default=FALSE, help="If set, keep non-covariate adjusted phenotypes")
args <- parser$parse_args()

## AUXILIARY FUNCTIONS

impute_na <- function(df, columns, method="mean") {
  for (column in columns) {
    if (method == "mean") {
      df[is.na(df[, column]), column] <- mean(df[, column], na.rm = TRUE)
    }
  }
  df
}

mean_across_visits <- function(df, columns, id_column="eid") {
  for (column in columns) {
    # new_colname <- glue::glue("X{ufc}")
    new_col <- df %>% select(starts_with(column)) %>% rowMeans(na.rm=TRUE)
    df[column] <- new_col
  }
  
  # remove columns for each individual visit to the assessment centre
  df <- df %>% select(id_column, columns)
 
  df
}

get_adjusted_phenotype <- function(df, phenotype, covariates) {
  # Generate formula: phenotype ~ covariate_1 + covariate_2 + ...
  formula_as_text <- "{phenotype} ~ {paste(covariates, collapse = \" + \")}"
  fit <- lm(formula=as.formula(glue::glue(formula_as_text)), data=df)
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
}

## Load covariates data
covariates_df <- read.csv(args$covariates_file)
covariates <- args$covariates

covariates_df <- mean_across_visits(covariates_df, covariates)

# impute columns with the population mean when missing
covariates_df <- impute_na(covariates_df, args$impute_with_mean_for)
  
# remove rows with missing data
covariates_df <- na.omit(covariates_df)
  
if (is.null(args$phenotypes)) {
  pheno_names <- colnames(pheno_df)
  pheno_names <- pheno_names[pheno_names != "ID"]
  print(pheno_names)
} else {
  pheno_names <- args$phenotypes
}

if (!is.null(args$phenotypes_black_list)) {
  pheno_names <- pheno_names[!pheno_names %in% args$phenotypes_black_list]
}

# Add covariates to the phenotype file to then perform linear regression
# ID | phenotype_1 | phenotype_2 | ... | covariate_1 | covariate_2 | ...
pheno_df <- pheno_df %>% left_join(covariates_df, by = c("ID"="eid"))

## Obtain covariate-adjusted phenotypes
new_colnames <- paste0(pheno_names, "_adj")

for (i in seq_along(pheno_names)) {
  pheno_df[, new_colnames[i]] <- NA
  pheno_df[, new_colnames[i]] <- get_adjusted_phenotype(pheno_df, pheno_names[i], covariates)
}

# Reorder columns
# ID | phenotype_1 | phenotype_2 | ... | covariate_1 | covariate_2 | ... | covariate_1_adj | covariate_2_adj | ...
pheno_df <- cbind(pheno_df %>% select(-covariates), pheno_df %>% select(covariates))

# Inverse-normalise
pheno_df_ <- pheno_df
pheno_df <- pheno_df_ %>% select(ID) %>% rename(IID=ID)

# Keep non-covariate adjusted phenotypes
if (args$keep_non_cov_adj) {
  inv_norm_non_cov_adj_df <- apply(
    pheno_df_ %>% select(pheno_names),
    MARGIN = 2,
    FUN = inverse_normalise
  )
  pheno_df <- cbind(pheno_df, as.data.frame(inv_norm_non_cov_adj_df))
  rm(inv_norm_non_cov_adj_df)
}

inv_norm_cov_adj_df <- apply(
  pheno_df_ %>% select(new_colnames),
  MARGIN = 2,
  FUN = inverse_normalise
)

pheno_df <- cbind(pheno_df, as.data.frame(inv_norm_cov_adj_df))
rm(inv_norm_cov_adj_df)

# Write output into file
write_delim(
  x = pheno_df,
  args$output_file, 
  col_names = TRUE, delim = ",", na = "NA"
)