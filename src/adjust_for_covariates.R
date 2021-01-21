# install.packages("argparse")
library(argparse)
library(yaml)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))
source("analysis/helperFunctions.R")

library(tidyverse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("-i", "--input_file", required=TRUE)
parser$add_argument("-o", "--output_file", required=TRUE)
parser$add_argument("-c", "--covariates_file", default="config_files/covariates/std_covariates_and_z5.yaml")
parser$add_argument("--covariates", nargs="+", default=c("X50", "X4079", "X4080", "X21001", "X21003", "X31"))
parser$add_argument("--impute_with_mean_for", nargs="+", default=c("X4079", "X4080", "z5"))
parser$add_argument("--samples_white_list", nargs="+", default="data/ids_list/cmr_british_ids.txt")
parser$add_argument("--samples_black_list", default=NULL)
parser$add_argument("--phenotypes", default=NULL, nargs="+")
parser$add_argument("--phenotypes_black_list", default=NULL, nargs="+")
parser$add_argument("--keep_non_cov_adj", action="store_true", default=FALSE, help="If set, keep non-covariate adjusted phenotypes in the output file.")
args <- parser$parse_args()

## AUXILIARY FUNCTIONS

get_sample_list <- function(white_lists=NULL, black_lists=NULL) {
  
  if (is.null(white_lists)) {
    wl <- read.delim("data/ids_list/cmr_british_ids.txt")[,1]
  } else if (length(white_lists) == 1) {
    wl <- read.delim(white_lists)[,1]
  } else {
    wl <- Reduce(intersect, sapply(white_lists, function(file) read.delim(file)[,1]))  
  }
  
  if (is.null(black_lists)) {
    bl <- character()
  } else if (length(black_lists)  == 1) {
    bl <- read.delim(black_lists)[,1]
  } else {
    bl <- Reduce(union, sapply(black_lists, function(file) read.delim(file)[,1]))
  }
  
  setdiff(wl, bl)
}

impute_na <- function(df, columns, method="mean") {
  # TODO: raise warning
  columns <- intersect(columns, colnames(df))
  for (column in columns) {
    if (method == "mean") {
      df[is.na(df[, column]), column] <- mean(df[, column], na.rm = TRUE)
    }
  }
  df
}

mean_across_visits <- function(df, columns, columns_to_reduce, id_column="ID") {
  columns_to_reduce <- sapply(strsplit(columns_to_reduce, "\\."), function(x) x[[1]])
  for (column_to_reduce in columns_to_reduce) {
    # new_colname <- glue::glue("X{ufc}")
    new_col <- df %>% select(starts_with(column_to_reduce)) %>% rowMeans(na.rm=TRUE)
    df[column_to_reduce] <- new_col
  }
  
  # remove columns for each individual visit to the assessment centre
  df <- df %>% select(id_column, columns)
 
  df
}

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

## Load covariates data
# covariates_df <- read.csv(args$covariates_file)
# covariates <- args$covariates
# 
# covariates_df <- mean_across_visits(covariates_df, covariates)
# 
# # impute columns with the population mean when missing
# covariates_df <- impute_na(covariates_df, args$impute_with_mean_for)
#   
# # remove rows with missing data
# covariates_df <- na.omit(covariates_df)
# covariates <- read_yaml("config_files/covariates/std_covariates_and_z5.yaml")
covariates <- read_yaml(args$covariates_file)

covariates_df <- data.frame()

covariate_names <- vector(length = 0)

# print(covariates)
for(covfile in names(covariates)) {
  df_ <- read.csv(covfile) %>% rename(ID=covariates[[covfile]][[1]][["id"]])
  df_$ID <- as.character(df_$ID)
  covariate_names <- c(covariate_names, unlist(covariates[[covfile]][2:length(covariates[[covfile]])]))
  print(covariate_names)
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

print(str(covariates_df))

covariates_df <- select(covariates_df, c("ID", covariate_names))

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
pheno_df$ID <- as.character(pheno_df$ID)
pheno_df <- pheno_df %>% left_join(covariates_df, by = "ID")

## Obtain covariate-adjusted phenotypes
new_colnames <- paste0(pheno_names, "_adj")

for (i in seq_along(pheno_names)) {
  pheno_df[, new_colnames[i]] <- NA
  new_col <- get_adjusted_phenotype(pheno_df, pheno_names[i], covariate_names)
  print(length(new_col))
  print(nrow(pheno_df))
  
  pheno_df[, new_colnames[i]] <- new_col
}

# Reorder columns
# ID | phenotype_1 | phenotype_2 | ... | covariate_1 | covariate_2 | ... | covariate_1_adj | covariate_2_adj | ...
pheno_df <- cbind(pheno_df %>% select(-covariate_names), pheno_df %>% select(covariate_names))

# Inverse-normalise
pheno_df_ <- pheno_df
pheno_df <- pheno_df_ %>% select(ID) %>% rename(FID=ID) %>% mutate(IID=FID)

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

dir.create(dirname(args$output_file), recursive = TRUE)

# Write output into file
write_delim(
  x = pheno_df,
  args$output_file, 
  col_names = TRUE, delim = "\t", na = "NA"
)