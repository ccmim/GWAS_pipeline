#install.packages("argparse")
library(argparse)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))
source("analysis/helperFunctions.R")

library(tidyverse)
library(stringr)

parser <- ArgumentParser()
parser$add_argument("-i", "--input_file")
parser$add_argument("-c", "--covariates_file")
parser$add_argument("-o", "--output_file")
parser$add_argument("--covariates")
parser$add_argument("--white_list", default=NULL)
parser$add_argument("--phenotypes", default=NULL, nargs="+")
parser$add_argument("--phenotypes_black_list", default=NULL, nargs="+")
parser$add_argument("--keep_non_inv_norm", action="store_true", default=FALSE)
args <- parser$parse_args()

field_codes <- c(
  "gender"="X31",
  "height"="X50",
  "BMI"="X21001",
  "weight"="X21002",
  "genetic_gender"="X22001",
  "ethnicity"="X21000",
  "age_at_recruitment"="X21003",
  "DBP"="X4079", #diastolic blood pressure
  "SBP"="X4080" #systolic blood pressure
)

covariates_df <- read.csv(args$covariates_file)

pheno_df <- read.table(args$input_file, sep = ",", header=TRUE) # %>% select(-subset) # subset has the partition (train, test) to which subject belongs. We don't need that now.

if (!is.null(args$phenotypes)) {
  pheno_names <- args$phenotypes
} else {
  pheno_names <- colnames(pheno_df)
}

if (!is.null(args$phenotypes_black_list)) {
  pheno_names <- pheno_names[!pheno_names %in% args$phenotypes_black_list]
}

pheno_df <- pheno_df %>% left_join(covariates_df, by = c("ID"="eid"))

# Remove X prefix
unique_field_codes <- unlist(sapply(colnames(pheno_df), strsplit, "\\.")) %>%
  .[grepl(pattern = "X", x = .)] %>%
  gsub(pattern = "X", replacement = "") %>% unique

for (ufc in unique_field_codes) {
  new_colname <- glue::glue("X{ufc}")
  new_col <- pheno_df %>% select(starts_with(glue::glue("X{ufc}."))) %>% rowMeans(na.rm=TRUE)
  pheno_df[new_colname] <- new_col
}

pheno_df <- pheno_df %>% select(-matches(".\\..\\."))
pheno_df[is.na(pheno_df[, "X4079"]), "X4079"] <- mean(pheno_df[, "X4079"], na.rm = TRUE)
pheno_df[is.na(pheno_df[, "X4080"]), "X4080"] <- mean(pheno_df[, "X4080"], na.rm = TRUE)
pheno_df <- na.omit(pheno_df)

# Obtain covariate-adjusted phenotypes
for (pheno in pheno_names) {
  new_colname <- paste0(pheno, "_adj")
  pheno_df[, new_colname] <- NA
  fit <- lm(formula=pheno_df[, pheno] ~ pheno_df$X50 + pheno_df$X4079 + pheno_df$X21001 + pheno_df$X21003 + pheno_df$X31)
  pheno_df[, new_colname] <- resid(fit)
}

# Reorder columns
pheno_df <- cbind(pheno_df %>% select(-starts_with("X")), pheno_df %>% select(starts_with("X")))
print(head(pheno_df))

# Retrieve adjusted phenotypes
adj_phenos <- colnames(pheno_df)[grepl(pattern = "adj", colnames(pheno_df))]

# Rank inverse-normalization
kk <- apply(pheno_df %>% select(-ID, -starts_with("X")), 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
pheno_df <- cbind(pheno_df %>% select(ID), as.data.frame(kk))
pheno_df <- pheno_df %>% rename(IID=ID)
print(head(pheno_df))
rm(kk)

filename <- args$output_file
write_delim(x = pheno_df, filename, col_names = TRUE, delim = ",", na = "NA")