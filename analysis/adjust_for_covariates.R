# install.packages("ukbtools")
# install.packages("devtools")
# library(devtools)
# library(ukbtools)

repo_dir <- "~/Leeds/doctorado/repos/GWAS_pipeline/analysis/"
setwd(repo_dir)

source("helperFunctions.R")
library(tidyverse)
library(stringr)

ONLY_WHITE <- TRUE
# ONLY_BRITISH <- TRUE

field_codes <- c(
  "gender"="X31",
  "height"="X50",
  "BMI"="X21001",
  "weight"="X21002",
  "genetic_gender"="X22001",
  "ethnicity"="X21000",
  "age_at_recruitment"="X21003",
  "DBP"="X4079",
  "SBP"="X4080"
)

# covariates_df <- read.table("../data/ukb30777_rahman_23092019.csv", sep = ",", header=TRUE)
# codes <- read.table("../data/code_mappings.txt", sep = "\t", header=FALSE)[,c(2,1)]

covariates_df <- read.csv("../data/ukb30777_rodrigo_12112019_selected.csv")
heart_indices <- read.table("../data/Cardiac_Function_Indexes_11350.tsv", sep = "\t", header=TRUE)
heart_indices <- heart_indices %>% left_join(covariates_df, by = c("Subject.ID"="eid"))

unique_field_codes <- unlist(sapply(colnames(heart_indices), strsplit, "\\.")) %>%
  .[grepl(pattern = "X", x = .)] %>%
  gsub(pattern = "X", replacement = "") %>% unique

for (ufc in unique_field_codes) {
  new_colname <- glue::glue("X{ufc}")
  new_col <- heart_indices %>% select(starts_with(glue::glue("X{ufc}."))) %>% rowMeans(na.rm=TRUE)
  heart_indices[new_colname] <- new_col
}

heart_indices <- heart_indices %>% select(-matches(".\\..\\."))
heart_indices[is.na(heart_indices[, "X4079"]), "X4079"] <- mean(heart_indices[, "X4079"], na.rm = TRUE)
heart_indices[is.na(heart_indices[, "X4080"]), "X4080"] <- mean(heart_indices[, "X4080"], na.rm = TRUE)

for (i in 2:19) {
  new_colname <- glue::glue("{colnames(heart_indices)[i]}_adj")
  heart_indices[, new_colname] <- NA
  
  fit <- lm(formula=heart_indices[,i] ~ heart_indices$X50 + heart_indices$X4079 + heart_indices$X21001 + heart_indices$X21003,
            subset=heart_indices$X31==0)
  heart_indices[heart_indices$X31==0 & !is.na(heart_indices[,i]), new_colname] <- resid(fit)
  
  fit <- lm(formula=heart_indices[,i] ~ heart_indices$X50 + heart_indices$X4079 + heart_indices$X21001 + heart_indices$X21003,
            subset=heart_indices$X31==1)
  heart_indices[heart_indices$X31==1 & !is.na(heart_indices[,i]), new_colname] <- resid(fit)
}

heart_indices <- cbind(heart_indices %>% select(-starts_with("X")), heart_indices %>% select(starts_with("X")))

adj_phenos <- colnames(heart_indices)[grepl(pattern = "adj", colnames(heart_indices))]

filename <- "../data/Cardiac_Adj_Function_Indexes_plus_covariates_11350.tsv"
write_delim(x = heart_indices, filename, col_names = TRUE, delim = "\t", na = "NA")

kk <- apply(heart_indices %>% select(-Subject.ID, -starts_with("X")), 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
kk <- cbind(heart_indices %>% select(Subject.ID), as.data.frame(kk))
filename <- "../data/Cardiac_InvNorm_Indexes_Adj_11350.tsv"
write_delim(x = kk, filename, col_names = TRUE, delim = "\t", na = "NA")