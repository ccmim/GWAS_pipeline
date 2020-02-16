library(tidyverse)

c_idx_df <- read.table("../data/Cardiac_Function_Indexes_11350.tsv", sep = "\t", header=T)
metadata_df <- read.table("../data/ukb30777_rahman_23092019.csv", sep = ",")

pheno_coding <- read.table("../data/code_mappings.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
names(pheno_coding) <- c("field","code")