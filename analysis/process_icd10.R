library(tidyverse)

setwd("~/data/PhD/UKBB/scripts/")

# df <- read.table("Extracted_rodrigo_p.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
df <- read.table("diagnoses.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
icd10 <- as.matrix(df[,2:ncol(df)])
rownames(icd10) <- df[,1]

icd_10_ix <- vector(mode = "list", length = nrow(icd10))
names(icd_10_ix) <- rownames(icd10)

for (i in 1:nrow(icd10)) {  
  icd_10_ix[[i]] <- unname(icd10[i,][grepl(pattern = "I", icd10[i,])])
}

pp <- unique(unlist(icd_10_ix))

indiv_per_disease <- cmr_indiv_per_disease <- vector(mode = "list", length = length(pp))

names(indiv_per_disease) <- pp
indiv_per_disease <- lapply(indiv_per_disease, function(x) unique(x))

names(cmr_indiv_per_disease) <- pp

for (i in 1:nrow(icd10)) {
  for (code in icd_10_ix[[i]])
    indiv_per_disease[[code]] <- c(indiv_per_disease[[code]], rownames(icd10)[i])
    # <- unname(icd10[i,][grepl(pattern = "I", icd10[i,])])
}

for (i in seq_along(indiv_per_disease)) {
  icd10 <- names(indiv_per_disease)[i]
  ids <- indiv_per_disease[[i]]
  output_file <- file.path("../subject_ids/icd10", paste0(icd10, ".txt"))
  writeLines(ids, output_file, sep='\n')  
}

# cmr_subj <- read_lines("~/data/PhD/UKBB/subject_ids/IDs_40k.csv")

# for (i in 1:length(cmr_indiv_per_disease)) {
#   cmr_indiv_per_disease[[i]] <- intersect(indiv_per_disease[[i]], cmr_subj)
# }

# code_mp_ <- read.table("~/data/PhD/UKBB/ICD10_code_mapping.tsv", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
# code_mp <- code_mp_[,2]
# names(code_mp) <- code_mp_[,1]