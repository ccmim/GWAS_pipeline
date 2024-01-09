  # all subjects for which we have genotypes
all_ids <- read.table("../data/calls/ukb11350_cal_chr9_v2_s488282.fam", header=FALSE)[,1]

cmr_ids <- read.table("../data/Cardiac_Function_Indexes_11350.tsv", header=TRUE, sep="\t") %>% .$`Subject.ID`
cmr_ids <- intersect(all_ids, cmr_ids)
  
cmr_manual_ids <- read.table("../data/Cardiac_Function_Indexes_11350.tsv", header=TRUE, sep="\t") %>% 
  filter(!is.na(LVEDV_manual)) %>% .$`Subject.ID`
cmr_manual_ids <- intersect(all_ids, cmr_manual_ids)

cov_df <- read.table("../data/ukb30777_rodrigo_12112019_selected.csv", header = TRUE, sep = ",")

all_ids <- cov_df$eid
cov_df <- cov_df %>% select(starts_with("X21000"))

white_ids <- all_ids[apply(cov_df == 1001 | cov_df == 1002 | cov_df == 1003, MARGIN = 1, FUN = all, na.rm=TRUE)]
white_ids <- intersect(all_ids, white_ids)

british_ids <- all_ids[apply(cov_df == 1001, MARGIN = 1, FUN = all, na.rm=TRUE)]
british_ids <- intersect(all_ids, british_ids)

cmr_white_ids <- intersect(cmr_ids, white_ids)
cmr_british_ids <- intersect(cmr_ids, british_ids)

write.table(data.frame("IID"=cmr_white_ids, "FID"=cmr_white_ids), file="../data/ids_list/cmr_white_ids.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame("IID"=cmr_british_ids, "FID"=cmr_british_ids), file="../data/ids_list/cmr_british_ids.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame("IID"=cmr_manual_ids, "FID"=cmr_manual_ids), file="../data/ids_list/cmr_manual_ids.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(data.frame("IID"=cmr_ids, "FID"=cmr_ids), file="../data/ids_list/cmr_ids.txt", quote = FALSE, sep = "\t", row.names = FALSE)