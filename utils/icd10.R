get_subjects_with_icd10 <- function(prefix) {
  dir <- "~/data/PhD/UKBB/subject_ids/icd10"
  file_list <- list.files(dir)
  file_list <- file_list[grepl(prefix, file_list)]
  file_list <- file.path(dir, file_list)
  unique(unlist(sapply(file_list, readLines)))
}