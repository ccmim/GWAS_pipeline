suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(glue)
})


get_include_list <- function(samples_to_include=NULL) {  
  
  if (is.null(samples_to_include)) {
    # TODO: change default
    # include_list <- read.delim("data/ids_list/british_ids.txt", sep = reader::get.delim("data/ids_list/british_ids.txt"))[,1]
    # include_list <- read.delim("data/ids_list/british_ids.txt")#, sep = reader::get.delim("data/ids_list/british_ids.txt"))[,1]
    include_list <- read.delim("data/ids_list/32580_british_from_35610.txt") #, sep = reader::get.delim("data/ids_list/british_ids.txt"))[,1]
  } else if (length(samples_to_include) == 1) {
    include_list <- read.delim(samples_to_include)#, sep = reader::get.delim(samples_to_include))[,1]
  } else {
    include_list <- Reduce(intersect, sapply(samples_to_include, function(file) read.delim(file, sep = reader::get.delim(file))[,1]))
  }
  
  include_list
}


get_exclude_list <- function(samples_to_exclude=NULL) {
  if (is.null(samples_to_exclude)) {
    exclude_list <- character()
  } else if (length(samples_to_exclude)  == 1) {
    exclude_list <- read.delim(samples_to_exclude, sep = reader::get.delim(samples_to_exclude))[,1]
  } else {
    exclude_list <- Reduce(union, sapply(samples_to_exclude, function(file) read.delim(file, sep = reader::get.delim(file))[,1]))
  }
  
  exclude_list
}


get_sample_list <- function(samples_to_include=NULL, samples_to_exclude=NULL) {
  include_list <- get_include_list(samples_to_include)$ID
  exclude_list <- get_exclude_list(samples_to_exclude)
  if (length(exclude_list) == 0){
    return(include_list)
  }
  setdiff(include_list, exclude_list)
}


exclude_samples <-  function(pheno_df, samples_to_include, samples_to_exclude, remove_rows=FALSE) {
  samples <- get_sample_list(samples_to_include, samples_to_exclude)
  
  if (remove_rows) {
    pheno_df <- pheno_df[samples,]  
  } else {
    pheno_df[samples,] <- NA
  }
  
  pheno_df %>% tibble::rownames_to_column("ID")
}
