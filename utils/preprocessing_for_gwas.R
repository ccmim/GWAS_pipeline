# install.packages("reader") # to use reader::get.delim (to infer table delimiter)
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})
 
get_include_list <- function(samples_to_include=NULL) {  
  if (is.null(samples_to_include)) {
    #TODO: change default
    wl <- read.delim("data/ids_list/cmr_british_ids.txt", sep = reader::get.delim("data/ids_list/cmr_british_ids.txt"))[,1]
  } else if (length(samples_to_include) == 1) {
    wl <- read.delim(samples_to_include, sep = reader::get.delim(samples_to_include))[,1]
  } else {
    wl <- Reduce(intersect, sapply(samples_to_include, function(file) read.delim(file, sep = reader::get.delim(file))[,1]))
  }
  wl
}

get_exclude_list <- function(samples_to_exclude=NULL) {
  if (is.null(samples_to_exclude)) {
    bl <- character()
  } else if (length(samples_to_exclude)  == 1) {
    bl <- read.delim(samples_to_exclude, sep = reader::get.delim(samples_to_exclude))[,1]
  } else {
    bl <- Reduce(union, sapply(samples_to_exclude, function(file) read.delim(file, sep = reader::get.delim(file))[,1]))
  }
  bl
}

get_sample_list <- function(samples_to_include=NULL, samples_to_exclude=NULL) {
  include_list <- get_include_list(samples_to_include)
  exclude_list <- get_exclude_list(samples_to_exclude)
  setdiff(include_list, exclude_list)
}


generate_covariates_df <- function(covariates_config_yaml, impute_with_mean_for=NULL) {
  
  # Example of config YAML for covariates:
  # data/covariates.csv:
  #   - id: "eid"
  #   - X50 # Height
  #   - X4079 # DBP
  # data/genomicPCs_unrelated_GBR.tsv:
  #   - id: "ID"
  #   - PC1
  #   - PC2
  
  covariates <- yaml::read_yaml(covariates_config_yaml)
  covariates_df <- data.frame()
  covariate_names <- character()
  
  for (covfile in names(covariates)) {

    df_ <- read.csv(covfile, sep=reader::get.delim(covfile)) 
    df_ <- df_ %>% rename(ID=covariates[[covfile]][[1]][["id"]]) %>% mutate(ID=as.character(ID))
    
    new_covariate_names <- unlist(covariates[[covfile]][2:length(covariates[[covfile]])])
    covariate_names <- c(covariate_names, new_covariate_names)
    
    if (nrow(covariates_df) == 0) {
      # TODO: add logging
      covariates_df <- df_  
    } else {
      covariates_df <- left_join(covariates_df, df_, by="ID")
    }
  }
  
  covariates_df <- mean_across_visits(covariates_df, covariate_names, colnames(covariates_df)[startsWith(colnames(covariates_df), "X")])
  
  if (!is.null(impute_with_mean_for))
    covariates_df <- covariates_df %>% impute_na(impute_with_mean_for)
  
  covariates_df <- covariates_df %>% na.omit
  covariates_df <- covariates_df %>% select(c("ID", all_of(covariate_names)))
  covariates_df
  
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
  df <- df %>% select(all_of(id_column), all_of(columns))
  df
}


read_raw_pheno <- function(pheno_file, pheno_names=NULL, exclude_columns=NULL) {
  
  # pheno_file 
  # pheno_names: vector of column names (phenotypes). If is.null(pheno_names), all columns are used except those in exclude_columns (if any)
  # exclude_columns: only used if is.null(pheno_names)
  
  if (!file.exists(pheno_file)) {
    #TODO: add logging
    quit(status = 1)
  } 
  
  delim <- reader::get.delim(pheno_file)
  pheno_df <- read.table(pheno_file, sep=delim, header=TRUE)
  rownames(pheno_df) <- as.character(pheno_df$ID)
  if (is.null(pheno_names)) {
    pheno_names <- colnames(pheno_df)
    if (!is.null(exclude_columns)) {
      pheno_names <- pheno_names[!pheno_names %in% exclude_columns]
    }
  }
  pheno_df[,pheno_names]
}


exclude_samples <-  function(pheno_df, samples_to_include, samples_to_exclude, remove_rows=FALSE) {
   samples <- get_sample_list(samples_to_include, samples_to_exclude)
   if (remove_rows) {
     pheno_df <- pheno_df[samples,]  
   } else {
     pheno_df[samples,] <- NA
   }
   pheno_df %>% rownames_to_column("ID")
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


create_adj_pheno_df <- function(raw_pheno_df, covariates_df) {

  # Add covariates to the phenotype file to then perform linear regression
  # ID | phenotype_1 | phenotype_2 | ... | covariate_1 | covariate_2 | ...
  pheno_names <- colnames(raw_pheno_df %>% select(-ID))
  covariate_names <- colnames(covariates_df %>% select(-ID))

  adj_pheno_df <- raw_pheno_df
  pheno_and_covar_df <- left_join(raw_pheno_df, covariates_df, by="ID")

  for (i in seq_along(pheno_names)) {
    adj_pheno <- get_adjusted_phenotype(pheno_and_covar_df, pheno_names[i], covariate_names)
    adj_pheno <- inverse_normalise(adj_pheno)
    adj_pheno_df[, pheno_names[i]] <- adj_pheno
  }
  
  adj_pheno_df
}


format_df_for_tool <- function(pheno_df, software="plink") {
  
  pheno_names <- colnames(pheno_df %>% select(-ID))

  if (tolower(software) == "plink") {
    # For compatibility with PLINK
    pheno_df <- pheno_df %>% rename(IID=ID) %>% mutate(FID=IID)
    pheno_df <- pheno_df[,c("FID", "IID", pheno_names)]
  } else if (tolower(software) == "bgenie") {
    pheno_df$ID <- NULL
    pheno_df <- pheno_df[samples, pheno_names]
  }
  pheno_df
}


generate_adj_pheno <- function(pheno_file, pheno_names, exclude_columns, samples_to_include, samples_to_exclude, covariates_config, gwas_software, output_file=NULL) {
  
  #TODO: Add logging
  raw_pheno_df <- read_raw_pheno(pheno_file, pheno_names, exclude_columns)
  raw_pheno_df <- raw_pheno_df %>% exclude_samples(samples_to_include, samples_to_exclude)
  # print(head(raw_pheno_df))
  covariates_df <- generate_covariates_df(covariates_config)
  # print(head(covariates_df))
  adj_pheno_df <- create_adj_pheno_df(raw_pheno_df, covariates_df)


  adj_pheno_df <- format_df_for_tool(adj_pheno_df, gwas_software)
  
  # Write output into file
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    write_delim(adj_pheno_df, output_file, col_names = TRUE, delim = "\t", na = "NA")
  } else {
    adj_pheno_df
  }
}
