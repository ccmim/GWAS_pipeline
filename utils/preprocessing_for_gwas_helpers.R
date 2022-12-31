options(error=recover)
# install.packages("reader") # to use reader::get.delim (to infer table delimiter)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(glue)
})

source("utils/samples_utils.R")

# Copied from ukbtools
ukb_gen_read_sample <- function (file, col.names = c("id_1", "id_2", "missing"), row.skip = 2) {
  sample <- readr::read_table(file, skip = row.skip, col_names = col.names)
  as.data.frame(sample)
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
    id_colname <- covariates[[covfile]][[1]][["id"]]
    df_ <- df_ %>% rename(ID=all_of(id_colname)) %>% mutate(ID=as.character(ID))
    
    new_covariate_names <- unlist(covariates[[covfile]][2:length(covariates[[covfile]])])
    covariate_names <- c(covariate_names, new_covariate_names)
 
    ifelse( nrow(covariates_df) == 0, covariates_df <- df_, covariates_df <- left_join(covariates_df, df_, by="ID") )   
    
    # if (nrow(covariates_df) == 0) {
    #   # TODO: add logging
    #    covariates_df <- df_  
    # } else {
    #    covariates_df <- left_join(covariates_df, df_, by="ID")
    # }
  }

  # covariates_df <- covariates_df %>% mean_across_visits(., covariate_names, colnames(.)[startsWith(colnames(.), "X")])

  if (!is.null(impute_with_mean_for))
    covariates_df <- covariates_df %>% impute_na(impute_with_mean_for)
  
  # print(head(covariates_df))
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

# Not currently used
mean_across_visits <- function(df, columns, columns_to_reduce, id_column="ID") {
  
  columns_to_reduce <- sapply(
    strsplit(columns_to_reduce, "\\."), function(x) x[[1]]
  )
  
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
    logging::logerror("File {pheno_file} was not found." %>% glue)
    quit(status = 1)
  } 
  
  delim <- reader::get.delim(pheno_file)
  pheno_df <- read.table(pheno_file, sep=delim, header=TRUE)
  if (is.null(pheno_df$ID)) {
    logging::logerror("Column ID seems to be absent in your phenotype file. Aborting...")
  }

  rownames(pheno_df) <- as.character(pheno_df$ID)
  if (is.null(pheno_names)) {
    pheno_names <- colnames(pheno_df)
    if (!is.null(exclude_columns)) {
      pheno_names <- pheno_names[!pheno_names %in% exclude_columns]
    }
  }
  pheno_df[,pheno_names]
}


fit_linear_model <- function(df, phenotype, covariates) {
  
  # Generate formula: phenotype ~ covariate_1 + covariate_2 + ...
  formula_as_text <- "{phenotype} ~ {paste(covariates, collapse = \" + \")}"
  formula_as_text <- glue::glue(formula_as_text)
  formula <- as.formula(formula_as_text)
  #TODO: add logging regarding number of excluded rows due to missing values
  fit <- lm(formula=formula, data=df, na.action = "na.exclude")
  fit
  # fit <- lm(formula=formula, data=df, na.action = "na.omit")
}


# Rank inverse-normalization
inverse_normalise <- function(x) 
    qnorm( (rank(x,na.last="keep")-0.5) / sum(!is.na(x)) )


adj_by_covariates <- function(raw_pheno_df, covariates_df) {

  # Add covariates to the phenotype file to then perform linear regression
  # ID | phenotype_1 | phenotype_2 | ... | covariate_1 | covariate_2 | ...
  pheno_names <- colnames(raw_pheno_df %>% select(-ID))
  covariate_names <- colnames(covariates_df %>% select(-ID))

  adj_pheno_df <- raw_pheno_df
  pheno_and_covar_df <- left_join(raw_pheno_df, covariates_df, by="ID")

  fit_summary_list <- list()
  
  for (i in seq_along(pheno_names)) {
    
    logging::loginfo(glue::glue("Processing phenotype {pheno_names[i]}..."))
    fit <- fit_linear_model(pheno_and_covar_df, pheno_names[i], covariate_names)

    adj_pheno = resid(fit)
    adj_pheno <- inverse_normalise(adj_pheno)
    adj_pheno_df[, pheno_names[i]] <- adj_pheno
    
    fit_summary <- summary(fit)
    
    fit_summary_list <- c(fit_summary_list, list(fit_summary$coefficients))
  }
  
  names(fit_summary_list) <- pheno_names
  
  list("adj_pheno_df"=adj_pheno_df, "fit_summaries"=fit_summary_list)
}


format_df_for_tool <- function(pheno_df, gwas_software="plink", ukb.sample=NULL) {
  
  pheno_names <- colnames(pheno_df %>% select(-ID))
  gwas_software <- tolower(gwas_software)
  
  if (gwas_software == "plink") {
    logging::loginfo("Formatting table for Plink...")
    pheno_df <- pheno_df %>% rename(IID=ID) %>% mutate(FID=IID)
    pheno_df <- pheno_df[, c("FID", "IID", pheno_names)]
  } else if (gwas_software == "bgenie") {
    logging::loginfo("Formatting table for BGENIE...")
    if (is.null(ukb.sample)) {
      logging::logerror("BGEN's sample file has not been provided and is required when running BGENIE. Aborting execution...")
      stop(2)
    }
    suppressMessages({
      sample_df <- ukb_gen_read_sample(file=ukb.sample)
    })
    names(sample_df)[1] <- "ID"
    sample_df <- mutate(sample_df, ID=as.character(ID))
    pheno_df <- sample_df %>% left_join(pheno_df, by = "ID") %>% select("ID", all_of(pheno_names)) # select(.dots = pheno_names)
    pheno_df <- pheno_df %>% arrange(ID)
  }
  
  pheno_df
  
}