library(tidyverse)

covariate_list <- c("31", "50", "4079", "21001", "21002", "21003")

covariate_df <- read_delim("../data/ukb30777_rodrigo_12112019_selected.csv", delim = ",")

formatted_cov_df <- covariate_df[,1] %>% rename(FID=eid) %>% mutate(IID=FID) %>% as.data.frame()

for (covariate in covariate_list) {
  new_colname <- as.character(glue::glue("X{covariate}")) # R adds "X" to the column names starting with a number
  new_col <- covariate_df %>% select(starts_with(glue::glue("{covariate}"))) %>% rowMeans(na.rm=TRUE)
  formatted_cov_df[new_colname] <- new_col
  
  # replace missing values for gender-specific averages
  if (covariate != "31") {
    gender_mean <- tapply(formatted_cov_df[,new_colname], factor(formatted_cov_df$X31), mean, na.rm=T)
    formatted_cov_df[formatted_cov_df$X31 == 0 & is.na(formatted_cov_df[,new_colname]), new_colname] <- gender_mean['0']
    formatted_cov_df[formatted_cov_df$X31 == 1 & is.na(formatted_cov_df[,new_colname]), new_colname] <- gender_mean['1']
  }
}

covariate_file <- "../data/covariates.tsv"
write_delim(x = formatted_cov_df, covariate_file, col_names = TRUE, delim = "\t", na = "-9")
