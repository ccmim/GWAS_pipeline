get_sample_list <- function(white_lists=NULL, black_lists=NULL) {
  
  if (is.null(white_lists)) {
    wl <- read.delim("data/ids_list/cmr_british_ids.txt")[,1]
  } else if (length(white_lists) == 1) {
    wl <- read.delim(white_lists)[,1]
  } else {
    wl <- Reduce(intersect, sapply(white_lists, function(file) read.delim(file)[,1]))
  }
  
  if (is.null(black_lists)) {
    bl <- character()
  } else if (length(black_lists)  == 1) {
    bl <- read.delim(black_lists)[,1]
  } else {
    bl <- Reduce(union, sapply(black_lists, function(file) read.delim(file)[,1]))
  }
  
  setdiff(wl, bl)
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
  df <- df %>% select(id_column, columns)
  df
}