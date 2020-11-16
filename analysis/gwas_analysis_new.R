library(tidyverse)
library(qqman)
library(glue)

library(argparse)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))

parser <- ArgumentParser()
parser$add_argument("--gwas_folder", nargs="+")
parser$add_argument("--phenotypes", default=sapply(0:15, function(x) paste0("z", x, "_adj")) )
parser$add_argument("--title", default=FALSE, action="store_true")
parser$add_argument("--cache_rds", action="store_true", default=FALSE)
parser$add_argument("--qqplot_pooled", action="store_true", default=FALSE)
args <- parser$parse_args()

run_ids <- args$gwas_folder
    
output_dir <- "output/coma/{run_id}"
gwas_fp <- file.path(output_dir, "GWAS__{pheno}.tsv")
gwas_fp_rds <- file.path(output_dir, "GWAS__{pheno}.rds")
manhattan_fp <- file.path(output_dir, "GWAS__{pheno}__manhattan.png")
qqplot_fp <- file.path(output_dir, "GWAS__{pheno}__QQ-plot.png")
qqplot_all_fp <- file.path(output_dir, "GWAS__all__QQ-plot.png")

color1 <- "mediumblue"
color2 <- "lightblue"

color1 <- "hotpink4"
color2 <- "palevioletred2"


# params_df <- read.csv("~/data/coma/run_parameters.csv", header=TRUE) %>% filter(run_id %in% run_ids)

for (run_id in run_ids) {
  
  print(glue("{run_id}..."))
  
  pvals <- vector(length = 0)
  
  # gwas_figs_dir <- file.path(run_id, "figs")
  # if (!dir.exists(gwas_figs_dir))
  #   dir.create(gwas_figs_dir)
  
  for (pheno in args$phenotypes) {
    
    if (!file.exists(glue(gwas_fp)) || file.info(glue(gwas_fp))$size == 0) {
      # print(glue(gwas_fp))
      next
    }
    
    if (file.exists(glue(gwas_fp_rds))) {
        gwas_f <- glue(gwas_fp_rds)
        gwas_df <- readRDS(gwas_f)
    } else if (args$cache_rds) {
        gwas_f <- glue(gwas_fp)
        gwas_df <- read_tsv(gwas_f, col_names = TRUE)
        saveRDS(gwas_df, glue(gwas_fp_rds))
    }
    
    print(glue("  {pheno}..."))
    
    gwas_df <- gwas_df %>% filter(!is.na(P) & P != 0)
    pvals <- c(pvals, gwas_df$P)

    plot_title <- ifelse(args$title, yes = glue("{run_id} - {pheno}"), no = "")
    
    # MANHATTAN PLOT
    png(glue(manhattan_fp), res=100, width = 3000, height = 1000)
    pp <- qqman::manhattan(gwas_df, main=plot_title, cex = 1.5, cex.lab = 2, cex.axis = 2, suggestiveline = F, col = c(color1, color2), ylim=c(0,15))
    print(pp)
    dev.off()
    
    # Q-Q PLOT
    png(glue(qqplot_fp), res=100, width = 1000, height = 1000)
    pp <- qqman::qq(gwas_df$P, main=plot_title, col = "blue4")
    print(pp)
    dev.off()
  }
  
  if (length(pvals) == 0)
    next
  
  # POOLED PHENOTYPES' Q-Q PLOT
  png(glue(qqplot_all_fp), res=100, width = 800, height = 800)
  pp <- qqman::qq(pvals, main=glue("{run_id}"), col = "blue4")
  print(pp)
  dev.off()
  
}

# names(gwas_lst) <- phenotypes
# df <- data.frame(P=gwas_df$P, P_adj=gwas_adj_df$P)
# if (!file.exists(filename_rds)) {
#   gwas_df <- read.table(file = filename, header=TRUE) %>% filter(!is.na(P))
#   saveRDS(gwas_df, filename_rds)
# } else {
#   readRDS(filename_rds)
# }

# filename <- file.path(params$gwas_results_dir, glue(params$gwas_results_adj_fp))
# filename_rds <- file.path(params$gwas_results_dir, glue(params$gwas_results_adj_fp_rds))
# if (!file.exists(filename_rds)) {
#   gwas_adj_df <- read.table(file = filename, header=TRUE) %>% filter(!is.na(P))
#   saveRDS(gwas_adj_df, filename_rds)
# } else {
#   readRDS(filename_rds)
# }