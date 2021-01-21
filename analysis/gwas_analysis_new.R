library(tidyverse)
library(qqman)
library(glue)

library(argparse)

setwd(system("git rev-parse --show-toplevel", intern = TRUE))

parser <- ArgumentParser()
parser$add_argument("--output_folder", default="output/coma")
parser$add_argument("--gwas_folder", nargs="+")
parser$add_argument("--gwas_pattern", default="GWAS__{phenotype}", help="File pattern including the \"{phenotype}\" field.")
parser$add_argument("--phenotypes", default=sapply(0:15, function(x) paste0("z", x)))
# parser$add_argument("--phenotypes", nargs="+", default=sapply(0:15, function(x) paste0("z", x, "_adj")) )

parser$add_argument("--title", default=FALSE, action="store_true")
parser$add_argument("--cache_rds", action="store_true", default=FALSE)
parser$add_argument("--qqplot_pooled", action="store_true", default=FALSE)
parser$add_argument("--color_odd_chr", default="mediumblue")
parser$add_argument("--color_even_chr", default="lightblue")
args <- parser$parse_args()

# gwas_fp <- file.path(output_dir, "GWAS__{pheno}_GBR.tsv")
output_dir <- file.path(args$output_folder, "{args$gwas_folder}")
# text files
gwas_fp <- file.path(output_dir, paste0(args$gwas_pattern, ".tsv"))
gwas_fp_rds <- file.path(output_dir, paste0(args$gwas_pattern, ".rds"))
gwas_summary_fp <- file.path(output_dir, paste0(args$gwas_pattern, "__regionwise_summary.tsv"))

print(gwas_fp)
print(gwas_fp_rds)
print(gwas_summary_fp)

# figures
figs_dir <- file.path(output_dir, "figures")
manhattan_fp <- file.path(figs_dir, paste0(args$gwas_pattern, "__manhattan.png"))
qqplot_fp <- file.path(figs_dir, paste0(args$gwas_pattern, "figures", "__QQ-plot.png"))
qqplot_all_fp <- file.path(figs_dir, "GWAS__all__QQ-plot.png")

color1 <- args$color_odd_chr
color2 <- args$color_even_chr

# color1 <- "hotpink4"
# color2 <- "palevioletred2"

# params_df <- read.csv("~/data/coma/run_parameters.csv", header=TRUE) %>% filter(run_id %in% run_ids)

# Process LD-independent genomic regions
regions <- read.delim("data/ld_indep_regions/fourier_ls-all_EUR_hg19.bed", stringsAsFactors = F)
regions <- regions %>% group_by(chr) %>% mutate(id = paste0(row_number()))
regions$chr <-  sub("\\s+$", "", regions$chr)
regions <- regions %>% mutate(id=paste(chr, id, sep = "_")) %>% ungroup()

gwas_files <- character()

for (run_id in args$gwas_folder) {
  
  print(glue("{run_id}..."))
  dir.create(glue::glue(figs_dir))
  pvals <- vector(length = 0)
  
  # gwas_figs_dir <- file.path(run_id, "figs")
  # if (!dir.exists(gwas_figs_dir))
  #   dir.create(gwas_figs_dir)
  
  for (phenotype in args$phenotypes) {
    
    if (!file.exists(glue(gwas_fp)) || file.info(glue(gwas_fp))$size == 0) {
      # print(glue(gwas_fp))
      next
    }
    
    gwas_files <- c(gwas_files, glue(gwas_fp))
    
    if (file.exists(glue(gwas_fp_rds))) {
        gwas_f <- glue(gwas_fp_rds)
        gwas_df <- readRDS(gwas_f)
    } else if (args$cache_rds) {
        gwas_f <- glue(gwas_fp)
        gwas_df <- read_tsv(gwas_f, col_names = TRUE)
        saveRDS(gwas_df, glue(gwas_fp_rds))
    }
    
    print(glue("  {phenotype}..."))
    
    gwas_df <- gwas_df %>% filter(!is.na(P) & P != 0)
    pvals <- c(pvals, gwas_df$P)

    plot_title <- ifelse(args$title, yes = glue("{run_id} - {phenotype}"), no = "")
    
    # MANHATTAN PLOT
    png(glue(manhattan_fp), res=100, width = 3000, height = 1000)
    pp <- qqman::manhattan(gwas_df, main=plot_title, cex = 1.5, cex.lab = 2, cex.axis = 2, suggestiveline = F, col = c(color1, color2), ylim=c(0,20))
    print(pp)
    dev.off()
    
    # Q-Q PLOT
    png(glue(qqplot_fp), res=100, width = 1000, height = 1000)
    pp <- qqman::qq(gwas_df$P, main=plot_title, cex.axis=2, col = "blue4")
    print(pp)
    dev.off()
    
    gwas_list <- list()
    
    # Produce region-wise summary (best p-value per region)
    for(chr_ in 1:22) {
      reduced_regions <- regions %>% filter(chr==paste0("chr", chr_))
      reduced_gwas <- gwas_df %>% filter(CHR==chr_)
      reduced_gwas <- reduced_gwas %>% mutate(region=cut(BP, breaks=reduced_regions$start, labels=head(reduced_regions$id, -1)))
      gwas_list <- c(gwas_list, list(reduced_gwas))
    }
    
    best_p_per_region <- bind_rows(gwas_list) %>% filter(!is.na(region)) %>% group_by(region) %>% slice(which.min(P)) %>% ungroup()
    write.csv(best_p_per_region, file = glue(gwas_summary_fp), sep = "\t", quote = FALSE, row.names = FALSE)
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