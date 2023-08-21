#install.packages("gprofiler2")
library(gprofiler2)
library(tidyverse)
library(biomaRt)

source("~/01_repos/GWAS_pipeline/src/postprocessing/retrieve_gene_positions.R")

SNPS_LOCATIONS = "~/01_repos/CardiacGWAS/results/snps_for_biomart__one_per_region.txt"
WINDOW_SIZE <- 2e6

get_gene_positions <- function() {
  
  selected_attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol", "gene_biotype")
  ensembl <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl", GRCh=37)
  dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  
  pp <- getBM(attributes = selected_attributes, mart = ensembl)
  
  pp <- pp %>% filter(chromosome_name %in% as.character(1:22))
  
  pp <- pp %>% filter(gene_biotype == "protein_coding")
  pp
  
}


get_genes_around_snps <- function(genes_df, snps_df, window_size=2e6) {
  genes_around_gwas_hits <- merge(genes_df, snps_locations, by.x = "chromosome_name", by.y = "CHR")
  genes_around_gwas_hits <- genes_around_gwas_hits %>% mutate(distance_from_tss=abs(BP-start_position))
  genes_around_gwas_hits <- genes_around_gwas_hits %>% filter(abs(BP-start_position) < window_size/2)
  genes_around_gwas_hits
}

get_genes_for_term <- function(term_name) {
  ensgs <- gprofiler2::gost(term_name)$meta$genes_metadata$query$query_1$ensgs
  genes_df %>% filter(ensembl_gene_id %in% ensgs)
}

snps_df <- read.csv(SNPS_LOCATIONS)
genes_df <- get_gene_positions()
genes_around_gwas_hits <- get_genes_around_snps(genes_df, snps_df, WINDOW_SIZE)

TERM_REGEX <- "cardiac|heart|muscle|sarco|calcium|contractile|myo|Z disc|I band|atri|cardi"

go_results <- gprofiler2::gost(genes_around_gwas_hits$hgnc_symbol)

relevant_term_ids <- go_results$result %>% .[grepl(TERM_REGEX, .$term_name), "term_id"]

################################################################################################################

genes_for_terms <- lapply(relevant_term_ids, get_genes_for_term)
names(genes_for_terms) <- relevant_term_ids

# get_genes_for_term(term_name) %>% filter(ensembl_gene_id %in% genes_around_gwas_hits$ensembl_gene_id)

genes_in_term <- sapply(
  genes_for_terms, 
  function(x) {
    intersect(x$ensembl_gene_id, genes_around_gwas_hits$ensembl_gene_id)
  }
)

candidate_genes <- genes_df %>% 
  filter(ensembl_gene_id %in% unlist(genes_in_term)) %>% 
  mutate(chromosome_name=as.numeric(chromosome_name)) 

# include distance to SNP
merge(
  candidate_genes %>% dplyr::select(-gene_biotype), 
  genes_around_gwas_hits %>% dplyr::select(ensembl_gene_id, SNP, BP, region, distance_from_tss), 
  by="ensembl_gene_id"
) %>% 
  arrange(chromosome_name, start_position) %>% 
  mutate("is_inside"=(sign(BP-start_position)*sign(BP-end_position) < 0)) %>%
  filter(is_inside)

heart_development_go_term <- "GO:0007507"
genes_around_gwas_hits %>% filter(ensembl_gene_id %in% get_genes_for_term(heart_development_go_term)$ensembl_gene_id) %>% 
  mutate("is_inside"=(sign(BP-start_position)*sign(BP-end_position) < 0)) %>%
  filter(is_inside) %>% 
  arrange(as.numeric(chromosome_name), start_position)