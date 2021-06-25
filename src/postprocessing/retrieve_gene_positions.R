library(biomaRt)

gene_list <- c("TTN", "MTSS1", "BAG3", "GRK5", "NMB", "CDKN1A")

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
pp <- getBM(attributes = c("start_position", "chromosome_name"), filters = c("hgnc_symbol"), values = gene_list, mart = ensembl)