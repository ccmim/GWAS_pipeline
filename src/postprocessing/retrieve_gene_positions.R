#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("biomaRt")

library(tidyverse)
library(biomaRt)

get_gene_positions <- function() {
  
  selected_attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol", "gene_biotype")
  ensembl <- useEnsembl(mirror = "useast", "ensembl",dataset="hsapiens_gene_ensembl", GRCh=37)
  dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  
  pp <- getBM(
    attributes = selected_attributes, 
    #filters = c("hgnc_symbol"), 
    values = gene_list, mart = ensembl
    
  )
  
  pp <- pp %>% filter(chromosome_name %in% as.character(1:22))
  
  pp <- pp %>% filter(gene_biotype == "protein_coding")
  pp
  
}

# gene_list <- c("TTN", "MTSS1", "BAG3", "GRK5", "NMB", "CDKN1A")

LD_INDEP_REGIONS = "~/01_repos/GWAS_pipeline/data/ld_indep_regions/fourier_ls-all_EUR_hg19.bed"
SNPS_LOCATIONS = "~/01_repos/CardiacGWAS/results/snps_for_biomart__one_per_region.txt"

snps_locations = read.csv(SNPS_LOCATIONS)

# pp[grepl("WNT", pp$hgnc_symbol),] %>% .[order(as.numeric(.$chromosome_name)),]


regions <- read.delim(LD_INDEP_REGIONS, sep = "\t")
regions <- regions %>% group_by(chr) %>% mutate(id = paste0(row_number())) %>% as.data.frame
regions$chr <-  sub("\\s+$", "", regions$chr)
regions <- regions %>% mutate(id=paste(chr, id, sep = "_")) %>% ungroup() 
regions$chr <- regions$chr %>% gsub("chr", "", .) %>% as.numeric()

# rownames(regions) <- regions$id
# regions <- regions %>% select(-id)


# Merge by region
# kk <- merge(pp, regions, by.x = "chromosome_name", by.y = "chr")
# kk <- kk %>% filter(start < start_position) %>% filter(start_position < stop)
# kk <- kk[kk$id %in% regions,]

# Merge by SNP

# "CD" -> gene; kk %>% filter(grepl(gene, hgnc_symbol)) # %>% .[,c("id", "hgnc_symbol")]

genes <- c(
"WNT16",
"TTN",
"TP73",
"TOMM70",
"TBX5",
"SKI",
"SHOX2",
"SHC1",
"SH3PXD2B",
"RYR2",
"PTPN11",
"PLN",
"PITX2",
"PIM1",
"PI16",
"PARP2",
"NKX2-5",
"NFATC4",
"NDUFV2",
"MYL2",
"MYH7",
"MYH6",
"MAPK14",
"LOX",
"LMNA",
"HSPB7",
"HEY2",
"GRK2",
"GJC1",
"GJB6",
"GATA6",
"FZD2",
"EFNA1",
"CPLANE2",
"CDKN1A",
"CDC42",
"CCDC103",
"ANK2",
"ADAMTS6",
"ADAM15",
"ACTN2",
"ABI3BP",
"ACTN2",
"LMNA",
"RYR2",
"CDC42",
"TP73",
"SKI",
"MYL2",
"TBX5",
"GJB6",
"MYH6",
"MYH7",
"PARP2",
"GJC1",
"NDUFV2",
"GATA6",
"TTN",
"TOMM70",
"SHOX2",
"NKX2-5",
"MAPK14",
"PIM1",
"HEY2",
"PI16",
"PLN",
"LMNA", 
"RYR2", 
"INPP5F", 
"ATP2A2", 
"MYH6", 
"MYH7", 
"PARP2", 
"GATA6", 
"TTN", 
"TOMM70", 
"PDE5A", 
"HEY2", 
"PI16",
"DPM3", 
"ACTN2",
"GABRD",
"PRKCZ",
"LMNA",
"PRDM16",
"HSPG2",
"SPEN",
"RYR2",
"SKI",
"CPT2",
"BAG3",
"VCL",
"KAT6B",
"MYL2",
"MYH6",
"MYH7",
"KCNJ2",
"TTN",
"TMEM43",
"SLC6A6",
"PLN",
"ACTN2",
"LMNA",
"RYR2",
"CPT2",
"MYL2",
"TBX5",
"PTPN11",
"MMP14",
"MYH6",
"MYH7",
"KCNJ2",
"GATA6",
"TTN",
"TMEM43",
"MYOZ2",
"NKX2-5"
)

# region_names = c(
#   "chr6_78", 
#   "chr17_27", 
#   "chr10_74", 
#   "chr1_11", 
#   "chr2_108", 
#   "chr2_23", 
#   "chr5_103",
#   "chr10_74",
#   "chr10_49",
#   "chr10_48",
#   "chr12_67",
#   "chr3_98",
#   "chr12_19",
#   "chr14_3",
#   "chr4_72",
#   "chr2_69",
#   "chr3_63",
#   "chr4_77",
#   "chr6_29",
#   "chr16_41",
#   "chr13_2",
#   "chr1_78",
#   "chr1_15",
#   "chr1_124",
#   "chr4_78",
#   "chr12_47",
#   "chr17_40",
#   "chr11_19",
#   "chr18_7",
#   "chr16_4",
#   "chr4_9",
#   "chr5_40",
#   "chr1_2",
#   "chr7_34",
#   "chr11_37",
#   "chr12_69",
#   "chr1_11",
#   "chr7_72",
#   "chr14_1",
#   "chr3_10",
#   "chr6_84",
#   "chr2_87",
#   "chr1_32",
#   "chr17_14",
#   "chr16_44",
#   "chr17_26",
#   "chr5_74",
#   "chr18_11",
#   "chr1_33",
#   "chr8_72"
# )# 