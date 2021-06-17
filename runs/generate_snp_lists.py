import pandas as pd
import os

def generate_snp_lists(chromosomes, hwe_pval=1e-4, maf=1e-2, missingness=0.1, info=0.3):
    
    df_lst = []
    snp_stats_file_p = "ukb_imp_chr{}_40k__snp-stats.txt" 

    for chromosome in chromosomes:
        print(chromosome)
        snp_stats_file = snp_stats_file_p.format(chromosome)
        snp_stats = pd.read_csv(snp_stats_file, sep=" ", comment="#")
        df_lst.append(snp_stats)
        
    df = pd.concat(df_lst)
    
    snps_hwe = df.rsid[df.HW_exact_p_value < hwe_pval]
    snps_maf = df.rsid[df.minor_allele_frequency < maf]
    snps_info = df.rsid[df.impute_info < info]
    
    return list(snps_hwe), list(snps_maf), list(snps_info)

snps_hwe, snps_maf, snps_info = generate_snp_lists(chromosomes=range(1,23))

with open("hwe_pval_lt_1e-5.txt", "wt") as ff:
    ff.write("\n".join(snps_hwe))

with open("maf_lt_1e-2.txt", "wt") as ff:
    ff.write("\n".join(snps_maf))

with open("info_lt_0.3.txt", "wt") as ff:
    ff.write("\n".join(snps_info))
