import os
import pandas as pd
import re
from subprocess	import call


def run(args):

    '''
    Provide path for a file containing snps and chromosome/position (position is actually not necessary)
    '''
    
    SNPID_COL=1
    CHROMOSOME_COL=0

    df = pd.read_csv(args.snp_file, sep='\t')
    snp_dict = {}

    for row in df.iterrows():
        snp_dict[row[1][CHROMOSOME_COL]] = snp_dict.get(row[1][CHROMOSOME_COL], [])
        if not row[1][SNPID_COL].startswith("Affx"):
            snp_dict[row[1][CHROMOSOME_COL]].append(row[1][SNPID_COL])
    snp_dict

    # from IPython import embed; embed()

    for chromosome in snp_dict.keys():
        bfile = args.bfile_pattern.format(chromosome)
        ofile = args.ofile_pattern.format(chromosome)
        snp_list = ",".join(snp_dict[chromosome])
        command = "plink -bfile {} -snps {} -recodeA --out {}".format(bfile, snp_list, ofile)
        command = command.split()
        call(command)
        # print(command)


if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile_pattern")
    parser.add_argument("--ofile_pattern")
    parser.add_argument("--snp_file")

    # parser.add_argument("--bim_file_pattern")
    # parser.add_argument("--bed_file_pattern")
    
    args = parser.parse_args()
    run(args)