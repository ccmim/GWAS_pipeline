import pandas as pd
from subprocess import call, check_output, check_call
import shlex
from split_bgen_by_region_qctool import build_qctool_command
from PATHS import *

def get_QCed_snps(bgen_file, snps_stats_file=None, snps_list_file=None, bgen_ofile=None):

  if snps_stats_file is None:
      snps_stats_file = bgen_file[:-5] + "_snps_stats.txt"
      snps_stats_file = str(snps_stats_file)

  if snps_list_file is None:
      snps_list_file = bgen_file[:-5] + "_snps_passing_QC.txt"

  if bgen_ofile is None:
      bgen_ofile = bgen_file[:-5] + "_snps_passing_QC.bgen"

  qctool_stats_command = "qctool -g {bgen_file} -snp-stats -osnp {snps_stats_file}".format(bgen_file=bgen_file, snps_stats_file=snps_stats_file)
  check_call(shlex.split(qctool_stats_command))
 
  snps_stats_df = pd.read_csv(snps_stats_file, sep=" ", comment="#")
  snps_stats_df = snps_stats_df[["rsid", "HW_exact_p_value", "impute_info", "missing_proportion", "minor_allele_frequency"]]                                                                                                                   
  
  qc_condition = (snps_stats_df.HW_exact_p_value > 1e-5) & (snps_stats_df.impute_info > 0.3) & (snps_stats_df.minor_allele_frequency > 0.01) 
  rsids = snps_stats_df[qc_condition].reset_index().rsid
  rsids.to_csv(snps_list_file, index=False)

  qctool_filter_bgen = build_qctool_command(input_genotype_file=bgen_file, output_genotype_file=bgen_ofile, incl_rsid_files=[snps_list_file]) 
  check_call(shlex.split(qctool_filter_bgen))


if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromosome")
    parser.add_argument("--start_pos")
    parser.add_argument("--end_pos")
    
    args = parser.parse_args()    
    region = args.__dict__
    region["chromosome_adjusted"] = region["chromosome"]

    bgen_file = OUTPUT_GENOTYPE_FILE_PATTERN.format(**region) 
    snps_stats_file = SNPS_STATS_FILE_PATTERN.format(**region)
    snps_list_file = SNPS_LIST_FILE_PATTERN.format(**region)
    bgen_ofile = FINAL_BGEN_FILE_PATTERN.format(**region)

    get_QCed_snps(bgen_file, snps_stats_file, snps_list_file, bgen_ofile)

