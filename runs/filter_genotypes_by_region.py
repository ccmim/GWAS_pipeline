import os
import pandas as pd
from subprocess import call
from SGE_utils import *

home = os.getenv("HOME")
df = pd.read_csv(os.path.join(home, "data/genetics/fourier_ls-all_EUR_hg19.bed"), sep="\t")
df.columns = [x.strip() for x in df.columns]

excl_rsid_files = [ os.path.join("../snp_files", x) for x in ["hwe_pval_lt_1e-5.txt", "maf_lt_1e-2.txt", "info_lt_0.3.txt"]]
incl_sample = "cmr_GBR_ids.txt"

suffix = "hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"
input_genotype_file_pattern = "../maf_gt_1e-3/ukb_imp_chr{}_40k_maf_gt_1e-3.bgen"
output_genotype_file_pattern = "../{suffix}/ukb_imp_chr{chromosome_adjusted}_{start_pos}-{end_pos}__{suffix}.bgen"

for region in df.iterrows():
    
    region = region[1] # element 0 is the row index
    chromosome = region.chr[3:].strip() # chrN ---> N
    chromosome_adjusted = "0"+chromosome if len(chromosome) == 1 else chromosome
    
    region = {
      "chromosome_adjusted": chromosome_adjusted, 
      "start_pos": region.start, 
      "end_pos": region.stop
    }
    
    input_genotype_file = input_genotype_file_pattern.format(chromosome)
    output_genotype_file = output_genotype_file_pattern.format(**region, suffix=suffix)
    
    if not os.path.exists(input_genotype_file):
        continue
    
    incl_range = "{chromosome_adjusted}:{start_pos}-{end_pos}".format(**region)
    
    path_extension = "PATH=/home/home01/scrb/bin:$PATH"
    
    command = [
      "qctool", 
      "-g", input_genotype_file, 
      "-og", output_genotype_file, 
      "-excl-rsids"] + excl_rsid_files + [
      "-incl-range", incl_range, 
      "-incl-samples", incl_sample
    ]
    
    command = " ".join(command)
    
    jobname = "chr{chromosome_adjusted}_{start_pos}-{end_pos}".format(**region)
    submit_sge_job(jobname, commands=[path_extension, command], memory_limit="8G", walltime="01:00:00", dry_run=False)
