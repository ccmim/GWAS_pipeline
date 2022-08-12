from split_bgen_by_region_qctool import get_regions_df, get_region_data
from PATHS import *
import shlex
from subprocess import call, check_output

import os, sys
os.chdir(check_output(shlex.split("git rev-parse --show-toplevel")).strip())
print(os.getcwd())
import pandas as pd
sys.path.append(".")
from utils.ARC_helpers.SGE_utils import *

if __name__ == "__main__":

    df = get_regions_df()

    chromosome = 1

    #for region in df.head().iterrows():

    for region in df.iterrows():

        region = get_region_data(region)
        
        new_chromosome = region["chromosome"].lstrip("0")
        if chromosome != new_chromosome:
            try:
                print("Created {count} batch job scripts for chromosome {chromosome}... ".format(count=count, chromosome=chromosome))
            except:
                pass
            count = 0
            chromosome = new_chromosome

        path_extension = "PATH=/home/home01/scrb/bin:$PATH"
 
        args = "--chromosome {chromosome} --start_pos {start_pos} --end_pos {end_pos}".format(**region)
        
        filter_bgen_py_command = "python data/code/filter_bgens_for_qc_criteria.py {args}".format(args=args)

        jobname = "chr{chromosome_adjusted}_{start_pos}-{end_pos}".format(**region)
        count += 1
        submit_sge_job(
            jobname,
            commands=[path_extension, "module load anaconda", "conda activate base", filter_bgen_py_command],
            memory_limit="8G", walltime="01:00:00", dry_run=False
        )

    print(f"Created {count} batch job scripts for chromosome {chromosome}... ")
