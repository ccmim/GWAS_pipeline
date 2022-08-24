import shlex
import re
import argparse
from subprocess import call, check_output

import os, sys
os.chdir(check_output(shlex.split("git rev-parse --show-toplevel")).strip())
import pandas as pd
sys.path.append(".")

from data.code.PATHS import *
from data.code.split_bgen_by_region_qctool import get_regions_df, get_region_data
from utils.ARC_helpers.SGE_utils import *

import itertools


# From: https://stackoverflow.com/questions/8991506/iterate-an-iterator-by-chunks-of-n-in-python
def grouper(n, iterable):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk
        

#def grouper(n, iterable, fillvalue=None):
#    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
#    args = [iter(iterable)] * n
#    return izip_longest(fillvalue=fillvalue, *args)


def run_bgenie(regions, args, previous_chromosome, job_count):
        
    chromosome = regions[0]["chromosome"].lstrip("0")
    
    if int(chromosome) != int(previous_chromosome):
        try:
            print("Created {job_count} batch job scripts for chromosome {chromosome}... ".format(job_count=job_count, chromosome=previous_chromosome))
        except:
            pass
        job_count = 0
        previous_chromosome = chromosome

    path_extension = "PATH=/home/home01/scrb/bin:$PATH"
 
    pheno_file = args.pheno_file
   
    bgenie_commands = []
      
    for region in regions:
    
        bgen_file = OUTPUT_GENOTYPE_FILE_PATTERN.format(**region)
        region_name = "chr{chromosome}_{start_pos}-{end_pos}.txt".format(**region)
        gwas_file_per_region = os.path.join(args.gwas_folder, region_name)
        bgenie_command = f"bgenie --bgen {bgen_file} --pheno {pheno_file} --out {gwas_file_per_region} --pvals --thread 8"
        bgenie_commands.append(bgenie_command)
    
    first_region = regions[0]
    last_region = region
    job_start_pos = "chr{chromosome}_{start_pos}".format(**first_region)
    job_end_pos = "chr{chromosome}_{start_pos}".format(**last_region)
    
    jobname = f"bgenie_{job_start_pos}-{job_end_pos}"

    job_count += 1     
    # add_job_as_pending = "echo ${JOB_ID} >> " + f"{args.gwas_folder}/pending_jobs"
    pending_jobs_file = os.path.join(os.path.dirname(args.gwas_folder), "pending_jobs")
    remove_job_from_pending = "sed -i '/'${JOB_ID}'/d' " +  f"{pending_jobs_file}"

    commands = [path_extension, "module load anaconda", "conda activate gwas", *bgenie_commands, remove_job_from_pending]
    return jobname, commands, chromosome, job_count


def rm_pending_jobs_file(job_ids_file):

    if os.path.exists(job_ids_file):
      print(f"Removing pre-existing pending jobs file: {job_ids_file}")
      os.remove(job_ids_file)


def jobid_from_sge_msg(sge_msg):
    
    jobid = (re.search("([0-9]+)", str(sge_message)).group(1))
    return jobid


def add_job_as_pending(jobid, job_ids_file):

    open(job_ids_file, "at").write("{}\n".format(jobid))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
     
    parser.add_argument("--gwas_folder")
    parser.add_argument("--pheno_file")
    parser.add_argument("--n_regions_per_job", default=8)
    parser.add_argument("--dry-run", "--dry_run", default=False, action="store_true")
    args = parser.parse_args()
    
    pheno_file = args.pheno_file
 
    previous_chromosome = 1
    job_count = 0

    job_ids_file = os.path.join(os.path.dirname(args.gwas_folder), "pending_jobs")
    rm_pending_jobs_file(job_ids_file)
    
    df = get_regions_df()
    
    # for region_row in df.head(100).iterrows():
    for region_rows in grouper(args.n_regions_per_job, df.iterrows()):
    
        regions_data = [get_region_data(region_row) for region_row in region_rows]        
        
        jobname, commands, previous_chromosome, job_count = run_bgenie(regions_data, args, previous_chromosome, job_count)
        
        # print(jobname)
        
        job_config = {
          "jobname": jobname,
          "memory_limit": "8G", 
          "walltime": "01:00:00",
          "open_mp_nodes": 8,           
          "stderr": f"runs/logs/bgenie/{jobname}.err", 
          "stdout": f"runs/logs/bgenie/{jobname}.out"
        }
       
        jobfile, sge_message = submit_sge_job(
            commands=commands,
            **job_config,
            dry_run=args.dry_run,
        )        
        
        print(jobfile)
        
        if not args.dry_run:
            jobid = jobid_from_sge_msg(sge_message)
            add_job_as_pending(jobid, job_ids_file)
        
        
    print(f"Created {job_count} batch job scripts for chromosome {previous_chromosome}... ")
