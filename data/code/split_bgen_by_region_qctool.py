# from IPython import embed
from data.code.PATHS import *
import shlex
from subprocess import call, check_output

import os, sys
os.chdir(check_output(shlex.split("git rev-parse --show-toplevel")).strip())
print(os.getcwd())
import pandas as pd
sys.path.append(".")
from utils.ARC_helpers.SGE_utils import *


def get_regions_df():
    df = pd.read_csv(REGIONS_FILE, sep = "\t")
    df.columns = [x.strip() for x in df.columns]
    return df


def generate_incl_snps_file(full_df, region, replace_previous=False):

    # df.columns = [x.strip() for x in df.columns]
    filename = REDUCED_SNP_FILE_PATTERN.format(chromosome=chromosome, start_pos=region["start_pos"], end_pos=region["end_pos"])
    if not os.path.exists(filename) or replace_previous:
        df_filtered = full_df[(region["start_pos"] < full_df.iloc[:,1]) & (region["end_pos"] > full_df.iloc[:,1])]
        included_snps_df = df_filtered.iloc[:,[0]]
        included_snps_df.to_csv(filename, index=False, header=False)
    return filename



def get_region_data(row):

    region = row[1]  # element 0 is the row index
    chromosome = region.chr[3:].strip()  # chrN ---> N

    chromosome_adjusted = "0" + chromosome if len(chromosome) == 1 else chromosome

    region = {
        "chromosome": chromosome_adjusted,
        "chromosome_adjusted": chromosome_adjusted,
        "start_pos": region.start,
        "end_pos": region.stop
    }

    return region


def build_bgenix_command(input_genotype_file, index_file, output_genotype_file, incl_range):
    command = ["bgenix", "-g", input_genotype_file, "-i", index_file, "-incl-range", incl_range, ">", output_genotype_file]
    command = " ".join(command)
    return command
      
 
def build_qctool_command(input_genotype_file, output_genotype_file, incl_rsid_files=None, excl_rsid_files=None, incl_range=None, incl_samples=None, samples=None):

    command = ["qctool", "-g", input_genotype_file, "-og", output_genotype_file]

    if samples is not None:
        command += ["-s", samples]

    if incl_rsid_files is not None:
        command += ["-incl-rsids"] + incl_rsid_files

    if excl_rsid_files is not None:
        command += ["-excl-rsids"] + excl_rsid_files

    if incl_range is not None:
        command += ["-incl-range", incl_range]

    if incl_samples is not None:
        command += ["-incl-samples", incl_samples]

    command = " ".join(command)
    return command


if __name__ == "__main__":

    df = get_regions_df()

    chromosome = 1
    #for region in df.head().iterrows():
    for region in df.iterrows():

        region = get_region_data(region)
        
        new_chromosome = region["chromosome"].lstrip("0")
        if chromosome != new_chromosome:
            try:
                print(f"Created {count} batch job scripts for chromosome {chromosome}... ")
            except:
                pass
            count = 0
            chromosome = new_chromosome
            # print(f"Processing chromosome {chromosome}")
            input_genotype_file = INPUT_GENOTYPE_FILE_PATTERN.format(chromosome=chromosome)
            index_file = INDEX_FILE_PATTERN.format(chromosome=chromosome)
            if not os.path.exists(input_genotype_file):
                print(f"File {input_genotype_file} does not exist. Skipping...")
                continue
            incl_rsid_file = INCLUDE_SNP_FILES[0].format(chromosome=chromosome)
            full_df = pd.read_csv(incl_rsid_file, sep=" ", header=None)
        reduced_incl_snp_file = generate_incl_snps_file(full_df, region)


        incl_range = "{chromosome_adjusted}:{start_pos}-{end_pos}".format(**region)
        temp_genotype_file = TEMP_GENOTYPE_FILE_PATTERN.format(**region)
        output_genotype_file = OUTPUT_GENOTYPE_FILE_PATTERN.format(**region)

        path_extension = "PATH=/home/home01/scrb/bin:$PATH"
        bgenix_command = build_bgenix_command(input_genotype_file, index_file, temp_genotype_file, incl_range)
        qctool_command = build_qctool_command(temp_genotype_file, output_genotype_file, incl_rsid_files=[reduced_incl_snp_file], incl_range=None, incl_samples=INCLUDE_SAMPLES, samples=SAMPLES_FILE)

        jobname = "chr{chromosome_adjusted}_{start_pos}-{end_pos}".format(**region)

        count += 1
        submit_sge_job(
            jobname,
            commands=[path_extension, bgenix_command, qctool_command, f"rm {temp_genotype_file}"],
            memory_limit="8G", walltime="01:00:00", dry_run=False
        )

    print(f"Created {count} batch job scripts for chromosome {chromosome}... ")
