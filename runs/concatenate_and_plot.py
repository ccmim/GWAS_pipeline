import os
import pandas as pd
from subprocess import call
from SGE_utils import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--experiments", nargs="+")
parser.add_argument("--z", nargs="+")
args = parser.parse_args()

for experiment in args.experiments:
  for z in args.z:  
    
    path_extension = "PATH=${HOME}/.bin:${PATH}"

    conda_env = "module load anaconda; source activate gwas"
    
    concat_command = [
      "Rscript", "concatenate_gwas.R",
      "--experiments", experiment,
      "--z", z
    ]
    
    concat_command = " ".join(concat_command)

    plot_command = [
      "Rscript", "../analysis/gwas_analysis.R", "--output_folder", "output/coma", 
      "--gwas_folder", "{}/BGEN".format(experiment),
      "--gwas_pattern", "GWAS__{phenotype}__std_covariates__GBR__BGEN__qc",
      "--phenotypes", z
    ]

    if experiment == "2020-09-30_12-36-48":
      plot_command += ["--color_odd_chr", "salmon4", "--color_even_chr", "salmon2"]
      # plot_command += ["--color_odd_chr", "hotpink4", "--color_even_chr", "red2"]

    plot_command = " ".join(plot_command)

    jobname = "experiment_{}__{}".format(experiment, z)
    # submit_sge_job(jobname, commands=[path_extension, conda_env, concat_command, plot_command], memory_limit="16G", walltime="01:00:00", dry_run=False)
    submit_sge_job(jobname, commands=[path_extension, conda_env, plot_command], memory_limit="16G", walltime="01:00:00", dry_run=False)
