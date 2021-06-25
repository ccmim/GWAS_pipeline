import os
import copy
from subprocess import call
import ARC_helpers as ARC
from ARC.SGE_utils import *
odir = "../hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"

experiment = "2020-09-11_02-13-41"
experiment = "2020-09-30_12-36-48"
experiment = "2020-09-30_12-36-43"

def adjust_for_covariates():
   
    command  = ["Rscript", "src/adjust_for_covariates.R"]
    command += ["-i", pre_phenofile]
    command += ["-o", pheno_file]
    command += ["--samples_white_list", samples]
    command += ["--covariates_file", covariates]
    command += ["--phenotypes_black_list", "ID", "subset"]
    command += ["--gwas_software", "bgenie"]
    print(" ".join(command))
    call(command)

def main():

  pre_phenofile = "data/coma_output/{}/latent_space.csv".format(experiment) 
  samples = "data/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3/samples.txt"
  covariates = "config_files/covariates/std_covariates.yaml"
  pheno_file= "data/tmp/{}/latent_space_adjusted.csv".format(experiment)
  adjust_for_covariates()
  
  bgen_dir = "data/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"
  odir = "output/coma/{}/BGEN".format(experiment)
  os.makedirs(odir, exist_ok=True)
  
  commands = []
  commands.append("PATH=/home/home01/scrb/bin:$PATH")

  for bgen_file in [ x for x in os.listdir(bgen_dir) if x.endswith("bgen") ]:
        
    region = bgen_file.split("__")[0].split("imp_")[1]
    bgen_file = os.path.join(bgen_dir, bgen_file)
    ofile = os.path.join(odir, "{}.txt".format(region))

    bgen_command = ["bgenie", "--bgen", bgen_file, "--pheno", pheno_file, "--pvals", "--out", ofile]
    bgen_command = " ".join(bgen_command)

    commands_ = copy.copy(commands)
    commands_.append(bgen_command)

    jobname = region
    
    submit_sge_job(
      jobname, commands=commands_, folder="runs/jobs", 
      memory_limit="4G", walltime="00:10:00", dry_run=False
    )

if __name__ == "__main__":
  main()
