import os
from subprocess import call
from SGE_utils import *
odir = "../hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"

experiment = "2020-09-11_02-13-41"
experiment = "2020-09-30_12-36-48"
experiment = "2020-09-30_12-36-43"
pre_phenofile = "data/coma_output/{}/latent_space.csv".format(experiment) 
samples = "data/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3/samples.txt"
covariates = "config_files/covariates/std_covariates.yaml"
pheno_file= "data/tmp/{}/latent_space_adjusted.csv".format(experiment)

def adjust_for_covariates():
    
    # config = yaml.load(open(os.path.join("config_files/coma", config_file)))
    # experiment = config["filename_patterns"]["gwas"].split("__")[0]

    command  = ["Rscript", "src/adjust_for_covariates.R"]
    command += ["-i", pre_phenofile]
    command += ["-o", pheno_file]
    command += ["--samples_white_list", samples]
    command += ["--covariates_file", covariates]
    command += ["--phenotypes_black_list", "ID", "subset"]
    command += ["--gwas_software", "bgenie"]
    print(" ".join(command))
    call(command)

adjust_for_covariates()
# for bgen_file in [x for x in os.listdir(odir) if x.endswith("bgen")]:
#  print(bgen_file)
bgen_dir = "data/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"
odir = "output/coma/{}/BGEN".format(experiment)
os.makedirs(odir, exist_ok=True)
path_extension = "PATH=/home/home01/scrb/bin:$PATH"

for bgen_file in [ x for x in os.listdir(bgen_dir) if x.endswith("bgen") ]:
  # bgen_file = "data/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3/ukb_imp_chr01_100826405-102041016__hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3.bgen"
  region = bgen_file.split("__")[0].split("imp_")[1]
  bgen_file = os.path.join(bgen_dir, bgen_file)
  ofile = os.path.join(odir, "{}.txt".format(region))
  bgen_command = ["bgenie", "--bgen", bgen_file, "--pheno", pheno_file, "--pvals", "--out", ofile]
  bgen_command = " ".join(bgen_command)
  jobname = region

  submit_sge_job(jobname, commands=[path_extension, bgen_command], folder="runs/jobs", memory_limit="4G", walltime="00:10:00", dry_run=False)
