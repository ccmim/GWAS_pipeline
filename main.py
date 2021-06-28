import os
import shlex
import copy
from subprocess import call

# odir = "../hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"

# experiment = "2020-09-11_02-13-41"

import os, sys, shlex
from subprocess import call, check_output
repo_rootdir = check_output(shlex.split("git rev-parse --show-toplevel")).strip().decode('ascii')
os.chdir(repo_rootdir)
sys.path.append(os.getcwd())

sys.path.append("utils/ARC-helpers")
from SGE_utils import *

# Import modules
import pandas as pd
import yaml
import re
from pprint import pprint
from string import Formatter
from copy import deepcopy

import src
from src.auxiliary import unfold_config
from src.run_gwas import GWAS_Run
import warnings

############################################################################################################
############################################################################################################


def extract_formatter_tokens(pattern):
    # From https://stackoverflow.com/questions/25996937/how-can-i-extract-keywords-from-a-python-format-string
    # return [fname for _, fname, _, _ in Formatter().parse(pattern) if fname]
    
    regex = re.compile('{[A-Za-z_]+}')
    tokens = {}
    for x in regex.findall(pattern):
        token = x[1:-1]
        tokens[token] = None
    return tokens


def prepare_config(args):

    # TODO: solve this warning. 
    # This is just a workaround to prevent too many warning messages to show up.
    # Warning message is: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")    
      config = yaml.load(open(args.yaml_config_file))
      ### Generate suffix
      name_rules = yaml.load(open(args.name_rules))

    if args.covariates is not None:
        config["covariates_config"] = args.covariates

    if args.sample_white_lists is not None:
        config["sample_white_lists"] = args.sample_white_lists

    if args.sample_black_lists is not None:
        config["sample_black_lists"] = args.sample_black_lists    
    
    if args.quality_control is not None:
        config["quality_control"] = args.quality_control

    # produce the actual suffix (config["suffix"]) from suffix template (args.suffix)
    tokens = extract_formatter_tokens(args.suffix)

    for token in tokens:
        if token in config.keys():
            if isinstance(config[token], list):
            # need to cast to tuple because lists cannot be dict keys
                option_value = tuple(config[token])
            else:
                option_value = config[token]            
            tokens[token] = name_rules[token][option_value]

    suffix = args.suffix.format(**tokens)    
    
    ###############

    # def overwrite_config_items(config, args):
    #     for attr, value in args.__dict__.items():
    #         if attr in config.keys() and value is not None:
    #             config[attr] = value

    config = unfold_config(args.yaml_config_file, no_unfolding_for=["covariates_config"])
    
    config["suffix"] = suffix

    config["filenames"] = {}

    if args.phenotype_file is not None:
        config["filename_patterns"]["phenotype"] = args.phenotype_file

    if args.gwas_file is not None:
        config["filename_patterns"]["gwas"] = args.gwas_file

    if args.coma_experiment is not None:
        config["experiment"] = args.coma_experiment
    
    for _fp in ["phenotype", "phenotype_intermediate", "tmpdir", "gwas"]:
        fp = config["filename_patterns"][_fp]
        tokens = extract_formatter_tokens(fp)
        filename = fp.format(**{token: config.get(token, None) for token in tokens})
        
        if _fp != "gwas":
            config["filename_patterns"][_fp] = None
            config["filenames"][_fp] = filename
        else:
            config["filename_patterns"][_fp] = filename

    if args.covariates is not None:
        config["covariates_config"] = args.covariates
    
    config["sample_white_lists"] = args.sample_white_lists
    config["sample_black_lists"] = args.sample_black_lists
    config["phenotype_list"] = args.phenotypes

    if args.chromosomes is not None:
        config["chromosomes"] = args.chromosomes

    config["gwas_software"] = args.gwas_software
    config["bgen_sample_file"] = args.bgen_sample_file
    
    return config


def adjust_for_covariates(config):
    
    # config = yaml.load(open(os.path.join("config_files/coma", config_file)))
    # experiment = config["filename_patterns"]["gwas"].split("__")[0]

    command  = "Rscript src/preprocess_files_for_GWAS.R\n"
    command += "--phenotype_file {}\n".format(config["filenames"]["phenotype"])
    if config["phenotype_list"] is not None:
      command += "--phenotypes {}\n".format(" ".join(config["phenotype_list"]))
    else:
      pass # default behaviour, i.e. use all phenotypes (all columns except those excluded in the following line)
    
    # TOFIX: this should not be hardcoded! It will limit the applicability to phenotype files with other formats
    command += "--columns_to_exclude ID subset\n"    
    command += "--covariates_config_yaml {}\n".format(config["covariates_config"])

    if config["sample_white_lists"] is not None:
      command += "--samples_to_include {}\n".format(" ".join(list(config["sample_white_lists"])))
    if config["sample_black_lists"] is not None:
      command += "--samples_to_exclude {}\n".format(" ".join(list(config["sample_black_lists"])))

    command += "--output_file {}\n".format(config["filenames"]["phenotype_intermediate"])
    command += "--gwas_software {}\n".format(config["gwas_software"])
    if config["bgen_sample_file"] is not None:
      command += "--bgen_sample_file {}\n".format(config["bgen_sample_file"])
    # command += "--overwrite_output\n"
    
    print("\nPreprocessing the phenotype file to perform GWAS on {}.".format(config["gwas_software"]))    
    print(command)
    # print("\n")
    
    call(shlex.split(command))
    print("\n")

def generate_summary_and_figures(config):
    
    ### TODO: Use a configuration file for the file name patterns

    gwas_folder = os.path.dirname(config["filename_patterns"]["gwas"])

    command  = "Rscript src/postprocessing/gwas_analysis.R\n"
    command += "--output_folder .\n"
    command += "--gwas_folder {}\n".format(gwas_folder) # "output/traditional_indices" + "/" + config["suffix"],        
    command += "--gwas_pattern GWAS__{{phenotype}}__{}\n".format(config["suffix"])
    if config["phenotype_list"] is not None:
      command += "--phenotypes {}\n".format(" ".join(config["phenotype_list"])) 
    command += "--qqplot_pooled\n"
    command += "--cache_rds\n"

    print("\nCreating Manhattan plots, Q-Q plots and region-wise summaries.")    
    print(command)
    print("\n")
    
    call(shlex.split(command))

def concatenate_gwas(config, submit_to_hpc_queue=False):
    
    ### TODO: Use a configuration file for the file name patterns

    gwas_folder = os.path.dirname(config["filename_patterns"]["gwas"])
    output_fp = "GWAS__{{z}}__{suffix}".format(suffix=config["suffix"])

    concat_command  = "Rscript src/postprocessing/concatenate_gwas.R\n"
    concat_command += "--experiment {}\n".format(config["experiment"])
    concat_command += "--z {z}\n"
    concat_command += "--input_results_dir {}\n".format(gwas_folder) # "output/traditional_indices" + "/" + config["suffix"],        
    concat_command += "--output_filename_pattern {}".format(output_fp)
    
    path_extension = "PATH=${HOME}/.bin:${PATH}"
    conda_env = "module load anaconda; source activate gwas"

    print("\nConcatenate region-wise GWAS summary statistics\n")
                
    #TOFIX: replace this, please
    for z in ["z"+str(i) for i in range(8)]:
      if submit_to_hpc_queue:
        concat_command_ = concat_command.replace("\n", " \\\n").format(z=z)
        jobname = "concat__{experiment}__{z}".format(experiment=config["experiment"], z=z)
        submit_sge_job(jobname, 
          commands=[path_extension, conda_env, concat_command_], 
          memory_limit="16G", walltime="01:00:00", dry_run=False
        )
      else:
        concat_command_ = concat_command.format(z=z)
        print(concat_command_)
        call(shlex.split(concat_command_))
        print("\n")


def main(config):
        
    adjust_for_covariates(config)
    BGENIE_Run(config).run()
    yaml.dump(config, open(os.path.join(os.path.dirname(config["filename_patterns"]["gwas"]), "config.yaml"), "w"))
    concatenate_gwas(config, submit_to_hpc_queue=True)
    generate_summary_and_figures(config)


class BGENIE_Run(object):
    
  def __init__(self, config):

    experiment = config["experiment"]
    self.bgen_dir = "data/genotypes/imputed/hwe_pval_gt_1e-5__maf_gt_1e-2__info_gt_0.3"
    self.odir = os.path.dirname(config["filename_patterns"]["gwas"]) # "output/coma/{}/BGEN".format(experiment)
    self.pheno_file = config["filenames"]["phenotype_intermediate"]
    os.makedirs(self.odir, exist_ok=True)
  
        
  def run(self, run_locally=False):
    
    commands = []
    commands.append("export PATH=/home/home01/scrb/.bin:$PATH")

    bgen_files = [ x for x in os.listdir(self.bgen_dir) if x.endswith("bgen") ]

    for bgen_file in bgen_files:
        
      region = bgen_file.split("__")[0].split("imp_")[1]
      bgen_file = os.path.join(self.bgen_dir, bgen_file)
      ofile = os.path.join(self.odir, "{}.txt".format(region))

      bgen_command = "bgenie\n"
      bgen_command += "--bgen {}\n".format(bgen_file)
      bgen_command += "--pheno {}\n".format(self.pheno_file)
      bgen_command += "--pvals\n"
      bgen_command += "--out {}\n".format(ofile)
      # bgen_command = " ".join(bgen_command)


      if run_locally:
        os.environ["PATH"] += ":"+os.path.join(os.environ["HOME"],".bin")
        for command_ in commands_:
          call(shlex.split(command_))
      else:
        jobname = region
        commands_ = copy.copy(commands)
        bgen_command = bgen_command.replace("\n", " \\\n")
        commands_.append(bgen_command)
        submit_sge_job(
          jobname, commands=commands_, folder="runs/jobs", 
          memory_limit="4G", walltime="00:05:00", dry_run=False
        )


if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Harmonise data, run GWAS and generate descriptive plots.")

    parser.add_argument("--yaml_config_file", "-c", default="config_files/ref_config.yaml", help="Reference configuration file. The rest of the command line arguments overwrite the configuration items.")
    parser.add_argument("--phenotype_file", default=None)
    parser.add_argument("--phenotypes", "-ph", nargs="+", default=None, help="List of phenotypes to perform GWAS on.")
    parser.add_argument("--covariates", default=None, help="YAML file with the configuration specifying the covariates to adjust for.")
    parser.add_argument("--phenotype_intermediate", default=None, help="File pattern for the intermediate phenotype file, with covariate-adjusted scores formatted according to the GWAS tool to be used.")
    parser.add_argument("--gwas_file", default=None, help="File pattern for the output GWAS file")
    parser.add_argument("--gwas_software", default="plink", help="GWAS tool. Currently only plink and BGENIE are supported")
    parser.add_argument("--sample_white_lists", nargs="+", default=None)
    parser.add_argument("--sample_black_lists", nargs="+", default=None)
    #TODO: add bgen_sample_file to the configuration file    
    parser.add_argument("--bgen_sample_file", default=None, help="Only required if --gwas_software option is BGENIE. It's the sample file linked to the BGEN file.")
    parser.add_argument("--coma_experiment", "-e", default=None, help="If the file patterns contain the {experiment} token, this is replaced by this argument. Meant to be used with CoMA experiment, hence its name.")
    parser.add_argument("--quality_control", "-qc", default=None)
    parser.add_argument("--chromosomes", "-chr", default=None, help="Chromosomes as a list of comma-separated ranges, e.g. \"1-4,6,10-15\"")    
    parser.add_argument("--name_rules", default="config_files/filename_rules/filename_rules.yaml")
    parser.add_argument("--suffix", default="{covariates}__{sample_white_lists}__{sample_black_lists}__{quality_control}", help="Subfolder ")
    
    args = parser.parse_args()

    config = prepare_config(args)

    main(config)
