import os, sys, shlex
from subprocess import call, check_output
repo_rootdir = check_output(shlex.split("git rev-parse --show-toplevel")).strip().decode('ascii')
os.chdir(repo_rootdir)
sys.path.append(os.getcwd())

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
import time
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
      config = yaml.safe_load(open(args.yaml_config_file))
      ### Generate suffix
      name_rules = yaml.full_load(open(args.name_rules))

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
    config["run_id"] = args.run_id
    config["experiment_id"] = args.experiment_id
    for token in tokens:
        if token in config.keys():
            if isinstance(config[token], list):
            # need to cast to tuple because lists cannot be dict keys
                option_value = tuple(config[token])
            else:
                option_value = config[token]            
            try:
                tokens[token] = name_rules[token][option_value]
            except:
                tokens[token] = option_value

    suffix = args.suffix.format(**tokens)    
    
    ###############

    config = unfold_config(args.yaml_config_file, no_unfolding_for=["covariates_config"])
    
    config["suffix"] = suffix

    config["filenames"] = {}

    if args.phenotype_file is not None:
        config["filename_patterns"]["phenotype"] = args.phenotype_file

    if args.gwas_file is not None:
        config["filename_patterns"]["gwas"] = args.gwas_file

    config["run_id"] = args.run_id
    config["experiment_id"] = args.experiment_id

    for _fp in ["phenotype_file", "phenotype_intermediate", "gwas"]:
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
    
    if args.bgen_sample_file is None:
       config["bgen_sample_file"] = config["filename_patterns"]["genotype"]["bgen_sample_file"]
       del config["filename_patterns"]["genotype"]["bgen_sample_file"]

    return config


def adjust_for_covariates(config):
    
    # config = yaml.load(open(os.path.join("config_files/coma", config_file)))
    
    command  = "Rscript src/preprocess_files_for_GWAS.R\n"
    command += "--phenotype_file {}\n".format(config["filenames"]["phenotype_file"])
    if config["phenotype_list"] is not None:
      command += "--phenotypes {}\n".format(config["phenotype_list"])
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
    
    command += "--overwrite\n"
    
    return command


def build_bgenie_command(config):
    
    gwas_folder = os.path.dirname(config["filename_patterns"]["gwas"])
    gwas_folder_by_region =gwas_folder + "/by_region"
    
    pheno_file = config["filenames"]["phenotype_intermediate"]
    os.makedirs(gwas_folder_by_region, exist_ok=True)    

    command = "python src/bgenie_per_region.py\n"
    command += f"--gwas_folder {gwas_folder_by_region}\n"
    command += f"--pheno_file {pheno_file}\n"    
    # command += f"--dry-run"
    return command


def postprocess_gwas_by_region(config):
    
    experiment_id = config["experiment_id"] 
    run_id = config["run_id"]
    gwas_file_pattern = config["filename_patterns"]["gwas"]
    gwas_folder = os.path.dirname(gwas_file_pattern)

    df = pd.read_csv(config["filenames"]["phenotype_file"])
    
    # phenotypes = list(df.columns[df.columns.str.startswith("z")])
    phenotypes = list(df.columns[df.columns != "ID"])
    phenotypes = " ".join(phenotypes)
        
    command = "python src/postprocessing/postprocess_gwas_by_region.py\n" 
    command += f"--experiment_id {experiment_id} \n"
    command += f"--run_id {run_id} \n"
    command += f"--input_results_dir {gwas_folder} \n"
    command += f"--output_filename_pattern {gwas_file_pattern}\n"
    command += f"--phenotypes {phenotypes}\n"
    
    return command


if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Harmonise data, run GWAS and generate descriptive plots.")

    parser.add_argument("--yaml_config_file", "-c", default="config_files/ref_config.yaml", help="Reference configuration file. The rest of the command line arguments overwrite the configuration items.")
    parser.add_argument("--phenotype_file", default=None)
    parser.add_argument("--phenotypes", "-ph", nargs="+", default=None, help="List of phenotypes to perform GWAS on.")
    parser.add_argument("--covariates", default=None, help="YAML file with the configuration specifying the covariates to adjust for.")
    parser.add_argument("--phenotype_intermediate", default=None, help="File pattern for the intermediate phenotype file, with covariate-adjusted scores formatted according to the GWAS tool to be used.")
    parser.add_argument("--gwas_file", default=None, help="File pattern for the output GWAS file")
    parser.add_argument("--gwas_software", default="bgenie", help="GWAS tool. Currently only plink and BGENIE are supported")
    parser.add_argument("--sample_white_lists", nargs="+", default=None)
    parser.add_argument("--sample_black_lists", nargs="+", default=None)
    #TODO: add bgen_sample_file to the configuration file    
    parser.add_argument("--bgen_sample_file", default=None, help="Only required if --gwas_software option is BGENIE. It's the sample file linked to the BGEN file.")
    parser.add_argument("--quality_control", "-qc", default=None)
    parser.add_argument("--chromosomes", "-chr", default=None, help="Chromosomes as a list of comma-separated ranges, e.g. \"1-4,6,10-15\"")    
    parser.add_argument("--name_rules", default="config_files/filename_rules/filename_rules.yaml")
    parser.add_argument("--suffix", default="{experiment_id}_{run_id}", help="")
    
    parser.add_argument("--experiment_id", help="MLflow experiment ID")
    parser.add_argument("--run_id", help="MLflow run ID")
    
    parser.add_argument("--no_print_commands", action="store_true", help="If set, commands will not be printed out.")
    parser.add_argument("--dry-run", "--dry_run", "--dryrun", action="store_true", help="If set, it will only print the commands that are to be run.")
    parser.add_argument("--print_config", action="store_true", help="If set, it will only print the commands that are to be run.")
    parser.add_argument("--steps_to_run", nargs="+", default=["1", "2", "3"], help="a sublist of [1,2,3], indicating which steps are to be performed.")
    
    args = parser.parse_args()
    config = prepare_config(args)

    if args.print_config:
        pprint(config)
        print("\n")

    adj_command = adjust_for_covariates(config)
    gwas_command = build_bgenie_command(config)    
    postproc_command = postprocess_gwas_by_region(config)
    
    commands = {
      "1": adj_command,
      "2": gwas_command,
      "3": postproc_command
    }  
    
    messages = {
      "1": "\nPreprocessing the phenotype file to perform GWAS on {}\n.".format(config["gwas_software"]),
      "2": "\nSubmitting GWAS jobs to the queue\n",
      "3": "\nConcatenating per-region GWAS files, creating Manhattan plots, Q-Q plots and region-wise summaries\n"
    }   
 
    for k, command in commands.items():

        if k in args.steps_to_run: 
            if not args.no_print_commands:
                print(command)
            if not args.dry_run:
                print(messages[k])
                call(shlex.split(command))
                
    yaml.dump(config, open(os.path.join(os.path.dirname(config["filename_patterns"]["gwas"]), "config.yaml"), "w"))
