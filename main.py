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

    command  = ["Rscript", "src/preprocess_files_for_GWAS.R"]    
    command += ["--phenotype_file", config["filenames"]["phenotype"]]    
    if config["phenotype_list"] is not None:
      command += ["--phenotypes"] + config["phenotype_list"]
    else:
      pass # default behaviour, i.e. use all phenotypes (all columns except those excluded in the following line)
    
    # TOFIX: this should not be hardcoded! It will limit the applicability to phenotype files with other formats
    command += ["--columns_to_exclude", "ID", "subset"]    
    command += ["--covariates_config_yaml"] + [config["covariates_config"]]    

    if config["sample_white_lists"] is not None:
      command += ["--samples_to_include"] + list(config["sample_white_lists"])
    if config["sample_black_lists"] is not None:
      command += ["--samples_to_exclude"] + list(config["sample_black_lists"])    

    command += ["--output_file", config["filenames"]["phenotype_intermediate"]]
    command += ["--gwas_software", config["gwas_software"]]
    if config["bgen_sample_file"] is not None:
      command += ["--bgen_sample_file", config["bgen_sample_file"]]
    
    print("\nPreprocessing the phenotype file to perform GWAS on {}.".format(config["gwas_software"]))    
    print(" ".join(command))
    print("\n")
    
    call(command)
    print("\n")

def generate_summary_and_figures(config):
    
    ### TODO: Use a configuration file for the file name patterns

    gwas_folder = os.path.dirname(config["filename_patterns"]["gwas"])
    command  = ["Rscript", "src/postprocessing/gwas_analysis.R"]
    command += ["--output_folder", "."]
    command += ["--gwas_folder", gwas_folder] # "output/traditional_indices" + "/" + config["suffix"],        
    command += ["--gwas_pattern", "GWAS__{phenotype}__" + config["suffix"]]
    if config["phenotype_list"] is not None:
      command += ["--phenotypes"] + config["phenotype_list"]
    command += ["--qqplot_pooled", "--cache_rds"]

    print("\nCreating Manhattan plots, Q-Q plots and region-wise summaries.")    
    print(" ".join(command))
    print("\n")
    
    call(command)


def main(config):
        
    adjust_for_covariates(config)
    # GWAS_Run(config).run()
    yaml.dump(config, open(os.path.join(os.path.dirname(config["filename_patterns"]["gwas"]), "config.yaml"), "w"))
    generate_summary_and_figures(config)


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
