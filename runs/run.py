import os, shlex
from subprocess import call, check_output
repo_rootdir = check_output(shlex.split("git rev-parse --show-toplevel")).strip().decode('ascii')
os.chdir(repo_rootdir)

import sys
sys.path.append(os.getcwd())

import src
from src.auxiliary import unfold_config
from src.run_gwas import GWAS_Run

# Import modules
import pandas as pd
import yaml
from copy import deepcopy
import re

from pprint import pprint

from string import Formatter

def adjust_for_covariates(config):
    
    # config = yaml.load(open(os.path.join("config_files/coma", config_file)))
    # experiment = config["filename_patterns"]["gwas"].split("__")[0]

    command  = ["Rscript", "src/adjust_for_covariates.R"]
    command += ["-i", config["filenames"]["phenotype"]]
    command += ["-o", config["filenames"]["phenotype_intermediate"]]
    command += ["--samples_white_list"] + list(config["sample_white_lists"])
    command += ["--covariates_file"] + [config["covariates"]]
    command += ["--phenotypes"] + [x.split("_adj")[0] for x in config["phenotype_list"]]
    command += ["--phenotypes_black_list", "ID", "subset"]
    command += ["--keep_non_cov_adj"]

    print(" ".join(command))
    call(command)


def generate_summary_and_figures(config):
    
    command  = ["Rscript", "analysis/gwas_analysis_new.R"]
    command += ["--output_folder", "."]
    command += ["--gwas_folder", os.path.dirname(config["filename_patterns"]["gwas"])] # "output/traditional_indices" + "/" + config["suffix"],        
    command += ["--gwas_pattern", "GWAS__{phenotype}__" + config["suffix"]]
    command += ["--phenotypes"] + config["phenotype_list"]
    command += ["--qqplot_pooled", "--cache_rds"]

    print(" ".join(command))
    call(command)


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
    
    config = yaml.load(open(args.yaml_config_file))
    
    ### Generate suffix
    name_rules = yaml.load(open(args.name_rules))

    if args.covariates is not None:
        config["covariates"] = args.covariates

    if args.sample_white_lists is not None:
        config["sample_white_lists"] = args.sample_white_lists

    if args.sample_black_lists is not None:
        config["sample_black_lists"] = args.sample_black_lists
    
    # produce suffix (config["suffix"]) from suffix template (args.suffix)
    tokens = extract_formatter_tokens(args.suffix)

    for token in tokens:
        print(token)
        if token in config.keys():
            if isinstance(config[token], list):
            # need to cast to tuple because lists cannot be dict keys
                option_value = tuple(config[token])
            else:
                option_value = config[token]
            tokens[token] = name_rules[token][option_value]

    suffix = args.suffix.format(**tokens)    
    
    # print(suffix)

    ###############

    # def overwrite_config_items(config, args):
    #     for attr, value in args.__dict__.items():
    #         if attr in config.keys() and value is not None:
    #             config[attr] = value

    config = unfold_config(args.yaml_config_file)
    
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
        config["covariates"] = args.covariates

    if args.sample_white_lists is not None:
        config["sample_white_lists"] = args.sample_white_lists

    if args.sample_black_lists is not None:
        config["sample_black_lists"] = args.sample_black_lists

    if args.phenotypes is not None:
        config["phenotype_list"] = args.phenotypes

    if args.chromosomes is not None:
        config["chromosomes"] = args.chromosomes

    # pprint(config)
    return config


def main(config):
    
    adjust_for_covariates(config)
    GWAS_Run(config).run()
    generate_summary_and_figures(config)


if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Harmonise data, run GWAS and generate descriptive plots.")

    parser.add_argument("--yaml_config_file", "-c", default="config_files/ref_config.yaml")
    parser.add_argument("--phenotype_file", default=None)
    parser.add_argument("--covariates", default=None)
    parser.add_argument("--intermediate_phenotype_file", default=None)
    parser.add_argument("--gwas_file", default=None)
    parser.add_argument("--sample_white_lists", nargs="+", default=None)
    parser.add_argument("--sample_black_lists", nargs="+", default=None)
    parser.add_argument("--coma_experiment", default=None)
    parser.add_argument("--chromosomes", "-chr", default=None, help="Chromosomes as a list of comma-separated ranges, e.g. \"1-4,6,10-15\"")
    parser.add_argument("--phenotypes", "-ph", nargs="+", default=None, help="List of phenotypes to perform GWAS on.")
    parser.add_argument("--name_rules", default="config_files/filename_rules/filename_rules.yaml")
    parser.add_argument("--suffix", default="{covariates}__{sample_white_lists}__{quality_control}", help="Subfolder ")
    
    # args = parser.parse_known_args()
    args = parser.parse_args()

    config = prepare_config(args)

    pprint(config)
    main(config)