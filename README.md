# GWAS pipeline

This repository contains code to execute GWAS on UK Biobank data, using the Plink and BGENIE tools, and perform downstream analysis on the results.

The folder `code/` contains scripts to perform pre-processing of the data and GWAS execution.
The folder `analysis/` contains scripts to perform statistical analysis on the GWAS results and generate figures (like Manhattan plots or Q-Q plot).
The folder `download_data/` contains scripts to download data from the UK Biobank and filter the genotype files.

## Requirements
This code has been tested on `Python 3.6.3` and `R 3.6`.
For pre-processing data, it requires the `tidyverse` package.
For performing GWAS, it requires installing the tools `plink 1.9` (for `bed/bim/fam` files) and/or `BGENIE 1.4.1` (for `bgen` files).
For performing downstream analysis, it requires the `qqman` package.

### Conda environment
_TO DO_: provide file with the requirements to build the environment.

## Usage
The pipeline consists of scripts for:
1) **Fetching data**: download genotype data from the UK Biobank and image-derived phenotypes from our AWS S3 buckets.
2) **Pre-processing data**: adjusting the phenotypes for a set of covariates and performing inverse-normalization.
3) **Executing GWAS**: self-explanatory. 
4) **Analyze results**: explore the output of the GWAS by generating Manhattan plots and Q-Q plots. Integrating data from other sources in order to interpret the results.   

Steps 2 to 4 rely on a single `yaml` configuration file.

#### Fetching data
Instructions on how to download each kind of genetic data can be found [in this link](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/ukbgene_instruct.html).

#### Pre-processing data
The script that performs this task is `src/preprocess_files_for_GWAS.R`.
Example of usage (GWAS on left-ventricular end-diastolic volume, run on unrelated British subjects using Plink):

```
Rscript src/preprocess_files_for_GWAS \
  --phenotype_file data/phenotypes/cardiac_phenotypes/lvedv.csv
  --phenotypes LVEDV
  --columns_to_exclude id 
  --samples_to_include data/ids_list/british_subjects.txt
  --samples_to_exclude data/ids_list/related_british_subjects.txt
  --covariates_config_yaml config_files/standard_covariates.yml
  --output_file output/lvedv_adjusted_british.tsv
  --gwas_software plink
  --overwrite_output
```

#### Executing GWAS
One needs to execute the command

`python main.py --yaml_config_file <YAML_FILE>`

The complexity is located in the YAML configuration file.

##### Configuration file

    chromosomes: 20-22  
    data_dir: <>  
    output_dir: <>  
    individuals: <>  
    filename_patterns: {
        genotype: {
            bed: <>,
            bim: <>, 
            fam: <>
        },
        phenotype: {
            phenotype_file: <>,
            phenotype_file_tmp: <>,
            phenotype_list: <>,   
            covariates: <>
        },
        gwas: "gwas/{phenotype}/gwas__{phenotype}{suffix}"
    }
    exec:
        plink: "plink"
    suffix: "{TOKEN1}__{TOKEN2}"
    options: {
        TOKEN1: VALUE1_i
        TOKEN2: VALUE2_j            
    }
    suffix_tokens: {
      TOKEN1: {
        VALUE1_1: NAME1_1,  
        VALUE1_2: NAME1_2 
      },
      TOKEN2: {
        VALUE2_1: NAME2_1,  
        VALUE2_2: NAME2_2
      }
    }

- `chromosomes` (optional): it consists of comma-separated ranges, where a range is a single number (e.g. `1`) or a proper range (e.g. `3-7`).
- `data_dir` (optional): directory where the input data is stored. If provided and the genotype and phenotype file paths are relative paths, they are interpreted as hanging from `data_dir`.
- `output_dir` (optional): directory where the output data will be stored.
- `individuals` (optional): path to a file containing the subset of individuals on which to run GWAS, one ID per line.
- `filename_patterns`: rules that determine the paths of the input and output files.
  - `genotype`:
    - `bed` (required): file name pattern containing that may contain the substring `{chromosome}` in which case it is replaced by the actual chromosome number. 
    - `bim` (required): idem previous.
    - `fam` (required): idem previous.
  - `phenotype`:
    - `phenotype_file` (required): file containing the phenotypes to run GWAS on.
    - `phenotype_file_delim` (optional, default: `,`): file delimiter.  
    - `phenotype_tmp_file` (optional): name of a temporary file that will be used as input to the GWAS software if the previous does not match the expected format.
    - `covariates` (not used yet).
  - `gwas` (required): name pattern for output files. It can contain the fields `{phenotype}` and `{suffix}`.
  - `gwas_suffix`: if `{suffix}` is present in the above field, this field should contain a string containing different .
- `options`: dictionary containing the names of the options as keys and values the names and values of the options, respectively. 
- `suffix_tokens`: nested dictionary, detailing the string that will be added to the output name for each of the options.
- `exec`: path of the executables.
  - `bgenie` (optional, default: `bgenie`): path to the BGENIE executable.
  - `plink` (optional, default: `plink`): path to the BGENIE executable.
  
#### Analyze results
The code for this task is contained in the folder `analysis/`.
