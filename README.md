# GWAS pipeline

This repository contains code to execute GWAS on UK Biobank data, using the Plink and BGENIE tools, and perform downstream analysis on the results.

The folder `src/` contains scripts to perform pre-processing of the data and GWAS execution.
The folder `analysis/` contains scripts to perform statistical analysis on the GWAS results and generate figures (like Manhattan plots or Q-Q plot).
The folder `download_data/` contains scripts to download data from the UK Biobank and filter the genotype files.

## Requirements
### Software environment
A Dockerfile is provided for running most of the code in this repository. You can also download the [corresponding Docker image](https://hub.docker.com/r/rbonazzola/gwas_pipeline) from DockerHub.
Alternatively, you can use the Dockerfile as a guide to build your own environment without Docker.

The image is based on Ubuntu 22.04 and contains the following tools:
- R 4.3.1
- Python 3.11
- qctool 2.2.0
- bgenie 1.3
- plink 1.9
- GreedyRelated (to remove related subjects)

If your are on a platform on which Singularity is supported (but Docker isn't), you should be able to convert the Docker image into a Singularity SIF file directly from DockerHub, by using the following command: 

```bash
singularity build <SING_IMAGE_NAME>.sif docker-daemon://<DOCKER_IMAGE_NAME>:<TAG>
```

For computing genetic PCs, a different Docker image is used, namely the one provided by the author of `flashpca`: see https://github.com/gabraham/flashpca/blob/master/docker.md.

For LocusZoom plots, another Docker image will be provided soon.

## Usage
The pipeline consists of scripts for:

1) **Fetching the data**: download genotype and covariate data from the UK Biobank.
2) **Data pre-processing**: 
	- 2a. Genetic data: filtering out subjects and genetic variants and, optionally, splitting the genome into small regions to ease with subsequent parallel processing.
	- 2b. Generate genetic PCs.
	- 2c. Phenotypic data: adjusting the phenotypes for a set of covariates and performing inverse-normalization on the phenotypic scores (custom R script)
	- 2d. Filtering for related subjects and other characteristics: scripts to produce the final files that will be the input to the GWAS.
3) **Executing GWAS**: self-explanatory. Currently, we support Plink and BGENIE.
4) **Compile results and generate figures**: compile the GWAS results and generate Manhattan plots and Q-Q plots. You can also generate LocusZoom plots.
5) **Downstream analysis**: Integrate data from other sources in order to interpret the results. We currently support: proximity analysis using `biomaRt`, gene ontology term enrichment with `g:Profiler`, transcriptome-wide association studies with `S-PrediXcan` and pleiotropy analysis using the IEU GWAS Open Project.

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

`python main.py -c <YAML_FILE>`

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
