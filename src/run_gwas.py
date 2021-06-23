import os
import re
import pandas as pd
import yaml
import subprocess
from .auxiliary import *
import shutil
import glob
# import logging


class GWAS_Run:

    def __init__(self, yaml_config_file):
        self.config = GWASConfig(yaml_config_file)
        # self.filename_rules = config["suffix_tokens"]

    def __str__(self):
        return ("\n".join([
            "bed file pattern: %s" % self.config.bed_fp,
            "bim file pattern: %s" % self.config.bim_fp,
            "fam file pattern: %s" % self.config.fam_fp,
            "Phenotype file name: %s" % self.pheno_f if self.pheno_f is not None else self.pheno_fp,
            "GWAS file name: %s" % self.gwas_fp # if self.gwas_f is not None else self.gwas_fp
        ]))

    def execute_one_gwas(self, chromosome, phenotype):
        # TODO: add log messages (redirect error output?)

        # output_basename = self.gwas_fp.format(phenotype=phenotype, suffix="%s__chr%s" % (self.output_suffix, chromosome))
        # print("    Generating file {}.qassoc...".format(output_basename))

        print("  Processing chromosome {}...".format(chromosome))
        command = [
            self.config.plink_exec,
            "--bed", self.config.bed_fp.format(chromosome=chromosome),
            "--bim", self.config.bim_fp.format(chromosome=chromosome),
            "--fam", self.config.fam_fp.format(chromosome=chromosome),
            "--pheno", self.config.pheno_f,
            "--pheno-name", phenotype,
            "--out", self.output_file,
            "--threads", "8"
        ]

        if self.config.covariates_f is None:
            command += ["--assoc"]
            # self.output_suffix = ".qassoc"
        else:
            command += [
                "--linear",
                "--covar", self.covariates_f,
                "--parameters", ",".join(self.cov_index_list)
            ]
        

        if self.config.indiv_f is not None:
            command += ["--keep", self.config.indiv_f]

        # Quality controls    
        if self.config.qc["hwe_pval_thres"] is not None:
            command += ["--hwe", str(self.config.qc["hwe_pval_thres"])]
        if self.config.qc["snp_missing_rate_thres"] is not None:
            command += ["--geno", str(self.config.qc["snp_missing_rate_thres"])]
        if self.config.qc["sample_missing_rate_thres"] is not None:
            command += ["--mind", str(self.config.qc["sample_missing_rate_thres"])]
        if self.config.qc["maf_thres"] is not None:
            command += ["--maf", str(self.config.qc["maf_thres"])]

        #TODO: add option to output PLINK messages to the console
        subprocess.call(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    #####################################################################

    def gwas_f_chr_suffix(self):
        return self.config.gwas_fp + "__chr{chr}"

    def create_gwas_dir(self, phenotype):        
        gwas_dir = os.path.dirname(self.gwas_fp_).format(phenotype=phenotype)        
        os.makedirs(gwas_dir, exist_ok=True)
        os.makedirs(os.path.join(gwas_dir, "logs"), exist_ok=True)
        return gwas_dir

    def sorted_chromosomes(self):
        return [ str(x) for x in sorted([int(y) for y in self.config.chromosomes]) ]

    def stop_because_existent(self):
        return os.path.exists(self.output_file) and not self.config.overwrite_output

    def merge_chromosomes(self, phenotype):
        self.merged_output_file = self.config.gwas_fp.format(phenotype=phenotype) + ".tsv"
        merged_fh = open(self.merged_output_file, "w")
        # TODO: suffix "qassoc" is added ad hoc, it would be better to include as attribute in the class
        for i, file_chr in enumerate(self.output_files_by_pheno):
            chr_fh = open(file_chr)
            for j, line in enumerate(chr_fh):
                if (i == 0 and j == 0) or (j != 0):
                    merged_fh.write("%s\n" % "\t".join(line.strip().split()))
            chr_fh.close()
            os.remove(file_chr)
        merged_fh.close()


    #### EXECUTE ALL THE GWAS ####
    def run(self):

        import datetime

        for phenotype in self.config.phenotypes:

            self.timestamp = datetime.datetime

            print("Processing {}...".format(phenotype))
            self.gwas_fp_ = self.gwas_f_chr_suffix()
            gwas_dir = self.create_gwas_dir(phenotype)
            self.output_files_by_pheno = []

            for chromosome in self.sorted_chromosomes():
                self.output_file = self.gwas_fp_.format(phenotype=phenotype, chr=chromosome)
                if self.stop_because_existent():
                    print("Output file already exists, if a new run is desired please delete the previous file.")
                    continue
                else:
                    self.execute_one_gwas(chromosome, phenotype)
                    # This is specific for PLINK, need to modify it in order to use BGEN
                    self.output_files_by_pheno.append(self.output_file + ".qassoc")

            if self.config.merge_chromosomes_flag:
                self.merge_chromosomes(phenotype)

            for log_file in glob.glob(os.path.join(gwas_dir, "*log")):                
                log_file_basename = os.path.basename(log_file)
                shutil.move(log_file, os.path.join(gwas_dir, "logs", log_file_basename))

            for irem_file in glob.glob(os.path.join(gwas_dir, "*irem")):                
                irem_file_basename = os.path.basename(irem_file)
                shutil.move(irem_file, os.path.join(gwas_dir, "logs", irem_file_basename))

        if self.config.delete_temp_flag:
            os.remove(self.config.pheno_f)

        


class GWASConfig:

    def __init__(self, config):
        
        '''
          Parameters (dict, string): configuration parameters for GWAS. It can be:
            1. a dictionary or 
            2. the path to a YAML file
        '''

        # if `config` is a path to a YAML, store the path in an attribute
        if is_yaml_file(config):
            self.yaml_config_file = config
            config = yaml.load(open(config))

        # TODO: when the dict key is a tuple for which element order is irrelevant, use frozenset instead.            
        # self.suffix = config["suffix"].format(**tokens)
        
        # End of actions that need to be performed prior to unnesting the nested yaml files

        config = unfold_config(config)

        self.data_dir = config.get("data_dir", None)
        
        # self.pheno_f_tmp = self.get_pheno_tmp_file(config) # temporal file with phenotype values

        self.pheno_f = self.get_pheno_file(config)
        # self.pheno_f_tmp = os.path.join(self.tmpdir, "phenotypes.tsv") # temporal file with phenotype values
        self.phenotypes = self.get_phenotypes(config)
        self.tmpdir = config["filenames"]["tmpdir"]
        os.makedirs(self.tmpdir, exist_ok=True)

        self.indiv_f = self.get_indiv_f(config)
        
        self.bed_fp, self.bim_fp, self.fam_fp = self.get_genotype_files(config)        
        self.chromosomes = self.get_chromosomes(config)
        self.gwas_fp = config["filename_patterns"]["gwas"]
        # self.get_gwas_f(config)
        
        self.plink_exec = self.get_executable_paths(config)

        self.qc = self.get_qc(config)

        # other options (not required)
        self.delete_temp_flag = config.get("delete_temp", True)
        self.merge_chromosomes_flag = config.get("merge_chromosomes", True)
        self.overwrite_output_flag = config.get("overwrite_output", True)
        
        #TODO: Implement this
        self.covariates_f = None

        #TODO: Implement TABIX indexing
        # self.generate_tabix = config.get("generate_tabix", False)

    #TODO: use @property decorator    
    def get_pheno_file(self, config):
        pheno_f = config["filenames"]["phenotype_intermediate"]
        if os.path.isabs(pheno_f) or self.data_dir is None or \
           pheno_f.startswith(self.data_dir) or os.path.exists(pheno_f):
           return pheno_f
        else:
            return os.path.join(self.data_dir, pheno_f)

    def get_phenotypes(self, config):

        repo_rootdir = get_repo_rootdir()
        yaml_dir = os.path.join(repo_rootdir, "config_files")
        phenotype_list = config.get("phenotype_list", None)
        # if None, use all the phenotypes
        if phenotype_list is None:
            #TODO: by default, guess the column separator
            phenotypes = open(self.pheno_f).readline().strip().split("\t")
            #TODO: This column name might vary across datasets, it's better to provide it as part of the phenotype configuration
            # print(phenotypes)
            phenotypes.remove("IID")
            phenotypes.remove("FID")
        elif isinstance(phenotype_list, str):
            phenotype_file = os.path.join(yaml_dir, phenotype_list)
            if os.path.exists(phenotype_file):
              phenotypes = [x.strip() for x in open(phenotype_file)]
        elif isinstance(phenotype_list, list):
            phenotypes = phenotype_list
        return phenotypes

    def get_indiv_f(self, config):

        """
          Create a temporal file containing the ID's of the subjects to be included
          and return its path
        """

        # TODO: establish default behaviour for when white list is not provided.
        sample_white_lists = config.get("sample_white_lists", None)
        sample_black_lists = config.get("sample_black_lists", None)

        wl = []

        if sample_white_lists is not None:
            for file in sample_white_lists:
                wl.append(set(pd.read_csv(file, sep="\t").iloc[:,0]))
            wl = set.intersection(*wl)
        else:
            return None

        if sample_black_lists is not None:
            for file in sample_black_lists:
                bl.append(set(pd.read_csv(file, sep="\t").iloc[:,0]))
            bl = set.union(*bl)
        else:
            bl = set()
        
        wl = wl - bl

        indiv_f = os.path.join(self.tmpdir, "subjects.txt")
        with open(indiv_f, "w") as indiv_fh:
          indiv_fh.write("\n".join([str(x) + "\t" + str(x) for x in wl]))

        return indiv_f # config["filename_patterns"].get("individuals", None)

    def get_genotype_files(self, config):
        """
        I am assuming that the genotype file name patterns contain "{data_dir}"
        :param config:
        :return:
        """
        bed_fp = config["filename_patterns"]["genotype"]["bed"] #.format(data_dir=config["data_dir"])
        bim_fp = config["filename_patterns"]["genotype"]["bim"] #.format(data_dir=config["data_dir"])
        fam_fp = config["filename_patterns"]["genotype"]["fam"] #.format(data_dir=config["data_dir"])
        return bed_fp, bim_fp, fam_fp

    def get_gwas_f(self, config):
        # output_dir = config.get("output_dir", "output")
        # gwas_fp = os.path.join(output_dir, config["filename_patterns"]["gwas"])
        #TODO: finish this!!!
        gwas_fp = gwas_fp.format(**self.__suffix_tokens)
        return gwas_fp

    def get_chromosomes(self, config):
        chromosomes = parseIntSet(str(config["chromosomes"]))
        chromosomes = [x for x in chromosomes if is_number(x)]
        chromosomes += [x for x in chromosomes if x in ("X", "Y")]
        return chromosomes

    def get_executable_paths(self, config):
        ## paths to executable files
        try:
            plink_exec = config["exec"]["plink"]
        except:
            plink_exec = "plink"
        return plink_exec

    def generate_output_name(self):
        pass

    def get_covariates(self, config):

        # Not fully implemented yet.
        # Currently adjustment for covariates is performed separately in R.
        # TODO: finish this!

        self.covariates_f = config["filename_patterns"].get("covariates", None)

        def get_covariate_column_indices(self):
            header_comps = open(self.covariates_f).readline().strip().split()
            self.cov_index_list = [str(header_comps.index(cov) - 1) for cov in self.covariate_list]

        ## covariates list
        if self.covariates_f is not None:
            self.covariate_list = [x.strip() for x in config["covariate_list"].split(",")]
            self.get_covariate_column_indices()

    def get_options(self, config):
        options = config["options"]
        options_suffix_tokens = {}
        if options is not None:
            for k in options.keys():
                options_suffix_tokens[k] = self.filename_rules[k][options[k]]
                self.output_suffix = config["suffix"].format(**options_suffix_tokens)
        else:
            self.output_suffix = config["suffix"]

    def get_qc(self, config):
        qc_params = ( "hwe_pval_thres", "snp_missing_rate_thres", "sample_missing_rate_thres", "maf_thres")
        return { x: config["quality_control"].get(x, None) for x in qc_params }
