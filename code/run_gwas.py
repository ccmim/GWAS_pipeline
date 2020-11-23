import os
import pandas as pd
import yaml
from subprocess import call
from code.auxiliary import *
# import logging


class GWAS_Run:

    def __init__(self, yaml_config_file):
        self.config = GWASConfig(yaml_config_file)
        # self.filename_rules = config["suffix_tokens"]

        # Generate
        self.generate_phenotype_file()

    def __str__(self):
        return ("\n".join([
            "bed file pattern: %s" % self.config.bed_fp,
            "bim file pattern: %s" % self.config.bim_fp,
            "fam file pattern: %s" % self.config.fam_fp,
            "Phenotype file name: %s" % self.pheno_f if self.pheno_f is not None else self.pheno_fp,
            "GWAS file name: %s" % self.gwas_fp # if self.gwas_f is not None else self.gwas_fp
        ]))

    def generate_phenotype_file(self, input_id_column="IID"):
        
        '''
        Create a phenotype file with the format as required by plink
        IID|FID|phenotype1|phenotype2|...
        '''

        df = pd.read_table(self.config.pheno_f, sep=',')

        df.rename(
            columns={input_id_column: "IID"},
            inplace=True
        )

        df["FID"] = df["IID"]  # duplicate ID column in FID column

        cols = df.columns.to_list()
        cols = [cols[-1]] + [cols[0]] + cols[1:-1] # reorder columns: IID|FID|phenotypes
        cols = [x for x in cols if x in (["IID", "FID"] + self.config.phenotypes)]
        df = df[cols]

        # TODO: raise warning if `not_found` is not empty
        # not_found = {x for x in self.phenotypes if x not in cols}

        df.to_csv(
            self.config.pheno_f_tmp,
            sep="\t", na_rep="NA", index=False
        )


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
            "--pheno", self.config.pheno_f_tmp,
            "--pheno-name", phenotype,
            "--out", self.output_file
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
            # self.output_suffix = ".assoc.linear"

        if self.config.indiv_f is not None:
            command += ["--keep", self.indiv_f]

        call(command)

    #####################################################################
    # def extract_metadata(fields): pass
    # def adjust_phenotype(covariates): pass
    # def manhattan_plot(gwas_results): pass

    def gwas_f_chr_suffix(self):
        return self.config.gwas_fp + "__chr{chr}"

    def create_gwas_dir(self, phenotype):
        os.makedirs(os.path.dirname(self.gwas_fp_).format(phenotype=phenotype), exist_ok=True)

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
            self.create_gwas_dir(phenotype)
            self.output_files_by_pheno = []

            for chromosome in self.sorted_chromosomes():
                self.output_file = self.gwas_fp_.format(phenotype=phenotype, chr=chromosome)
                if self.stop_because_existent():
                    print("Output file already exists, if a new run is desired please delete the previous file.")
                    continue
                else:
                    self.execute_one_gwas(chromosome, phenotype)
                    self.output_files_by_pheno.append(self.output_file + ".qassoc")

            if self.config.merge_chromosomes_flag:
                self.merge_chromosomes(phenotype)

        if self.config.delete_temp_flag:
            os.remove(self.config.pheno_f_tmp)


class GWASConfig:

    def __init__(self, config):
        
        if is_yaml_file(config):
          self.yaml_config_file = yaml_config_file
          config = unfold_config(yaml_config_file)

        self.data_dir = self.get_data_dir(config)
        self.pheno_f = self.get_pheno_file(config)
        self.pheno_f_tmp = self.get_pheno_tmp_file(config) # temporal file with phenotype values
        self.phenotypes = self.get_phenotypes(config)

        self.indiv_f = self.get_indiv_f(config)
        self.bed_fp, self.bim_fp, self.fam_fp = self.get_genotype_files(config)
        self.__suffix_tokens = config.get("suffix_tokens", {})
        self.gwas_fp = self.get_gwas_f(config)

        # other options (not required)
        self.delete_temp_flag = False # self.get_delete_tmp_flag(config)
        self.merge_chromosomes_flag = self.get_merge_chromosomes_flag(config)
        self.overwrite_output_flag = self.get_overwrite_output_flag(config)
        self.chromosomes = self.get_chromosomes(config)

        self.plink_exec = self.get_executable_paths(config)

        #TODO: Implement this
        self.covariates_f = None

        #TODO: Implement TABIX indexing
        # self.generate_tabix = config.get("generate_tabix", False)


    #TODO: use @property decorator
    def get_data_dir(self, config):
        return config.get("data_dir", None)

    def get_pheno_file(self, config):
        pheno_f = config["filename_patterns"]["phenotype"]["phenotype_file"]
        if os.path.isabs(pheno_f) or self.data_dir is None or os.path.exists(pheno_f):
            return pheno_f
        else:
            return os.path.join(self.data_dir, pheno_f)

    def get_pheno_tmp_file(self, config):
        pheno_f_tmp = config["filename_patterns"]["phenotype"]["phenotype_file_tmp"]
        if os.path.isabs(pheno_f_tmp) or self.data_dir is None:
            return pheno_f_tmp
        else:
            return os.path.join(self.data_dir, pheno_f_tmp)

    def get_phenotypes(self, config):

        repo_rootdir = get_repo_rootdir()
        yaml_dir = os.path.join(repo_rootdir, "config_files")
        phenotype_list = config.get("phenotype_list", None)
        # if None, use all the phenotypes
        if phenotype_list is None:
            phenotypes = open(self.pheno_f).readline().strip().split(",")
            print(phenotypes)
            #TODO: This column name might vary across datasets!
            phenotypes.remove("IID")
        elif isinstance(phenotype_list, str):
            phenotype_file = os.path.join(yaml_dir, phenotype_list)
            if os.path.exists(phenotype_file):
              phenotypes = [x.strip() for x in open(phenotype_file)]
        elif isinstance(phenotype_list, list):
            phenotypes = phenotype_list
        return phenotypes

    def get_indiv_f(self, config):
        return config["filename_patterns"].get("individuals", None)

    def get_genotype_files(self, config):
        """
        I am assuming that the genotype file name patterns contain "{data_dir}"
        :param config:
        :return:
        """
        bed_fp = config["filename_patterns"]["genotype"]["bed"].format(data_dir=config["data_dir"])
        bim_fp = config["filename_patterns"]["genotype"]["bim"].format(data_dir=config["data_dir"])
        fam_fp = config["filename_patterns"]["genotype"]["fam"].format(data_dir=config["data_dir"])
        return bed_fp, bim_fp, fam_fp

    def get_gwas_f(self, config):
        output_dir = config.get("output_dir", "output")
        gwas_fp = os.path.join(output_dir, config["filename_patterns"]["gwas"])
        #TODO: finish this!!!
        gwas_fp = gwas_fp.format(**self.__suffix_tokens)
        return gwas_fp

    def get_delete_tmp_flag(self, config):
        return config.get("delete_temp", True)

    def get_merge_chromosomes_flag(self, config):
        return config.get("merge_chromosomes", True)

    def get_overwrite_output_flag(self, config):
        return config.get("overwrite_output", True)

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

        # Not implemented yet.
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
