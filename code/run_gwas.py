import os
import pandas as pd
import yaml
from subprocess import call
from .auxiliary import is_yaml_file, parseIntSet, is_number
# import logging


class GWAS_Run:

    def __init__(self, yaml_config_file):

        self.yaml_config_file = yaml_config_file
        config = self.unfold_config(yaml_config_file)

        # dir stands for directory
        # f stands for file
        # fp stands for file pattern

        # phenotype_list_f = config.get("phenotype_list", None)
        # if phenotype_list_f is not None:
        #     self.phenotypes = [x.strip() for x in open(phenotype_list_f)]
        # else:
        #     self.phenotypes = config.get("phenotypes_list").strip().split()

        phenotype_list = config["phenotype_list"]
        if os.path.exists(os.path.join(os.path.dirname(self.yaml_config_file), phenotype_list)):
            phenotype_list = os.path.join(os.path.dirname(self.yaml_config_file), phenotype_list)
            self.phenotypes_all = [x.strip() for x in open(phenotype_list)]
            self.phenotypes = self.phenotypes_all
        else:
            self.phenotypes = phenotype_list

        ## file name patterns
        fp_config = config["filename_patterns"]

        # self.indiv_f = fp_config.get("individuals", None)

        self.bed_fp = fp_config["genotype"]["bed"].format(data_dir=config["data_dir"])
        self.bim_fp = fp_config["genotype"]["bim"].format(data_dir=config["data_dir"])
        self.fam_fp = fp_config["genotype"]["fam"].format(data_dir=config["data_dir"])

        # file with phenotype values
        self.pheno_f = os.path.join(config["data_dir"], fp_config["phenotype"]["phenotype_file"])
        # temporal file with phenotype values
        self.pheno_f_tmp = os.path.join(config["data_dir"], fp_config["phenotype"]["phenotype_file_tmp"])

        self.gwas_fp = os.path.join(config["output_dir"], fp_config["gwas"])

        options = config["options"]
        self.filename_rules = config["suffix_tokens"]

        options_suffix_tokens = {}
        for k in options.keys():
            options_suffix_tokens[k] = self.filename_rules[k][options[k]]

        self.output_suffix = config["suffix"].format(**options_suffix_tokens)

        ## covariates list
        self.covariates_f = fp_config.get("covariates", None)
        if self.covariates_f is not None:
            self.covariate_list = [x.strip() for x in config["covariate_list"].split(",")]
            self.get_covariate_column_indices()

        ## paths to executable files
        ## self.sbat_exec = config["exec"]["sbat"]
        self.plink_exec = config["exec"]["plink"]

        # Generate
        self.generate_phenotype_file()

        chromosomes = parseIntSet(str(config["chromosomes"]))
        self.chromosomes = [x for x in chromosomes if is_number(x)]
        self.chromosomes += [x for x in chromosomes if x in ("X", "Y")]

        ################################

        # other options
        self.delete_temp = config.get("delete_temp", True)
        self.merge_chromosomes = config.get("merge_chromosomes", True)
        self.ids = config.get("individuals", None).format(data_dir=config["data_dir"])

        self.overwrite_output = config.get("overwrite_output", True)
        # self.generate_tabix = config.get("generate_tabix", False)
        # self.adjust_for_covariates =

        ## TODO: FIX THIS TEMPORARY PATCH
        self.config = config

    def get_covariate_column_indices(self):
        header_comps = open(self.covariates_f).readline().strip().split()
        self.cov_index_list = [str(header_comps.index(cov)-1) for cov in self.covariate_list]

    def __str__(self):
        return ("\n".join([
            "bed file pattern: %s" % self.bed_fp,
            "bim file pattern: %s" % self.bim_fp,
            "fam file pattern: %s" % self.fam_fp,
            "Phenotype file name: %s" % self.pheno_f if self.pheno_f is not None else self.pheno_fp,
            "GWAS file name: %s" % self.gwas_fp # if self.gwas_f is not None else self.gwas_fp
        ]))

    def generate_phenotype_file(self):
        '''
        Create a phenotype file with the format as required by plink
        IID|FID|phenotype1|phenotype2|...
        '''

        df = pd.read_table(self.pheno_f)
        df.rename(columns={"Subject ID": "IID"}, inplace=True)

        df["FID"] = df["IID"]  # duplicate ID column in FID column
        # print(df.head())

        cols = df.columns.to_list()
        cols = [cols[-1]] + [cols[0]] + cols[1:-1]

        # ID|FID|pheno1|pheno2|...

        cols = [x for x in cols if x in (["IID", "FID"] + self.phenotypes)]

        df = df[cols]

        # TODO: raise warning if `not_found` is not empty
        # not_found = {x for x in self.phenotypes if x not in cols}

        df.to_csv(self.pheno_f_tmp, sep="\t", na_rep="NA", index=False)


    def unfold_config(self, token):
        '''
        This function reads a yaml configuration file
        If some the fields are paths to other yaml files,
        it will load their contents as values of the corresponding keys
        '''
        if is_yaml_file(token):
            try:
                token = yaml.safe_load(open(os.path.join(os.path.dirname(self.yaml_config_file),token)))
            except:
                token = yaml.safe_load(open(token))
        if isinstance(token, dict):
            for k, v in token.items():
                # print("{0}: {1}".format(k, v))
                token[k] = self.unfold_config(v)
        return token

    def format_genotype_file(origin="bed", destiny="bgen"):

        '''
        run shell commands in order to convert from the origin format to the destiny format (the one required by the GWAS tool, e.g. plink or BGENIE)
        '''

        pass

        # TODO: add log messages (redirect error output?)
        if (origin, destiny) == ("bed", "bgen"):
            output_f = genotype_f[:-3] + "bgen"
            command = [plink, "--pfile", genotype_f, "--export", "bgen-1.2", "--out", output_f]
        else:
            return

        call(command)

    def filter_individuals(ids=None, conditions=None, fields=None, tmp=True):
        '''
        filter by ids, or a set of conditions imposed on a field or fields (e.g. gender, ethnicity and/or age).
        if tmp == True, the files generated herein will be erased in the end
        '''
        pass

    def execute(self, chromosome, phenotype):
        # TODO: add log messages (redirect error output?)
        print("  Processing chromosome {}...".format(chromosome))
        output_basename = self.gwas_fp.format(phenotype=phenotype, suffix="%s__chr%s" % (self.output_suffix, chromosome))
        # print("    Generating file {}.qassoc...".format(output_basename))
        command = [self.plink_exec,
                   "--bed", self.bed_fp.format(chromosome=chromosome),
                   "--bim", self.bim_fp.format(chromosome=chromosome),
                   "--fam", self.fam_fp.format(chromosome=chromosome),
                   "--pheno", self.pheno_f_tmp,
                   "--pheno-name", phenotype,
                   "--out", output_basename]

        if self.covariates_f is None:
            command += ["--assoc"]
            # self.output_suffix = ".qassoc"
        else:
            command += ["--linear",
                        "--covar", self.covariates_f,
                        "--parameters", ",".join(self.cov_index_list)]
            # self.output_suffix = ".assoc.linear"

        if self.ids is not None:
            command += ["--keep", self.ids]

        call(command)


    def extract_metadata(fields):
        command = ["ukbconv", ""]
        pass

    def adjust_phenotype(covariates):
        pass

    def manhattan_plot(gwas_results):
        pass

    def run(self):

        for phenotype in self.phenotypes:
            print("Processing {}...".format(phenotype))
            output_filename = self.gwas_fp.format(phenotype=phenotype, suffix=self.output_suffix)
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)
            if os.path.exists(output_filename) and not self.overwrite_output:
                print("Output file already exists, if a new run is desired please delete the previous file.")
                continue
            else:
                for chr in sorted(int(i) for i in self.chromosomes):
                    self.execute(str(chr), phenotype)

            if self.merge_chromosomes:

                with open(output_filename, "w") as merged_fh:
                    ## TODO: suffix "qassoc" is added ad hoc, it would be better to include as attribute in the class
                    partitioned_gwas_files = [self.gwas_fp.format(phenotype=phenotype, suffix="%s__chr%s.qassoc" % (self.output_suffix, chr)) for chr in self.chromosomes]
                    for i, file_chr in enumerate(partitioned_gwas_files):
                        try:
                            with open(file_chr) as chr_fh:
                                for j, line in enumerate(chr_fh):
                                    if (i == 0 and j == 0) or (j != 0):
                                        merged_fh.write("%s\n" % "\t".join(line.strip().split()))
                                os.remove(file_chr)
                        except:
                           from IPython import embed; embed()

        if self.delete_temp:
            os.remove(self.pheno_f_tmp)


def main(args):
    GWAS_Run(args.yaml_config_file).run()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Harmonize data, run GWAS and generate descriptive plots.")
    parser.add_argument("--yaml_config_file", default="config.yaml")
    args = parser.parse_args()

    main(args)
