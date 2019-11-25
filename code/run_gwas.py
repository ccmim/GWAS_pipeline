import os
import pandas as pd
import yaml
from subprocess import call
# import logging


# copied from http://thoughtsbyclayg.blogspot.com/2008/10/parsing-list-of-numbers-in-python.html
# with a few modifications
def parseIntSet(nputstr=""):
    selection = set()
    invalid = set()
    # tokens are comma separated values

    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        try:
            # autosomal chromosomes are integers between 1 and 22
            if int(i) <= 22 and int(i) >= 1:
                selection.add(i)
            else:
                invalid.add(i)
        except:
            # if the token is not a number, then it might be a range, like 5-9
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items seperated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token) - 1]
                    for x in range(first, last + 1):
                        selection.add(str(x))
            except:
                if i == "X" or i == "Y":
                    selection.add(i)
                else:
                    # if not an int nor a range nor X/Y...
                    invalid.add(i)

    # Report invalid tokens before returning valid selection
    # print("Invalid set: " + str(invalid))

    return selection

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class GWAS_Run:

    def __init__(self, yaml_config_file):

        config = yaml.safe_load(open(yaml_config_file))

        # dir stands for directory
        # f stands for file
        # fp stands for file pattern

        ## paths to data
        data_dir = config["data_dirs"]
        genotype_dir = data_dir["genotype"]
        pheno_dir = data_dir["phenotype"]
        output_dir = data_dir["output"]
        tmp_dir = data_dir["tmp"]

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        ## file name patterns
        fp_config = config["filename_pattern"]
        phenotype_list_f = fp_config.get("phenotype_list", None)
        # self.indiv_f = fp_config.get("individuals", None)
        self.genotype_fp = os.path.join(genotype_dir, fp_config["genotype"])
        self.bim_fp = os.path.join(genotype_dir, fp_config["bim"])
        self.fam_fp = os.path.join(genotype_dir, fp_config["fam"])
        self.pheno_f = os.path.join(pheno_dir, fp_config["phenotype"])
        self.gwas_fp = os.path.join(output_dir, fp_config["gwas"])
        self.pheno_f_tmp = os.path.join(tmp_dir, fp_config["phenotype_tmp"])


        self.covariates_f = fp_config.get("covariates", None)
        if self.covariates_f is not None:
            self.covariate_list = [x.strip() for x in config["covariate_list"].split(",")]
            self.get_covariate_columns()

        # paths to executable files
        # self.sbat_exec = config["exec"]["sbat"]
        self.plink_exec = config["exec"]["plink"]

        # other parameters of the run
        if phenotype_list_f is not None:
            self.phenotypes = [x.strip() for x in open(phenotype_list_f)]
        else:
            self.phenotypes =  fp.config_get("phenotypes").strip().split()
        self.generate_phenotype_file(self.phenotypes)

        chromosomes = parseIntSet(str(config["chromosomes"]))
        self.chromosomes = [x for x in chromosomes if is_number(x)] + [x for x in chromosomes if x in ("X", "Y")]

        ################################

        # other options
        self.delete_temp = config.get("delete_temp", True)
        self.merge_chromosomes = config.get("merge_chromosomes", True)
        # self.overwrite_output = config.get("overwrite_output", False)
        # self.generate_tabix = config.get("generate_tabix", False)
        self.ids = config.get("individuals", None)

    def get_covariate_columns(self):
        header_comps = open(self.covariates_f).readline().strip().split()
        print(header_comps)
        self.cov_index_list = [str(header_comps.index(cov)-1) for cov in self.covariate_list]

    def __str__(self):
        return ("\n".join([
            "Genotype file pattern: %s" % self.genotype_fp,
            "bim file pattern: %s" % self.bim_fp,
            "fam file pattern: %s" % self.fam_fp,
            "Phenotype file name: %s" % self.pheno_f,
            "GWAS file name: %s" % self.gwas_f
        ]))

    def generate_phenotype_file(self, phenotypes, tmp=True):
        '''
        Create a phenotype file with the format as required by plink
        IID|FID|phenotype1|phenotype2|...
        '''

        df = pd.read_table(self.pheno_f)
        df.rename(columns={"Subject.ID": "IID"}, inplace=True)
        df["FID"] = df["IID"]  # duplicate ID column in FID column
        # print(df.head())

        cols = df.columns.to_list()
        cols = [cols[-1]] + [cols[0]] + cols[1:-1]

        # ID|FID|pheno1|pheno2|...
        # os.path.isdir()
        # filter
        cols = [x for x in cols if x in (["IID", "FID"] + phenotypes)]
        not_found = {x for x in phenotypes if x not in cols}
        df = df[cols]
        # TODO: raise warning if `not_found` is not empty

        df.to_csv(self.pheno_f_tmp, sep="\t", na_rep="NA", index=False)

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
        output_basename = self.gwas_fp.format(phenotype=phenotype, suffix="__chr%s" % chromosome)
        command = [self.plink_exec,
                   "--bed", self.genotype_fp.format(chromosome=chromosome),
                   "--bim", self.bim_fp.format(chromosome=chromosome),
                   "--fam", self.fam_fp.format(chromosome=chromosome),
                   "--pheno", self.pheno_f_tmp,
                   "--pheno-name", phenotype,
                   "--out", output_basename]
        if self.covariates_f is None:
            command += ["--assoc"]
            self.output_suffix = ".qassoc"
        else:
            command += ["--linear",
                       "--covar", self.covariates_f,
                       "--parameters", ",".join(self.cov_index_list)]
            self.output_suffix = ".assoc.linear"

        if self.ids is not None:
            command += ["--keep", self.ids]

        # return(" ".join(command))
        # command = [plink, "--bed", bed_file, "--bim", bim_file, "--fam", fam_file, "--linear", "assoc", "--pheno", phenotype_f, "--pheno-name", phenotype, "--out", gwas_f]

        print("command")
        call(command)

    def extract_metadata(fields):
        command = ["ukbconv", ""]
        pass

    def adjust_phenotype(covariates):
        pass

    def manhattan_plot(gwas_results):
        pass


def run(args):

    gwas_run = GWAS_Run(args.yaml_config_file)

    for phenotype in gwas_run.phenotypes:

        for chr in gwas_run.chromosomes:
            gwas_run.execute(chr, phenotype)

        if gwas_run.merge_chromosomes:
            with open(gwas_run.gwas_fp.format(phenotype=phenotype, suffix=gwas_run.output_suffix), "w") as merged_fh:
                partitioned_gwas_files = [gwas_run.gwas_fp.format(phenotype=phenotype, suffix="__chr%s%s" % (chr, gwas_run.output_suffix)) for chr in gwas_run.chromosomes]
                for i, file_chr in enumerate(partitioned_gwas_files):
                     with open(file_chr) as chr_fh:
                         for j, line in enumerate(chr_fh):
                              if (i == 0 and j == 0) or (j != 0):
                                  merged_fh.write("%s\n" % "\t".join(line.strip().split()))
                         os.remove(file_chr)

    if gwas_run.delete_temp:
        os.remove(gwas_run.pheno_f_tmp)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Harmonize data, run GWAS and generate descriptive plots.")
    parser.add_argument("--yaml_config_file", default="config.yaml")
    args = parser.parse_args()

    run(args)