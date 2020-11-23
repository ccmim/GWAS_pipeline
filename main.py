from code.run_gwas import GWAS_Run

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Harmonise data, run GWAS and generate descriptive plots.")
    parser.add_argument("--yaml_config_file", "-c", default="config.yaml")
    parser.add_argument("--phenotype_file", default=None)
    parser.add_argument("--gwas_file", default=None)
    parser.add_argument("--ids_file", default=None)
    args = parser.parse_args()

    from code.auxiliary import unfold_config
    config = unfold_config(args.yaml_config_file)
    
    def overwrite_config_items(config, args):
        for attr, value in args.__dict__.items():
            if attr in config.keys() and value is not None:
                config[attr] = value

    print(config)

    if args.phenotype_file is not None:
        config["filename_patterns"]["phenotype"]["phenotype_file"] = args.phenotype_file
        config["filename_patterns"]["phenotype"]["phenotype_file"]


    if args.gwas_file is not None:
        config["filename_patterns"]["gwas"] = args.gwas_file

    if args.ids_file is not None:
        config["individuals"] = args.ids_file
        
    
    gwas_run = GWAS_Run(config)

    gwas_run.run()