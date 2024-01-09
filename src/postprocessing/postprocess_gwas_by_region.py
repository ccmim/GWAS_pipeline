import argparse
import os, sys
from subprocess import call, check_output
import shlex

repo_rootdir = check_output(shlex.split("git rev-parse --show-toplevel")).strip().decode('ascii')
os.chdir(repo_rootdir)
sys.path.append(os.getcwd())

from utils.ARC_helpers.SGE_utils import *

def build_concatenate_command(args, phenotypes):
     
    experiment_id = args.experiment_id 
    run_id = args.run_id
        
    phenotypes = " ".join(phenotypes)
    command = "Rscript src/postprocessing/concatenate_gwas.R\n" 
    command += f"--experiment_ids {experiment_id} \n"
    command += f"--run_ids {run_id} \n"
    command += f"--input_results_dir {args.input_results_dir}/by_region \n"
    command += f"--output_filename_pattern {args.output_filename_pattern}\n"
    command += f"--phenotypes {phenotypes}\n"
    
    command = " ".join(shlex.split(command))
    return command
    

def generate_summary_and_figures(args, phenotypes):
    
    ### TODO: Use a configuration file for the file name patterns

    gwas_pattern = args.output_filename_pattern
    gwas_folder = os.path.dirname(gwas_pattern)
    gwas_file_pattern = os.path.basename(gwas_pattern)   
    
    phenotypes = " ".join(phenotypes)
    command  = "Rscript src/postprocessing/gwas_analysis.R\n"
    # command += "--output_folder .\n"
    command += f"--gwas_folder {gwas_folder}\n"        
    command += f"--gwas_pattern {gwas_file_pattern}\n"
    command += f"--phenotypes {phenotypes}\n"
    command += "--cache_rds\n"
    
    command = " ".join(shlex.split(command))
    
    return command
    
    #if config["phenotype_list"] is not None:
    # command += "--qqplot_pooled\n"    
    
    
def sort_gwas_file_command(config, phenotypes):
    
    gwas_filename_pattern = args.output_filename_pattern    
    
    commands = []
    for phenotype in phenotypes:
        command = "bash src/postprocessing/sort_by_chr_and_pos.sh " + gwas_filename_pattern.format(phenotype=phenotype) + ".tsv"
        commands.append(command)
        
    commands = " ".join(commands)
    return commands                
    

def build_commands(args, phenotype):
            
    phenotype = [phenotype]
    
    path_extension = "PATH=/home/home01/scrb/bin:$PATH"
    conda_command = "module load anaconda; source activate gwas"
    
    concat_command = build_concatenate_command(args, phenotype)
    postproc_command = generate_summary_and_figures(args, phenotype)
    sort_command = sort_gwas_file_command(args, phenotype)
    
    commands = [
      path_extension, 
      conda_command, 
      #concat_command, 
      #sort_command, 
      postproc_command
    ]
        
        
    jobname = f"job_{args.run_id}_{args.experiment_id}_{phenotype[0]}"
    
    return jobname, commands
    



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
     
    parser.add_argument("--experiment_id", help="MLflow experiment ID")
    parser.add_argument("--run_id", help="MLflow run ID")
    parser.add_argument("--phenotypes", nargs="+", default=None)
    parser.add_argument("--input_results_dir", default=None)
    parser.add_argument("--output_filename_pattern", default=None)    
    parser.add_argument("--dry-run", "--dry_run", default=False)    
           
    args = parser.parse_args()
    
    for phenotype in args.phenotypes:
    
        jobname, commands = build_commands(args, phenotype)
        
        #print(jobname)
        #print(commands)
        
        job_config = {
          "jobname": jobname,
          "memory_limit": "16G", 
          "walltime": "01:00:00",           
          "folder": f"runs/jobs/gwas_postproc/{args.experiment_id}/{args.run_id}/",
          "stderr": f"runs/logs/gwas_postproc/{args.experiment_id}/{args.run_id}/{jobname}.err", 
          "stdout": f"runs/logs/gwas_postproc/{args.experiment_id}/{args.run_id}/{jobname}.out"
        }
        
        jobfile, sge_message = submit_sge_job(
            commands=commands,
            **job_config,
            dry_run=args.dry_run,
        )        
        
        print(jobfile)
        
        #if not args.dry_run:
        #    jobid = jobid_from_sge_msg(sge_message)
        #    add_job_as_pending(jobid, job_ids_file)