import os 
from subprocess import call
def submit_sge_job(jobname, commands, dry_run=True, folder="jobs", stdout=None, stderr=None, memory_limit="16G", cwd=True, walltime=None):
  
  SGE_header = []

  # Run with current environment (-V) and in the current directory (-cwd)

  if cwd:
    SGE_header.append("#$ -cwd")
  
  if walltime is None:
    walltime = "48:00:00"
    
  SGE_header.append("#$ -l h_rt={}".format(walltime))
  SGE_header.append("#$ -l h_vmem={}".format(memory_limit))
  SGE_header.append("#$ -N {}".format(jobname))

  if stdout is None:
    stdout = "logs/{}.out".format(jobname)
          
  if stderr is None:
    stderr = "logs/{}.err".format(jobname)

  SGE_header.append("#$ -o {}".format(stdout))
  SGE_header.append("#$ -e {}".format(stderr))
    
  SGE_header.extend(commands)
    
  jobfile = os.path.join(folder, "{}.sh".format(jobname))

  if not os.path.exists(os.path.dirname(jobfile)):
    os.makedirs(os.path.dirname(jobfile), exist_ok=True)
    
  # I will assume stdout and stderr are written to the same folder
  if not os.path.exists(os.path.dirname(stdout)):
    os.makedirs(os.path.dirname(stdout), exist_ok=True)

  with open(jobfile, "wt") as jf:
    jf.write("\n".join(SGE_header))
        
  if not dry_run:
    call(["qsub", jobfile])

