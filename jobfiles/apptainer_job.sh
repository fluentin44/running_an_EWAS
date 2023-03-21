#!/bin/bash

####################################
#  Iridis 5 slurm script template  
#                                  
#  Submit script: sbatch filename  
#                                  
####################################
#SBATCH --job-name=align_test               # Name
#SBATCH --ntasks=64                         # Number of processor cores (i.e. tasks)
#SBATCH --nodes=1                           # Number of nodes requested
#SBATCH --ntasks-per-node=64                # Tasks per node
#SBATCH --cpus-per-task=1                   # Threads per task
#SBATCH --time=60:00:00                     # walltime
#SBATCH --output=output_files/slurm-%j.out
#SBATCH --error=error_files/slurm-%j.err
#SBATCH --mail-type=ALL                     # mail alert at start, end and abortion of execution
#SBATCH --mail-user=<USERNAME>@soton.ac.uk     # send mail to this address

cd /scratch/<USERNAME>/<PROJECT_FOLDER>/

module load apptainer

apptainer exec jobfiles/methylation_2.0.sif r src/<script.r>
