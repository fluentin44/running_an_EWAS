#!/bin/bash

#SBATCH --job-name=<NAME OF JOB>            # Name of job - can be anything
#SBATCH --ntasks-per-node=1                 # Tasks per node
#SBATCH --mem=30G                           # Memory required
#SBATCH --nodes=1                           # Number of nodes requested (always 1)
#SBATCH --mail-type=ALL                     # mail alert at start, end and abortion of execution
#SBATCH --mail-user=moh1r17@soton.ac.uk     # send mail to this address
#SBATCH --time=60:00:00                     # walltime - time needed for job - 60:00:00 = 60 hours, max time allowed
#SBATCH --output=output_files/slurm-%j.out
#SBATCH --error=error_files/slurm-%j.err

module load gcc/10.3.0
module load R/4.1.1

cd "/mainfs/scratch/<USERNAME>/<PATH TO FOLDER CONTAINING JOBFILE and SCRIPT FILE>"

R CMD BATCH ./<NAME OF SCRIPT FILE YOU WANT TO RUN>
