#!/bin/bash -l

# usage sbatch GSB_per_build_size.sh

# Define job name
#SBATCH -J GSB
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 60:00:00
# Define array length 
#SBATCH --array=1-1 
# Define cores  
#SBATCH -c 1

module load python

# List all modules
module list
#
echo =========================================================
echo SLURM job: submitted date = $(date)
date_start=$(date +%s)
hostname
echo Current directory: $(pwd)
echo "Print the following environmental variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"

#############################################

cd /home/hleboswe/hc-storage/

python3 /home/hleboswe/hc-storage/yeast_GSB_per_build_size/02_GSB_counts_per_build_size.py file_list_v2.txt NA NA NA 2 1 rev_ CONTAM_ phospho:STY:A