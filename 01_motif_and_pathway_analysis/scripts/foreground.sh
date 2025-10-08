#!/bin/bash -l

# usage sbatch foreground.sh

# Define job name
#SBATCH -J foreground
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 60:00:00
# Define array length - number of files to search in samples.txt
#SBATCH --array=1-1 # !!!!!!! Set the second number to match exactly the count of raw files (lines) in samples.txt !!!!!
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

echo python3 /mnt/hc-storage/users/hleboswe/yeast_motif/03_motif_foreground.py
python3 /mnt/hc-storage/users/hleboswe/yeast_motif/03_motif_foreground.py