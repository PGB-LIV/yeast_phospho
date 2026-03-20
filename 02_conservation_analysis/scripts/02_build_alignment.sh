#!/bin/bash -l

# usage sbatch 02_build_alignment.sh

# Define job name
#SBATCH -J 02_build_alignment.sh
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 160:00:00
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

echo python3 02_Build_yeast_ptm_alignments_syntenic.py
python3 02_Build_yeast_ptm_alignments_syntenic.py
