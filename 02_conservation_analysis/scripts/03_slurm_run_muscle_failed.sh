#!/bin/bash -l


# Define job name
#SBATCH -J muscle
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 20:00:00
# Define array length - number of files to search in samples.txt
#SBATCH --array=1-22
# Define cores
#SBATCH -c 10

echo =========================================================
echo SLURM job: submitted date = $(date)
date_start=$(date +%s)
hostname
echo Current directory: $(pwd)
echo "Print the following environmental variables:"
echo "Job name                     : $SLURM_JOB_NAME"
echo "Job ID                       : $SLURM_JOB_ID"
echo "Job array index              : $SLURM_ARRAY_TASK_ID"


FILE=`sed -n ${SLURM_ARRAY_TASK_ID}p samples_failed.txt`

OUTPUT=`basename -s _notaligned.fasta ${FILE}`.afa
muscle/muscle5.1.linux_intel64 -align ${FILE} -output ${OUTPUT} -threads 10

