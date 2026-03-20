#!/bin/bash -l

# usage sbatch background_tyrosine.sh

# Define job name
#SBATCH -J background
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o slurm_%J.%N.out
# Define a standard error file
#SBATCH -e slurm_%J.%N.err
# Define time limit
#SBATCH -t 60:00:00
# Define array length - number of files 
#SBATCH --array=1-48 # !!!!!!! Set the second number to match exactly the count of  files (lines) in .txt file!!!!!
# Define cores  
#SBATCH -c 1

module load python/3.10.0

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
# change directory (also contains files.txt)
cd /home/hleboswe/hc-storage/

# output directory for files per experiment
OUTPUT_DIR="/home/hleboswe/hc-storage/yeast_motif/background_per_file_tyrosine"

# file_list_v2.txt contains list of FDR files
# make this in linux as windows and linux have different end of line characters!
# for n in array, get 15mers
FILE=`sed -n ${SLURM_ARRAY_TASK_ID}p file_list_v2.txt`

# PXD as basename
BASE_NAME=$(echo "$FILE" | tr '/' '_')
# output file: background_PXD.csv
OUTPUT_FILE="$OUTPUT_DIR/tyr_background_${BASE_NAME}.csv"
echo "Generating background: $FILE -> $OUTPUT_FILE"

# Run the Python script
python3 /mnt/hc-storage/users/hleboswe/yeast_motif/06_15mer_background_tyrosine.py "$FILE" "$OUTPUT_FILE"