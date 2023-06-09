#!/bin/bash

# ----------------------------------------------------------
# Script to run MultiQC on in-house IVF ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J multiqc
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err

#  Parallel environment settings 
#  For more information on these please see the documentation 
#  Allowed parameters: 
#   -c, --cpus-per-task 
#   -N, --nodes 
#   -n, --ntasks 

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"


# Define variables
IN=$1  # input_dir
OUT=$1  # output_dir (same as input)
MODULES=$2
REPORT_NAME=multiqc_report_"$MODULES".html 

# Load modules
module load MultiQC/1.9-foss-2019b-Python-3.7.4

# Run MultiQC
multiqc -i "$IN" -n "$REPORT_NAME" -o "$OUT"


echo "###########################################################"
echo "MultiQC finished."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
