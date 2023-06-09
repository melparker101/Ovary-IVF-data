#!/bin/bash

# ----------------------------------------------------------
# Script to run quality check on in-house IVF ovary data (trimmed reads)
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J qc_trimmed_reads
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err
#SBATCH -a 1-15

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


IN=$1  # trimmed_reads
OUT=$2  # qc_trimmed_results

# Merge index is a file with list of 8 digit numbers (part of the file name)
# This means input file is the task_id'th line of that list
INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' $IN/index.txt)

module load FastQC/0.11.9-Java-11

echo $INPUT_FILE
fastqc $IN/"$INPUT_FILE"_R1_001_trimmed_U.fastq.gz -o $OUT
fastqc $IN/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz -o $OUT
fastqc $IN/"$INPUT_FILE"_R2_001_trimmed_U.fastq.gz -o $OUT
fastqc $IN/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz -o $OUT


echo "###########################################################"
echo "QC finished."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
