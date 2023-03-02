#!/bin/bash

# ----------------------------------------------------------
# Script to create perform quality control checks on raw reads
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J qc_raw_reads
#SBATCH -o output.out
#SBATCH -e error.err
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


IN=$1  # raw_reads
OUT=$2  # qc_raw_results

# The input file is the name on the task_id'th line of the index list
INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' $IN/index.txt)

module load FastQC/0.11.9-Java-11

echo $INPUT_FILE

fastqc $IN/"$INPUT_FILE"_R1_001.fastq.gz -o $OUT &
fastqc $IN/"$INPUT_FILE"_R2_001.fastq.gz -o $OUT &
wait


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
