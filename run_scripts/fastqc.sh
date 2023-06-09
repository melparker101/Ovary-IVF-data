#!/bin/bash

# ----------------------------------------------------------
# Script to create perform quality control checks on reads using fastqc
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J qc_reads
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


# Define variables
IN=$1  # Input dir
OUT=$2  # Output dir
DATA_TYPE=$3  # File type, e.g. fastq, bam
INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $IN/index.txt)

# Load modules
module load FastQC/0.11.9-Java-11

# Create output directory
mkdir -p $OUT

# Run FastQC
if [[ $DATA_TYPE == fastq ]]; then

  echo "Data type: $DATA_TYPE."
  fastqc $IN/"$INPUT_FILE"_R1_001*.fastq.gz -o $OUT &
  fastqc $IN/"$INPUT_FILE"_R2_001*.fastq.gz -o $OUT &
  wait
  
elif [[ $DATA_TYPE == bam ]]; then

  echo "Data type: $DATA_TYPE."
  fastqc $IN/"$INPUT_FILE"*.bam -o $OUT
  
else 
  echo "Data type not recognised. Please use \"fastq\" or \"bam\" as an argument."
  
fi


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
