#!/bin/bash

# ----------------------------------------------------------
# Script to make an index file for in-house IVF RNA-seq ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J make_index
#SBATCH -o output.out
#SBATCH -e error.err

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
IN=$1  # fastq files/reads
DATA_TYPE=$2  # fastq or bam

# Create an index file in the fastq directory containing a list of fastq sample names
if [[ $DATA_TYPE == fastq ]]; then

  for f in "$IN"/IVF*R1*.fastq.gz; do echo $(basename${f%_R1_001*}) >> "$IN"/index.txt ; done
  
elif [[ $DATA_TYPE == bam ]]; then 

  for f in "$IN"/IVF*.bam; do echo $(basename ${f%Aligned*}) >> "$IN"/index.txt ; done
  
else
  echo "Data type not recognised. Please use \"fastq\" or \"bam\"."
  
fi


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
