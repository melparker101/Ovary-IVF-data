#!/bin/bash

#############################################################
## Index BAM files
## melodyjparker14@gmail.com - Dec 22
## This code create an index file for each BAM file using samtools
#############################################################

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J index_bam_files
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


# Define variables
IN=$1
INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' $IN/index.txt)

# Load modules
module load samtools/1.8-gcc5.4.0

# Index BAM files
samtools index "$INPUT_FILE" "$INPUT_FILE".bai


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
