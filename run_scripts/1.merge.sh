#!/bin/bash

# ----------------------------------------------------------
# Script to create merge RNA reads
# Edit the read names in the script based on the run ids
# Based on a script from Saskia Reibe
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------


#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J merge_reads
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


IN=$1  # index_dir
IN1=$2  # fastq_lane1
IN2=$3  # fastq_lane2
OUT=$4  # merged_fastq


INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' $IN/index.txt)


cat $IN1/IVF_123456_"$INPUT_FILE"_1.fastq.gz $IN2/IVF_987654_"$INPUT_FILE"_1.fastq.gz > $OUT/"$INPUT_FILE"_R1.fastq.gz

cat $IN/IVF_123456_"$INPUT_FILE"_2.fastq.gz $IN2/IVF_987654_"$INPUT_FILE"_2.fastq.gz > $OUT/"$INPUT_FILE"_R2.fastq.gz


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
