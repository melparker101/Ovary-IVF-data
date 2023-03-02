#!/bin/bash

# ----------------------------------------------------------
# Script to trim Nextera adapters from in-house IVF ovary RNA-seq data 
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J remove_adapters
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


IN=$1  # raw_reads
OUT=$2  # trimmed_reads

INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' IN/index.txt)

module load Trimmomatic/0.39-Java-11

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE \
$IN/"$INPUT_FILE"_R1_001.fastq.gz \
$IN/"$INPUT_FILE"_R2_001.fastq.gz \
$OUT/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz \
$OUT/"$INPUT_FILE"_R1_001_trimmed_U.fastq.gz \
$OUT/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz \
$OUT/"$INPUT_FILE"_R2_001_trimmed_U.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10 \
MINLEN:20


echo "###########################################################"
echo "Adapter removed."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
