#!/bin/bash

# ----------------------------------------------------------
# Script to trim adapters from in-house IVF ovary RNA-seq data 
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N trim_adapters

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o logs/output.log
#$ -e logs/error.log

# Parallel environemnt settings
#  For more information on these please see the wiki
#  Allowed settings:
#   shmem
#   mpi
#   node_mpi
#   ramdisk
#$ -pe shmem 4
#$ -t 1:15

fastq=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/trimmed_reads2

# index is a file with list of 8 digit numbers (part of the file name)
# this means input file is the task_id'th line of that list
INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' fastq/index.txt)

module load Trimmomatic/0.39-Java-11

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE \
$fastq/"$INPUT_FILE"_R1_001.fastq.gz \
$fastq/"$INPUT_FILE"_R2_001.fastq.gz \
$OUT/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz \
$OUT/"$INPUT_FILE"_R1_001_trimmed_U.fastq.gz \
$OUT/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz \
$OUT/"$INPUT_FILE"_R2_001_trimmed_U.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10 \
MINLEN:20

echo "Adapter removed."
