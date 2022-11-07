#!/bin/bash

# ----------------------------------------------------------
# Script to run quality check on in-house IVF ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N qc

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o logs/
#$ -e logs/

# Parallel environemnt settings
#  For more information on these please see the wiki
#  Allowed settings:
#   shmem
#   mpi
#   node_mpi
#   ramdisk
#$ -pe shmem 4
#$ -t 1-15

fastq=//well/lindgren/users/mzf347/alignment/ivf_cumulus/trimmed_reads
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/qc_trimmed_results

# merge index is a file with list of 8 digit numbers (part of the file name)
# this means input file is the task_id'th line of that list
INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $fastq/../merge_index)

module load FastQC/0.11.9-Java-11

echo $INPUT_FILE
fastqc $fastq/IVF00"$INPUT_FILE"_R1_001_trimmed_U.fastq.gz -o $OUT
fastqc $fastq/IVF00"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz -o $OUT
fastqc $fastq/IVF00"$INPUT_FILE"_R2_001_trimmed_U.fastq.gz -o $OUT
fastqc $fastq/IVF00"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz -o $OUT

echo "QC finished."

