#!/bin/bash

# ----------------------------------------------------------
# Script to run quality check on in-house IVF ovary data (trimmed reads)
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N qc_trimmed_reads

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

IN=$1  # trimmed_reads
OUT=$2  # qc_trimmed_results

# Merge index is a file with list of 8 digit numbers (part of the file name)
# This means input file is the task_id'th line of that list
INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $IN/index.txt)


module load FastQC/0.11.9-Java-11


echo $INPUT_FILE
fastqc $IN/"$INPUT_FILE"_R1_001_trimmed_U.fastq.gz -o $OUT
fastqc $IN/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz -o $OUT
fastqc $IN/"$INPUT_FILE"_R2_001_trimmed_U.fastq.gz -o $OUT
fastqc $IN/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz -o $OUT


echo "QC finished."


# End of script job
