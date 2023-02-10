#!/bin/bash

# ----------------------------------------------------------
# Script to create merge RNA reads
# Based on a script from Saskia Reibe
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------

# Specify a job name
#$ -N merge_RNA_reads

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe 

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
###$ -o logs/output.log
###$ -e logs/error.log

# Parallel environemnt settings
#  For more information on these please see the wiki
#  Allowed settings:
#   shmem
#   mpi
#   node_mpi
#   ramdisk
#$ -pe shmem 4
#$ -t 1-15

# Some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Task ID: $SGE_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

IN=index_dir
IN1=fastq_lane1
IN2=fastq_lane2
OUT=merged_fastq


INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $IN/index.txt)


cat $IN1/IVF_123456_"$INPUT_FILE"_1.fastq.gz $IN2/IVF_987654_"$INPUT_FILE"_1.fastq.gz > $OUT/"$INPUT_FILE"_R1.fastq.gz

cat $IN/IVF_987654_"$INPUT_FILE"_2.fastq.gz $IN2/IVF_987654_"$INPUT_FILE"_2.fastq.gz > $OUT/"$INPUT_FILE"_R2.fastq.gz


# End of job script
