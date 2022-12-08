#!/bin/bash

# ----------------------------------------------------------
# Script to run quality check on in-house IVF ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N fastqc

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


IN=$1  # Input dir
OUT=$2  # Output dir


# For BAM files:
# for f in *bam; do echo $f; done > index.txt


INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $IN/index.txt)


module load FastQC/0.11.9-Java-11

if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi

echo $INPUT_FILE
fastqc $IN/"$INPUT_FILE" -o $OUT &
fastqc $IN/"$INPUT_FILE -o $OUT &
wait

echo "FastQC finished."


# End of job script
