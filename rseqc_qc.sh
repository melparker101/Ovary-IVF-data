#!/bin/bash

IN=$1  # Input dir
OUT=rseq
OUT=$2  # Output dir

# ----------------------------------------------------------
# Script to run RSeQC quality checks on in-house IVF ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N fastqc

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe

# Run the job in the current working directory
#$ -wd $IN 

# Instead of separate output and error log files, send the error output to the regular output log file
#$ -j y

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


# Create an index of file names
if [ ! -f file_index.txt ]; then
  for f in *bam; do echo $f; done > file_index.txt
fi


INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $IN/file_index.txt)


module load RSeQC/3.0.0-foss-2018b-Python-3.6.6


if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi


# See if we can run these in parallel too
inner_distance.py -i $INPUT_FILE -r ref/gencode.v42.primary_assembly.bed12 -o rseqc/$INPUT_FILE

echo "inner distance done"

junction_annotation.py -i $INPUT_FILE -r ref/gencode.v42.primary_assembly.bed12 -o rseqc/$INPUT_FILE

echo "junction annotation done"

junction_saturation.py -i $INPUT_FILE -r ref/gencode.v42.primary_assembly.bed12 -o rseqc/$INPUT_FILE

echo "junction saturation done"

read_distribution.py  -i $INPUT_FILE -r ref/gencode.v42.primary_assembly.bed12 > rseqc/$INPUT_FILE.read_dist.txt

echo "read dist done"

RNA_fragment_size.py -i $INPUT_FILE -r ref/gencode.v42.primary_assembly.bed12 > rseqc/$INPUT_FILE.frag_size.txt

echo "RNA frag size done"

bam_stat.py -i $INPUT_FILE > $INPUT_FILE.bam_stat.txt

echo "BAM stat done"


echo "RSeQC QC finished."


# End of job script
