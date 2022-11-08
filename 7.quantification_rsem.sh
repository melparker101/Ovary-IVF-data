#!/bin/bash

# ----------------------------------------------------------
# Based on Saskia Reibe's code
# Script to Perform gene and isoform quantification on
# in-house IVF ovary RNA-seq data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

# Specify a job name
#$ -N rna-seq_rsem

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe 

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
###$ -o output.log
###$ -e error.log

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

# Begin writing your script here
module purge
module load RSEM/1.3.2-foss-2018b

fastq=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads
REF=//well/lindgren/users/mzf347/ref_genomes/homo_sapiens/gencode/GRCh38.p13
IN=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/rsem

if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi

INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $fastq/index.txt)


rsem-calculate-expression --bam --paired-end -p 4 $IN/"$INPUT_FILE"_Aligned.toTranscriptome.out.bam $REF $OUT/rsem_"$INPUT_FILE"


# End of job script
