#!/bin/bash

# ----------------------------------------------------------
# Script to generate a genome index for in-house IVF ovary RNA-seq data using STAR
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

# Specify a job name
#$ -N rna-seq_star

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

# Some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Task ID: $SGE_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# Load modules
module load STAR/2.7.9a-GCC-11.2.0

REF_GENOME=//well/lindgren/users/mzf347/ref_genomes/homo_sapiens/gencode/GRCh38.p13
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/ref

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $OUT \
--genomeFastaFiles $REF_GENOME/GRCh38.primary_assembly.genome.fa.gz \
--sjdbGTFfile $REF_GENOME/gencode.v42.primary_assembly.annotation.gtf.gz \
--sjdbOverhang 149  # read-length - 1

# End of job script
