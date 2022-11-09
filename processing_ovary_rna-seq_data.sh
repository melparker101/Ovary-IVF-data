#!/bin/bash

#############################################################
## Processing Lindgren In-House Ovary Data
## melodyjparker14@gmail.com - Nov 22
## This script takes raw, paired fastq files, removes adapters, quality checks, then maps them to the hg38 108 reference genome. 
## From this, a count matrix is produced which then enables downstream analysis.
## Based on a workflow provided by Saskia Reibe.
#############################################################

###################################
# 0 - WORKFLOW
###################################
# 0. make index for raw fastq file names if it doesn't exist
# 1. make star genome index
# 2. make rsem ref
# wait till this is done
# 3. run script2.sh
### 1. fastq
### 2. trim adapters
### 3. fastq again
### 4. map reads (star)
### 5. quantification
# make count matrix
# 4. run script3.r
### 1. downstream analysis

###################################
# 1 - JOB SETUP
###################################

# Specify a job name
#$ -N create_rsem_ref

# Project name and target queue
#$ -P lindgren.prjc
#$ -q short.qe 

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
###$ -o logs/output.log
###$ -e logs/error.log

# Required slots
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

###################################
# 2 - DEFINE PATHS AND FILE NAMES
###################################
# Paths:
FASTQ=$1  # Path of the directory containing the raw reads
REF_GENOME=$2  # Path of the directory containing the reference genome filesref_genome
STAR_INDEX=star_index
RSEM=rsem

# FASTQ=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads
# REF_GENOME=//well/lindgren/users/mzf347/ref_genomes/homo_sapiens/gencode/GRCh38.p13
# STAR_INDEX=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star_index   
# RSEM=//well/lindgren/users/mzf347/alignment/ivf_cumulus/rsem

# File names:
ref_genome_fasta=GRCh38.primary_assembly.genome.fa
ref_genome_gtf=gencode.v42.primary_assembly.annotation.gtf

###################################
# 2 - LOAD MODULES
###################################
module load STAR/2.7.9a-GCC-11.2.0
module load RSEM/1.3.2-foss-2018b
module load MultiQC/1.9-foss-2019b-Python-3.7.4

###################################
# 3 - START MAIN CODE
###################################
# Make an index for the raw reads file names
if [ ! -f $FASTQ/index.txt ]
then
  for f in $fastq/IVF*R1*; do g="${f%_R1*}" ; echo ${g##*/} >> $FASTQ/index.txt ; done
  echo "Index file created."
else
  echo "Index file already exists."
fi


# Make a genome index using STAR
if [ ! -d "$STAR_INDEX" ]
then
  mkdir -p $STAR_INDEX
  echo "STAR index directory created."
else
  echo "STAR index directory already exists."
fi

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $STAR_INDEX \
--genomeFastaFiles $REF_GENOME/ref_genome_fasta \
--sjdbGTFfile $REF_GENOME/ref_genome_gtf \
--sjdbOverhang 149  # read-length - 1

# Build a transcript reference using RSEM.

if [ ! -d $RSEM_REF ]
then
  mkdir -p $RSEM_REF
  echo "RSEM reference directory created."
fi
  echo "RSEM reference directory already exists."

rsem-prepare-reference --gtf $REF_GENOME/ref_genome_gtf \
   				     $REF_GENOME/ref_genome_fasta $RSEM_REF/human
              
# The next steps will be run in parallel so require a separate script.

# Run this script to perform QC, trim adapters, run QC again, if required.
# sh quality_control.sh
multiqc qc_raw_results -o qc_raw_results
multiqc qc_trimmed_results -o qc_trimmed_results

# map reads to reference genome, then quantify genes and isoforms.
sh mapping_and_quantifying.sh

# Generate a count matrix from the quantification results.
rsem-generate-data-matrix $RSEM/*.genes.results >> full.count.genes.matrix.txt

# Run downstream analysis
sh downstream_analysis.r


