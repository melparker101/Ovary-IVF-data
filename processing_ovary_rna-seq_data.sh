#!/bin/bash

#############################################################
## Processing Lindgren In-House Ovary Data
## melodyjparker14@gmail.com - Nov 22
## This script processes our in-house, ovary, RNA-seq data.
## It inputs raw reads and a reference genome.
## It outputs a count matrix.
## It uses the hg38 108 reference genome.
## Based on a workflow provided by Saskia Reibe.
#############################################################

###################################
# 0 - WORKFLOW
###################################
# - make index file for raw fastq file names
# - make STAR genome index
# - make RSEM transcript ref
# - perform quality control
# - run trim_reads.sh
# - use multiqc on fastqc results
# - run mapping_and_quantification.sh
# - generate count matrix
# - run downstream_analysis.sh

###################################
# 1 - JOB SETUP
###################################
# Specify a job name
#$ -N ovary_rna-seq_processing

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
# 2 - CONSTANTS
###################################
# Inputs:
RAW_READS=$1  # Path of the directory containing the raw reads
REF_GENOME=$2  # Path of the directory containing the reference genome files
ref_genome_fasta=GRCh38.primary_assembly.genome.fa
ref_genome_gtf=gencode.v42.primary_assembly.annotation.gtf

# Output directory names:
TRIMMED_READS=trimmed_reads
STAR_INDEX=star_index
STAR_OUT=star
RSEM_REF=rsem_ref
RSEM_OUT=rsem
QC=qc_results


###################################
# 3 - LOAD MODULES
###################################
module load STAR/2.7.9a-GCC-11.2.0
module load RSEM/1.3.2-foss-2018b
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4

###################################
# 4 - START MAIN CODE
###################################
# Have script stop if there is an error
set -e

# Make an index for the raw reads file names
if [ ! -f $RAW_READS/index.txt ]
then
  for f in $RAW_READS/IVF*R1*; do g="${f%_R1*}" ; echo ${g##*/} >> $RAW_READS/index.txt ; done
  echo "Index file created."
fi


# Make a genome index using STAR
if [ ! -p "$STAR_INDEX" ]
then
  mkdir -p $STAR_INDEX
fi

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $STAR_INDEX \
--genomeFastaFiles $REF_GENOME/ref_genome_fasta \
--sjdbGTFfile $REF_GENOME/ref_genome_gtf \
--sjdbOverhang 149  # read-length - 1

# Build a transcript reference using RSEM.
if [ ! -p $RSEM_REF ]
then
  mkdir -p $RSEM_REF
fi

rsem-prepare-reference --gtf $REF_GENOME/ref_genome_gtf \
   				     $REF_GENOME/ref_genome_fasta $RSEM_REF/human
              

# Perform QC on raw reads
mkdir $QC/qc_raw_results
fastqc $RAW_READS/*.fq -o $QC/qc_raw_results
multiqc $QC/qc_raw_results -o $QC/qc_raw_results
echo "QC on raw reads is finished."


# Trim reads
# This will be run in parallel so requires a separate script
if [ ! -p $TRIMMED_READS ]
then
  sh trim_reads.sh $RAW_READS $TRIMMED_READS
fi


# Perform QC on trimmed reads
mkdir $QC/qc_trimmed_results
fastqc $TRIMMED_READS/*.fq -o $QC/qc_trimmed_results
multiqc $QC/qc_trimmed_results -o $QC/qc_trimmed_results
echo "QC on trimmed reads is finished."


# map reads to reference genome, then quantify genes and isoforms.
sh mapping_and_quantification.sh -raw $RAW_READS -trimmed $TRIMMED_READS -ref $REF_GENOME \
                                 -star_index $STAR_INDEX -star $STAR_OUT \
                                 -rsem_ref $RSEM_REF -rsem $RSEM_OUT

# Generate a count matrix from the quantification results.
rsem-generate-data-matrix $RSEM_OUT/*.genes.results >> full.count.genes.matrix.txt

# Run downstream analysis
# sh downstream_analysis.r
