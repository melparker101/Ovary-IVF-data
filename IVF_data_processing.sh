#!/bin/bash

#############################################################
## Processing Lindgren In-House Ovary Data
## melodyjparker14@gmail.com - Nov 22
## This script processes our in-house, ovary, RNA-seq data.
## It inputs raw reads and a reference genome.
## It outputs a count matrix.
## It uses the hg38 108 reference genome.
## Although it can be run quickly as a single pipeline, it is recommended to run run-scripts separately, checking QC at each stage and adjusting accordingly
## Based on a workflow provided by Saskia Reibe.
#############################################################

###################################
# 0 - JOB SETUP
###################################
#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J main_ovary_pipeline
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err

# check that the nested scripts called from this script are allocated more than 1 cpu when this only has 1 cpu
echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

###################################
# 1 - WORKFLOW
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
# 2 - VARIABLES
###################################
# Inputs:
RAW_READS=$1  # Path of the directory containing the raw reads
REF_GENOME=$2  # Path of the directory containing the reference genome files  # ref_genomes/homo_sapiens/gencode/GRCh38.p13
ref_genome_fasta=GRCh38.primary_assembly.genome.fa  # Name of reference genome fasta file
ref_genome_gtf=gencode.v42.primary_assembly.annotation.gtf  # Name of corresponding annotations file

# Output directory names:
TRIMMED_READS=trimmed_reads
STAR_INDEX=star_index
STAR_OUT=star
RSEM_REF=rsem_ref
RSEM_OUT=rsem
QC=qc_results

###################################
# 3 - CALL SCRIPTS
###################################
# Have script stop if there is an error
set -e

# A. Make an index file
# $1 = Raw reads dir.
JOBA_ID=$(sbatch --parsable -p short make_index.sh $RAW_READS) 

# B. Merge fastq files
# This is not necessary for our samples, but add arguments if running
JOBB_ID=$(sbatch --parsable -p short -d afterok:$JOBA_ID merge.sh)

# C. QC the raw reads
# $1 = Raw reads dir.
# $2 = QC raw reads results dir.
# JOBC_ID=$(sbatch --parsable -p short -d afterok:$JOBB_ID qc_reads.sh $RAW_READS $QC)  # Use this line instead if the fastq merge step was used
JOBC_ID=$(sbatch --parsable -p short -d afterok:$JOBA_ID qc_reads.sh $RAW_READS $QC/qc_raw_results)  # Use this line if fastq merging step was skipped

# D. Aggregate the QC results and visualise in a report
# $1 = Input/output dir
# $2 = Report name
JOBD_ID=$(sbatch --parsable -p short -d afterok:$JOBC_ID multiqc.sh $QC/qc_raw_results "fastqc_reads")

# E. Trim off the adapters from the reads
# $1 = Raw reads dir.
# $2 = Trimmed reads dir.
JOBE_IE=$(sbatch --parsable -p short -d afterok:$JOBD_ID trim_adapters.sh $RAW_READS $TRIMMED_READS)

# F. Make an index file for the trimmed reads
# $1 = Trimmed reads dir.
JOBE_IF=$(sbatch --parsable -p short -d afterok:$JOBE_ID make_index.sh $TRIMMED_READS)

# G. QC the trimmed reads
# $1 = Trimmed reads dir.
# $2 = QC trimmed reads results dir.
JOBG_ID=$(sbatch --parsable -p short -d afterok:$JOBF_ID qc_reads.sh $TRIMMED_READS $QC/qc_trimmed_results)  # Use this line if fastq merging step was skipped

# H. Aggregate the QC results and visualise in a report
# $1 = Input/output dir
# $2 = Report name
JOBH_ID=$(sbatch --parsable -p short -d afterok:$JOBG_ID multiqc.sh $QC/qc_trimmed_results "fastqc_reads")

# I. Generate a genome index for STAR
# $1 = Reference genome directory
# $2 = STAR index directory # star_index
JOBI_ID=$(sbatch --parsable -p short -d afterok:$JOBH_ID generate_genome_index.sh $REF_GENOME $STAR_INDEX)

# J. Align reads to reference genome using STAR
# $1 = Trimmed reads dir.
# $2 = STAR output directory
# $3 = STAR index directory
# $4 = Reference genome dir
# $5 = Reference genome annotations filename
JOBJ_ID=$(sbatch --parsable -p short -d afterok:$JOBI_ID alignment_star.sh $TRIMMED_READS $STAR $STAR_INDEX $REF_GENOME $ref_genome_gtf)

# K. Perform AlignmentQC
# This is annother nested script
JOBK_ID=$(sbatch --parsable -p short -d afterok:$JOBJ_ID qc_alignment.sh)

# L. Generate RSEM reference
# $1 = RSEM reference
# $2 = Reference genome annotations file: gencode.v42.primary_assembly.annotation.gtf
# $3 = Reference genome fasta file: GRCh38.primary_assembly.genome.fa
# Do not wait for alignment QC to finish because it takes ages - start running asap
JOBL_ID=$(sbatch --parsable -p short create_rsem_ref.sh $RSEM_REF $ref_genome_gtf $ref_genome_fasta)

# M. quantification_rsem.sh

