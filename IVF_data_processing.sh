#!/bin/bash

#############################################################
## Processing Lindgren In-House Ovary Data
## melodyjparker14@gmail.com - Nov 22
## This script processes our in-house, ovary, RNA-seq data.
## It inputs raw reads and a reference genome.
## It outputs a count matrix.
## It uses the hg38 108 reference genome.
## Although it can be run quickly as a single pipeline, it is recommended to run run-scripts separately, checking QC at each stage and adjusting accordingly
## The downstream analysis requires the metadata file "IVF_phenotypic_data_revised.csv".
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

# Check that the nested scripts called from this script are allocated more than 1 cpu when this only has 1 cpu
# Check script names have been changed  x
# Add runscripts dir before all script names  x 
# Check that only paired reads are being used after trimming
# Sort out the alignment script

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
# - create an index for the trimmed reads
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
STAR=star
RSEM_REF=rsem_ref
RSEM=rsem
QC=qc_results

###################################
# 3 - CALL SCRIPTS
###################################
# Have script stop if there is an error
set -e


### A. Make an index file
#
# $1 = Raw reads dir.
JOBA_ID=$(sbatch --parsable -p short run_scripts/make_index.sh $RAW_READS) 
echo "Creating sample index for the raw reads. Job ID: $JOBA_ID."


### B. Merge fastq files
# Dependent on job A
# This is not necessary for our samples, but edit if required
#
JOBB_ID=$(sbatch --parsable -p short -d afterok:$JOBA_ID run_scripts/merge.sh)
echo "Merging fastq files. Job ID: $JOBB_ID."


### C. QC the raw reads
# Dependent on job A (or B if fastq files need merging - edit if required)
#
# $1 = Raw reads dir.
# $2 = QC raw reads results dir.
JOBC_ID=$(sbatch --parsable -p short -d afterok:$JOBA_ID run_scripts/qc_reads.sh $RAW_READS $QC/qc_raw_results)  # Use this line if fastq merging step was skipped
echo "Running FastQC on raw reads. Job ID: $JOBC_ID."


### D. Aggregate the QC results and visualise in a report
# Dependent on job C
#
# $1 = Input/output dir
# $2 = Report name
JOBD_ID=$(sbatch --parsable -p short -d afterok:$JOBC_ID run_scripts/multiqc.sh $QC/qc_raw_results "fastqc_raw_reads")
echo "Running MultiQC to aggregate QC results. Job ID: $JOBD_ID."


### E. Trim off the adapters from the reads
# Dependent on job A (or B if fastq files need merging - edit if required)
#
# $1 = Raw reads dir.
# $2 = Trimmed reads dir.
JOBE_IE=$(sbatch --parsable -p short -d afterok:$JOBA_ID run_scripts/trim_adapters.sh $RAW_READS $TRIMMED_READS)
echo "Trimming adapters from raw reads. Job ID: $JOBE_ID."


### F. Make an index file for the trimmed reads
# Dependent on job E
#
# $1 = Trimmed reads dir.
JOBE_IF=$(sbatch --parsable -p short -d afterok:$JOBE_ID run_scripts/make_index.sh $TRIMMED_READS)
echo "Creating a sample index for the trimmed reads. Job ID: $JOBF_ID."


### G. QC the trimmed reads
# Dependent on job F
#
# $1 = Trimmed reads dir.
# $2 = QC trimmed reads results dir.
JOBG_ID=$(sbatch --parsable -p short -d afterok:$JOBF_ID run_scripts/qc_reads.sh $TRIMMED_READS $QC/qc_trimmed_results)  # Use this line if fastq merging step was skipped
echo "Running FastQC on trimmed reads. Job ID: $JOBG_ID."


### H. Aggregate the QC results and visualise in a report
# Dependent on job G
#
# $1 = Input/output dir
# $2 = Report name
JOBH_ID=$(sbatch --parsable -p short -d afterok:$JOBG_ID run_scripts/multiqc.sh $QC/qc_trimmed_results "fastqc_trimmed_reads")
echo "Running MultiQC to aggregate QC results. Job ID: $JOBH_ID."


### I. Generate a genome index for STAR
# Not dependent on anything
#
# $1 = Reference genome directory
# $2 = STAR index directory # star_index
JOBI_ID=$(sbatch --parsable -p short run_scripts/generate_genome_index.sh $REF_GENOME $STAR_INDEX)
echo "Generating STAR genome index. Job ID: $JOBI_ID."


### J. Align reads to reference genome using STAR
# Dependent on jobs F and I
#
# $1 = Trimmed reads dir.
# $2 = STAR output directory
# $3 = STAR index directory
# $4 = Reference genome dir
# $5 = Reference genome annotations filename
JOBJ_ID=$(sbatch --parsable -p short -d afterok:$JOBF_ID:afterok:$JOBI_ID run_scripts/alignment_star.sh $TRIMMED_READS $STAR $STAR_INDEX $REF_GENOME $ref_genome_gtf)
echo "Aligning reads to reference genome using STAR. Job ID: $JOBJ_ID."


### K. Perform Alignment QC
# Dependent on job J
# This is annother nested script
# Do we need arguments?
#
JOBK_ID=$(sbatch --parsable -p short -d afterok:$JOBJ_ID run_scripts/qc_alignment.sh)
echo "Running alignment QC script. Job ID: $JOBK_ID."


### L. Generate RSEM reference
# Not dependent on anything
#
# $1 = RSEM reference
# $2 = Reference genome annotations file: gencode.v42.primary_assembly.annotation.gtf
# $3 = Reference genome fasta file: GRCh38.primary_assembly.genome.fa
JOBL_ID=$(sbatch --parsable -p short run_scripts/create_rsem_ref.sh $RSEM_REF $ref_genome_gtf $ref_genome_fasta)
echo "Generating RSEM reference. Job ID: $JOBL_ID."


### M. Quantify genes and isoforms
# Dependent on jobs J and L
#
# $1 = Input: STAR output dir.
# $2 = Output: RSEM output dir.
# $3 = RSEM reference
JOBM_ID=$(sbatch --parsable -p short -d afterok:$JOBJ_ID:afterok:$JOBL_ID run_scripts/quantification_rsem.sh $STAR $RSEM $RSEM_REF)
echo "Quantifying genes and isoforms using RSEM. Job ID: $JOBM_ID."


### N. Generate count matrix
# Dependent on job M
# 
# $1 = RSEM output dir.
JOBN_ID=$(sbatch --parsable -p short -d afterok:$JOBM_ID run_scripts/generate_count_matrix.sh $RSEM)
echo "Generating count matrix. Job ID: $JOBN_ID."


### O. Run downstream analysis using DESeq2
# Dependent on job N
#
JOBO_ID=$(sbatch --parsable -p short -d afterok:$JOBN_ID RScript run_scripts/downstream_analysis_DESeq2.R)
echo "Running downstream analysis using DESeq2. Job ID: $JOBO_ID."


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
