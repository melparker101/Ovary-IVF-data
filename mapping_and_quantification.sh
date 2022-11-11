#!/bin/bash

#############################################################
## Processing Lindgren In-House Ovary Data
## melodyjparker14@gmail.com - Nov 22
## This script maps reads to a reference genome 
## and quantifies genes and isoforms.
## It inputs raw reads and a reference genome.
## It outputs a count matrix.
## It uses the hg38 108 reference genome.
## Based on a workflow provided by Saskia Reibe.
#############################################################

###################################
# 0 - WORKFLOW
###################################
# - map reads to reference genome
# - quantify genes and isoforms

###################################
# 1 - JOB SETUP
###################################
# Specify a job name
#$ -N mapping_and_quantification.sh

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

###################################
# 2 - CONSTANTS
###################################
while getopts raw:trimmed:ref:star_index:star:rsem_ref:rsem: flag
do
    case "${flag}" in
        raw) RAW=${OPTARG};;
        trimmed) TRIMMED=${OPTARG};;
        ref) REF=${OPTARG};;
        star_index) STAR_INDEX=${OPTARG};;
        star) STAR=${OPTARG};;
        rsem_ref) RSEM_REF=${OPTARG};;
        rsem) RSEM=${OPTARG};;
    esac
done

# RAW_READS=$1  # path of raw reads
# TRIMMED_READS=$2  # path of trimmed reads
# REF=$4  # path of reference genome
# STAR_INDEX=$3  # path of STAR genome index 
# STAR_OUT=$4  # path of output directory
# RSEM_REF=$5
# RSEM=$6

# RAW_READS=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads
# IN=//well/lindgren/users/mzf347/alignment/ivf_cumulus/trimmed_reads
# REF=//well/lindgren/users/mzf347/ref_genomes/homo_sapiens/gencode/GRCh38.p13
# STAR_INDEX=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star_index
# STAR_OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star

INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' RAW_READS/index.txt)

###################################
# 3 - LOAD MODULES
###################################
module load STAR/2.7.9a-GCC-11.2.0
module load RSEM/1.3.2-foss-2018b

###################################
# 4 - START MAIN CODE
###################################
# Map reads to reference genome
mkdir -p $STAR

STAR --runThreadN 6 \
--readFilesIn $TRIMMED/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz $IN/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix $STAR_OUT/$INPUT_FILE \
--genomeDir $STAR_INDEX \
--sjdbGTFfile $REF/gencode.v42.primary_assembly.annotation.gtf --outSJfilterReads Unique --sjdbOverhang 149 \
--outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --chimSegmentMin 20 \
--outSAMtype BAM SortedByCoordinate --outSAMattributes All \
--quantMode TranscriptomeSAM  # output alignments translated into transcript coordinates

echo "Mapping finished."


# Quantify genes and isoforms
mkdir -p $RSEM
rsem-calculate-expression --bam --paired-end -p 4 $STAR_OUT/"$INPUT_FILE"_Aligned.toTranscriptome.out.bam $RSEM_REF/human $RSEM/rsem_"$INPUT_FILE"

echo "Quantification finished."


# End of job script
