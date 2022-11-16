#!/bin/bash

#############################################################
## Processing Lindgren In-House Ovary Data - Trimming Reads
## melodyjparker14@gmail.com - Nov 22
## This script trims raw reads/removes adapters from raw reads.
## It inputs raw reads and outputs trimmed reads.
## Default is for Nextera PE adapters.
## It requires an index file in the input directory.
#############################################################

###################################
# 1 - JOB SETUP
###################################
# Specify a job name
#$ -N trim_reads.sh

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
IN=$1  # raw reads directory
OUT=$2  # trimmed reads directory
ADAPTER=NexteraPE-PE.fa

INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $IN/index.txt)

###################################
# 3 - LOAD MODULES
###################################
module load Trimmomatic/0.39-Java-11

###################################
# 4 - START MAIN CODE
###################################
# make output directory
mkdir -p $OUT

# Trim adapters from raw reads
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE \
$IN/"$INPUT_FILE"_R1_001.fastq.gz \
$IN/"$INPUT_FILE"_R2_001.fastq.gz \
$OUT/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz \
$OUT/"$INPUT_FILE"_R1_001_trimmed_U.fastq.gz \
$OUT/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz \
$OUT/"$INPUT_FILE"_R2_001_trimmed_U.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/"$ADAPTER":2:30:10 \
MINLEN:20


echo "Reads trimmed: adapter removed."


# Job finished
