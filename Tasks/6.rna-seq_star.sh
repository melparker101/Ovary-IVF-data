#!/bin/bash

# ----------------------------------------------------------
# Based on Saskia Reibe's code
# Script to map in-house IVF ovary RNA-seq data to genome
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
###$ -o logs/output.log
###$ -e logs/error.log

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
# module purge 
# module use -a /apps/eb/dev/{skylake,ivybridge}/modules/all
module load STAR/2.7.9a-GCC-11.2.0

fastq=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads
STAR_INDEX=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star_index
REF=//well/lindgren/users/mzf347/ref_genomes/homo_sapiens/gencode/GRCh38.p13
IN=//well/lindgren/users/mzf347/alignment/ivf_cumulus/trimmed_reads
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star


INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' /index.txt)


STAR --runThreadN 6 \
--readFilesIn $IN/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz $IN/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix $OUT/$INPUT_FILE \
--genomeDir $STAR_INDEX \
--sjdbGTFfile $REF/gencode.v42.primary_assembly.annotation.gtf --outSJfilterReads Unique --sjdbOverhang 149 \
--outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --chimSegmentMin 20 \
--outSAMtype BAM SortedByCoordinate --outSAMattributes All \
--quantMode TranscriptomeSAM  # output alignments translated into transcript coordinates

# End of job script
