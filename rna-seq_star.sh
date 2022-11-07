#!/bin/bash

# ----------------------------------------------------------
# Based on Saskia Reibe's code
# Script to run quality check on in-house IVF ovary data
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
module load STAR/2.7.3a-GCC-9.3.0

REF=//well/lindgren/users/mzf347/alignment/ivf_cumulus/ref_genome
IN==//well/lindgren/users/mzf347/alignment/ivf_cumulus/qc_trimmed_results
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star


INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' $fastq/index.txt)


STAR --runThreadN 6 \
--readFilesIn $IN/"$INPUT_FILE"_1_val_1.fq.gz $IN/"$INPUT_FILE"_2_val_2.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix $OUT/$INPUT_FILE \
--genomeDir $REF/star \
--sjdbGTFfile $REF/Homo_sapiens.GRCh38.99.gtf --outSJfilterReads Unique --sjdbOverhang 100 \
--outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 \
--chimSegmentMin 20 --outSAMattributes All \
--quantMode TranscriptomeSAM 

# End of job script
