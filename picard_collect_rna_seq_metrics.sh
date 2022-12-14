#!/bin/bash

# ----------------------------------------------------------
# A script to run Picard's CollectRnaSeqMetrics on BAM files for alignment QC in parallel.
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------

# Specify a job name
#$ -N CollectRnaSeqMetrics

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


if [ ! -f file_index.txt ]; then
  for f in *.bam; do echo $f >> file_index.txt; done
fi


INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' file_index.txt)
OUT=picard
REF=ref
GENE_PRED=gencode.v42.primary_assembly.ref_flat.txt  # genePred
RIB_INT=ref_ribosome.interval_list  # Interval list file for the ribosomal sequence location in the reference


if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi


module load picard/2.23.0-Java-11


java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics I=$INPUT_FILE O=$OUT/"$INPUT_FILE".RNA_Metrics REF_FLAT=$REF/$GENE_PRED STRAND=NONE RIBOSOMAL_INTERVALS=$REF/$RIB_INT


echo "CollectRnaSeqMetrics complete."


# End of job script

