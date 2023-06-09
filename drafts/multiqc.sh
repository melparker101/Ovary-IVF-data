#!/bin/bash

# ----------------------------------------------------------
# Script to run multiqc directly after fastqc has finished running
# The default creates a multiqc report in the pwd.
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N multiqc

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o logs/
#$ -e logs/

# Parallel environemnt settings
#  For more information on these please see the wiki
#  Allowed settings:
#   shmem
#   mpi
#   node_mpi
#   ramdisk
#$ -pe shmem 1


IN=$1  # Input dir
OUT=IN
MODULES=$2
# FQCJOBID=$2  # FastQC job ID
# more efficient way so that we're not wasting nodes


module load MultiQC/1.9-foss-2019b-Python-3.7.4

if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi

# while
# qstat | grep "$FQCJOBID" > /dev/null; do sleep 60; echo "sleeping" ; done

echo "FastQC complete."


multiqc $IN -o $OUT -n multiqc_report_"$MODULES".html 

# multiqc .
# python3 -m multiqc ./  # Use a specific python interpreter


echo "Multiqc complete."

duration=$(( SECONDS - start ))

# End of job script
