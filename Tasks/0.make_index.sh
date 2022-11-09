#!/bin/bash

# ----------------------------------------------------------
# Script to make an index file for in-house IVF RNA-seq ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

# Specify a job name
#$ -N make_index

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
#$ -pe shmem 1

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
fastq=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads


for f in $fastq/IVF*R1*; do g="${f%_R1*}" ; echo ${g##*/} >> $fastq/index.txt ; done


# End of job script
