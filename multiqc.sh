#!/bin/bash

# ----------------------------------------------------------
# Script to run MultiQC on in-house IVF ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

SECONDS=0

# Specify a job name
#$ -N multiqc

# Project name and target queue choose short or long
#$ -P lindgren.prjc
#$ -q short.qe

# Run the job in the current working directory
#$ -wd //well/lindgren/users/mzf347/alignment/ivf_cumulus/qc_trimmed_results

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
#$ -pe shmem 4

module load MultiQC/1.9-foss-2019b-Python-3.7.4
multiqc .
