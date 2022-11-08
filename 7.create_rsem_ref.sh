#!/bin/bash

# ----------------------------------------------------------
# Script to create RSEM reference
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

# Specify a job name
#$ -N create_rsem_ref

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
# module load STAR/2.7.9a-GCC-11.2.0
module load RSEM/1.3.2-foss-2018b

fastq=//well/lindgren/users/mzf347/alignment/ivf_cumulus/raw_reads
STAR_INDEX=//well/lindgren/users/mzf347/alignment/ivf_cumulus/star_index
STAR_PATH=/gpfs3/apps/eb/2020b/skylake/software/STAR/2.7.9a-GCC-11.2.0/bin/STAR
REF_GENOME=//well/lindgren/users/mzf347/ref_genomes/homo_sapiens/gencode/GRCh38.p13
IN=//well/lindgren/users/mzf347/alignment/ivf_cumulus/trimmed_reads
OUT=//well/lindgren/users/mzf347/alignment/ivf_cumulus/rsem_ref

if [ ! -d $OUT]; then
  mkdir rsem_ref
fi


rsem-prepare-reference --gtf $REF_GENOME/Mus_musculus.GRCm38.82.chr.gtf \
   				     $REF_GENOME/GRCh38.primary_assembly.genome.fa \
               $OUT                    


# End of job script 
 
 
 
