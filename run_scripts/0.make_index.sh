#!/bin/bash

# ----------------------------------------------------------
# Script to make an index file for in-house IVF RNA-seq ovary data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J make_index
#SBATCH -o output.out
#SBATCH -e error.err

#  Parallel environment settings 
#  For more information on these please see the documentation 
#  Allowed parameters: 
#   -c, --cpus-per-task 
#   -N, --nodes 
#   -n, --ntasks 

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"


fastq=$1  # fastq files/reads 

<<comment
for f in $fastq/IVF*R1*; do g="${f%_R1*}" ; echo ${g##*/} >> $fastq/index.txt ; done
comment

if [ ! -f index.txt ]; then
for f in IVF*R1*.fastq.gz; do echo ${f%_R1_001*} >> index.txt ; done
fi


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
