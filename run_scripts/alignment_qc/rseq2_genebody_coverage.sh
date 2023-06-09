#!/bin/bash

# ----------------------------------------------------------
# Script to run gene body coverage
# melodyjparker14@gmail.com - June 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J create_rsem_ref
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err

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

IN=$1
OUT=$2
REF_GENOME=ref

geneBody_coverage.py -i $IN -r "$REF_GENOME"/gencode.v42.primary_assembly.bed12 -o "$OUT"/rseqc

echo "###########################################################"
echo "RSEM reference generated."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
