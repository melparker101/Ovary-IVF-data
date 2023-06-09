#!/bin/bash

# ----------------------------------------------------------
# Script to generate a genome index for in-house IVF ovary RNA-seq data using STAR
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J star_index
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


# Load modules
module load STAR/2.7.9a-GCC-11.2.0

REF_GENOME=$1  # ref_genomes/homo_sapiens/gencode/GRCh38.p13
OUT=$2  # star_index

# Make an output directory
mkdir -p $OUT

# Generate STAR index
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir $OUT \
--genomeFastaFiles $REF_GENOME/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile $REF_GENOME/gencode.v42.primary_assembly.annotation.gtf \
--sjdbOverhang 149  # read-length - 1


echo "###########################################################"
echo "Star index generated."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
