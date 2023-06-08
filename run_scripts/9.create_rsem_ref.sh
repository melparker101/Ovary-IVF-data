#!/bin/bash

# ----------------------------------------------------------
# Script to create RSEM reference
# melodyjparker14@gmail.com - Nov 22
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

# Load modules
module load RSEM/1.3.2-foss-2018b

# Define variables
OUT=$6  # rsem_ref
ref_genome_gtf=$  # gencode.v42.primary_assembly.annotation.gtf
ref_genome_fasta  # GRCh38.primary_assembly.genome.fa

# Create output directory for RSEM reference
mkdir -p $OUT

# Create RSEM reference
rsem-prepare-reference --gtf $REF_GENOME/$ref_genome_gtf \
   				     $REF_GENOME/$ref_genome_fasta $OUT/human                    


echo "###########################################################"
echo "RSEM reference generated."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
