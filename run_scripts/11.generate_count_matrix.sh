#!/bin/bash

# ----------------------------------------------------------
# Script to merge RSEM results into a count matrix
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J generate_count_matrix
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

# Generate count matrix
rsem-generate-data-matrix *.genes.results >> full.count.genes.matrix.txt

# Trim the version number off the ensembl IDs
sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' full.count.genes.matrix.txt > count_matrix.txt


echo "###########################################################"
echo "Count matrix generated."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
