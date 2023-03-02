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

# Begin writing your script here
module load RSEM/1.3.2-foss-2018b

fastq=$1  # raw_reads
STAR_INDEX=$2  # star_index
STAR_PATH=$3  # software/STAR/2.7.9a-GCC-11.2.0/bin/STAR
REF_GENOME=$4  # ref_genomes/homo_sapiens/gencode/GRCh38.p13
IN=$5  # ivf_cumulus/trimmed_reads
OUT=$6  # rsem_ref

if [ ! -d $OUT ]; then
  mkdir -p $OUT
fi


rsem-prepare-reference --gtf $REF_GENOME/gencode.v42.primary_assembly.annotation.gtf \
   				     $REF_GENOME/GRCh38.primary_assembly.genome.fa $OUT/human                    


echo "###########################################################"
echo "RSEM reference generated."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
