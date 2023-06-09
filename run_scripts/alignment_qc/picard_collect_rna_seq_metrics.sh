#!/bin/bash

# ----------------------------------------------------------
# A script to run Picard's CollectRnaSeqMetrics on BAM files for alignment QC in parallel.
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J CollectRnaSeqMetrics
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err
#SBATCH -a 1-15

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


# Make a picard output directory
mkdir -p 

IN=$1  # directory containin the bam files
OUT=$2
REF=$3  # ref
GENE_PRED=gencode.v42.primary_assembly.ref_flat.txt  # genePred
RIB_INT=ref_ribosome.interval_list  # interval list file for the ribosomal sequence location in the reference

INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' "$IN"/index.txt)

if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi

module load picard/2.23.0-Java-11

java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics I="$INPUT_FILE"*.bam O=$OUT/"$INPUT_FILE".RNA_Metrics REF_FLAT=$REF/$GENE_PRED STRAND=NONE RIBOSOMAL_INTERVALS=$REF/$RIB_INT


echo "###########################################################"
echo "CollectRnaSeqMetrics complete."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
