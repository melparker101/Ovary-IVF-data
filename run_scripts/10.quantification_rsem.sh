#!/bin/bash

# ----------------------------------------------------------
# Script to perform gene and isoform quantification on
# in-house IVF ovary RNA-seq data
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------


#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J quantification_rsem
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


IN=$1  # trimmed_reads
OUT=$2  # qc_trimmed_results

# Begin writing your script here
module purge
module load RSEM/1.3.2-foss-2018b

fastq=raw_reads
RSEM_REF=rsem_ref
IN=star
OUT=rsem

if [ ! -d "$OUT" ]; then
  mkdir -p $OUT
fi

# Merge index is a file with list of 8 digit numbers (part of the file name)
# This means input file is the task_id'th line of that list
INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' $IN/index.txt)


rsem-calculate-expression --bam --paired-end -p 4 $IN/"$INPUT_FILE"_Aligned.toTranscriptome.out.bam $RSEM_REF/human $OUT/rsem_"$INPUT_FILE"


echo "###########################################################"
echo "Genes and isoforms quantified."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
