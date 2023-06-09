#!/bin/bash

# ----------------------------------------------------------
# An array script to run Picard's function MarkDuplicates
# melodyjparker14@gmail.com - Dec 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J MarkDuplicates
#SBATCH -o logs/
#SBATCH -e logs/
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


# Define variables
IN=$1
OUT=bam_dedup
INPUT_FILE=$(sed "$SGE_TASK_ID"'q;d' index.txt)

# Load modules
module load picard/2.23.0-Java-11

# Make output directory
mkdir -p $OUT

# Mark Duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT="$INPUT_FILE"*.bam OUTPUT="$OUT"/dedup_"$INPUT_FILE"_Aligned.sortedByCoord.out.bam M="$OUT"/"$INPUT_FILE"_metrics.txt; done 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
