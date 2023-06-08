#!/bin/bash

# ----------------------------------------------------------
# Script to map in-house IVF ovary RNA-seq data to a reference genome
# melodyjparker14@gmail.com - Nov 22
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J alignment_star
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


# module purge 
# module use -a /apps/eb/dev/{skylake,ivybridge}/modules/all
module load STAR/2.7.9a-GCC-11.2.0

IN=$1  # trimmed_reads
OUT=$2  # star
STAR_INDEX=$3  # star_index
REF_GENOME=$4  # ref_genomes/homo_sapiens/gencode/GRCh38.p13
ref_genome_gtf=$5  # gencode.v42.primary_assembly.annotation.gtf

INPUT_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' /index.txt)

mkdir -p $OUT

STAR --runThreadN 6 \
--readFilesIn "$IN"/"$INPUT_FILE"_R1_001_trimmed_P.fastq.gz "$IN"/"$INPUT_FILE"_R2_001_trimmed_P.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix "$OUT"/"$INPUT_FILE"_ \
--genomeDir $STAR_INDEX \
--sjdbGTFfile "$REF_GENOME"/"$ref_genome_gtf" --outSJfilterReads Unique --sjdbOverhang 149 \
--outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --chimSegmentMin 20 \
--outSAMtype BAM SortedByCoordinate --outSAMattributes All \
--quantMode TranscriptomeSAM  # output alignments translated into transcript coordinates


echo "###########################################################"
echo "Read mapping finished."
echo "Finished at: "`date`
echo "###########################################################"
exit 0
