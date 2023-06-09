#!/bin/bash

#############################################################
## Alignment QC
## melodyjparker14@gmail.com - Dec 22
## This code performs quality control checks on aligned RNA-seq data (BAM files) from in-house IVF ovary samples.
## Tutorial followed:
## https://rnabio.org/module-02-alignment/0002/06/01/Alignment_QC/
#############################################################

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J alignment_qc
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

###################################
# WORKFLOW
###################################
# 1. Samtools
# 2. FastQC
# 3. Picard
# 4. RSeQC
# 5. MultiQC

###################################
# LOAD MODULES
###################################
source $HOME/.bashrc

# conda create -n ucsc
conda activate ucsc
# conda install -c bioconda ucsc-gtftogenepred
# conda install -c bioconda ucsc-genepredtobed

echo "conda ucasc env activated."

module load RSeQC/3.0.0-foss-2018b-Python-3.6.6
# module load samtools/1.8-gcc5.4.0 
module load BEDOPS/2.4.35-foss-2018b
module load picard/2.23.0-Java-11
module load MultiQC/1.7-foss-2018b-Python-3.6.6

echo "modules loaded"

###################################
# Variables
###################################
STAR=$1  # star
REF_GENOME=$2  # ref_gc
BAM=bam_files
ALIGNMENT_QC=run_scripts/alignment_qc
QC=qc_results

###################################
# SETUP
###################################
# Make a directory for BAM files
mkdir -p $BAM

# Copy across the sorted bam files from the STAR output folder
cp $STAR/*sortedByCoord.out.bam $BAM/

# Copy gencode reference genome folder across
# $ref_gc contains the latest gencode reference fa file and annotations gtf file
cp -R $REF_GENOME ref

# Create an index file containing a list of all the BAM files
JOBA_ID=$(sbatch --parsable -p short run_scripts/make_index.sh $BAM bam) 
echo "Creating an index file of sample names for BAM files. Job ID: $JOBA_ID."

# Index BAM files
# module load parallel
# parallel  samtools index ::: *.bam
# Loop:
# for f in *.bam; do samtools index "$f" "$f".bai ; done

# Send off script to do this in paralell
# Dependent on job A
JOBB_ID=$(sbatch --parsable -p short -d afterok:$JOBA_ID "$ALIGNMENT_QC"/index_bam_files.sh $BAM) 
echo "Indexing BAM files. Job ID: $JOBB_ID."

###################################
# 1 - FASTQC
###################################
# Run FastQC in parallel - send off a script
# Dependent on job A
JOBC_ID=$(sbatch --parsable -p short -d afterok:$JOBA_ID run_scripts/fastqc.sh $BAM $QC/qc_bam bam)
echo "Running FastQC on BAM files. Job ID: $JOBC_ID."

# Run Multiqc once FastQC has finished
JOBD_ID=$(sbatch --parsable -p short -d afterok:$JOBC_ID run_scripts/multiqc.sh $QC/$QC/qc_bam fastqc_bam)
echo "Running MultiQC to aggregate QC results. Job ID: $JOBD_ID."

###################################
# 2 - PICARD
###################################
cd ref

# Create a dictionary from the reference fasta file
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict

# Create a bed file of the location of ribosomal sequences in our reference
grep --color=none -i -P "rrna" gencode.v42.primary_assembly.annotation.gtf > ref_ribosome.gtf
gff2bed < ref_ribosome.gtf > ref_ribosome.bed

# Create an interval list file for the ribosomal sequence location in the reference
java -jar $EBROOTPICARD/picard.jar BedToIntervalList I=ref_ribosome.bed O=ref_ribosome.interval_list SD=GRCh38.primary_assembly.genome.dict

# Create a genePred file for the reference transcriptome
gtfToGenePred -genePredExt gencode.v42.primary_assembly.annotation.gtf gencode.v42.primary_assembly.ref_flat.txt

# Reformat this genePred file
cat gencode.v42.primary_assembly.ref_flat.txt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > tmp.txt
mv tmp.txt gencode.v42.primary_assembly.ref_flat.txt

cd ..

# find *sortedByCoord.out.bam -exec echo java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics I={} O=picard/{}.RNA_Metrics REF_FLAT=ref/gencode.v42.primary_assembly.ref_flat.txt STRAND=NONE RIBOSOMAL_INTERVALS=ref/ref_ribosome.interval_list \; | sh
# This step takes a while---run this in parallel
# Requires indexed BAMs
sbatch -p short-d afterok:$JOBB_ID picard_collect_rna_seq_metrics.sh bam_files "$QC"/picard ref
echo "Running Picard CollectRnaSeqMetrics."

# Mark duplicates
# This requires BAM index files, so depends on job B. Collect job ID
# $1 = input dir, $2 = output dir, $3 part of the report name
JOBE_ID=$(sbatch --parsable -p short -d afterok:$JOBB_ID "$QC"/picard_mark_duplicates.sh $BAM)
echo "Running Picard MarkDuplicates. Job ID: $JOBE_ID."

###################################
# 2 - SAMTOOLS
###################################
# Perform flagstat on the bam files which have duplicates marked

# Make output dir
mkdir "$QC"/flagstat

# Wait for MarkDuplicates to finish running
squeue --job $job_id
while squeue --job $JOBE_ID >/dev/null 2>&1; do
  sleep 1
done

# Run samtools flagstat
for f in bam_dedup/*.bam; do samtools flagstat "$f" > "$QC"/flagstat/$(basename "${f%.out*}").flagstat; done
echo "Flagstats created."

# Make an index for bam_dedup
JOBF_ID=$(sbatch --parsable run_scripts/make_index.sh bam_dedup)
echo "Making an index for deduplicated bam files"

###################################
# 4 - RSEQC
###################################
# conda create -n ucsc
# conda activate ucsc
# conda install -c bioconda ucsc-gtftogenepred
# conda install -c bioconda ucsc-genepredtobed
# module load RSeQC/3.0.0-foss-2018b-Python-3.6.6
# module load samtools/1.8-gcc5.4.0  # for indexing bam files
# module load BEDOPS/2.4.35-foss-2018b
# module load picard/2.23.0-Java-11

# Create the required reference files
if [[ -f ref/gencode.v42.primary_assembly.genePred ]]; then
	# Convert Gtf to genePred
	gtfToGenePred ref/gencode.v42.primary_assembly.annotation.gtf ref/gencode.v42.primary_assembly.genePred
fi

if [[ -f ref/gencode.v42.primary_assembly.bed12 ]]; then
	# Convert genPred to bed12
	genePredToBed ref/gencode.v42.primary_assembly.genePred ref/gencode.v42.primary_assembly.bed12
fi

# Create the output file
mkdir -p "$QC"/rseqc

# This step takes a while
# Send a separate script off
# bam_files=$(ls *bam | paste -sd "," -)

# geneBody_coverage.py -i $bam_files -r ref/gencode.v42.primary_assembly.bed12 -o rseqc

sbatch 
JOBF_ID=$(sbatch --parsable "$QC"/rseq2_genebody_coverage.sh $BAM "$QC")
echo "Running geneBody_coverage by rseq2."

echo "Gene converage script sent off"

# Send off script for the rest of the RSeQC analysis


###################################
# 5 - MULTIQC
###################################
# Wait until scripts are finished
JOBH_ID=$(sbatch --parsable -p short -d afterok:$JOBB_ID run_scripts/picard_mark_duplicates.sh $BAM)
echo "Running MultiQC. Job ID: $JOBH_ID."


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
