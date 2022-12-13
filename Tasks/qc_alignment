#############################################################
## Alignment QC
## melodyjparker14@gmail.com - Dec 22
## This code performs quality control checks on aligned RNA-seq data (BAM files) from in-house IVF ovary samples.
## Tutorial followed:
## https://rnabio.org/module-02-alignment/0002/06/01/Alignment_QC/
#############################################################

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
module load samtools/1.8-gcc5.4.0  # for indexing bam files
module load BEDOPS/2.4.35-foss-2018b
module load picard/2.23.0-Java-11
module load MultiQC/1.7-foss-2018b-Python-3.6.6

echo "modules loaded"
###################################
# SETUP
###################################
# Set wd
# mkdir sortedByCoord
# cp *sortedByCoord.out.bam sortedByCoord/
cd sortedByCoord

# source $HOME/.bashrc

# Copy gencode reference genome folder across
# $ref_gc contains the latest gencode reference fa file and annotations gtf file
cp -R $ref_gc ref

# Index BAM files
# Consider using GNU parallel
# module load parallel
# parallel  samtools index ::: *.bam
# The 'parallel' module on the cluster clashes with the GCC version of other modules
# Use loop
for f in *.bam; do samtools index "$f" "$f".bai ; done
echo "BAM file indexing complete."

###################################
# FUNCTIONS
###################################
# Capture the output of a command
# https://stackoverflow.com/questions/24283097/reusing-output-from-last-command-in-bash
cap () { tee /tmp/capture.out; }  # Capture
ret () { cat /tmp/capture.out; }  # Retrieve



###################################
# 1 - FASTQC
###################################
# Run FastQC in parallel - send off a script
qsub fastqc.sh pwd fastqc bam | cap

# Retrieve job ID
fastqc_job_id=$(ret | awk -v RS='[0-9]+' '$0=RT' | head -1)  

# Run Multiqc once FastQC has finished
qsub multiqc.sh fastqc fastqc_job_id

###################################
# 2 - PICARD
###################################
cd ref

# Create a dictionary from the reference fasta file
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=GRCh38.primary_assembly.genome.fa O=GRCh38.primary_assembly.genome.dict

# Create a bed file of the location of ribosomal sequences in out reference
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
mkdir picard
# find *sortedByCoord.out.bam -exec echo java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics I={} O=picard/{}.RNA_Metrics REF_FLAT=ref/gencode.v42.primary_assembly.ref_flat.txt STRAND=NONE RIBOSOMAL_INTERVALS=ref/ref_ribosome.interval_list \; | sh
# This step takes a while---run this in parallel
JOBIDA=$(qsub -terse picard_collect_rna_seq_metrics.sh)

# Mark duplicates
# for f in *bam; do java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT="$f" OUTPUT=dedup/dedup_"$f" M=dedup/"${f%.sorted*}"_metrics.txt; done
# This is also slow, so also run in parallel
JOBIDB=$(qsub -terse picard_mark_duplicates.sh)

echo "Picard checks complete."

###################################
# 2 - SAMTOOLS
###################################
# Perform flagstat on the bam files which have duplicates marked
mkdir flagstat
for f in dedup/*.bam; do samtools flagstat "$f" > flagstat/"$f".flagstat; done
echo "Flagstats created."

cd dedup
ls *bam > index.txt
cd ..

# Make an index file
(cd dedup && ls *bam) > dedup/index.txt

# Loop through the above index file
while IFS= read -r file
do
	samtools flagstat dedup/"$file" > flagstat/"${file%.out*}".flagstat
done < "dedup/index.txt"


while IFS= read -r file; do	samtools flagstat dedup/"$file" > flagstat/"${file%.out*}".flagstat; done < "dedup/index.txt"

# Alternatively:
# find *sortedByCoord*.bam -exec echo samtools flagstat {} \> flagstat/{}.flagstat \; | sh

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

cd ref


# Convert Gtf to genePred
# gtfToGenePred gencode.v42.primary_assembly.annotation.gtf gencode.v42.primary_assembly.genePred

# Convert genPred to bed12
# genePredToBed gencode.v42.primary_assembly.genePred gencode.v42.primary_assembly.bed12

cd ..
mkdir rseqc

# This step takes a while
# Send a separate script off
# bam_files=$(ls *bam | paste -sd "," -)
# geneBody_coverage.py -i $bam_files -r ref/gencode.v42.primary_assembly.bed12 -o rseqc

echo "Gene converage script sent off"

# Send off script for the rest of the RSeQC analysis



###################################
# 5 - MULTIQC
###################################
# Wait until scripts are finished
qsub -q short.qc -hold_jid $JOBA_ID,$JOBIDB multiqc.sh $PWD  # check if you can wait on more than one



"MultiQC script sent off."

conda deactivate
