# Ovary-IVF-data

A workflow for analysing ovary bulk RNA-seq data. The samples were from patients who had undergone IVF, and were of cumulus cells or folicular fluid. This pipline uses as input the raw read data that were aquired, after external sequencing. The reads were aligned to a reference genome, genes were quantified, then downstream analysis was performed. The steps of the pipeline are as follows:

#### 1. Merge files
This was not necessary as the files were already merged.
#### 2. Perform QC on raw reads
QC was performed on the data using FastQC, then MultiQC was run to aggregate and visualise the results. All files failed on adapter content and also on per base sequence content.
- https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help
#### 3. Trim files
The adapters were removed using trimmomatic. No other trimming was necessary (https://sequencing.qcfail.com/articles/positional-sequence-bias-in-random-primed-libraries/).
```
IN=<raw_reads_dir>
OUT=<trimmed_reads_dir>

qsub trim_reads.sh $IN $OUT
```
#### 4. Perform QC on trimmed files
Step 2 was repeated on the trimmed files to ensure the quality of these is sufficient.
#### 5. Create an index file for STAR aligner
STAR was used to create an index file. The STAR manual recommends using the genome PRI assembly sequencing files. <br />
https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf <br />
The current version (Nov-22) is Release 42 (GRCh38.p13). 
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz
gzip -d GRCh38.primary_assembly.genome.fa.gz gencode.v42.primary_assembly.annotation.gtf.gz
```
An index file was created.
#### 6. Map reads to the genome
Star was used for read mapping. Alternatives include HISAT2, Bowtie, and TopHat2. 
#### 7. Create RSEM reference
#### 8. Perform gene and isoform quantification
Gene and isoform expression of rna-seq was estimated using RSEM. Alternatives include eXpress and Salmon.
#### 9. Create counts matrix
RSEM was used to create a count matrix.
#### 10. Perform downstream analysis
DESeq2 was used to perform downstream analysis. The analysis included finding the most highly expressed genes, generating PCA plots and performing differential gene expression. <br />
https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
