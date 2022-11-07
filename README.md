# Ovary-IVF-data
Workflow for aligning raw ovary reads. Fastq files were taken and a bam file was produced.

### 1. Merge files
This was not necessary as the files were already merged.
### 2. Perform QC
QC was performed on the data using FastQC, then MultiQC was run to visualise the results together. All files failed on adapter content and also on per base sequence content.
### 3. Trim files
The adapters were removed using trimmomatic. No other trimming was necessary (https://sequencing.qcfail.com/articles/positional-sequence-bias-in-random-primed-libraries/).
### 4. Perform QC on trimmed files
Step 2 was repeated on the trimmed files to ensure the quality of these is sufficient.
### 5. Use STAR to create an index file
An annotation file was downloaded from ensemble http://www.ensembl.org/info/data/ftp/index.html/
```
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
```
