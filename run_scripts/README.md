### These scripts are being altered to run with slurm, but have not yet been tested.

### Pipeline overview
1. Create index file for raw fastq file names
2. Perform quality control (QC) on raw reads
3. Aggregate QC results for raw reads
4. Trim adapters
5. Create an index for the trimmed reads
6. Perform QC on trimmed reads
7. Aggregate QC results for trimmed reads
8. Generate STAR genome index
9. Align genes to reference genome
10. Generate RSEM transcript reference
11. Quantify genes and isoforms
12. Generate count matrix
13. Run downstream analysis
