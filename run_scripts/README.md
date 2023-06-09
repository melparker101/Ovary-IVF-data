### These scripts are being altered to run with slurm, but have not yet been tested.

### Pipeline overview
1. Create index file for raw fastq file names: [make_index.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/make_index.sh)
2. Merge fastq files: [merge.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/merge.sh)
3. Perform quality control (QC) on raw reads: [qc_reads.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/qc_reads.sh)
4. Aggregate QC results for raw reads: [multiqc.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/multiqc.sh)
5. Trim adapters: [trim_adapters.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/trim_adapters.sh)
6. Create an index for the trimmed reads: [make_index.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/make_index.sh)
7. Perform QC on trimmed reads: [qc_reads.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/qc_reads.sh)
8. Aggregate QC results for trimmed reads: [multiqc.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/multiqc.sh)
9. Generate STAR genome index: [generate_genome_index.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/generate_genome_index.sh)
11. Align genes to reference genome: [alignment_star.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/alignment_star.sh)
12. Perform alignment QC: [qc_alignment.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/qc_alignment.sh)
13. Generate RSEM transcript reference: [create_rsem_ref.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/create_rsem_ref.sh)
14. Quantify genes and isoforms: [quantification_rsem.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/quantification_rsem.sh)
15. Generate count matrix: [generate_count_matrix.sh](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_script/generate_count_matrix.sh)
16. Run downstream analysis: [downstream_analysis_DESeq2.R](https://github.com/melparker101/Ovary-IVF-data/blob/main/run_scripts/downstream_analysis_DESeq2.R)
