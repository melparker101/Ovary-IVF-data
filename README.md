# Process Ovary RNA-Seq data

### 1. Clone the respository
```
git clone https://github.com/melparker101/Ovary-IVF-data
```

### 2. Process the data
Send a job to the cluster.

```
raw_reads=<raw_reads_dir>
ref=<ref_genome_dir>

qsub processing_ovary_rna-seq_data.sh $raw_reads $ref
```

# Trim data only
```

```

# Map and quantify only
```
# Input directory names/paths
TRIMMED_READS=<trimmed_reads>
STAR_INDEX=<star_index>
RSEM_REF=<rsem_ref>

# Output directory names/paths:
STAR_OUT=<star_output>
RSEM_OUT=<rsem_output>

sh mapping_and_quantification.sh -raw $RAW_READS -trimmed $TRIMMED_READS -ref $REF_GENOME \
                                 -star_index $STAR_INDEX -star $STAR_OUT \
                                 -rsem_ref $RSEM_REF -rsem $RSEM_OUT
```
