# Process Ovary RNA-Seq data

### 1. Clone the respository
```
git clone https://github.com/melparker101/Ovary-IVF-data

cd Ovary-IVF-data
```

### 2. Download reference genome files
```
mkdir <ref_genome_dir>
cd <ref_genome_dir>

# Download files
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz

# Decompress files
for f in *.gz ; do gunzip -c "$f" > "${f%.*}" ; done

cd ..
```

### 3. Process the data
Send a job to the cluster.

```
# Input directories
raw_reads=<raw_reads_dir>
ref=<ref_genome_dir>

qsub processing_ovary_rna-seq_data.sh $raw_reads $ref
```

# Trim data only
```
IN=<raw_reads_dir>
OUT=<trimmed_reads_dir>

qsub trim_reads.sh $IN $OUT
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

qsub mapping_and_quantification.sh -raw $RAW_READS -trimmed $TRIMMED_READS -ref $REF_GENOME \
                                 -star_index $STAR_INDEX -star $STAR_OUT \
                                 -rsem_ref $RSEM_REF -rsem $RSEM_OUT
```
