# Process Ovary RNA-Seq data
The main script is [IVF_data_processing.sh]([https://github.com/melparker101/Ovary-IVF-data/tree/main](https://github.com/melparker101/Ovary-IVF-data/blob/main/IVF_data_processing.sh)). This runs the whole pipeline altogether. However, it is recommended to run the run-scripts separately to assess the QC reports along the way and adjust accordingly.

### 1. Clone the respository
```
git clone https://github.com/melparker101/Ovary-IVF-data

cd Ovary-IVF-data
```

### 2. Download reference genome files
```
# Make a directory for the reference genome
mkdir ref

# Alternatively, download gencode files (this is recommended in the STAR manual)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz -P ref
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz -P ref

# Decompress files
# Gunzip invokes gzip -d
for f in ref/*.gz ; do gunzip -c "$f" > "${f%.*}" ; done
```

```
# Alternatively, download ensembl files
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
```

GMT files can be downloaded from:
https://wikipathways-data.wmcloud.org/current/gmt/

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
