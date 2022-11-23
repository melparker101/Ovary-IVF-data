#############################################################
## Downstream Analysis
## melodyjparker14@gmail.com - Nov 22
## Based off Saskia Reibe's code and following the DESeq2 Tutorial
## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#############################################################

##############################
# 0 - Workflow
##############################
## - Create a DESeqDataSet object
## - 
##

##############################
# 1 - Load librairies
##############################
# library(tidyverse)
# library(DESeq2)
# library("AnnotationDbi")
# library("org.Hs.eg.db")
# library(geneplotter)
# library(tidyverse)
# library("pheatmap")
# library("RColorBrewer")
# library(ggbiplot)
# library(gplots)

############################## 
# 2 - Set working directory
##############################
setwd("ivf_cumulus")
dir()

############################## 
# 3 - Source file
##############################
# input_file <- "full.count.genes.matrix.txt"
in_counts_file <- "count_matrix2.txt"
in_phenotype_file <- "IVF_phenotypic_data_repox_15nov2022_trimmed.csv"
out_norm_counts <- "normalized_counts.txt"

############################## 
# 4 - Start Code
##############################
# Read in data
count_data <- read.delim(in_counts_file, header = T, sep="\t", row.names = 1,stringsAsFactors=FALSE)
pheno_data <- read.csv(in_phenotype_file)

samples <- []
condition <- []

# Trim column names to match the IDs in the phenotype data
samples <- colnames(count_data)
samples <- gsub("\\..*","",samples)  # Removes the end after .
samples <- gsub(".*_","",samples)  # Removed the start up to _
colnames(count_data) <- samples

pheno_data <- subset(pheno_data, ID %in% samples)
dim(pheno_data)

# Check that they match
all(colnames(count_data) %in% pheno_data$ID)
all(colnames(count_data) == pheno_data$ID)
# Use match() to rearrange if needed

# View just CaseControl and ReasonforIVF columns
pheno_data[,c("ID","CaseControl","ReasonForIVF")]

# Group male factor together
pheno_data$ReasonForIVF[grep("Male|Sperm",pheno_data$ReasonForIVF)]<-"Male Factor"

condition <- pheno_data$ReasonForIVF
coldata <- as.data.frame(cbind(samples,condition))

# Create a matrix
count_mat <- as.matrix(count_data)
storage.mode(count_mat) = "integer"

# Quick QC
# This isn't actually necessary for our data
qc <- complete.cases(count_mat)
summary(qc)
countdata <- na.omit(count_mat)

# Create a DESeqDataSet object
ddsMat <- DESeqDataSetFromMatrix(countData = count_mat,
                                 colData = coldata,
                                 design = ~condition)
                                 
head(assay(ddsMat))  

# View in the script editor
View(counts(dds))

# Filter out empty rows
# Strict filtering to increase power is automatically applied later
keep <- rowSums(counts(ddsMat)) > 0
ddsMat <- ddsMat[keep,]

# We can also filter like this
# keep <- rowSums(counts(ddsMat)) >= 10
# ddsMat <- ddsMat[keep,]

# Add control features for estimating size factors
ddsMat <- estimateSizeFactors(ddsMat, controlGenes=ctrlGenes)
ddsMat <- DESeq(ddsMat)

# View normalisation factor applied and compare non-normalised data with normalised data
sizeFactors(dds)
head(counts(ddsMat, normalized = F))
head(counts(ddsMat, normalized = T))

# Write normalised table to file
normalized_counts <- counts(ddsMat, normalized=TRUE)
write.table(normalized_counts, file=out_norm_counts, sep="\t", quote=F, col.names=NA)

# This tells us if the rows contain all non-zero values
GeneCounts <- counts(ddsMat)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)

# Plots
# These didn't work
# multidensity( counts(ddsMat, normalized = T)[idx.nz ,],    
#							xlab="mean counts", xlim=c(0, 1000),main="Normalised data")
# multidensity( counts(ddsMat, normalized = F)[idx.nz ,],         
#							xlab="mean counts", xlim=c(0, 1000),main="Not normalised data")

# Variance distribution
vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

rld <- rlog(ddsMat, blind = FALSE)
head(assay(rld), 3)
colData(rld)                                
                                 
                                 
                                 
                                 
                                 
                                 
                                 

# Add symbol, entrez, and uniprot columns to the dataframe
count_data$symbol <- mapIds(org.Hs.eg.db,
                            keys=row.names(count_data),                                                      
                            column="SYMBOL",
                            keytype="ENSEMBL",          
                            multiVals="first")

count_data$entrez <- mapIds(org.Hs.eg.db,
                            keys=row.names(count_data),
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")

count_data$uniprot <- mapIds(org.Hs.eg.db,
                             keys=row.names(count_data),
                             column="UNIPROT",
                             keytype="ENSEMBL",
                             multiVals="first")            








