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
## - Perform exploratory analysis and visualization
## - Perform differential expression analysis

##############################
# 1 - Load librairies
##############################
# library(tidyverse)
# library("AnnotationDbi")
# library("org.Hs.eg.db")
# library(geneplotter)
# library(tidyverse)
# library(ggbiplot)
# library(gplots)

library(DESeq2)  # For the DESeqDataSet
library("dplyr")
library("ggplot2")   # For variance distribution plots
library("glmpca")  # For GLM-PCA
library("pheatmap")  # For the heatmap
library("RColorBrewer")  # For the heatmap colours
library("PoiClaClu")

############################## 
# 2 - Set working directory
##############################
setwd("ivf_cumulus")
dir()

############################## 
# 3 - Source file
##############################
# input_file <- "full.count.genes.matrix.txt"
in_counts_file <- "count_matrix2.txt"  # Counts should NOT be normalised
in_phenotype_file <- "IVF_phenotypic_data_repox_15nov2022_trimmed.csv"
out_norm_counts <- "normalized_counts.txt"

############################## 
# 4 - Start Code
##############################
### CREATE DESEQDATASET OBJECT ###

# Read in data
count_data <- read.delim(in_counts_file, header = T, sep="\t", row.names = 1,stringsAsFactors=FALSE)
pheno_data <- read.csv(in_phenotype_file)

samples <- c()
condition <- c()

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

### PERFORM EXPLORATORY ANALYSIS AND VISUALISATION ###

# Variance stabilizing transformation
vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

# Regularized-logarithm transformation
rld <- rlog(ddsMat, blind = FALSE)
head(assay(rld), 3)
colData(rld)      

# Compare variance distribution for log2+1, vst and rlog transformations
df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
  
# From the graphs, it appears that rlog gives the best variance distribution
# Use rlog for PCA later

# Sample distance
sampleDists <- dist(t(assay(rld)))
sampleDists

# Visualize the distances in a heatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, vsd$samples, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
         
# Heatmap of sample-to-sample distances using the variance stabilising transformed values
poisd <- PoissonDistance(t(counts(ddsMat)))

# Poisson Distance
# This measure of dissimilarity between counts also takes the inherent variance structure of counts 
# into consideration when calculating the distances between samples
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$condition, dds$samples, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# Basic PCA
plotPCA(rld)

# Annotated PCA
plotPCA(rld, intgroup = c("condition", "samples"))

# Make a PCA table from the rld data
pcaData <- plotPCA(rld, intgroup = c("condition", "samples"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

# Make a PCA plot using shapes and colours
ggplot(pcaData, aes(x = PC1, y = PC2, color = samples, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rlog data")

# RNA counts, generally, are never normally distributed    
ggplot(count_data) +
  geom_histogram(aes(x = IVF0003), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
                                 
# GLM-PCA
gpca <- glmpca(counts(ddsMat), L=2)
gpca.dat <- gpca$factors
gpca.dat$condition <- dds$condition
gpca.dat$samples <- dds$samples

# Plot GLM-PCA
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = samples, shape = condition)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# Plot just the samples
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = samples)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
     
# MSD plot
# Use a matrix of distances
sampleDists <- dist(t(assay(rld)))
sampleDists

mds <- as.data.frame(colData(rld))  %>%
         cbind(cmdscale(sampleDistMatrix))
         
ggplot(mds, aes(x = `1`, y = `2`, color = samples, shape = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with rld data")
  
### DE ANALYSIS ###

res <- results(ddsMat)
res
mcols(res, use.names = TRUE)
summary(res)

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
                            
# Plot counts per condition for top gene
# Save to pdf
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("condition"))     

pdf(file="draft_counts_plot.pdf", width=30)
{par(lwd = 2)
plotCounts(dds, gene = topGene, intgroup=c("condition"))
}
dev.off()

# Volcano plot
# https://lashlock.github.io/compbio/R_presentation.html
# Reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add coloured points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


                                 



        








