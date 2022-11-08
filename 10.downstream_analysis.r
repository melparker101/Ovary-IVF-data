##example but try working alongside tutorial


library(tidyverse)

library(DESeq2)
library("AnnotationDbi")


library("org.Hs.eg.db")
library(geneplotter)
library(tidyverse)

library("pheatmap")

library("RColorBrewer")

library(ggbiplot)
library(gplots)

setwd("/Volumes/saskiar/BDI/RNA-seq/adipo/rna_seq")

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/work/pcos")

dir()


countdata0 <- read.delim("GSE193123_gene_count_pcos_healthy.txt", header = T, sep="\t", row.names = 1,stringsAsFactors=FALSE)

head(countdata0)

dim(countdata0)

countdata <- as.matrix(countdata0)



storage.mode(countdata) = "integer"





condition <- c(rep("healthy",3), rep("pcos",3))

samples <-colnames(countdata)





qc <- complete.cases(countdata)


summary(qc)

countdata <- na.omit(countdata)



countdata0$symbol <- mapIds(org.Hs.eg.db,
                            keys=row.names(countdata0),
                            
                            
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            
                            
                            
                            multiVals="first")




countdata0$entrez <- mapIds(org.Hs.eg.db,
                            
                            
                            keys=row.names(countdata0),
                            
                            column="ENTREZID",
                            
                            keytype="ENSEMBL",
                            
                            
                            multiVals="first")




countdata0$uniprot <- mapIds(org.Hs.eg.db,
                             
                             
                             keys=row.names(countdata0),
                             
                             column="UNIPROT",
                             
                             keytype="ENSEMBL",
                             
                             
                             multiVals="first")

head(countdata0)



coldata <- as.data.frame(cbind(samples,condition))



dim(countdata) #dimension of the counttable (identified genes and samples)

head(countdata,5) #first 5 rows of the counttable




ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 
                                 design = ~condition)


head(assay(ddsMat))

#normalised vs not normalised

ddsMat <- estimateSizeFactors(ddsMat)

sizeFactors(ddsMat)

GeneCounts <- counts(ddsMat)


idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})


sum(idx.nz)



multidensity( counts(ddsMat, normalized = F)[idx.nz ,],    
              
              
              xlab="mean counts", xlim=c(0, 1000),main="Not normalised data")


head(counts(ddsMat, normalized = F))

multidensity( counts(ddsMat, normalized = T)[idx.nz ,],    
              
              xlab="mean counts", xlim=c(0, 1000),main="Normalised data")



dds <- ddsMat


nrow(dds)    #number of rows in raw data


dds <- dds[ rowSums(counts(dds)) > 1, ]


nrow(dds)  #number of rows after filtering




rld <- vst(dds)


#PCA

plotPCA(rld)

head(assay(rld))



dds <- DESeq(dds)

#dds1 <-dds

#save(dds1,file="hWAT_dds.Rda")

par(mar=c(12,5,2,2))

boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)



#res <- results(dds, cooksCutoff=FALSE,contrast=c("condition","chow","hfd"))


res <- results(dds, contrast=c("condition","pcos","healthy"))





plotMA(res)



res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     
                     
                     column="ENTREZID",
                     
                     keytype="ENSEMBL",
                     
                     
                     multiVals="first")



res$uniprot <- mapIds(org.Hs.eg.db,
                      
                      
                      keys=row.names(res),
                      
                      
                      column="UNIPROT",
                      
                      
                      keytype="ENSEMBL",
                      
                      multiVals="first")





resOrdered <- res[order(res$padj),]


summary(res)


final <- subset(resOrdered, padj<0.05 & abs(log2FoldChange)>1)

dim(final)

head(final)

write.table(final,"pcos_healjty_fc1_p05.txt", sep = '\t', col.names = NA)


##
topGene <- rownames(res)[which.min(res$padj)]

topGene <-"ENSG00000135218"



#gene_name <- res[topGene,7]


(gene_name <- countdata0[topGene,7])


#head(countdata0)



data <- plotCounts(dds, gene=topGene, intgroup=c("condition"), returnData=TRUE)




#data$condition  <- factor(data$condition,levels = c("hWAT_D0","hWAT_D3","hWAT_D8","hWAT_D24","hBAT_D0","hBAT_D3","hBAT_D8","hBAT_D14","SGBS_D0","SGBS_D3","SGBS_D8","SGBS_D14"))

col_palette <- colorRampPalette(brewer.pal(9, "Set1")) 

plot <- ggplot(data, aes(x=condition, y=count, color = condition)) +
  #scale_y_log10() + 
  
  
  geom_point(position=position_jitter(width=.1,height=0), size=3)+ 
  
  ggtitle(gene_name)    +
  
 # ylim(0,30) +
  
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +

 
  
  theme(plot.title = element_text(lineheight=4, face="bold", size=12)) +
  
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=10),axis.text.x  = element_text(angle=90, vjust=0.5, size=10),axis.text.y  = element_text(angle=0, vjust=0.5, size=10)) +
  
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(legend.position="none") 





plot +  scale_colour_manual(values = col_palette(12)) +
  geom_vline(xintercept = 4.5, linetype="dotted", 
             color = "blue", size=1.5) 



ggsave("USP19_cell_line_expression.png",width = 20, height = 20, units = "cm")
#heatmaps




library(clusterProfiler)
library(org.Hs.eg.db)
library(ggridges)

xx <- as.data.frame(final)
gene <- xx[,8]

head(gene)

hbat <- gene
hWAT <- gene
sgbs <- gene

data(gcSample)
str(gcSample) 

cell_3 <- list(hbat,hWAT, sgbs)
str(cell_3)



cell_3 <- unique(cell_3)
names(cell_3)<-c("hBAT", "hWAT", "SGBS")
ck <- compareCluster(geneCluster = cell_3, fun = enrichKEGG,organism="hsa", pvalueCutoff=0.05)

ck <- compareCluster(geneCluster = cell_3, fun = enrichKEGG)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 

library(ReactomePA)
ck <- compareCluster(geneCluster = cell_3, fun = enrichPathway)

ck <- compareCluster(geneCluster = cell_3, fun = enrichWP, organism = org.Hs.eg.db)

head(ck)
#fun = One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .
dotplot(ck,showCategory = NULL, title = "KEGG pathway enrichment")

ppi <- 150
png("cell_lines_kegg_compare_all.png", width=10*ppi, height=12*ppi, res=ppi)
dotplot(ck, showCategory = NULL,title = "KEGG pathway enrichment")
dev.off()

ego <- enrichGO(gene          = gene,
                #universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

barplot(ego, showCategory=20,title = "GO term enrichment Biological Pathways")


kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk,20)

barplot(kk)

browseKEGG(kk, 'hsa04923')


geneList = xx[,2]
names(geneList) = as.character(xx[,8])
geneList = sort(geneList, decreasing = TRUE)
head(geneList)



ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

head(ego3,25)

ridgeplot(ego3)


#library(magrittr)



wp2gene <- read.gmt("~/Library/Mobile Documents/com~apple~CloudDocs/work/resources/wikipathways-20220310-gmt-Homo_sapiens.gmt")

head(wp2gene)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


head(wpid2gene)

ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp)


dotplot(ewp)



ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)

head(ewp2)

ridgeplot(ewp2)



p= ridgeplot(ewp2)
p + labs(title="GSEA of hWAT cell line D3 vs D8 (Wikipathways)")
ggsave("gsea_hWAT_d3_d8.png",width = 40, height = 20, units = "cm")
