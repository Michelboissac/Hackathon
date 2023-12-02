  #1)creation environnement R avec 3 librairires
#conda create -n r_env r-essentials r-base
#conda activate r_env

  #2)puis ouvrir R en ligne de commande en tapant R
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install.packages("pheatmap")
#install.packages("RColorBrewer")


library("DESeq2")
library("pheatmap")
library("RColorBrewer")


#1#################Creation des tableau cts et coldata necessaire pour deseq2
tab1=read.table("output_dir/SRR10379721_trimmed.bam.count.txt",head = TRUE)
tab2=read.table("output_dir/SRR10379722_trimmed.bam.count.txt",head = TRUE)
tab3=read.table("output_dir/SRR10379723_trimmed.bam.count.txt",head = TRUE)
tab4=read.table("output_dir/SRR10379724_trimmed.bam.count.txt",head = TRUE)
tab5=read.table("output_dir/SRR10379725_trimmed.bam.count.txt",head = TRUE)
tab6=read.table("output_dir/SRR10379726_trimmed.bam.count.txt",head = TRUE)

cts=data.frame(tab1$SRR10379721_trimmed.bam,
               tab2$SRR10379722_trimmed.bam,
               tab3$SRR10379723_trimmed.bam,
               tab4$SRR10379724_trimmed.bam,
               tab5$SRR10379725_trimmed.bam,
               tab6$SRR10379726_trimmed.bam)

rownames(cts)=tab1$Name
x=c("tab1.SRR10379721_trimmed.bam",
    "tab2.SRR10379722_trimmed.bam",
    "tab3.SRR10379723_trimmed.bam",
    "tab4.SRR10379724_trimmed.bam",
    "tab5.SRR10379725_trimmed.bam",
    "tab6.SRR10379726_trimmed.bam"
)

coldata=data.frame(condition = x)
rownames(coldata)=coldata$condition
coldata$condition=c("treated","untreated")
coldata$type=c("1","2","3","1","2","3")
colnames(coldata)=c("condition","type")

coldata$condition = factor(coldata$condition)
coldata$type = factor(coldata$type)


#2#################################Analyse Deseq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ type + condition)


dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
summary(res)
png(file = "plotMA.png", width = 800, height = 700)
plotMA(res, ylim=c(-2,2))
dev.off()


keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
# Spécifier le facteur de niveau
dds$condition <- relevel(dds$condition, ref = "untreated")
# Exécuter DEsq
dds <- DESeq(dds)
res <- results(dds)
res05 <-results(dds, alpha = 0.05)
# Explore results
summary(res05)
# MA plot

png(file = "plotres05.png", width = 800, height = 700)
plotMA(res05)
dev.off()


png(file = "plotCounts.png", width = 800, height = 700)
plotCounts(dds = dds,
           gene = which.min(res$padj),
           intgroup = "condition")
dev.off()

ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(file = "pheatmap.png", width = 800, height = 700)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


png(file = "plotPCA.png", width = 800, height = 700)
plotPCA(vsd, intgroup=c("condition"))
dev.off()

#############################################################################
