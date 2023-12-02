# Packages et librairies nécessaires
 (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # Installation de readr si besoin 
  
BiocManager::install(version = "3.10")
BiocManager::install("EnrichmentBrowser")
BiocManager::install('KEGGgraph')
BiocManager::install()
if (!require("readr", quietly = TRUE))
  install.packages("readr")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("plyr", quietly = TRUE))
  install.packages("plyr")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("downloader", quietly = TRUE))
  install.packages("downloader")

if (!require("jsonlite", quietly = TRUE))
  install.packages("jsonlite")

if (!require("RJSONIO", quietly = TRUE))
install.packages("RJSONIO")

BiocManager::install("DESeq2")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")


BiocManager::install("DESeq2")
if (!require("tibble", quietly = TRUE))
  install.packages("tibble")
install.packages("ggmaplot")
if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr")
install.packages("nloptr")
install.packages("lme4")
install.packages("pbkrtest")
install.packages("car")

install.packages("systemfonts")
# Appel des librairies installé
library("RJSONIO")
library("EnrichmentBrowser")
library("KEGGgraph")
library("DESeq2")
library("tibble")
library("readr")
library("dplyr")
library("plyr")
library("pheatmap")
library("RColorBrewer")
library("EnrichmentBrowser")
library("utils")
library("ggpubr")
library("downloader")
library("jsonlite")
library("tidyverse")
############################## graphic test #######################


#################################### Construction de la matrice de comptage:  ##################
# Créer une chaine de caractères avec les fichier à extension .txt selon le repertoire de travail

count_matrix <- function(dir){
  #set working direcroty 
  setwd(dir = dir)
  # # Créer une chaine de caractères avec les fichier à extension .txt selon le repertoire de travail
  files_txt <- list.files( pattern = "\\.txt$")
  ## stocker les dataframe de samples dans une liste en gardant que le Geneid et les featurecounts
  data_list <- lapply(files_txt, function(x)read.table(x, skip = 1)[, c("V1", "V8")])
  # Fusionner par le Geneid
  df <- plyr::join_all(data_list, by= "V1", type = "inner") 
  # Remplacer la première colonnes par la première ligne
  names(df) <- as.matrix(df[1, ]) # selectionner la première colonne
  df <- df[-1, ] # supprimer celle-ci
  # Convertir les colonnes dans leur type appropriés
  df[] <- lapply(df, function(x) type.convert(as.character(x))) 
  rownames(df) <- NULL
  # Construction de la matrice de comptage 
  cts <- data.frame(df, row.names = 1) # Remplacer l'index par l'identifiants des gênes pour obtenir la matrice de comptage
  cts <- with(data = cts, expr = cts[,order(c(4,5,6,1,2,3))])
  colnames(cts) <- c("ctrl4", "ctrl5","ctrl6", "IP1", "IP2", "IP3")
  return(cts)
}
Deseq2_object <- function(cts, thrsd){
  # Le facteur de design
  smple_trimmed <- colnames(cts)
  sample_infos = data.frame(annot=smple_trimmed,
                            Staphylococcus_aureus=as.factor(c( "treated", 
                                                               "treated", 
                                                               "treated", 
                                                               "untreated",
                                                               "untreated", 
                                                               "untreated")),
                            row.names = 1
                            ) 
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = sample_infos,
                                design = ~ Staphylococcus_aureus
                                )
  dds$Staphylococcus_aureus <- relevel(dds$Staphylococcus_aureus, ref = "untreated")
  # Exécuter DEsq
  dds <- DESeq(dds)
  #Explorer les résultats
  res <-results(dds, alpha = thrsd)
  return(res,dds)
  
}
visualisation <- function(){
  
}
url = "https://www.genome.jp/kegg-bin/download_htext?htext=sao00001.keg&format=json"
destfile = "/home/kngadi/Documents"
download.file(url = url, destfile=destfile)
res.05 <- Deseq2_object(cts = cts, thrsd = 0.05)
summary(res.05)

cts <- count_matrix(dir = "/home/kngadi/Documents/outputfeaturecounts_inputDeseq2(-s 1)")


trnsl_genid = read.csv(file = "Gene_Names_1col_true.csv") # ID  de gènes impliqués dans la tradutions
trnsl_genid = paste0('gene-', trnsl_genid$Name)
trnsl_cts = cts[c(trnsl_genid),]


dds_trnsl$Staphylococcus_aureus <- relevel(dds_trnsl$Staphylococcus_aureus,
                                           ref = "untreated")
# Exécuter DEsq

dds_trnsl <- DESeq(dds_trnsl)
res_trnsl05 <- results(dds_trnsl, alpha = 0.05)
res <- results(dds)
res05 <-results(dds, alpha = 0.05)
# Explore results
summary(res05)
summary(object = res_trnsl05)
# MA plot
# Default plot
ggmaplot(res05,
         fdr = 0.05, size = 0.4,
         genenames = NULL,
         select.top.method = c("padj"),
         ggtheme = ggplot2::theme_minimal(),
         top = 0
         )
ggmaplot(res.05,
         fdr = 0.05, fc = 1.5, size = 0.4,
         genenames = NULL,
         ggtheme = ggplot2::theme_minimal(),
         top = 0
)

plotMA(res05, alpha= 0.05, colSig="red")
plotMA(res_trnsl05, alpha= 0.05, colSig="red")
plotCounts(dds = dds,
           gene = which.min(res$padj),
           intgroup = "Staphylococcus_aureus")
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Staphylococcus_aureus, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Staphylococcus_aureus"))
