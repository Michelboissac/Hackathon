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
library("ggpubr")



#1#################Creation des tableau cts et coldata necessaire pour deseq2
#dir = "/home/kngadi/Documents/outputfeaturecounts_inputDeseq2(-s 1)/"
#setwd(dir = dir)

count_matrix <- function(){
  #set working direcroty 
  #setwd(dir = dir)
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
  #Ordonner les sampples selon le tableau de reference pour comparaison
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
  cts <- na.omit(cts)
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = sample_infos,
                                design = ~ Staphylococcus_aureus
  )
  dds$Staphylococcus_aureus <- relevel(dds$Staphylococcus_aureus, ref = "untreated")
  # Exécuter DEsq
  dds <- DESeq(dds)
  #Explorer les résultats
  res <-results(dds, alpha = thrsd)
  return(list(res=res,dds=dds))
  
}

dds <- dsq_objt$dds
res<- dsq_objt$res
  # Plot Ma-plot
png(file = "plotMA.png", width = 800, height = 700)
plotMA(res, alpha=0.05, colNonSig="black", colSig="red")
dev.off()

ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))


sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# clustrering hierachique
png(file = "pheatmap.png", width = 800, height = 700)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

png(file = "plotPCA.png", width = 800, height = 700)
plotPCA(vsd, intgroup=c("Staphylococcus_aureus"))
dev.off()
