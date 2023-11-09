# Packages et librairies nécessaires
 (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # Installation de readr si besoin 
  
BiocManager::install(version = "3.10")

if (!require("readr", quietly = TRUE))
  install.packages("readr")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("plyr", quietly = TRUE))
  install.packages("plyr")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

BiocManager::install("DESeq2")

BiocManager::install("DESeq2")
if (!require("tibble", quietly = TRUE))
  install.packages("tibble")
# Appel des librairies installé
library("ggplot2")
library("DESeq2")
library("tibble")
library("pasilla")
library("tidyverse")
library("readr")
library("dplyr")
library("plyr")

#################################### Our Data: FeatureCounts output ##################

files_txt <- list.files(pattern = "\\.txt$") # Créer une chaine de caractères avec les fichier à extrension .txt selon le repertoire de travail
data_list <- lapply(files_txt, function(x)read.table(x, skip = 1)[, c("V1", "V8")])# stocker les dataframe de samples dans une liste en gardant que le Geneid et les featurecounts
df <- plyr::join_all(data_list, by= "V1", type = "inner") # Fusionner par le Geneid

#head(df)
# Remplacer la première colonnes par la première ligne
names(df) <- as.matrix(df[1, ]) # selectionner la première colonne
df <- df[-1, ] # supprimer celle-ci
df[] <- lapply(df, function(x) type.convert(as.character(x))) # Convertir les colonnes dans leur type appropriés
rownames(df) <- NULL
# Construction de la matrice de comptage 
cts <- data.frame(df, row.names = 1) # Remplacer l'index par l'identifiants des gênes pour obtenir la matrice de comptage
smple_trimmed <- colnames(cts)
sample_infos = data.frame(annot=smple_trimmed,
                          Staphylococcus_aureus=as.factor(c("untreated","untreated", "untreated", "treated", "treated", "treated")),
                          row.names = 1) # Le factuer de design
# Construire l'objet DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sample_infos,
                              design = ~ Staphylococcus_aureus)
# Supprimer les lignes avec de faibles nombre de gênes
# Garder les lignes qui ont au moins 10 lectures au total
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
# Spécifier le facteur de niveau
dds$Staphylococcus_aureus <- relevel(dds$Staphylococcus_aureus, ref = "untreated")
# Exécuter DEsq
dds <- DESeq(dds)
res <- results(dds)
res05 <-results(dds, alpha = 0.05)
# Explore results
summary(res)
# MA plot
plotMA(res)
