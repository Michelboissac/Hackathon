
params.publishDir = './output_dir2'
output_dir1 = file("output_dir1")

process createDir1 {

    script:
    """
    mkdir -p ${output_dir1}
    """
}

        //PROCESSUS

process trimgalore {
    input:
    path x
    output:
    path "${x.baseName}_trimmed.fq"

    script:
    """
    trim_galore -q 20 --phred33 --length 25 ${x}
    """
}

process genome_index {
    input:
    file f
    file g
    output:
    path "*"

    script:
    """
    bowtie-build ${f} ${g}
    """
}
process mapping {
    input:
    file trimmed
    file ref
    output:
    path "*bam"
    script:
    """
    bowtie -p 16 -S reference.gff <(cat ${trimmed}) | \\
    samtools sort -@ 16 -o ${trimmed.baseName}.bam -
    samtools index ${trimmed.baseName}.bam
    """
}
process featurecount {
        input :
        file x
        file ref
        output:
        path "*.txt"
        script:
        """
        featureCounts --extraAttributes Name -t gene -g ID -F GTF -T 16 -a ${ref} -o ${x}.count.txt ${x}

        cp *.txt ${output_dir1}
         """
}


process deseq2 {

    publishDir(
        path: "${params.publishDir}/",
        mode: 'copy',
        overwrite: false
    )



    input:

        file x


        output:

        path "*.png"
    script:
    """


#!Rscript
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
tab1=read.table("SRR10379721_trimmed.bam.count.txt",head = TRUE)
tab2=read.table("SRR10379722_trimmed.bam.count.txt",head = TRUE)
tab3=read.table("SRR10379723_trimmed.bam.count.txt",head = TRUE)
tab4=read.table("SRR10379724_trimmed.bam.count.txt",head = TRUE)
tab5=read.table("SRR10379725_trimmed.bam.count.txt",head = TRUE)
tab6=read.table("SRR10379726_trimmed.bam.count.txt",head = TRUE)

# Affichage de la data frame


cts=data.frame(
        tab1[["SRR10379721_trimmed.bam"]],
        tab2[["SRR10379722_trimmed.bam"]],
        tab3[["SRR10379723_trimmed.bam"]],
        tab4[["SRR10379724_trimmed.bam"]],
        tab5[["SRR10379725_trimmed.bam"]],
        tab6[["SRR10379726_trimmed.bam"]])


x=c("tab1.SRR10379721_trimmed.bam",
    "tab2.SRR10379722_trimmed.bam",
    "tab3.SRR10379723_trimmed.bam",
    "tab4.SRR10379724_trimmed.bam",
    "tab5.SRR10379725_trimmed.bam",
    "tab6.SRR10379726_trimmed.bam"
)


rownames(cts)=tab1[["Name"]]
colnames(cts)=x

coldata <- data.frame(
  condition=c("treated","treated","treated","untreated","untreated","untreated"),
  type = c("1", "2", "3", "1", "2", "3")
)
rownames(coldata) <- x





write.csv(coldata, "mini_dataframe.csv", row.names = TRUE)

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
dds[["condition"]] <- relevel(dds[["condition"]], ref = "untreated")
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
           gene = which.min(res[["padj"]]),
           intgroup = "condition")
dev.off()

        """


}

workflow {
        //CREATION REPERTOIRE OUTPUT
        createDir1()

        geno = genome_index(file("reference.fasta"),file("reference.gff"))



        def fastq = Channel.fromPath('*.fastq')
        trimq = trimgalore(fastq)
        mapped = mapping(trimq,geno)
        count = featurecount(mapped,file("reference.gff"))

        count_grouped = count.collect()
        Deseq2 = deseq2(count_grouped)
        Deseq2.view { "Result: ${it}" }



}



