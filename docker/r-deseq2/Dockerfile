### start from bioconductor/bioconductor_docker
FROM bioconductor/bioconductor_docker

### system libraries
### Try to only install system libraries you actually need
### Package Manager is a good resource to help discover system deps
### RUN apt-get update && apt-get install -y

### install R packages required
### Change the packages list to suit your needs


RUN Rscript --vanilla -e 'install.packages("BiocManager",\
version = "3.16",repos="https://cloud.r-project.org/")'

RUN R -e 'BiocManager::install("DESeq2", type = "source", force = TRUE)'