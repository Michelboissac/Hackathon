FROM bioconductor/bioconductor_docker:latest

RUN apt-get -y update \
    && apt-get install -y software-properties-common \
    && apt-get -y update \
    && apt-get -y update --fix-missing \
    && apt-get install -y \
        build-essential \
        libssh2-1-dev \
		libxml2-dev \
		libcairo2-dev \
		libxt-dev \
		libssl-dev \
		libv8-dev \
        libmagick++-dev \
        libz-dev \
        libcurl4-openssl-dev \
        r-base \
        r-base-dev \
        r-cran-rcpp \
        openjdk-11-jdk \
        littler \
        python3 \
        python3-dev \
        python3-matplotlib \
        make \
    && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
    && apt-get clean

RUN install.r \
        devtools \
        matrixStats\
        methods \
        stats4 \
        locfit \
        rcpp \
        dplyr\
        plyr\
        docopt \
        lattice \
        RCurl \
        openssl \
        Hmisc \
        tidyr \
        readr \
        stringr \
        ggplot2 \
        EnrichmentBrowser \
        tibble \
        pheatmap \
        ggplot2 \
        RColorBrewer \
        BiocManager\ 
        && rm -rf /var/lib/apt/lists/*

#install other required R packages from github
RUN Rscript -e "\
BiocManager::install('DESeq2'); \
BiocManager::install('GenomicFeatures'); \
BiocManager::install('AnnotationDbi'); \
BiocManager::install('pheatmap'); \
BiocManager::install('RColorBrewer'); \
BiocManager::install('MatrixGenerics'); \
BiocManager::install('Biobase'); \
BiocManager::install('BiocParallel')"