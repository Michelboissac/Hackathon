FROM ubuntu:22.04 

RUN apt-get update --fix-missing \
&& apt-get install -y wget tar gzip build-essential unzip make gcc libbz2-dev zlib1g zlib1g-dev liblzma5 liblzma-dev libncurses5 libncurses5-dev bzip2


RUN cd /usr/local/ \
&& wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
&& tar -xjvf samtools.tar.bz2 \
&& rm -rf samtools.tar.bz2 

RUN cd /usr/local/samtools-1.9 \
&& ./configure \
&& make \
&& make install \
&& cd /usr/local \
&& rm -rf /usr/local/samtools-1.9




RUN wget --no-check-certificate -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download \
&& unzip bowtie.zip -d /opt/ \
&& ln -s /opt/bowtie-0.12.7/ /opt/bowtie \
&& rm bowtie.zip 


ENV PATH $PATH:/opt/bowtie-0.12.7/