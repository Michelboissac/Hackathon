FROM ubuntu:22.04 

RUN apt-get update --fix-missing \
&& apt-get install -y wget unzip

RUN wget --no-check-certificate -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download \
&& unzip bowtie.zip -d /opt/ \
&& ln -s /opt/bowtie-0.12.7/ /opt/bowtie \
&& rm bowtie.zip 

ENV PATH $PATH:/opt/bowtie-0.12.7/  