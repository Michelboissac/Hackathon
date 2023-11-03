FROM ubuntu:22.04 
WORKDIR /opt/

RUN apt-get update --fix-missing \
&& apt-get install -y wget

RUN wget --no-check-certificate https://sourceforge.net/projects/subread/files/subread-1.4.6/subread-1.4.6-Linux-x86_64.tar.gz/download 

RUN tar -xzf download && \
    ln -s /opt/subread-1.4.6-Linux-x86_64/bin/featureCounts /usr/local/bin/

RUN rm download
