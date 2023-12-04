FROM ubuntu:22.04 
WORKDIR /home

RUN apt-get update -y \
	&& apt-get install -y --no-install-recommends \
	build-essential \
	software-properties-common \
	apt-utils \
  unzip \
  wget \
	&& wget -O subread.tar.gz https://downloads.sourceforge.net/project/subread/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz \
	&& tar zxvf subread.tar.gz \
	&& rm subread.tar.gz \
	&& ln -sf /home/subread-2.0.0-Linux-x86_64/bin/* /usr/local/bin/

