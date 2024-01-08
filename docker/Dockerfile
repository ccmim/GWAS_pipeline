FROM ubuntu:20.04

RUN apt-get update &&  apt-get install -y \
    curl \
    ca-certificates \
    vim \
    sudo \
    git \
    bzip2 \
    libx11-6 \
    wget \
    zip tar unzip \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update
RUN apt-get install -y gcc
RUN apt-get install -y tmux
RUN apt-get install -y make
RUN apt-get install -y --reinstall build-essential
RUN apt-get install -y libboost-dev

RUN TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update -y && \
    apt-get install -y build-essential libfuse-dev libcurl4-openssl-dev libxml2-dev pkg-config libssl-dev mime-support automake libtool 
RUN apt-get install lsb-release -y

RUN apt-get update
RUN apt-get install ffmpeg libsm6 libxext6  -y


RUN apt-get update 
RUN apt-get install r-base -y

# Create a working directory.
RUN mkdir /app
WORKDIR /app

# Create a non-root user and switch to it.
RUN adduser --disabled-password --gecos '' --shell /bin/bash user \
 && chown -R user:user /app
RUN echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-user
USER user

# All users can use /home/user as their home directory.
ENV HOME=/home/user
RUN chmod 777 /home/user

# install Miniconda
RUN curl -so ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
 && chmod +x ~/miniconda.sh \
 && ~/miniconda.sh -b -p ~/miniconda \
 && rm ~/miniconda.sh
ENV PATH=/home/user/miniconda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false

# Create a Python 3.11 environment.
RUN /home/user/miniconda/bin/conda install conda-build
RUN /home/user/miniconda/bin/conda create -y --name gwas python=3.11
RUN /home/user/miniconda/bin/conda clean -ya

ENV CONDA_DEFAULT_ENV=gwas
ENV CONDA_PREFIX=/home/user/miniconda/envs/$CONDA_DEFAULT_ENV
ENV PATH=$CONDA_PREFIX/bin:$PATH

# RUN conda install pytorch=${TORCH} torchvision torchaudio cudatoolkit=11.3 -c pytorch
# Install python
RUN conda install -c conda-forge tqdm pandas
RUN conda install jupyter -c anaconda
RUN conda install -c anaconda pyyaml

USER root

# Install qctool
COPY qctool /app/software/qctool
WORKDIR /app/software/qctool
RUN ls -l
RUN ./waf configure
RUN ./waf

RUN ln -s /app/software/build/release/apps/qctool_v2.2.0 /bin/qctool

# Install GreedyRelated
WORKDIR /app/software
RUN git clone https://gitlab.com/choishingwan/GreedyRelated.git
RUN apt-get install cmake -y
WORKDIR /app/software/GreedyRelated
RUN cmake .
RUN make
RUN ln -s /app/software/GreedyRelated/bin/GreedyRelated /bin/GreedyRelated
RUN GreedyRelated -h

# Install Plink
WORKDIR /app/software/plink
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
RUN unzip plink_linux_x86_64_20231211.zip
RUN ln -s /app/software/plink/plink /bin/plink

# Install BGEN (bgenix, cat-bgen, etc.)
WORKDIR /app/software
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
RUN tar xzvf bgen.tgz
RUN mv bgen.tgz bgen
WORKDIR /app/software/bgen
RUN ./waf configure
RUN ./waf
RUN ln -s /app/software/bgen/build/apps/bgenix /bin/bgenix
RUN ln -s /app/software/bgen/build/apps/cat-bgen /bin/cat-bgen

# Install bgenie
COPY bgenie/bgenie_v1.3_static1 /app/software/bgenie
RUN ln -s /app/software/bgenie/bgenie_v1.3_static1 /bin/bgenie

# Install R and required packages
RUN conda install -c conda-forge r=4
RUN conda install -c r r-tidyverse
RUN conda install -c r r-argparse
RUN conda install -c conda-forge r-rcpp
RUN conda install -c conda-forge r-gprofiler2
RUN conda install -c bioconda bioconductor-biomart
RUN conda install -c bioconda r-qqman

RUN conda install -c conda-forge r-logging
RUN conda install -c conda-forge r-reader
