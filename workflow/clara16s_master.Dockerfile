# ================================================
# DATMA16S Integrated Dockerfile (Ubuntu 22.04)
# ================================================
FROM ubuntu:22.04 AS build

ARG DEBIAN_FRONTEND=noninteractive

# --- Base system packages ---
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        git-all \
        make \
        gcc \
        g++ \
        wget \
        libssl3 \
        libssl-dev \
        libcurl4 \
        ca-certificates \
        libc6-dev \
        zlib1g-dev \
        xz-utils \
        zip \
        unzip \
        bc \
        libbz2-dev \
        liblzma-dev \
        curl \
        python3-dev \
        python3-pip \
        npm \
        cargo \
        && rm -rf /var/lib/apt/lists/*

# --- Python dependencies (global) ---
RUN pip3 install --no-cache-dir \
    numpy \
    joblib \
    HTSeq \
    pysam \
    scikit-learn \
    scipy \
    six \
    pandas \
    biopython

RUN apt update
RUN apt install \
    fastqc \
    cutadapt \
    unzip \
    default-jre    
# --- Frontend dependencies ---
#ADD DATMA16s_frontend /DATMA16s_frontend
#WORKDIR /DATMA16s_frontend
#RUN npm install --force


# --- Miniconda ---
WORKDIR /DATMA16s_tools
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    chmod +x /tmp/miniconda.sh && \
    /tmp/miniconda.sh -b -p /miniconda && \
    rm /tmp/miniconda.sh
ENV PATH="/miniconda/bin:$PATH"

# --- QIIME2 ---
ENV CONDA_ACCEPT_CHANNELS=true
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda env create -n qiime2-env \
   --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml && \
    conda clean --all -y
RUN mkdir -p /.config /.cache && chmod -R a+rw /.config /.cache

# RUN conda config --add channels conda-forge
# RUN conda config --add channels bioconda
# RUN conda config --add channels defaults
# RUN conda create -n clara_report \
#     fastqc \
#     multiqc \
#     cutadapt \
#     pigz \
#     unzip \
#     python=3.11 \
#     -c conda-forge \
#     -c bioconda

#WORKDIR /krona_tool
#RUN git clone https://github.com/marbl/Krona.git \
#    && cd Krona/KronaTools \
#    && ./install.pl --prefix /opt/krona
#ENV PATH="/opt/krona/bin:${PATH}"

# ENV CONDA_ACCEPT_CHANNELS=true
# RUN conda config --add channels conda-forge
# RUN conda config --add channels bioconda
# RUN conda config --add channels defaults
# RUN conda create -n clara_report \
#     r-base \
#     r-tidyverse \
#     r-plotly \
#     r-dt \
#     r-vegan \
#     r-dplyr \
#     r-tidyr \
#     r-ggplot2 \
#     r-vegan  \
#     quarto \
#     bioconductor-phyloseq \
#     -c conda-forge \
#     -c bioconda 
    
RUN pip3 install biopython
RUN pip3 install pandas

# --- Copy project resources ---
ADD resources/*.py /usr/bin/
ADD resources/*.sh /usr/bin/

WORKDIR /data

