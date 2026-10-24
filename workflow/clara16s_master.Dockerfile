# ================================================
# CLARA16s Master Dockerfile
# ARM64 / Raspberry Pi 5
# Ubuntu 22.04
# ================================================

FROM ubuntu:22.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
#ARG TARGETARCH

# --- Base system packages ---
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        git-all \
        make \
        gcc \
        g++ \
        wget \
        curl \
        ca-certificates \
        libssl3 \
        libssl-dev \
        libcurl4 \
        libc6-dev \
        zlib1g-dev \
        xz-utils \
        zip \
        unzip \
        bc \
        libbz2-dev \
        liblzma-dev \
        python3-dev \
        python3-pip \
        npm \
        fastqc \
        cutadapt \
        default-jre \
        && \
    rm -rf /var/lib/apt/lists/*


# --- Python dependencies ---
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


# =================================================
# Miniforge
# =================================================

WORKDIR /CLARA16s_tools

RUN wget  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
        -O /tmp/miniforge.sh && \
    chmod +x /tmp/miniforge.sh && \
    /tmp/miniforge.sh -b -p /miniforge && \
    rm /tmp/miniforge.sh

ENV PATH="/miniforge/bin:${PATH}"

RUN conda --version && \
    python --version && \
    conda info

# =================================================
# Conda configuration
# =================================================

ENV CONDA_ACCEPT_CHANNELS=true
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda config --system --add channels conda-forge && \
    conda config --system --add channels bioconda && \
    conda clean --all --yes


# =================================================
# Krona
# =================================================

WORKDIR /krona_tool

RUN git clone https://github.com/marbl/Krona.git && \
    cd Krona/KronaTools && \
    ./install.pl --prefix /opt/krona

ENV PATH="/opt/krona/bin:${PATH}"


# =================================================
# CLARA16s report environment
# =================================================

RUN conda create -y \
        -n clara_report \
        r-base \
        r-tidyverse \
        r-plotly \
        r-dt \
        r-vegan \
        r-dplyr \
        r-tidyr \
        r-ggplot2 \
        quarto \
        bioconductor-phyloseq \
        -c conda-forge \
        -c bioconda && \
    conda clean --all --yes


# =================================================
# Additional Python tools
# =================================================

RUN pip3 install --no-cache-dir \
        biopython \
        pandas \
        multiqc


# =================================================
# Copy CLARA16s resources
# =================================================

ADD resources/*.py /usr/bin/
ADD resources/*.sh /usr/bin/

WORKDIR /data
