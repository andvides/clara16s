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

# --- Frontend dependencies ---
#ADD DATMA16s_frontend /DATMA16s_frontend
#WORKDIR /DATMA16s_frontend
#RUN npm install --force

# --- HyperEx ---
RUN cargo install hyperex --locked && cp /root/.cargo/bin/hyperex /usr/bin/

# --- selectFasta ---
WORKDIR /DATMA16s_tools
RUN wget https://github.com/andvides/selectFasta/releases/download/v3.1/selectFasta_v3.1.zip && \
    unzip selectFasta_v3.1.zip && \
    make -C selectFasta_v3.1 && \
    cp selectFasta_v3.1/selectFasta /usr/bin/

# --- seqkit ---
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz && \
    tar -zxvf seqkit_linux_amd64.tar.gz && \
    cp seqkit /usr/bin/

# --- ART simulator ---
#WORKDIR /
#ADD resources/art_bin_MountRainier /art
#ENV PATH="/art:${PATH}"

# --- Miniconda ---
WORKDIR /DATMA16s_tools
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    chmod +x /tmp/miniconda.sh && \
    /tmp/miniconda.sh -b -p /miniconda && \
    rm /tmp/miniconda.sh
ENV PATH="/miniconda/bin:$PATH"
#ENV LD_LIBRARY_PATH="/miniconda/lib:/miniconda/pkgs/libcurl-8.11.1-hc9e6f67_0/lib:$LD_LIBRARY_PATH"

# --- QIIME2 ---
ENV CONDA_ACCEPT_CHANNELS=true
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda env create -n qiime2-env \
   --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml && \
    conda clean --all -y
RUN mkdir -p /.config /.cache && chmod -R a+rw /.config /.cache


# --- MHASS ---
#RUN apt-get update && apt-get install -y --no-install-recommends \
#    libp11-kit0 libffi8 && \
#    rm -rf /var/lib/apt/lists/*

WORKDIR /DATMA16s_tools
#RUN git clone https://github.com/rhowardstone/MHASS.git
#WORKDIR /DATMA16s_tools/MHASS
#RUN bash install_dependencies.sh 
#ENV PATH="/DATMA16s_tools/MHASS:$PATH"

#RUN conda run -n mhass Rscript -e '\
#    packages_needed <- c("optparse", "BiocManager"); \
#    packages_to_install <- packages_needed[!sapply(packages_needed, requireNamespace, quietly = TRUE)]; \
#    if(length(packages_to_install) > 0) install.packages(packages_to_install, repos="https://cloud.r-project.org/")'

#RUN conda install -n mhass -c bioconda -y samtools

RUN pip3 install biopython
RUN pip3 install pandas

# --- Copy project resources ---
ADD resources/*.py /usr/bin/
ADD resources/*.sh /usr/bin/

WORKDIR /data
